"""
Zenodo QC - quality check of the marine datasets discovered by the Zenodo Inventory (Task 1).

Task 2 of the "Zenodo Inventory" project (EMODnet Biology WP4, Task 7.2).

The checks mirror the EMODnet Biology data quality criteria, as implemented by the
LifeWatch & EMODnet Biology QC tool and by the EMODnetBiocheck R package:

  Level A - metadata QC, for every inventory record (Zenodo landing page + API):
            open access, reusable licence, persistent identifier, contact, temporal
            coverage, tabular data files, inventory status, documentation.
  Level B - data QC, for open records with tabular files: the files are downloaded,
            read (CSV / TSV / TXT / XLSX / XLS / ODS / ZIP / Darwin Core Archive),
            their columns mapped to Darwin Core terms and the content checked:
            required fields, coordinates (range, region, on land via OBIS xylookup),
            dates (ISO 8601), depth vs bathymetry, taxonomy vs WoRMS, controlled
            vocabularies, unique identifiers, measurement units.

Usage:
    python zenodo_qc.py --inventory data.json --out output
    python zenodo_qc.py --ids 8094440,7654321 --force
    python zenodo_qc.py --skip-data          # metadata QC only
"""

import argparse
import csv
import io
import json
import os
import re
import threading
import time
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import warnings

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

TOOL_VERSION = 'zenodo_qc v1.0'

ZENODO_API = 'https://zenodo.org/api/records/'
XYLOOKUP_API = 'https://api.obis.org/xylookup'
WORMS_API = 'https://www.marinespecies.org/rest/AphiaRecordsByMatchNames'

ZENODO_TOKEN = os.environ.get('ZENODO_TOKEN', '')

warnings.filterwarnings('ignore', category=UserWarning)

MAX_FILE_MB = 30            # files above this size are skipped
MAX_FILES_PER_RECORD = 8    # tabular files read per record
MAX_ROWS = 300_000          # rows read per table
MAX_WORMS_NAMES = 200       # unique names sent to WoRMS per dataset (sampled above)
MAX_XY_POINTS = 3000        # unique coordinates sent to xylookup per dataset

# Documentation / helper files are read last (data files first).
HELPER_FILES = (r'readme|licen|citation|md5|checksum|glossary|mapping|workflow|dictionary|codebook|'
                r'sessioninfo|metadata|description|filter')
TABULAR_EXT = {'.csv', '.tsv', '.txt', '.tab', '.xlsx', '.xlsm', '.xls', '.ods'}
ARCHIVE_EXT = {'.zip'}

# ---------------------------------------------------------------------------
# Criteria
# ---------------------------------------------------------------------------
CRITERIA = [
    # Level A - metadata
    ('M01', 'Open access', 'metadata', 'error',
     'Files are openly accessible on Zenodo (not embargoed, restricted or closed).'),
    ('M02', 'Reusable licence', 'metadata', 'error',
     'A licence is declared. CC0 / CC-BY pass; NC, ND, SA or non-data licences are flagged.'),
    ('M03', 'Persistent identifier', 'metadata', 'error',
     'The dataset has a DOI that can be cited.'),
    ('M04', 'Contact / provider', 'metadata', 'warning',
     'A data provider can be contacted (e-mail in description, or creator with affiliation).'),
    ('M05', 'Temporal coverage', 'metadata', 'warning',
     'The sampling period can be derived from the metadata.'),
    ('M06', 'Tabular data files', 'metadata', 'error',
     'At least one machine-readable tabular file (CSV, TSV, TXT, XLS(X), ODS, ZIP/DwC-A).'),
    ('M07', 'Inventory relevance', 'metadata', 'warning',
     'Record was "Included" by the Task 1 marine-relevance screening (not "Manual review").'),
    ('M08', 'Documentation', 'metadata', 'warning',
     'Description of at least 300 characters and keywords present.'),
    # Level B - data
    ('D01', 'Readable data', 'data', 'error',
     'At least one tabular file could be downloaded and parsed into a table.'),
    ('D02', 'Core fields: when, where, what', 'data', 'error',
     'eventDate (or year), decimalLatitude, decimalLongitude and scientificName can be mapped.'),
    ('D03', 'Darwin Core identifiers & terms', 'data', 'warning',
     'occurrenceID, scientificNameID, occurrenceStatus and basisOfRecord are present.'),
    ('D04', 'Coordinate validity', 'data', 'error',
     'Latitude in [-90, 90], longitude in [-180, 180], no (0, 0), decimal degrees.'),
    ('D05', 'Coordinates in declared region', 'data', 'warning',
     'Points fall within the bounding box of the regions declared in the inventory.'),
    ('D06', 'Points at sea', 'data', 'error',
     'Points are not on land (OBIS xylookup shoredistance; >25 % on land fails).'),
    ('D07', 'Date validity', 'data', 'error',
     'Dates parse, are ISO 8601 (YYYY-MM-DD) and are not in the future.'),
    ('D08', 'Depth consistency', 'data', 'warning',
     'Depths are numeric, non-negative, min <= max and not deeper than the bathymetry.'),
    ('D09', 'Taxon names in WoRMS', 'data', 'error',
     'Scientific names match the World Register of Marine Species (>= 90 % pass, >= 60 % warn).'),
    ('D10', 'Marine / brackish taxa', 'data', 'warning',
     'Matched taxa are flagged marine or brackish in WoRMS.'),
    ('D11', 'scientificNameID (WoRMS LSID)', 'data', 'warning',
     'scientificNameID values are valid urn:lsid:marinespecies.org:taxname: LSIDs.'),
    ('D12', 'Controlled vocabularies', 'data', 'warning',
     'basisOfRecord and occurrenceStatus use the Darwin Core vocabularies.'),
    ('D13', 'Unique occurrenceID', 'data', 'error',
     'occurrenceID values are unique and not empty.'),
    ('D14', 'Measurements described', 'data', 'warning',
     'Measured values (abundance, biomass, ...) have a type and a unit (eMoF-ready).'),
    ('D15', 'Long (tidy) layout', 'data', 'warning',
     'One record per row; species-by-sample matrices ("wide") need reshaping before DwC.'),
]
CRIT = {c[0]: dict(id=c[0], name=c[1], level=c[2], severity=c[3], description=c[4]) for c in CRITERIA}

# Approximate bounding boxes (lon_min, lat_min, lon_max, lat_max) of the inventory regions.
REGION_BBOX = {
    'Arctic': [(-180, 60, 180, 90)],
    'Atlantic': [(-100, -60, 25, 75)],
    'Baltic Sea': [(7.5, 53, 31, 66.5)],
    'Black Sea': [(26, 40, 42.5, 48)],
    'Mediterranean Sea': [(-6.5, 30, 37, 46.5)],
    'North Sea': [(-6, 48, 13, 62.5)],
}

# ---------------------------------------------------------------------------
# Column mapping: header synonyms -> Darwin Core terms
# ---------------------------------------------------------------------------
def _norm(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


SYNONYMS = {
    'decimalLatitude': ['decimallatitude', 'latitude', 'lat', 'latdec', 'latitudedecimal', 'latdd',
                        'ycoord', 'latitudedeg', 'latitudedegrees', 'latituden', 'latn', 'startlatitude', 'latstart', 'latbegin', 'latitudebegin',
                        'latitudestart', 'lat1', 'decimallat', 'latwgs84', 'latitudewgs84', 'y'],
    'decimalLongitude': ['decimallongitude', 'longitude', 'lon', 'long', 'lng', 'londec', 'longitudedecimal',
                         'londd', 'xcoord', 'longitudedeg', 'longitudedegrees', 'longitudee', 'lone',
                         'startlongitude', 'lonstart', 'lonbegin', 'longitudebegin', 'longitudestart', 'lon1', 'decimallon', 'lonwgs84',
                         'longitudewgs84', 'x'],
    'eventDate': ['eventdate', 'date', 'sampledate', 'samplingdate', 'datetime', 'collectiondate',
                  'dateofsampling', 'surveydate', 'observationdate', 'datecollected', 'eventdatetime',
                  'datum', 'fecha', 'data', 'dateutc', 'startdate', 'datestart', 'datebegin', 'sightingdate', 'timestamp'],
    'year': ['year', 'yr', 'sampleyear', 'samplingyear', 'yyyy', 'surveyyear', 'ano', 'jahr'],
    'month': ['month', 'mon', 'mm', 'samplemonth'],
    'day': ['day', 'dd', 'sampleday'],
    'scientificName': ['scientificname', 'species', 'speciesname', 'taxon', 'taxonname', 'taxa', 'scientific',
                       'latinname', 'sciname', 'speciesscientificname', 'validname', 'acceptedname',
                       'organism', 'organismname', 'name', 'taxonomicname', 'specie', 'especie', 'art',
                       'aphianame', 'identification', 'speciesid'],
    'scientificNameID': ['scientificnameid', 'aphiaid', 'lsid', 'wormsid', 'wormsaphiaid', 'taxonid', 'aphia'],
    'occurrenceID': ['occurrenceid'],
    'eventID': ['eventid', 'sampleid', 'sample', 'samplecode', 'samplename', 'haulid', 'stationid'],
    'locality': ['locality', 'station', 'site', 'stationname', 'sitename', 'location', 'locationname',
                 'stn', 'area', 'transect'],
    'basisOfRecord': ['basisofrecord'],
    'occurrenceStatus': ['occurrencestatus', 'presence', 'presenceabsence', 'pa'],
    'minimumDepthInMeters': ['minimumdepthinmeters', 'mindepth', 'depthmin', 'minimumdepth', 'depthfrom',
                             'upperdepth', 'depthupper', 'startdepth'],
    'maximumDepthInMeters': ['maximumdepthinmeters', 'maxdepth', 'depthmax', 'maximumdepth', 'depthto',
                             'lowerdepth', 'depthlower', 'enddepth'],
    'depth': ['depth', 'depthm', 'depthinmeters', 'waterdepth', 'sampledepth', 'samplingdepth', 'depthmeters',
              'bottomdepth', 'z'],
    'individualCount': ['individualcount', 'count', 'counts', 'number', 'numberofindividuals', 'n', 'nind',
                        'abundancecount', 'total', 'quantity'],
    'organismQuantity': ['organismquantity'],
    'measurementType': ['measurementtype'],
    'measurementValue': ['measurementvalue'],
    'measurementUnit': ['measurementunit', 'unit', 'units'],
}
PARTIAL_TERMS = [
    (r'^decimallatitude', 'decimalLatitude'), (r'^decimallongitude', 'decimalLongitude'),
    (r'^(accepted|verbatim|valid)?scientificname$', 'scientificName'), (r'^eventdate', 'eventDate'),
    (r'^(latitude|lat)(start|begin|from|1)', 'decimalLatitude'),
    (r'^(longitude|lon|long)(start|begin|from|1)', 'decimalLongitude'),
    (r'^(accepted)?aphiaid$', 'scientificNameID'),
]
SYN_LOOKUP = {}
for term, syns in SYNONYMS.items():
    for s in syns:
        SYN_LOOKUP.setdefault(s, term)

# Measured-value columns (prefix match on the normalised header).
MEASURE_PREFIXES = ('abundance', 'density', 'biomass', 'count', 'individualcount', 'cover', 'coverage',
                    'wetweight', 'dryweight', 'weight', 'cpue', 'frequency', 'percentcover', 'numberper',
                    'indm', 'cellsl', 'cells', 'nindm', 'organismquantity', 'measurementvalue')
UNIT_HINT = re.compile(r'(\(|\[|/|per\b|_m2|_m3|m-2|m-3|m\^|ind|cells|mg|µg|ug\b|\bg\b|kg|%|l-1|/l|percent)',
                       re.IGNORECASE)

BASIS_OF_RECORD = {'humanobservation', 'machineobservation', 'materialsample', 'preservedspecimen',
                   'livingspecimen', 'fossilspecimen', 'materialcitation', 'occurrence', 'observation',
                   'materialentity'}
OCCURRENCE_STATUS = {'present', 'absent'}

LSID_RE = re.compile(r'^urn:lsid:marinespecies\.org:taxname:\d+$')
ISO_DATE_RE = re.compile(
    r'^\d{4}(-\d{2}(-\d{2}([T ]\d{2}(:\d{2}(:\d{2}(\.\d+)?)?)?(Z|[+-]\d{2}(:?\d{2})?)?)?)?)?'
    r'(/\d{4}(-\d{2}(-\d{2}([T ]\d{2}(:\d{2}(:\d{2}(\.\d+)?)?)?(Z|[+-]\d{2}(:?\d{2})?)?)?)?)?)?$')
BINOMIAL_RE = re.compile(r'^[A-Z][a-z]+(\s|_|\.)+(\(?[A-Z]?[a-z]+\)?(\s|_|\.)*)?[a-z]+\.?$')
GENUS_RE = re.compile(r'^[A-Z][a-z]{3,}(\s+(sp|spp|cf|aff)\.?)?$')

# ---------------------------------------------------------------------------
# HTTP utilities (polite, thread-safe rate limiting)
# ---------------------------------------------------------------------------
_thread_local = threading.local()


def session():
    s = getattr(_thread_local, 'session', None)
    if s is None:
        s = requests.Session()
        retries = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
                        allowed_methods=['GET', 'POST'], raise_on_status=False,
                        respect_retry_after_header=True)
        s.mount('https://', HTTPAdapter(max_retries=retries))
        s.headers['User-Agent'] = 'zenodo_qc/1.0 (EMODnet Biology, NIMRD)'
        _thread_local.session = s
    return s


class RateLimiter:
    def __init__(self, interval):
        self.interval = interval
        self.lock = threading.Lock()
        self.last = 0.0

    def wait(self):
        with self.lock:
            delay = self.last + self.interval - time.time()
            if delay > 0:
                time.sleep(delay)
            self.last = time.time()


ZENODO_LIMIT = RateLimiter(1.0)
WORMS_LIMIT = RateLimiter(0.4)
OBIS_LIMIT = RateLimiter(0.4)


def zenodo_get(url, **kw):
    params = kw.pop('params', {}) or {}
    if ZENODO_TOKEN:
        params['access_token'] = ZENODO_TOKEN
    for attempt in range(4):
        ZENODO_LIMIT.wait()
        r = session().get(url, params=params, timeout=kw.pop('timeout', 120), **kw)
        if r.status_code == 429:
            time.sleep(15 * (attempt + 1))
            continue
        return r
    return r


# ---------------------------------------------------------------------------
# Zenodo record + file cache
# ---------------------------------------------------------------------------
def safe_name(name):
    return re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', name)[:180]


def get_record(rec_id, cache_dir):
    folder = os.path.join(cache_dir, rec_id)
    path = os.path.join(folder, '_record.json')
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    r = zenodo_get(ZENODO_API + rec_id, timeout=60)
    if r.status_code != 200:
        return None
    data = r.json()
    os.makedirs(folder, exist_ok=True)
    tmp = path + f'.{threading.get_ident()}.part'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f)
    os.replace(tmp, path)
    return data


def prefetch_records(ids, cache_dir, log):
    """Fetch record metadata in batches of 25 through the search API (one request per batch
    instead of one per record - much faster when the Zenodo API is slow)."""
    todo = [i for i in ids if not os.path.exists(os.path.join(cache_dir, i, '_record.json'))]
    if not todo:
        return
    log(f'Prefetching metadata of {len(todo)} records in batches of 25...')
    for b in range(0, len(todo), 25):
        batch = todo[b:b + 25]
        try:
            r = zenodo_get(ZENODO_API.rstrip('/'), params={'q': 'recid:(' + ' OR '.join(batch) + ')',
                                                          'size': 25, 'all_versions': 'true'}, timeout=180)
            hits = r.json().get('hits', {}).get('hits', []) if r.status_code == 200 else []
        except Exception as ex:
            log(f'   batch {b // 25 + 1}: {ex}')
            continue
        for h in hits:
            rid = str(h.get('id'))
            folder = os.path.join(cache_dir, rid)
            os.makedirs(folder, exist_ok=True)
            path = os.path.join(folder, '_record.json')
            tmp = path + f'.{threading.get_ident()}.part'
            with open(tmp, 'w', encoding='utf-8') as f:
                json.dump(h, f)
            os.replace(tmp, path)
        log(f'   batch {b // 25 + 1}/{(len(todo) + 24) // 25}: {len(hits)} records')


def download(url, dest):
    """Download-if-missing; .part + rename so concurrent processes never read half files."""
    if os.path.exists(dest):
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + f'.{os.getpid()}.{threading.get_ident()}.part'
    try:
        ZENODO_LIMIT.wait()
        params = {'access_token': ZENODO_TOKEN} if ZENODO_TOKEN else {}
        with session().get(url, params=params, stream=True, timeout=300) as r:
            if r.status_code != 200:
                return False
            with open(tmp, 'wb') as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        if not os.path.exists(dest):
            os.replace(tmp, dest)
        return True
    except Exception:
        return False
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def ext_of(name):
    name = name.lower()
    if name.endswith('.csv.gz') or name.endswith('.txt.gz'):
        return '.gz'
    return os.path.splitext(name)[1]


# ---------------------------------------------------------------------------
# Reading tables
# ---------------------------------------------------------------------------
def _decode(raw):
    for enc in ('utf-8-sig', 'cp1252', 'latin-1'):
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return raw.decode('latin-1', errors='replace'), 'latin-1'


def _unique_columns(df):
    seen = {}
    cols = []
    for c in df.columns:
        c = str(c).strip()
        if c in seen:
            seen[c] += 1
            c = f'{c}.{seen[c]}'
        else:
            seen[c] = 0
        cols.append(c)
    df.columns = cols
    return df


def _fix_header(df):
    return _unique_columns(_find_header(df))


def _find_header(df):
    """If the real header is not on the first row, find the row with the most text cells."""
    unnamed = sum(1 for c in df.columns if str(c).startswith('Unnamed') or str(c).strip() == '')
    if unnamed <= len(df.columns) / 2 or len(df) < 2:
        return df
    best, best_score = None, 0
    for i in range(min(15, len(df))):
        row = df.iloc[i]
        score = sum(1 for v in row if isinstance(v, str) and v.strip() and not _is_number(v))
        if score > best_score:
            best, best_score = i, score
    if best is None or best_score < 2:
        return df
    new = df.iloc[best + 1:].copy()
    new.columns = [str(v).strip() if pd.notna(v) else f'col{j}' for j, v in enumerate(df.iloc[best])]
    return new.reset_index(drop=True)


def _is_number(v):
    try:
        float(str(v).replace(',', '.'))
        return True
    except ValueError:
        return False


def read_text_table(path):
    with open(path, 'rb') as f:
        raw = f.read()
    text, _ = _decode(raw)
    sample = text[:65536]
    lines = [l for l in sample.splitlines() if l.strip()][:50]
    if len(lines) < 2:
        return None
    try:
        sep = csv.Sniffer().sniff('\n'.join(lines), delimiters=',;\t|').delimiter
    except csv.Error:
        counts = {d: sum(l.count(d) for l in lines) for d in [',', ';', '\t', '|']}
        sep = max(counts, key=counts.get)
        if counts[sep] < len(lines):
            # whitespace-separated?
            if all(len(l.split()) == len(lines[0].split()) > 1 for l in lines[:10]):
                sep = r'\s+'
            else:
                return None
    try:
        df = pd.read_csv(io.StringIO(text), sep=sep, nrows=MAX_ROWS, dtype=str,
                         engine='python' if sep == r'\s+' else 'c', on_bad_lines='skip',
                         low_memory=False)
    except Exception:
        try:
            df = pd.read_csv(io.StringIO(text), sep=sep, nrows=MAX_ROWS, dtype=str, engine='python',
                             on_bad_lines='skip')
        except Exception:
            return None
    if df.shape[1] < 2 or df.shape[0] < 1:
        return None
    return _fix_header(df)


def read_excel_tables(path):
    tables = []
    try:
        sheets = pd.read_excel(path, sheet_name=None, nrows=MAX_ROWS, dtype=object)
    except Exception:
        return tables
    for name, df in sheets.items():
        df = df.dropna(how='all').dropna(axis=1, how='all')
        if df.shape[1] >= 2 and df.shape[0] >= 1:
            tables.append((name, _fix_header(df)))
    return tables


def parse_meta_xml(xml_bytes):
    """Return {'core': (file, rowType, {index: term}, sep, skip), 'extensions': [...]}."""
    root = ET.fromstring(xml_bytes)
    out = {'core': None, 'extensions': []}
    for node in root:
        tag = node.tag.split('}')[-1]
        if tag not in ('core', 'extension'):
            continue
        loc = node.find('.//{*}location')
        fields = {}
        for fld in node.findall('{*}field'):
            if fld.get('index') is not None:
                fields[int(fld.get('index'))] = fld.get('term', '').rsplit('/', 1)[-1]
        idnode = node.find('{*}id') if tag == 'core' else node.find('{*}coreid')
        if idnode is not None and idnode.get('index') is not None:
            fields.setdefault(int(idnode.get('index')), 'id' if tag == 'core' else 'coreid')
        sep = node.get('fieldsTerminatedBy', '\t').encode().decode('unicode_escape')
        skip = int(node.get('ignoreHeaderLines', '0') or 0)
        item = (loc.text.strip() if loc is not None else None, node.get('rowType', ''), fields, sep, skip)
        if tag == 'core':
            out['core'] = item
        else:
            out['extensions'].append(item)
    return out


def read_dwca(zf):
    """Read a Darwin Core Archive from an open ZipFile; return dict of DataFrames."""
    names = {os.path.basename(n): n for n in zf.namelist()}
    meta = parse_meta_xml(zf.read(names['meta.xml']))
    frames = {}
    for kind, item in [('core', meta['core'])] + [('ext', e) for e in meta['extensions']]:
        if not item or not item[0] or item[0] not in names:
            continue
        loc, row_type, fields, sep, skip = item
        info = zf.getinfo(names[loc])
        if info.file_size > MAX_FILE_MB * 4 * 1024 * 1024:
            continue
        raw = zf.read(names[loc])
        text, _ = _decode(raw)
        df = pd.read_csv(io.StringIO(text), sep=sep, header=None, skiprows=skip, dtype=str,
                         nrows=MAX_ROWS, quoting=csv.QUOTE_NONE, on_bad_lines='skip')
        df.columns = [fields.get(i, f'col{i}') for i in range(df.shape[1])]
        short = row_type.rsplit('/', 1)[-1] or loc
        frames[short] = df
    return frames, meta


def load_tables(rec_id, record, cache_dir, log):
    """Download and read the tabular files of a record.

    Returns (tables, files_checked, skipped) where tables = list of dict(name, df, layout, source).
    """
    files = record.get('files') or []
    cand = []
    for f in files:
        key = f.get('key', '')
        e = ext_of(key)
        if e in TABULAR_EXT or e in ARCHIVE_EXT:
            cand.append(f)
    # README-like text files last, small files first
    cand.sort(key=lambda f: (bool(re.search(HELPER_FILES, f['key'], re.I)),
                             f.get('size', 0)))
    tables, checked, skipped = [], [], []
    for f in cand:
        if len(checked) >= MAX_FILES_PER_RECORD:
            break
        key, size = f['key'], f.get('size', 0) or 0
        if size > MAX_FILE_MB * 1024 * 1024:
            skipped.append(f'{key} ({size / 1048576:.0f} MB > {MAX_FILE_MB} MB)')
            continue
        url = (f.get('links') or {}).get('self') or f'{ZENODO_API}{rec_id}/files/{key}/content'
        dest = os.path.join(cache_dir, rec_id, safe_name(key))
        if not download(url, dest):
            skipped.append(f'{key} (download failed)')
            continue
        e = ext_of(key)
        try:
            if e in ARCHIVE_EXT:
                got = read_zip(dest, key, log)
                if got:
                    for t in got:
                        tables.append(t)
                        checked.append(dict(name=t['name'], format=t['format'], rows=int(t['df'].shape[0]),
                                            columns=int(t['df'].shape[1]), layout=t['layout']))
                else:
                    skipped.append(f'{key} (no tabular content in archive)')
                continue
            if e in ('.xlsx', '.xlsm', '.xls', '.ods'):
                got = [(f'{key} [{sn}]', df) for sn, df in read_excel_tables(dest)]
            else:
                df = read_text_table(dest)
                got = [(key, df)] if df is not None else []
            if not got:
                skipped.append(f'{key} (not a table)')
                continue
            for name, df in got[:8]:
                layout = 'wide' if detect_wide(df)[0] else 'long'
                tables.append(dict(name=name, df=df, layout=layout, format=e.lstrip('.'), source=key))
                checked.append(dict(name=name, format=e.lstrip('.'), rows=int(df.shape[0]),
                                    columns=int(df.shape[1]), layout=layout))
        except Exception as ex:  # unreadable file: keep going with the others
            skipped.append(f'{key} (unreadable: {type(ex).__name__})')
    return tables, checked, skipped


def read_zip(path, key, log):
    out = []
    try:
        zf = zipfile.ZipFile(path)
    except zipfile.BadZipFile:
        return out
    with zf:
        names = [n for n in zf.namelist() if not n.endswith('/') and '__MACOSX' not in n]
        if any(os.path.basename(n) == 'meta.xml' for n in names):
            try:
                frames, _ = read_dwca(zf)
            except Exception as ex:
                log(f'   DwC-A parse error in {key}: {ex}')
                frames = {}
            if frames:
                merged = merge_dwca(frames)
                out.append(dict(name=f'{key} (DwC-A)', df=merged, layout='dwca', format='dwca', source=key,
                                emof=frames.get('MeasurementOrFact') if 'MeasurementOrFact' in frames
                                else frames.get('ExtendedMeasurementOrFact')))
                return out
        members = [n for n in names if ext_of(n) in TABULAR_EXT]
        members.sort(key=lambda n: (bool(re.search(HELPER_FILES, n, re.I)), zf.getinfo(n).file_size))
        extract_root = os.path.join(os.path.dirname(path), '_unzipped', safe_name(key))
        for n in members[:MAX_FILES_PER_RECORD]:
            info = zf.getinfo(n)
            if info.file_size > MAX_FILE_MB * 1024 * 1024:
                continue
            dest = os.path.join(extract_root, safe_name(n.replace('/', '__')))
            if not os.path.exists(dest):
                os.makedirs(extract_root, exist_ok=True)
                tmp = dest + f'.{os.getpid()}.{threading.get_ident()}.part'
                with open(tmp, 'wb') as f:
                    f.write(zf.read(n))
                os.replace(tmp, dest)
            e = ext_of(n)
            try:
                if e in ('.xlsx', '.xlsm', '.xls', '.ods'):
                    got = [(f'{key}/{n} [{sn}]', df) for sn, df in read_excel_tables(dest)]
                else:
                    df = read_text_table(dest)
                    got = [(f'{key}/{n}', df)] if df is not None else []
            except Exception:
                got = []
            for name, df in got[:4]:
                out.append(dict(name=name, df=df, layout='wide' if detect_wide(df)[0] else 'long',
                                format=e.lstrip('.'), source=key))
    return out


def merge_dwca(frames):
    """Flatten a DwC-A to one occurrence-level table (occurrence joined with event when present)."""
    occ = frames.get('Occurrence')
    ev = frames.get('Event')
    if occ is not None and ev is not None:
        key_occ = 'coreid' if 'coreid' in occ.columns else ('eventID' if 'eventID' in occ.columns else None)
        key_ev = 'id' if 'id' in ev.columns else ('eventID' if 'eventID' in ev.columns else None)
        if key_occ and key_ev:
            evc = ev.drop(columns=[c for c in ev.columns if c in occ.columns and c != key_ev])
            return occ.merge(evc, left_on=key_occ, right_on=key_ev, how='left')
        return occ
    if occ is not None:
        return occ
    if ev is not None:
        return ev
    return next(iter(frames.values()))


# ---------------------------------------------------------------------------
# Column mapping
# ---------------------------------------------------------------------------
def map_columns(df):
    """Return {dwcTerm: column} for a table, using exact DwC names first, then synonyms."""
    mapping = {}
    cols = list(df.columns)
    normed = {c: _norm(c) for c in cols}
    # 1. exact match on the DwC term, then synonym dictionary
    for c in cols:
        n = normed[c]
        term = SYN_LOOKUP.get(n)
        if term is None:
            # strip common unit suffixes: "Latitude (N)", "depth_m", "Depth [m]"
            n2 = re.sub(r'(degrees?|deg|dd|wgs84|decimal|inmeters|meters|metres|m)$', '', n)
            term = SYN_LOOKUP.get(n2)
        if term is None:
            # Darwin Core variants: acceptedScientificName, decimalLatitudeStart, eventDateUTC ...
            for pat, t in PARTIAL_TERMS:
                if re.search(pat, n):
                    term = t
                    break
        if term and term not in mapping:
            mapping[term] = c
    # 2. weak names ('x', 'y', 'z', 'name', 'n', 'data', 'total') only when the content fits
    for term, col in list(mapping.items()):
        n = normed[col]
        if term in ('decimalLatitude', 'decimalLongitude') and n in ('x', 'y'):
            if numeric(df[col]).notna().mean() < 0.8:
                del mapping[term]
        if term == 'scientificName' and n in ('name', 'organism', 'art', 'identification', 'speciesid'):
            vals = df[col].dropna().astype(str).head(200)
            if len(vals) == 0 or vals.map(lambda v: bool(BINOMIAL_RE.match(v.strip()) or
                                                         GENUS_RE.match(v.strip()))).mean() < 0.3:
                del mapping[term]
        if term == 'eventDate' and n == 'data':
            if pd.to_datetime(df[col].dropna().astype(str).head(100), errors='coerce',
                              format='mixed').notna().mean() < 0.6:
                del mapping[term]
    # 3. a lone "depth" doubles as min/max depth
    if 'depth' in mapping:
        mapping.setdefault('minimumDepthInMeters', mapping['depth'])
        mapping.setdefault('maximumDepthInMeters', mapping['depth'])
    return mapping


def measure_columns(df):
    out = []
    for c in df.columns:
        n = _norm(c)
        if n.startswith(MEASURE_PREFIXES):
            if numeric(df[c]).notna().mean() >= 0.5:
                out.append(c)
    return out


def detect_wide(df):
    """Species-by-sample matrix: many columns named like taxa, holding numbers."""
    taxa_cols = []
    for c in df.columns:
        s = str(c).strip()
        if BINOMIAL_RE.match(s) or GENUS_RE.match(s) and not _norm(s) in SYN_LOOKUP:
            taxa_cols.append(c)
    if len(taxa_cols) >= 5 and len(taxa_cols) >= 0.3 * len(df.columns):
        num = sum(numeric(df[c]).notna().mean() > 0.5 for c in taxa_cols[:30])
        if num >= 0.6 * min(30, len(taxa_cols)):
            return True, taxa_cols
    # transposed matrix: first column holds taxon names, other columns numeric samples
    first = df.columns[0]
    vals = df[first].dropna().astype(str).head(200)
    if len(vals) >= 5 and df.shape[1] >= 5:
        taxa_like = vals.map(lambda v: bool(BINOMIAL_RE.match(v.strip()))).mean()
        num_cols = sum(numeric(df[c]).notna().mean() > 0.7 for c in df.columns[1:20])
        if taxa_like > 0.6 and num_cols >= 0.7 * min(19, df.shape[1] - 1) and \
                SYN_LOOKUP.get(_norm(first)) != 'scientificName':
            return True, ['__rows__']
    return False, []


def numeric(series):
    if series.dtype.kind in 'if':
        return series.astype(float)
    s = series.astype(str).str.strip().str.replace(',', '.', regex=False)
    s = s.str.replace(r'[^\d.\-eE+]', '', regex=True)
    return pd.to_numeric(s, errors='coerce')


def parse_dates(raw):
    """Parse a date column: datetimes, ISO 8601 (incl. intervals), dd/mm/yyyy or mm/dd/yyyy."""
    if raw.map(lambda v: isinstance(v, (datetime, pd.Timestamp))).mean() > 0.9:
        return pd.to_datetime(raw, errors='coerce', utc=True), 1.0
    s = raw.astype(str).str.strip()
    iso = s.map(lambda v: bool(ISO_DATE_RE.match(v))).mean()
    interval = s.str.match(r'^\d{4}[^/]*/\d{4}')
    s = s.where(~interval, s.str.split('/').str[0])
    first = pd.to_numeric(s.str.extract(r'^(\d{1,2})[./-]\d{1,2}[./-]\d{2,4}')[0], errors='coerce')
    dayfirst = bool((first > 12).any()) or not bool(first.notna().any())
    parsed = pd.to_datetime(s, errors='coerce', format='mixed', utc=True, dayfirst=dayfirst)
    return parsed, iso


def names_from_title(title, worms):
    """Single-species datasets often name the taxon only in the title."""
    cands = re.findall(r'\b([A-Z][a-z]{2,} [a-z]{3,})\b', title or '')
    cands = [c for c in dict.fromkeys(cands) if c.split()[1] not in
             ('and', 'the', 'from', 'data', 'for', 'with', 'sea', 'bay', 'coast', 'ocean', 'in', 'of')]
    if not cands:
        return []
    res = worms.match(cands[:10])
    return [n for n, v in res.items() if v]


def table_score(mapping):
    core = ['eventDate', 'decimalLatitude', 'decimalLongitude', 'scientificName']
    s = sum(3 for t in core if t in mapping)
    if 'eventDate' not in mapping and 'year' in mapping:
        s += 2
    s += sum(1 for t in ['occurrenceID', 'scientificNameID', 'basisOfRecord', 'occurrenceStatus',
                         'minimumDepthInMeters', 'eventID', 'locality'] if t in mapping)
    return s


# ---------------------------------------------------------------------------
# WoRMS + OBIS lookups (cached)
# ---------------------------------------------------------------------------
class WormsCache:
    def __init__(self, path):
        self.path = path
        self.lock = threading.Lock()
        self.data = {}
        if os.path.exists(path):
            try:
                with open(path, encoding='utf-8') as f:
                    self.data = json.load(f)
            except Exception:
                self.data = {}

    def save(self):
        with self.lock:
            tmp = self.path + '.part'
            with open(tmp, 'w', encoding='utf-8') as f:
                json.dump(self.data, f)
            os.replace(tmp, self.path)

    def match(self, names):
        names = [n for n in dict.fromkeys(names) if n]
        todo = [n for n in names if n not in self.data]
        for i in range(0, len(todo), 50):
            batch = todo[i:i + 50]
            params = [('scientificnames[]', n) for n in batch] + [('marine_only', 'false')]
            res = None
            for attempt in range(3):
                try:
                    WORMS_LIMIT.wait()
                    r = session().get(WORMS_API, params=params, timeout=90)
                    if r.status_code == 204:
                        res = [[] for _ in batch]
                    elif r.status_code == 200:
                        res = r.json()
                    if res is not None:
                        break
                except Exception:
                    time.sleep(3)
            if res is None:
                continue
            with self.lock:
                for n, hits in zip(batch, res):
                    hit = next((h for h in (hits or []) if h), None)
                    if hit:
                        self.data[n] = dict(aphia=hit.get('valid_AphiaID') or hit.get('AphiaID'),
                                            match=hit.get('match_type'),
                                            marine=bool(hit.get('isMarine')),
                                            brackish=bool(hit.get('isBrackish')),
                                            freshwater=bool(hit.get('isFreshwater')),
                                            terrestrial=bool(hit.get('isTerrestrial')))
                    else:
                        self.data[n] = None
        return {n: self.data.get(n) for n in names}


def clean_name(v):
    s = str(v).strip().replace('_', ' ')
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'\b(sp|spp|cf|aff|indet|var|subsp|ssp)\.?(\s|$).*', '', s, flags=re.I).strip()
    s = re.sub(r'\(.*?\)', '', s).strip()
    words = s.split(' ')
    if not words or not re.match(r'^[A-Z][a-z]+$', words[0]):
        return ''
    if len(words) >= 2 and re.match(r'^[a-z][a-z-]+$', words[1]):
        return f'{words[0]} {words[1]}'
    return words[0]


def xylookup(points):
    """points: list of (lon, lat) -> list of dict(shoredistance, bathymetry) or None."""
    out = []
    for i in range(0, len(points), 1000):
        batch = points[i:i + 1000]
        res = None
        for attempt in range(3):
            try:
                OBIS_LIMIT.wait()
                r = session().post(XYLOOKUP_API, json={'points': [list(p) for p in batch],
                                                        'areas': False, 'grids': True, 'shoredistance': True},
                                   timeout=120)
                if r.status_code == 200:
                    res = r.json()
                    break
            except Exception:
                time.sleep(3)
        if res is None or len(res) != len(batch):
            out.extend([None] * len(batch))
        else:
            out.extend(res)
    return out


# ---------------------------------------------------------------------------
# Level A - metadata QC
# ---------------------------------------------------------------------------
def check(cid, result, details=''):
    return dict(id=cid, result=result, details=details)


def metadata_checks(item, record):
    md = (record or {}).get('metadata', {})
    out = []
    access = item.get('Access conditions', 'N/A')
    access_api = md.get('access_right')
    if access_api:
        access = access_api.capitalize()
    out.append(check('M01', 'pass' if access == 'Open' else 'fail', f'Access: {access}'))

    lic = ((md.get('license') or {}).get('id') if record else None) or item.get('Data license', 'N/A')
    lic_l = (lic or '').lower()
    if lic_l in ('cc-zero', 'cc0-1.0', 'cc-by-4.0', 'cc-by-3.0', 'cc-by', 'odc-by-1.0', 'pddl-1.0', 'cc-by-2.0'):
        out.append(check('M02', 'pass', lic))
    elif lic_l in ('', 'n/a', 'notspecified', 'none', 'other-closed'):
        out.append(check('M02', 'fail', f'No licence declared ({lic})'))
    else:
        out.append(check('M02', 'warn', f'{lic}: restrictions (NC/ND/SA or non-data licence)'))

    doi = md.get('doi') or (record or {}).get('doi') or item.get('Data URL/DOI', '')
    out.append(check('M03', 'pass' if doi and str(doi).startswith('10.') else 'fail', str(doi or 'none')))

    email = item.get('Email of data provider', 'N/A')
    creators = md.get('creators') or []
    has_aff = any(c.get('affiliation') for c in creators)
    if email and email != 'N/A':
        out.append(check('M04', 'pass', f'E-mail: {email}'))
    elif has_aff or ('(' in item.get('Provider', '')):
        out.append(check('M04', 'warn', 'No e-mail; creator affiliation available'))
    else:
        out.append(check('M04', 'fail', 'No contact information'))

    tc = str(item.get('Temporal coverage', 'N/A'))
    pub_year = str(md.get('publication_date', ''))[:4]
    if tc in ('N/A', '', 'nan'):
        out.append(check('M05', 'fail', 'Unknown'))
    elif tc.strip() == pub_year:
        out.append(check('M05', 'warn', f'{tc} (publication year only)'))
    else:
        out.append(check('M05', 'pass', tc))

    if record and record.get('files') is not None:
        exts = {ext_of(f.get('key', '')) for f in record.get('files') or []}
    else:
        exts = set(item.get('File type', '').split(', '))
    tab = sorted(e for e in exts if e in TABULAR_EXT | ARCHIVE_EXT)
    if access != 'Open' and not (record or {}).get('files'):
        out.append(check('M06', 'na', 'Files not public'))
    else:
        out.append(check('M06', 'pass' if tab else 'fail', ', '.join(tab) if tab else
                         f'No tabular files ({", ".join(sorted(e for e in exts if e)) or "no files"})'))

    st = item.get('Status', '')
    out.append(check('M07', 'pass' if st == 'Included' else 'warn', st))

    desc = re.sub('<[^<]+?>', ' ', md.get('description', '') or '')
    kws = md.get('keywords') or []
    if len(desc) >= 300 and kws:
        out.append(check('M08', 'pass', f'{len(desc)} chars, {len(kws)} keywords'))
    elif record is None:
        out.append(check('M08', 'na', 'Record metadata unavailable'))
    else:
        out.append(check('M08', 'warn', f'{len(desc)} chars, {len(kws)} keywords'))
    return out, access, lic, doi, tab


def score(checks):
    """Weighted score: error-severity checks weigh 2, warnings 1; pass=1, warn=0.5, fail=0."""
    num = den = 0.0
    for c in checks:
        if c['result'] == 'na':
            continue
        w = 2.0 if CRIT[c['id']]['severity'] == 'error' else 1.0
        den += w
        num += w * {'pass': 1.0, 'warn': 0.5, 'fail': 0.0}[c['result']]
    return int(round(100 * num / den)) if den else None


# ---------------------------------------------------------------------------
# Level B - data QC
# ---------------------------------------------------------------------------
def pct(a, b):
    return round(100.0 * a / b, 1) if b else None


def data_checks(tables, regions, worms, title=''):
    checks = []
    issues = []
    stats = dict(n_records=None, n_taxa=None, worms_matched_pct=None, onland_pct=None)

    if not tables:
        checks.append(check('D01', 'fail', 'No file could be read as a table'))
        return checks, issues, stats
    checks.append(check('D01', 'pass', f'{len(tables)} table(s) read'))

    # primary table = best Darwin Core coverage
    for t in tables:
        t['map'] = map_columns(t['df'])
        t['score'] = table_score(t['map']) + (3 if t['layout'] == 'dwca' else 0)
    # the primary table is the one holding the taxa (occurrence level), then best DwC coverage
    tables.sort(key=lambda t: (not ('scientificName' in t['map'] or t['layout'] in ('wide', 'dwca')),
                               -t['score'], -t['df'].shape[0]))
    prim = tables[0]
    df, m = prim['df'], prim['map']
    wide, taxa_cols = detect_wide(df) if prim['layout'] != 'dwca' else (False, [])
    stats['n_records'] = int(df.shape[0])

    # coordinates may live in a separate station table
    coord_tab = prim if {'decimalLatitude', 'decimalLongitude'} <= set(m) else next(
        (t for t in tables[1:] if {'decimalLatitude', 'decimalLongitude'} <= set(t['map'])), None)

    # dates may live in a separate sample / station table as well
    has_d = lambda mp: 'eventDate' in mp or 'year' in mp
    date_tab = prim if has_d(m) else next((t for t in tables[1:] if has_d(t['map'])), None)

    # ---- D02 core fields
    has_date = date_tab is not None
    has_taxon = 'scientificName' in m or wide
    title_taxa = [] if has_taxon else names_from_title(title, worms)
    has_coord = coord_tab is not None
    missing = [n for n, ok in [('date', has_date), ('coordinates', has_coord),
                               ('scientific name', has_taxon or bool(title_taxa))] if not ok]
    notes = []
    if title_taxa:
        notes.append(f'no taxon column, taxon from title: {", ".join(title_taxa)}')
    if has_coord and coord_tab is not prim:
        notes.append(f'coordinates in separate table "{coord_tab["name"]}"')
    if has_date and date_tab is not prim:
        notes.append(f'dates in separate table "{date_tab["name"]}"')
    if not missing and notes:
        checks.append(check('D02', 'warn', 'Core fields found, but ' + '; '.join(notes) +
                            (' - join needed' if (coord_tab is not prim or date_tab is not prim) else '')))
    elif not missing:
        checks.append(check('D02', 'pass', f'Mapped from "{prim["name"]}": ' +
                            ', '.join(f'{k}<-{v}' for k, v in m.items() if k in
                                      ('eventDate', 'year', 'decimalLatitude', 'decimalLongitude',
                                       'scientificName')) + (' (taxa as columns)' if wide else '')))
    else:
        checks.append(check('D02', 'fail', 'Missing: ' + ', '.join(missing)))
        issues.append('Missing ' + ', '.join(missing))

    # ---- D03 DwC identifiers & terms
    dwc_terms = ['occurrenceID', 'scientificNameID', 'occurrenceStatus', 'basisOfRecord']
    present = [t for t in dwc_terms if t in m]
    if len(present) == 4:
        checks.append(check('D03', 'pass', 'All present'))
    else:
        checks.append(check('D03', 'warn' if present else 'fail',
                            f'Missing: {", ".join(t for t in dwc_terms if t not in present)} '
                            f'(to be generated at DwC conversion)'))

    # ---- D04 / D05 / D06 coordinates
    xy_info = {}
    if coord_tab is not None:
        cdf, cm = coord_tab['df'], coord_tab['map']
        lat = numeric(cdf[cm['decimalLatitude']])
        lon = numeric(cdf[cm['decimalLongitude']])
        n = int(lat.notna().sum())
        both = lat.notna() & lon.notna()
        bad_range = both & ((lat.abs() > 90) | (lon.abs() > 180))
        zero = both & (lat == 0) & (lon == 0)
        raw_str = cdf[cm['decimalLatitude']].dropna().astype(str).head(200)
        dms = raw_str.str.contains(r'[°\'"]|\d+\s+\d+\s+\d+', regex=True).mean() if len(raw_str) else 0
        nb = int(both.sum())
        if nb == 0:
            checks.append(check('D04', 'fail', 'Coordinate columns are empty or not numeric' +
                                (' (degrees-minutes-seconds?)' if dms > 0.3 else '')))
            issues.append('Coordinates not in decimal degrees')
        else:
            bad = int(bad_range.sum() + zero.sum())
            miss = len(cdf) - nb
            res = 'pass' if bad == 0 and miss <= 0.05 * len(cdf) else ('fail' if bad > 0.05 * nb else 'warn')
            checks.append(check('D04', res, f'{nb} coordinate pairs; {int(bad_range.sum())} out of range, '
                                f'{int(zero.sum())} at (0,0), {miss} missing'))
            if res != 'pass':
                issues.append('Invalid or missing coordinates')
            ok = both & ~bad_range & ~zero
            pts = pd.DataFrame({'lon': lon[ok].round(4), 'lat': lat[ok].round(4)}).drop_duplicates()
            # D05 region bbox
            boxes = [b for r in regions for b in REGION_BBOX.get(r, [])]
            if boxes and len(pts):
                inside = pd.Series(False, index=pts.index)
                for (x0, y0, x1, y1) in boxes:
                    inside |= (pts.lon >= x0) & (pts.lon <= x1) & (pts.lat >= y0) & (pts.lat <= y1)
                p = pct(int(inside.sum()), len(pts))
                checks.append(check('D05', 'pass' if p >= 90 else ('warn' if p >= 30 else 'fail'),
                                    f'{p}% of {len(pts)} unique positions inside {", ".join(regions)}'))
                if p < 30:
                    issues.append('Positions outside the declared region')
            else:
                checks.append(check('D05', 'na', 'No region bounding box'))
            # D06 on land
            sample = pts.sample(MAX_XY_POINTS, random_state=1) if len(pts) > MAX_XY_POINTS else pts
            plist = list(zip(sample.lon.tolist(), sample.lat.tolist()))
            res_xy = xylookup(plist) if plist else []
            valid = [(p, r) for p, r in zip(plist, res_xy) if r]
            if valid:
                land = sum(1 for _, r in valid if (r.get('shoredistance') or 0) < 0)
                p = pct(land, len(valid))
                stats['onland_pct'] = p
                res = 'pass' if p <= 5 else ('warn' if p <= 25 else 'fail')
                checks.append(check('D06', res, f'{land} of {len(valid)} unique positions on land ({p}%)'))
                if res != 'pass':
                    issues.append(f'{p}% of positions on land')
                xy_info = {p: (r.get('grids') or {}).get('bathymetry') for p, r in valid}
            else:
                checks.append(check('D06', 'na', 'OBIS xylookup unavailable'))
    else:
        checks.append(check('D04', 'na', 'No coordinates'))
        checks.append(check('D05', 'na', 'No coordinates'))
        checks.append(check('D06', 'na', 'No coordinates'))

    # ---- D07 dates
    ddf, dm = (date_tab['df'], date_tab['map']) if date_tab is not None else (df, {})
    if 'eventDate' in dm:
        raw = ddf[dm['eventDate']].dropna()
        n = len(raw)
        if n == 0:
            checks.append(check('D07', 'fail', 'Date column empty'))
        else:
            parsed, iso = parse_dates(raw)
            ok = parsed.notna().mean()
            future = int((parsed > pd.Timestamp.now(tz='UTC')).sum())
            early = int((parsed < pd.Timestamp('1700-01-01', tz='UTC')).sum())
            if ok < 0.8 or future > 0.05 * n:
                res = 'fail'
            elif iso < 0.95 or future or early or ok < 0.95:
                res = 'warn'
            else:
                res = 'pass'
            span = ''
            if parsed.notna().any():
                span = f'; range {parsed.min().date()} .. {parsed.max().date()}'
            checks.append(check('D07', res, f'{pct(ok * n, n)}% parseable, {round(iso * 100, 1)}% ISO 8601, '
                                f'{future} in future{span}'))
            if res != 'pass':
                issues.append('Dates not ISO 8601 / invalid' if ok >= 0.8 else 'Dates not parseable')
    elif 'year' in dm:
        y = numeric(ddf[dm['year']])
        ok = y.dropna().between(1700, datetime.now().year).mean() if y.notna().any() else 0
        miss = int(y.isna().sum())
        checks.append(check('D07', 'warn' if ok >= 0.9 else 'fail',
                            f'Only year/month/day columns ({round(ok * 100, 1)}% valid years, {miss} empty) - '
                            'eventDate to be built'))
        if ok >= 0.9:
            issues.append('eventDate to be built from year/month/day')
    else:
        checks.append(check('D07', 'na', 'No date'))

    # ---- D08 depth
    if 'minimumDepthInMeters' in m:
        dmin = numeric(df[m['minimumDepthInMeters']])
        dmax = numeric(df[m['maximumDepthInMeters']])
        n = int(dmin.notna().sum())
        if n == 0:
            checks.append(check('D08', 'warn', 'Depth column not numeric'))
        else:
            neg = int((dmin < 0).sum())
            inv = int((dmin > dmax).sum())
            deep = 0
            if xy_info and coord_tab is prim:
                lat = numeric(df[m['decimalLatitude']]).round(4)
                lon = numeric(df[m['decimalLongitude']]).round(4)
                bathy = pd.Series([xy_info.get((a, b)) for a, b in zip(lon, lat)], index=df.index, dtype=float)
                deep = int(((dmax > bathy * 1.1 + 50) & bathy.notna() & (bathy > 0)).sum())
            bad = neg + inv + deep
            res = 'pass' if bad == 0 else ('warn' if bad <= 0.2 * n else 'fail')
            checks.append(check('D08', res, f'{n} depth values; {neg} negative, {inv} min>max, '
                                f'{deep} deeper than bathymetry'))
            if res != 'pass':
                issues.append('Depth inconsistencies')
    else:
        checks.append(check('D08', 'na', 'No depth'))

    # ---- D09 / D10 taxonomy
    names = []
    if wide and taxa_cols and taxa_cols != ['__rows__']:
        names = [str(c) for c in taxa_cols]
    elif wide:
        names = df[df.columns[0]].dropna().astype(str).tolist()
    elif 'scientificName' in m:
        names = df[m['scientificName']].dropna().astype(str).tolist()
    elif title_taxa:
        names = title_taxa
    cleaned = list(dict.fromkeys(filter(None, (clean_name(n) for n in names))))
    stats['n_taxa'] = len(cleaned) if names else None
    if cleaned:
        sample = cleaned
        if len(cleaned) > MAX_WORMS_NAMES:
            sample = pd.Series(cleaned).sample(MAX_WORMS_NAMES, random_state=1).tolist()
        res_w = worms.match(sample)
        matched = [v for v in res_w.values() if v]
        p = pct(len(matched), len(sample))
        stats['worms_matched_pct'] = p
        res = 'pass' if p >= 90 else ('warn' if p >= 60 else 'fail')
        checks.append(check('D09', res, f'{len(matched)} of {len(sample)} unique names matched'
                            + (' (sample)' if len(sample) < len(cleaned) else '')))
        if res != 'pass':
            issues.append(f'Only {p}% of names found in WoRMS')
        if matched:
            mar = sum(1 for v in matched if v['marine'] or v['brackish'])
            p2 = pct(mar, len(matched))
            checks.append(check('D10', 'pass' if p2 >= 80 else ('warn' if p2 >= 50 else 'fail'),
                                f'{p2}% of matched taxa marine/brackish'))
            if p2 < 50:
                issues.append('Mostly non-marine taxa')
        else:
            checks.append(check('D10', 'na', 'No matched taxa'))
    elif names:
        checks.append(check('D09', 'fail', 'Taxon column holds no recognisable scientific names'))
        checks.append(check('D10', 'na', 'No matched taxa'))
        issues.append('No recognisable scientific names')
    else:
        checks.append(check('D09', 'na', 'No taxon column'))
        checks.append(check('D10', 'na', 'No taxon column'))

    # ---- D11 scientificNameID
    if 'scientificNameID' in m:
        v = df[m['scientificNameID']].dropna().astype(str).str.strip()
        lsid = v.map(lambda s: bool(LSID_RE.match(s))).mean() if len(v) else 0
        aphia_num = v.str.fullmatch(r'\d+').mean() if len(v) else 0
        if lsid >= 0.95:
            checks.append(check('D11', 'pass', f'{round(lsid * 100, 1)}% valid WoRMS LSIDs'))
        elif aphia_num >= 0.9:
            checks.append(check('D11', 'warn', 'Numeric AphiaIDs - convert to urn:lsid:marinespecies.org:taxname:'))
        else:
            checks.append(check('D11', 'fail', f'Only {round(lsid * 100, 1)}% valid WoRMS LSIDs'))
    else:
        checks.append(check('D11', 'na', 'Not present'))

    # ---- D12 vocabularies
    vocab = []
    for term, allowed in [('basisOfRecord', BASIS_OF_RECORD), ('occurrenceStatus', OCCURRENCE_STATUS)]:
        if term in m and _norm(m[term]) == _norm(term):
            v = df[m[term]].dropna().astype(str).map(_norm)
            if len(v):
                vocab.append((term, v.isin(allowed).mean()))
    if vocab:
        worst = min(p for _, p in vocab)
        checks.append(check('D12', 'pass' if worst >= 0.99 else ('warn' if worst >= 0.8 else 'fail'),
                            '; '.join(f'{t}: {round(p * 100, 1)}% valid' for t, p in vocab)))
    else:
        checks.append(check('D12', 'na', 'Terms not present'))

    # ---- D13 occurrenceID
    if 'occurrenceID' in m:
        v = df[m['occurrenceID']].astype(str).str.strip().str.replace(' ', ' ')
        empty = int((v.isin(['', 'nan', 'None'])).sum())
        dup = int(v[~v.isin(['', 'nan', 'None'])].duplicated().sum())
        checks.append(check('D13', 'pass' if empty == dup == 0 else 'fail', f'{dup} duplicates, {empty} empty'))
        if dup or empty:
            issues.append('Duplicate / empty occurrenceID')
    else:
        checks.append(check('D13', 'na', 'Not present'))

    # ---- D14 measurements
    emof = prim.get('emof')
    if emof is not None and len(emof):
        t_ok = emof['measurementType'].notna().mean() if 'measurementType' in emof else 0
        u_ok = emof['measurementUnit'].notna().mean() if 'measurementUnit' in emof else 0
        res = 'pass' if t_ok >= 0.99 and u_ok >= 0.6 else 'warn'
        checks.append(check('D14', res, f'eMoF: {len(emof)} rows, type {round(t_ok * 100)}%, '
                            f'unit {round(u_ok * 100)}%'))
    else:
        mcols = measure_columns(df)
        if wide:
            mcols = mcols or ['(matrix values)']
        if 'measurementValue' in m:
            ok = 'measurementType' in m and 'measurementUnit' in m
            checks.append(check('D14', 'pass' if ok else 'warn', 'measurementValue with ' +
                                ('type and unit' if ok else 'missing type/unit')))
        elif mcols:
            with_unit = [c for c in mcols if UNIT_HINT.search(str(c))]
            unit_col = 'measurementUnit' in m
            res = 'pass' if unit_col or len(with_unit) == len(mcols) else 'warn'
            checks.append(check('D14', res, f'{len(mcols)} measurement column(s): '
                                f'{", ".join(map(str, mcols[:5]))}; units ' +
                                ('given' if res == 'pass' else 'not stated in header')))
            if res != 'pass':
                issues.append('Measurement units not stated')
        else:
            checks.append(check('D14', 'na', 'No measurement columns'))

    # ---- D15 layout
    if prim['layout'] == 'dwca':
        checks.append(check('D15', 'pass', 'Darwin Core Archive'))
    elif wide:
        checks.append(check('D15', 'warn', 'Species-by-sample matrix - reshape to long format'))
        issues.append('Wide (matrix) layout')
    else:
        checks.append(check('D15', 'pass', 'Long format'))

    return checks, issues, stats


def decide(meta_checks, dchecks, data_checked, reason_not_checked, issues, layout_dwca):
    if not data_checked:
        return 'Not checked', f'Not checked: {reason_not_checked}'
    res = {c['id']: c['result'] for c in dchecks}
    errors = [c for c in dchecks if c['result'] == 'fail' and CRIT[c['id']]['severity'] == 'error']
    warns = [c for c in dchecks if c['result'] in ('warn', 'fail')]
    if errors:
        status = 'Fail'
    elif warns:
        status = 'Warnings'
    else:
        status = 'Pass'

    if res.get('D01') == 'fail':
        rec = 'Files could not be read - manual inspection needed'
    elif res.get('D02') == 'fail':
        miss = next(c['details'] for c in dchecks if c['id'] == 'D02').replace('Missing: ', '')
        rec = f'Needs {miss} before DwC conversion (contact provider / check paper)'
    elif status == 'Fail':
        bad = [CRIT[c['id']]['name'] for c in errors]
        rec = 'Fix before conversion: ' + ', '.join(bad)
    elif layout_dwca:
        rec = 'Already Darwin Core - candidate for direct EMODnet/IPT ingestion'
    elif status == 'Warnings':
        fixes = [i for i in issues if 'skipped' not in i]
        if res.get('D02') == 'warn':
            d02 = next(c['details'] for c in dchecks if c['id'] == 'D02')
            if 'taxon from title' in d02:
                fixes.insert(0, 'taxon taken from title')
            if 'separate table' in d02:
                fixes.insert(0, 'join sample/station tables')
        rec = ('Convertible to DwC after fixes: ' + ', '.join(dict.fromkeys(fixes[:3]))) if fixes else \
            'Ready for DwC conversion (identifiers/vocabularies to be generated)'
    else:
        rec = 'Ready for DwC conversion'
    return status, rec


# ---------------------------------------------------------------------------
# Per-record pipeline
# ---------------------------------------------------------------------------
def record_id_of(item):
    url = str(item.get('Metadata URL', ''))
    m = re.search(r'/records/(\d+)', url)
    return m.group(1) if m else None


def process(item, args, worms, log):
    rec_id = record_id_of(item)
    record = None
    try:
        record = get_record(rec_id, args.cache)
    except Exception as ex:
        log(f'   [{rec_id}] API error: {ex}')
    mchecks, access, lic, doi, tab = metadata_checks(item, record)
    regions = [r.strip() for r in str(item.get('Geographic coverage', '')).split(',') if r.strip()]

    out = dict(record_id=rec_id,
               title=item.get('Data title'),
               doi=str(doi) if doi else None,
               url=f'https://zenodo.org/records/{rec_id}',
               taxa=item.get('Taxa covered'),
               regions=item.get('Geographic coverage'),
               inventory_status=item.get('Status'),
               access=access,
               license=lic,
               metadata_score=score(mchecks),
               data_checked=False,
               files_checked=[],
               n_records=None, n_taxa=None, worms_matched_pct=None, onland_pct=None,
               data_score=None,
               qc_status='Not checked',
               checks=list(mchecks),
               issues=[],
               recommendation='')

    reason = None
    if record is None:
        reason = 'Zenodo record unavailable'
    elif access != 'Open':
        reason = f'{access} access'
    elif not tab:
        reason = 'no tabular data files'
    elif args.skip_data:
        reason = 'data QC skipped (--skip-data)'

    dchecks, issues, layout_dwca = [], [], False
    if reason is None:
        tables, checked, skipped = load_tables(rec_id, record, args.cache, log)
        out['files_checked'] = checked
        if not tables and not any('unreadable' in s for s in skipped):
            if skipped and all('MB >' in s or 'download failed' in s for s in skipped):
                reason = 'all tabular files too large / unavailable: ' + '; '.join(skipped[:3])
            else:
                reason = 'no tabular content found in files' + (f' ({"; ".join(skipped[:2])})' if skipped else '')
        if reason is None:
            dchecks, issues, stats = data_checks(tables, regions, worms, item.get('Data title', ''))
            out.update(stats)
            out['data_checked'] = True
            out['data_score'] = score(dchecks)
            layout_dwca = any(t['layout'] == 'dwca' for t in tables[:1])
            if skipped:
                issues.append(f'{len(skipped)} file(s) skipped: ' + '; '.join(skipped[:3]))
        for t in tables:
            t.pop('df', None)

    out['checks'] = mchecks + dchecks
    out['issues'] = list(dict.fromkeys(issues))
    out['qc_status'], out['recommendation'] = decide(mchecks, dchecks, out['data_checked'], reason,
                                                     out['issues'], layout_dwca)
    return out


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def build_summary(results):
    st = [r['qc_status'] for r in results]
    ms = [r['metadata_score'] for r in results if r['metadata_score'] is not None]
    ds = [r['data_score'] for r in results if r['data_score'] is not None]
    return dict(total=len(results),
                data_checked=sum(1 for r in results if r['data_checked']),
                **{'pass': st.count('Pass')},
                warnings=st.count('Warnings'),
                fail=st.count('Fail'),
                not_checked=st.count('Not checked'),
                avg_metadata_score=round(sum(ms) / len(ms), 1) if ms else 0.0,
                avg_data_score=round(sum(ds) / len(ds), 1) if ds else 0.0)


def write_outputs(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    summary = build_summary(results)
    doc = dict(generated=datetime.now(timezone.utc).astimezone().isoformat(timespec='seconds'),
               tool=TOOL_VERSION,
               summary=summary,
               criteria=[CRIT[c[0]] for c in CRITERIA],
               datasets=results)
    with open(os.path.join(out_dir, 'qc_results.json'), 'w', encoding='utf-8') as f:
        json.dump(doc, f, ensure_ascii=False, indent=1, default=str)

    rows = []
    for r in results:
        row = {k: r[k] for k in ('record_id', 'title', 'doi', 'url', 'taxa', 'regions', 'inventory_status',
                                 'access', 'license', 'qc_status', 'metadata_score', 'data_score',
                                 'data_checked', 'n_records', 'n_taxa', 'worms_matched_pct', 'onland_pct',
                                 'recommendation')}
        row['issues'] = '; '.join(r['issues'])
        row['files_checked'] = '; '.join(f['name'] for f in r['files_checked'])
        for c in r['checks']:
            row[c['id']] = c['result']
        rows.append(row)
    cols_order = list(rows[0].keys()) if rows else []
    for c in CRITERIA:
        if c[0] not in cols_order:
            cols_order.append(c[0])
    datasets = pd.DataFrame(rows).reindex(columns=cols_order)
    datasets.to_csv(os.path.join(out_dir, 'qc_results.csv'), index=False, encoding='utf-8-sig')

    long = pd.DataFrame([dict(record_id=r['record_id'], title=r['title'], check_id=c['id'],
                              check=CRIT[c['id']]['name'], level=CRIT[c['id']]['level'],
                              severity=CRIT[c['id']]['severity'], result=c['result'], details=c['details'])
                         for r in results for c in r['checks']])
    summ = pd.DataFrame([dict(metric=k, value=v) for k, v in summary.items()] +
                        [dict(metric='generated', value=doc['generated']), dict(metric='tool', value=TOOL_VERSION)])
    with pd.ExcelWriter(os.path.join(out_dir, 'QC_Report_Zenodo.xlsx'), engine='openpyxl') as w:
        summ.to_excel(w, sheet_name='Summary', index=False)
        datasets.to_excel(w, sheet_name='Datasets', index=False)
        long.to_excel(w, sheet_name='Checks', index=False)
        pd.DataFrame([CRIT[c[0]] for c in CRITERIA]).to_excel(w, sheet_name='Criteria', index=False)
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='EMODnet Biology quality check of Zenodo inventory datasets')
    ap.add_argument('--inventory', default='data.json', help='inventory JSON produced by Task 1')
    ap.add_argument('--cache', default='zenodo_cache', help='download cache folder (shared with Task 3)')
    ap.add_argument('--out', default='output', help='output folder')
    ap.add_argument('--limit', type=int, default=0, help='data-check at most N records')
    ap.add_argument('--ids', default='', help='comma-separated Zenodo record ids to process')
    ap.add_argument('--skip-data', action='store_true', help='metadata QC only')
    ap.add_argument('--force', action='store_true', help='ignore per-record results from earlier runs')
    ap.add_argument('--workers', type=int, default=6, help='parallel records')
    args = ap.parse_args()

    with open(args.inventory, encoding='utf-8') as f:
        inventory = json.load(f)
    if args.ids:
        wanted = set(args.ids.split(','))
        inventory = [i for i in inventory if record_id_of(i) in wanted]
    # Included first, then Manual review
    inventory.sort(key=lambda i: 0 if i.get('Status') == 'Included' else 1)

    os.makedirs(args.cache, exist_ok=True)
    partial_dir = os.path.join(args.out, '_records')
    os.makedirs(partial_dir, exist_ok=True)
    worms = WormsCache(os.path.join(args.out, 'worms_cache.json'))
    lock = threading.Lock()

    def log(msg):
        with lock:
            print(msg, flush=True)

    if not ZENODO_TOKEN:
        log('ZENODO_TOKEN not set - anonymous API access (lower rate limits).')

    t0 = time.time()
    results = {}
    data_budget = args.limit or None
    jobs = []
    for item in inventory:
        rid = record_id_of(item)
        if not rid:
            continue
        cached = os.path.join(partial_dir, f'{rid}.json')
        if os.path.exists(cached) and not args.force:
            with open(cached, encoding='utf-8') as f:
                results[rid] = json.load(f)
            continue
        jobs.append(item)

    prefetch_records([record_id_of(i) for i in jobs], args.cache, log)

    def run(item, skip):
        a = argparse.Namespace(**vars(args))
        a.skip_data = args.skip_data or skip
        r = process(item, a, worms, log)
        if not a.skip_data and r['recommendation'] != 'Not checked: Zenodo record unavailable':
            with open(os.path.join(partial_dir, f'{r["record_id"]}.json'), 'w', encoding='utf-8') as f:
                json.dump(r, f, ensure_ascii=False, default=str)
        return r

    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {}
        for idx, item in enumerate(jobs):
            skip = data_budget is not None and idx >= data_budget
            futs[ex.submit(run, item, skip)] = item
        for fut in as_completed(futs):
            done += 1
            try:
                r = fut.result()
            except Exception as e:
                rid = record_id_of(futs[fut])
                log(f'[{done}/{len(jobs)}] {rid} ERROR {type(e).__name__}: {e}')
                continue
            results[r['record_id']] = r
            log(f'[{done}/{len(jobs)}] {r["record_id"]} {r["qc_status"]:<11} meta={r["metadata_score"]} '
                f'data={r["data_score"]}  {str(r["title"])[:60]}')
            if done % 25 == 0:
                worms.save()

    worms.save()
    ordered = [results[record_id_of(i)] for i in inventory if record_id_of(i) in results]
    summary = write_outputs(ordered, args.out)
    log(f'\nDone in {(time.time() - t0) / 60:.1f} min')
    for k, v in summary.items():
        log(f'  {k:<20}: {v}')
    log(f'Outputs: {os.path.join(args.out, "qc_results.json")}, qc_results.csv, QC_Report_Zenodo.xlsx')


if __name__ == '__main__':
    main()
