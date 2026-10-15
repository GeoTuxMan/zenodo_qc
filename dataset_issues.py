"""
Errors and possible fixes for every dataset (Task 2 report over Tasks 2 and 3).

Collects, per Zenodo record, the problems found by
  * the Darwin Core conversion (Task 3, conversion_results.json): status, missing required terms,
    per-table issues;
  * the final QC of the archives (stage 2, qc_dwca_results.json): failed / warning checks D01-D15,
    M01-M08 and the EMODnetBiocheck / obistools messages;
  * the pre-screening of the raw files (stage 1, qc_results.json) for records without an archive;
  * the search for dates / coordinates in the linked papers (paper_values.json);
  * the data types of the files (data_types.json): files not read and why;
and attaches to every problem the ways to solve it.

    python dataset_issues.py
    python dataset_issues.py --dwc ../zenodo_dwc_repo/output --out output

Outputs: output/dataset_issues.json (consumed by the dashboard) and Dataset_Issues.xlsx
(sheets Summary, Datasets, Issues).
"""
import argparse
import json
import os
import re
from collections import Counter
from datetime import date

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# (regex on the message, category, severity or None = keep, fix)
RULES = [
    # ---- conversion (Task 3)
    (r'no taxon column / species matrix detected', 'No species data', 'error',
     'If the data are biological, give the converter a mapping override (<zenodo id>.json: column -> '
     'scientificName) or a species lookup table for the codes; if they are environmental, acoustic, model or '
     'tracking data, no Darwin Core Archive is needed - use the CSV / XLSX / JSON export (Data types).'),
    (r'no readable data tables found', 'No readable table', 'error',
     'The files are not tables (images, audio, rasters, scripts, PDF). Ask the authors for the data as CSV / '
     'Excel, or use the files as they are (listed under Data types).'),
    (r'^missing eventDate', 'Missing dates', 'error',
     'Take the sampling dates from the paper / supplement (Methods, station table) or ask the authors; add them '
     'as a column or a mapping override and re-run the converter. A Zenodo "collected" date also works.'),
    (r'^missing decimalLatitude|^missing decimalLongitude', 'Missing coordinates', 'error',
     'Take the station positions from the paper / supplement or ask the authors; a station table (station + '
     'lat + lon) joined by station code, or a place-name column resolvable in Marine Regions, is enough.'),
    (r'^missing scientificNameID', 'Taxa not matched in WoRMS', 'error',
     'Check the unmatched names with the WoRMS Taxon Match tool (typos, abbreviations, codes, authorship in '
     'the name) and correct them in a lookup table or mapping override.'),
    (r'names in the taxon column .* are taxa in WoRMS', 'Column mistaken for taxa', 'warning',
     'That column holds codes / sites / metadata rather than taxa; set the right taxon column in a mapping '
     'override.'),
    (r'taxon names not matched in WoRMS', 'Taxa not matched in WoRMS', 'warning',
     'Correct the names (WoRMS Taxon Match) or map vernacular names to scientific names.'),
    (r'scientificName .* taken from the dataset title', 'Taxon from title', 'info',
     'The species comes from the title (single-species dataset): check that every row is that species.'),
    (r'invalid coordinates', 'Invalid coordinates', 'warning',
     'Check the coordinate columns (degrees-minutes written as decimals, swapped lat / lon, missing sign '
     'for W / S, projected coordinates).'),
    (r'unparseable dates', 'Invalid dates', 'warning',
     'Write the dates in ISO 8601 (YYYY-MM-DD); check day / month order and text values.'),
    (r'only the first \d+ MB read|truncated at|zip too large|read remotely', 'Partly read (size)', 'info',
     'Large files were read only in part; run the converter locally with --large-files (and more disk / '
     'memory) to take everything.'),
    (r'derived from the sampling period stated in the paper', 'Date from paper', 'info',
     'eventDate taken from the paper (labelled in eventRemarks): check it; reject it with date_apply=false in '
     'paper_values.json.'),
    (r'derived from the study area of the paper', 'Coordinates from paper', 'info',
     'Coordinates = centre of the study area in the paper, with coordinateUncertaintyInMeters: replace them '
     'with station positions when available.'),
    (r'derived from the place name|derived from the collection period', 'Derived value', 'info',
     'Value derived from the place name (Marine Regions) or the Zenodo record, labelled in the archive: '
     'check it.'),
    (r'no quantitative measurements', 'Presence only', 'info',
     'Presence-only data: fine for EMODnet Biology; add abundance / biomass if the authors have it.'),
    (r'taxon reference list', 'Reference list skipped', 'info', 'Checklist / code list used as a lookup, not data.'),
    (r'columns joined from|took the date of the same sample|downloaded whole|taken from a station table',
     'Processing note', 'info', 'Done automatically (joins, station / sample look-ups); no action needed.'),
    (r'duplicate occurrence rows removed', 'Duplicate rows removed', 'info',
     'The same rows appeared in several tables / files and were kept once; check that the source files do not '
     'repeat the data.'),
    (r'conflicting occurrenceIDs made unique', 'occurrenceID made unique', 'warning',
     'The source IDs were not unique (several taxa / rows per ID) and were suffixed; give one ID per record '
     'in the source if possible.'),
    (r'event-level measurements dropped', 'Measurements dropped', 'warning',
     'Measurements of events that have no occurrence were dropped; keep them as a separate environmental '
     'table (Data types export).'),
    # ---- QC checks (stage 2)
    (r'^D01', 'No readable table', 'error', 'The archive could not be read; re-run the converter.'),
    (r'^D02', 'Missing core fields', None,
     'Complete date, coordinates and taxon for every record (see the missing-term fixes).'),
    (r'^D03', 'DwC terms', None, 'Add the missing Darwin Core identifiers / terms (occurrenceID, basisOfRecord...).'),
    (r'^D04', 'Invalid coordinates', None,
     'Coordinates out of range or (0,0): fix the values or the column mapping.'),
    (r'^D05', 'Outside declared region', None,
     'Positions outside the region of the inventory: check the sign of the longitude / swapped lat-lon, or '
     'correct the region.'),
    (r'^D06|located on land|Marine taxon located on land', 'Points on land', None,
     'Check sign of longitude (W must be negative), swapped lat / lon, degrees-minutes read as decimals; '
     'coastal / lagoon stations can be moved to the nearest water point with a larger '
     'coordinateUncertaintyInMeters.'),
    (r'^D07', 'Invalid dates', None, 'Fix unparseable or future dates (ISO 8601), check day / month order.'),
    (r'^D08|greater than the value found in the bathymetry', 'Depth vs bathymetry', None,
     'Check depth units (m vs cm), sign and min / max order; ROV / seamount stations can legitimately differ '
     'from the coarse bathymetry - document it in eventRemarks.'),
    (r'^D09', 'Taxa not matched in WoRMS', None,
     'Correct the names with the WoRMS Taxon Match tool.'),
    (r'^D10', 'Non-marine taxa', None,
     'Terrestrial / freshwater taxa: check identifications or remove them from the marine dataset.'),
    (r'^D11|scientificNameID does not resolve', 'scientificNameID', None,
     'Use WoRMS LSIDs (urn:lsid:marinespecies.org:taxname:<AphiaID>) for every record.'),
    (r'^D12', 'Vocabularies', None, 'Use the controlled values for basisOfRecord, occurrenceStatus, sex, lifeStage.'),
    (r'^D13', 'Duplicate occurrenceID', None, 'Make occurrenceID unique (dataset + event + taxon + row).'),
    (r'^D14|measurementTypeID is missing|measurementValueID|measurementUnitID', 'Measurement vocabularies', None,
     'Give measurementTypeID / measurementUnitID from the NERC vocabularies (P01 / P06) in the eMoF.'),
    (r'^D15', 'Wide layout', None, 'Wide tables were reshaped to long format; check the result.'),
    (r'^M01', 'Not open', None, 'Embargoed / restricted: wait for the embargo or ask the authors for access.'),
    (r'^M02', 'Licence', None,
     'NC / ND / no licence limits reuse: ask the authors to use CC-BY or CC0 for the data.'),
    (r'^M03', 'No DOI', None, 'Ask for a DOI (Zenodo provides one on publication).'),
    (r'^M04', 'No contact', None, 'Find a contact e-mail (paper corresponding author) for questions.'),
    (r'^M05', 'No temporal coverage', None, 'Add the sampling period to the metadata.'),
    (r'^M06', 'No tabular files', None, 'Ask the authors for tabular data (CSV / Excel).'),
    (r'^M08', 'Poor documentation', None, 'Ask for a README / data dictionary (columns, units, methods).'),
    # ---- EMODnetBiocheck
    (r'Empty value for required field (eventDate)', 'Missing dates', None,
     'Some events have no date: complete them (same sample / cruise date, paper) or drop them.'),
    (r'Empty value for required field (decimalLatitude|decimalLongitude)', 'Missing coordinates', None,
     'Some events have no position: complete them from the station list or drop them.'),
    (r'Required field (decimalLatitude|decimalLongitude) is missing', 'Missing coordinates', None,
     'No coordinates in the dataset: station positions from the paper / authors.'),
    (r'Required field eventDate is missing', 'Missing dates', None,
     'No dates in the dataset: sampling dates from the paper / authors.'),
    (r'(Empty value for required field|Required field) scientificNameID', 'Taxa not matched in WoRMS', None,
     'Names without a WoRMS match: correct them (WoRMS Taxon Match).'),
    (r'Required field institutionCode is missing', 'institutionCode', None,
     'Add the institution (EDMO code) in institutions.json of the converter.'),
    (r'parentEventID has no corresponding eventID', 'Event hierarchy', None,
     'Add the parent events (cruise / station) or remove parentEventID.'),
    (r'Potential duplicate record', 'Possible duplicates', None,
     'Identical records: check whether they are replicates (add a replicate / sample ID) or true duplicates.'),
    (r'Recommended field (minimum|maximum)DepthInMeters', 'Depth missing', None,
     'Add sampling depth (or depth range) if known.'),
    (r'No sampling (instrument|descriptors)', 'Sampling description', None,
     'Add samplingProtocol, gear and effort (sampleSizeValue / Unit) from the methods.'),
    (r'Biological value of 0 while occurrenceStatus is present', 'Zero abundance as present', None,
     'Rows with 0 must have occurrenceStatus = absent.'),
    (r'Duplicate measurementType', 'Duplicate measurements', None,
     'The same measurementType appears twice for one occurrence / event: distinguish them (units, method).'),
    (r'basisOfRecord', 'basisOfRecord', None, 'Use a valid basisOfRecord (HumanObservation, PreservedSpecimen...).'),
    # ---- papers / data types
    (r'^paper: closed', 'Paper not open access', 'info',
     'The linked paper is not open access: read it through a library subscription or ask the authors for '
     'the station list and sampling dates.'),
    (r'^file not read', 'Files not processed', 'info',
     'Larger than the size limit or inside a very large zip: process locally with higher limits '
     '(--large-files, --max-mb, --max-nc-mb).'),
]
RULES = [(re.compile(rx, re.I), cat, sev, fix) for rx, cat, sev, fix in RULES]
SEV_ORDER = {'error': 0, 'warning': 1, 'info': 2}


def classify(message, severity):
    for rx, cat, sev, fix in RULES:
        if rx.search(message):
            return cat, sev or severity, fix
    return 'Other', severity, 'Review the message and the source file.'


def load(path, default=None):
    try:
        with open(path, encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def record_issues(rid, conv, qc2, qc1, paper, files):
    out = []

    def add(stage, message, severity):
        cat, sev, fix = classify(message, severity)
        key = (stage, cat, message)
        if key not in {(i['stage'], i['category'], i['problem']) for i in out}:
            out.append({'stage': stage, 'severity': sev, 'category': cat, 'problem': message, 'fix': fix})

    if conv:
        for term in conv.get('missing_required') or []:
            add('Conversion', f'missing {term} (required term empty in the whole dataset)', 'error')
        for i in conv.get('issues') or []:
            add('Conversion', i, 'warning')
    if qc2:
        for c in qc2.get('checks') or []:
            if c.get('result') in ('fail', 'warn'):
                fail = c['result'] == 'fail' and CRIT_SEV.get(c['id']) == 'error'
                add('QC (archive)', f'{c["id"]} {CRIT.get(c["id"], "")}: {c.get("details", "")}',
                    'error' if fail else 'warning')
        bc = qc2.get('biocheck') or {}
        for m in bc.get('error_messages') or []:
            add('EMODnetBiocheck', m, 'error')
        for m in bc.get('warning_messages') or []:
            add('EMODnetBiocheck', m, 'warning')
    elif qc1:
        for c in qc1.get('checks') or []:
            if c.get('result') == 'fail':
                add('QC (raw files)', f'{c["id"]} {CRIT.get(c["id"], "")}: {c.get("details", "")}', 'warning')
    if paper and paper.get('papers'):
        closed = [p['doi'] for p in paper['papers'] if 'not available' in p.get('fulltext', '')]
        if closed and paper.get('status') != 'values found':
            add('Paper', f'paper: closed access ({", ".join(closed[:2])}) - dates / stations not extracted', 'info')
    skipped = [f for f in files if not f.get('processed') and re.search(r'larger|not downloaded|listed only|'
                                                                        r'zip larger', f.get('note', ''))]
    if skipped:
        add('Data types', f'file not read: {len(skipped)} file(s) over the size limits '
                          f'(e.g. {skipped[0]["file"][:80]})', 'info')
    out.sort(key=lambda i: (SEV_ORDER[i['severity']], i['stage'], i['category']))
    return out


CRIT = {}
CRIT_SEV = {}


def main(argv=None):
    ap = argparse.ArgumentParser(description='Errors and possible fixes for every dataset')
    ap.add_argument('--dwc', default=os.path.join(HERE, '..', 'zenodo_dwc_repo', 'output'),
                    help='output folder of the Darwin Core converter (Task 3)')
    ap.add_argument('--out', default=os.path.join(HERE, 'output'))
    args = ap.parse_args(argv)
    conv = {str(d['record_id']): d for d in load(os.path.join(args.dwc, 'conversion_results.json'))['datasets']}
    q2 = load(os.path.join(args.out, 'qc_dwca_results.json'), {'datasets': [], 'criteria': []})
    q1 = load(os.path.join(args.out, 'qc_results.json'), {'datasets': [], 'criteria': []})
    CRIT.update({c['id']: c['name'] for c in q2.get('criteria', []) + q1.get('criteria', [])})
    CRIT_SEV.update({c['id']: c['severity'] for c in q2.get('criteria', []) + q1.get('criteria', [])})
    qc2 = {str(d['record_id']): d for d in q2['datasets']}
    qc1 = {str(d['record_id']): d for d in q1['datasets']}
    papers = load(os.path.join(args.dwc, 'paper_values.json'), {})
    files = {}
    for f in (load(os.path.join(args.dwc, 'data_types.json'), {}) or {}).get('files', []):
        files.setdefault(str(f['record_id']), []).append(f)

    datasets, rows = [], []
    for rid, c in conv.items():
        issues = record_issues(rid, c, qc2.get(rid), qc1.get(rid), papers.get(rid), files.get(rid, []))
        n = Counter(i['severity'] for i in issues)
        d = {'record_id': rid, 'title': c.get('title', ''), 'doi': c.get('doi', ''), 'url': c.get('url', ''),
             'conversion_status': c.get('status', ''), 'qc_status': (qc2.get(rid) or {}).get('qc_status', ''),
             'errors': n['error'], 'warnings': n['warning'], 'infos': n['info'],
             'categories': sorted({i['category'] for i in issues if i['severity'] != 'info'}),
             'main_fix': next((i['fix'] for i in issues if i['severity'] == 'error'), ''),
             'issues': issues}
        datasets.append(d)
        rows += [{'record_id': rid, 'title': d['title'], **i} for i in issues]
    cats = Counter(i['category'] for d in datasets for i in d['issues'] if i['severity'] == 'error')
    by_dataset = Counter(c for d in datasets for c in {i['category'] for i in d['issues'] if i['severity'] == 'error'})
    summary = {'datasets': len(datasets), 'with_errors': sum(1 for d in datasets if d['errors']),
               'with_warnings_only': sum(1 for d in datasets if not d['errors'] and d['warnings']),
               'clean': sum(1 for d in datasets if not d['errors'] and not d['warnings']),
               'error_categories': dict(by_dataset.most_common()), 'error_messages': sum(cats.values())}
    payload = {'generated': date.today().isoformat(), 'summary': summary, 'datasets': datasets}
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, 'dataset_issues.json'), 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    with pd.ExcelWriter(os.path.join(args.out, 'Dataset_Issues.xlsx'), engine='xlsxwriter') as xw:
        pd.DataFrame([{'datasets with this error': v, 'error category': k} for k, v in by_dataset.most_common()]
                     ).to_excel(xw, sheet_name='Summary', index=False)
        pd.DataFrame([{k: (', '.join(v) if isinstance(v, list) else v) for k, v in d.items() if k != 'issues'}
                      for d in datasets]).to_excel(xw, sheet_name='Datasets', index=False)
        pd.DataFrame(rows).to_excel(xw, sheet_name='Issues', index=False)
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
