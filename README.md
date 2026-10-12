# Zenodo QC — Task 2

Quality check of the marine biology datasets discovered on Zenodo by the
[Zenodo Inventory](https://github.com/GeoTuxMan/zenodo_inventory), against the
**EMODnet Biology data quality criteria** (EMODnet Biology WP4, Task 7.2).

EMODnet Biology performs quality control with the **LifeWatch & EMODnet Biology QC tool** and the
**EMODnetBiocheck** R package. Datasets found on Zenodo are rarely in Darwin Core, so this tool
reproduces the same criteria in Python, on the raw files as published (CSV, Excel, ZIP, DwC-A),
and tells for every dataset whether it is fit — or what is missing — for mobilisation into EMODnet.

| Task | Repository | What it does |
|---|---|---|
| 1. Inventory | [`zenodo_inventory`](https://github.com/GeoTuxMan/zenodo_inventory) | Searches Zenodo, filters for marine relevance, writes Excel + JSON |
| 2. Quality Check | **this repo** | Checks the discovered datasets against the EMODnet Biology quality criteria |
| 3. Darwin Core conversion | [`zenodo_dwc_converter`](https://github.com/GeoTuxMan/zenodo_dwc_converter) | Converts the datasets to Darwin Core Archives (Event / Occurrence / eMoF) |

Results dashboard: <https://sdn2.rmri.ro/zenodo_inventory_new/>

---

## 1. Methodology

The input is the inventory JSON of Task 1 (Included + Manual review records). Every record goes
through two levels:

### Level A — metadata QC (all records)
Uses the inventory fields and the Zenodo record API (`https://zenodo.org/api/records/<id>`).

| ID | Criterion | Severity | Rule |
|---|---|---|---|
| M01 | Open access | error | files open (not embargoed / restricted / closed) |
| M02 | Reusable licence | error | CC0 / CC-BY pass; NC, ND, SA, software licences warn; none fails |
| M03 | Persistent identifier | error | DOI present |
| M04 | Contact / provider | warning | e-mail, or creator with affiliation |
| M05 | Temporal coverage | warning | sampling years known (publication year only → warn) |
| M06 | Tabular data files | error | CSV, TSV, TXT, XLS(X), ODS or ZIP/DwC-A present |
| M07 | Inventory relevance | warning | Task 1 status "Included" (Manual review → warn) |
| M08 | Documentation | warning | description ≥ 300 characters and keywords |

### Level B — data QC (open records with tabular files)
Files ≤ 30 MB are downloaded to a shared cache (`zenodo_cache/<record_id>/`), read (delimiter
sniffing, encoding fallback, header-row detection, every Excel sheet, ZIP members, Darwin Core
Archives via `meta.xml`) and their columns are mapped to Darwin Core terms through a synonym
dictionary (`lat`, `Latitude (N)` → `decimalLatitude`; `species`, `taxon` → `scientificName`;
`date`, `sampling_date` → `eventDate`; `depth_m` → depth; …). The table with the best Darwin Core
coverage is the primary table; station tables are used for coordinates when the primary table
has none. Species-by-sample matrices ("wide" layout) are detected.

| ID | Criterion | Severity | Rule (pass / warn / fail) | EMODnetBiocheck / obistools equivalent |
|---|---|---|---|---|
| D01 | Readable data | error | at least one file parsed | — |
| D02 | Core fields: when, where, what | error | date (or year), latitude, longitude, scientific name mapped; taxon from title or coordinates in a separate table → warn | `check_fields` (required) |
| D03 | Darwin Core identifiers & terms | warning | occurrenceID, scientificNameID, occurrenceStatus, basisOfRecord | `check_fields` |
| D04 | Coordinate validity | error | ranges, no (0,0), decimal degrees | `check_fields`, range checks |
| D05 | Coordinates in declared region | warning | ≥ 90 % inside region bounding box (≥ 30 % warn) | geographic consistency |
| D06 | Points at sea | error | OBIS xylookup shoredistance < 0 = land; ≤ 5 % pass, ≤ 25 % warn | `check_onland` |
| D07 | Date validity | error | ≥ 95 % parseable and ISO 8601, none in future | `check_eventdate` |
| D08 | Depth consistency | warning | numeric, ≥ 0, min ≤ max, ≤ bathymetry × 1.1 + 50 m | `check_depth` |
| D09 | Taxon names in WoRMS | error | WoRMS `AphiaRecordsByMatchNames`: ≥ 90 % pass, ≥ 60 % warn | `match_taxa` |
| D10 | Marine / brackish taxa | warning | ≥ 80 % of matched taxa marine/brackish | WoRMS environment flags |
| D11 | scientificNameID (WoRMS LSID) | warning | `urn:lsid:marinespecies.org:taxname:<id>` | Biocheck LSID check |
| D12 | Controlled vocabularies | warning | basisOfRecord, occurrenceStatus DwC values | Biocheck vocab checks |
| D13 | Unique occurrenceID | error | no duplicates / empties | Biocheck ID checks |
| D14 | Measurements described | warning | abundance/biomass/… columns have type + unit (eMoF-ready) | Biocheck eMoF checks |
| D15 | Long (tidy) layout | warning | matrices must be reshaped before DwC | — |

Taxon names are cleaned (authorship, `sp.`, `cf.` removed) and up to 200 unique names per
dataset are matched against WoRMS (random sample above that). Up to 3000 unique positions per
dataset are sent to the OBIS xylookup service. WoRMS answers are cached in
`output/worms_cache.json`.

### Scoring and status
* Each check is `pass`, `warn`, `fail` or `na`. **Score** (0–100) = weighted mean over the
  non-`na` checks (error-severity weight 2, warning 1; pass = 1, warn = 0.5, fail = 0), computed
  separately for metadata (`metadata_score`) and data (`data_score`).
* **QC status**
  * `Pass` — no warning or failure in the data checks;
  * `Warnings` — only warnings or failed warning-level checks: convertible after fixes;
  * `Fail` — at least one error-level data check failed (e.g. no coordinates);
  * `Not checked` — embargoed / restricted, no tabular files, files > 30 MB, record unavailable.
* A short **recommendation** is written for every dataset (e.g. *"Ready for DwC conversion"*,
  *"Needs date, coordinates before DwC conversion"*, *"Already Darwin Core – candidate for
  direct EMODnet/IPT ingestion"*).

Note: missing Darwin Core identifiers (occurrenceID, basisOfRecord, …) only produce warnings,
because Task 3 generates them during conversion. Missing *information* (where, when, what)
cannot be generated and is an error.

---

## 2. Installation

```bash
git clone https://github.com/GeoTuxMan/zenodo_qc.git
cd zenodo_qc
pip install -r requirements.txt
```

Optional: a Zenodo token for higher rate limits (never commit it):

```bash
export ZENODO_TOKEN=your_token          # PowerShell: $env:ZENODO_TOKEN = "your_token"
```

## 3. Usage

```bash
# full run (all inventory records; Included first, then Manual review)
python zenodo_qc.py --inventory ../zenodo_inventory/data.json --cache zenodo_cache --out output

# metadata QC only (fast, no downloads)
python zenodo_qc.py --inventory data.json --skip-data

# selected records / first N records with data QC
python zenodo_qc.py --inventory data.json --ids 8094440,19549677 --force
python zenodo_qc.py --inventory data.json --limit 50
```

| Option | Default | Meaning |
|---|---|---|
| `--inventory` | `data.json` | Task 1 inventory JSON |
| `--cache` | `zenodo_cache` | download cache (can be shared with Task 3) |
| `--out` | `output` | output folder |
| `--limit N` | all | data-check at most N records (the others get metadata QC only) |
| `--ids` | — | comma-separated Zenodo record ids |
| `--skip-data` | off | metadata QC only |
| `--force` | off | recompute records already processed (per-record results are kept in `output/_records/` so an interrupted run resumes) |
| `--workers` | 6 | records processed in parallel (API calls are rate-limited globally) |

### Outputs (`output/`)
* **`qc_results.json`** — consumed by the dashboard:
  `{generated, tool, summary{total, data_checked, pass, warnings, fail, not_checked, avg_metadata_score, avg_data_score}, criteria[], datasets[{record_id, title, doi, url, taxa, regions, inventory_status, access, license, metadata_score, data_checked, files_checked[], n_records, n_taxa, worms_matched_pct, onland_pct, data_score, qc_status, checks[{id, result, details}], issues[], recommendation}]}`
* **`qc_results.csv`** — one row per dataset, one column per check.
* **`QC_Report_Zenodo.xlsx`** — sheets `Summary`, `Datasets`, `Checks` (long format, one row
  per dataset × check with details), `Criteria`.

## 4. Optional R route — EMODnetBiocheck / obistools

Once a dataset is in Darwin Core (from Task 3, or already published on an IPT), the official
EMODnet Biology checks can be run with [`R/run_biocheck.R`](R/run_biocheck.R):

```bash
Rscript R/run_biocheck.R --dwca output_dwca/12345_dwca.zip --out output/biocheck_12345
Rscript R/run_biocheck.R --dir  path/to/unzipped_dwca
Rscript R/run_biocheck.R --ipt  "https://ipt.vliz.be/eurobis/resource?r=ro_black_sea_fpk_2006_2015"
```

The script installs `EMODnetBiocheck` from GitHub (`remotes::install_github("EMODnet/EMODnetBiocheck")`)
and runs `checkdataset()`; if the installation fails it still runs the `obistools` checks
(`check_fields`, `check_onland`, `check_depth`, `check_eventdate`, `match_taxa`) and writes
`obistools_issues.csv` and `summary.txt`.

## 5. Latest run (12 October 2026)

Input: 673 inventory records (436 Included + 237 Manual review), anonymous API access.

| QC status | Records |
|---|---|
| Pass | 0 |
| Warnings (convertible after fixes) | 23 |
| Fail | 398 |
| Not checked | 252 (109 no tabular files, 40 no tabular content, 40 embargoed, 26 restricted, 37 files > 30 MB) |

Data-checked: 421 records (1,543 tables). Average metadata score 84.9, average data score 59.7.
Main blocker: 370 of the 398 failures miss at least one core field (date, coordinates or
scientific name) — typical of paper supplements that report aggregated or sample-level results.
No dataset reaches *Pass* because none carries the full set of Darwin Core identifiers
(occurrenceID, scientificNameID, occurrenceStatus, basisOfRecord); these are generated by Task 3.

## 6. Limitations
* Column mapping is heuristic: unusual headers (other languages, codes) may be missed — such
  datasets appear as *Fail – Needs …* and should be checked by hand.
* Files > 30 MB and archives of non-tabular formats (NetCDF, GeoJSON, shapefiles, FASTA) are not
  analysed.
* Coordinates in degrees-minutes-seconds or projected systems (UTM) are reported as invalid.
* Region bounding boxes are approximate.

## 7. License & credits
Developed at NIMRD "Grigore Antipa", Constanța, Romania, for EMODnet Biology. Uses the
[Zenodo REST API](https://developers.zenodo.org/), [WoRMS REST](https://www.marinespecies.org/rest/)
and [OBIS xylookup](https://api.obis.org/). Each dataset keeps its own license.
