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

### Two stages: pre-screening and final QC

The official EMODnet Biology QC (LifeWatch QC tool, EMODnetBiocheck) is defined on **Darwin
Core** data: it checks `occurrenceID`, `scientificNameID`, `basisOfRecord`, ISO `eventDate`,
Event–Occurrence–eMoF links. Raw Zenodo files almost never carry these terms, so QC runs twice:

| Stage | Input | Purpose | Output |
|---|---|---|---|
| **1. Pre-screening** | raw files published on Zenodo | Which datasets *can* be converted (date, coordinates, taxon present?) and what is missing | `qc_results.*`, `QC_Report_Zenodo.xlsx` |
| **2. Final QC** | Darwin Core Archives produced by Task 3 | Authoritative Pass / Warnings / Fail against the EMODnet criteria + EMODnetBiocheck | `qc_dwca_results.*`, `QC_Report_DwC.xlsx` |

Workflow: Task 1 inventory → **stage 1** → Task 3 conversion of the convertible datasets →
**stage 2** final QC.

The input of stage 1 is the inventory JSON of Task 1 (Included + Manual review records). Every
record goes through two levels:

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
| D06 | Points at sea | error | OBIS xylookup: land = inside the coastline (shoredistance < 0) **and** > 5 m above sea level; inland points at sea level (lagoons, estuaries, deltas - left on land by the simplified xylookup coastline) are accepted and counted apart; ≤ 5 % pass, ≤ 25 % warn | `check_onland` |
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
  *"Needs date, coordinates before DwC conversion"*, *"Already Darwin Core"*).

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
| `--ids` | — | comma-separated Zenodo record ids (with `--dwca`, the earlier stage-2 results of the other archives are kept) |
| `--skip-data` | off | metadata QC only |
| `--force` | off | recompute records already processed (per-record results are kept in `output/_records/` so an interrupted run resumes) |
| `--workers` | 6 | records processed in parallel (API calls are rate-limited globally) |

### Stage 2 — final QC on the Darwin Core Archives (after Task 3)

```bash
# Python checks D01–D15 on every DwC-A listed in Task 3's conversion_results.json
python zenodo_qc.py --dwca ../zenodo_dwc_converter/output/conversion_results.json --out output

# ... plus EMODnetBiocheck / obistools on each archive (needs R)
python zenodo_qc.py --dwca ../zenodo_dwc_converter/output/conversion_results.json --run-biocheck
```

Metadata checks (M01–M08) are reused from the stage 1 results in the same `--out` folder.
Status rules in stage 2: any failed *error* criterion or any EMODnetBiocheck error → **Fail**;
warnings or EMODnetBiocheck warnings → **Warnings**; otherwise **Pass – meets all EMODnet
Biology criteria**.

EMODnetBiocheck's `general_issues` table has no severity column, so the messages are classified
(`classify_biocheck()` in `zenodo_qc.py`, same rule in `R/run_biocheck.R`):

| EMODnetBiocheck message | Severity |
|---|---|
| `Required field … is missing`, `Empty value for required field …` | error |
| `… has no corresponding eventID`, IDs not unique / not linked | error |
| `Recommended field … is missing`, `Empty value for recommended field …` | warning |
| `Potential duplicate record`, `No sampling instrument / descriptors present` | warning |
| `Coordinates are located on land`, depth deeper than bathymetry, marine taxa on land | warning |
| eMoF vocabulary issues (`mof_issues`), non-exclusive unit / type IDs | warning |

Tested on the published EurOBIS archive *Romanian Black Sea Phytoplankton 2006–2015*: data score
98, 0 Biocheck errors, 5 Biocheck warnings (`measurementUnitID does not resolve`) and 65 % of the
taxa flagged marine/brackish in WoRMS → *Warnings*.

### Outputs (`output/`)
* **`qc_results.json`** — consumed by the dashboard:
  `{generated, tool, summary{total, data_checked, pass, warnings, fail, not_checked, avg_metadata_score, avg_data_score}, criteria[], datasets[{record_id, title, doi, url, taxa, regions, inventory_status, access, license, metadata_score, data_checked, files_checked[], n_records, n_taxa, worms_matched_pct, onland_pct, data_score, qc_status, checks[{id, result, details}], issues[], recommendation}]}`
* **`qc_results.csv`** — one row per dataset, one column per check.
* **`QC_Report_Zenodo.xlsx`** — sheets `Summary`, `Datasets`, `Checks` (long format, one row
  per dataset × check with details), `Criteria`.
* Stage 2: **`qc_dwca_results.json`** (same schema plus `stage`, `conversion_status`,
  `raw_qc_status`, `dwca_file`, `biocheck{errors, warnings, messages…}`), `qc_dwca_results.csv`,
  `QC_Report_DwC.xlsx`; Biocheck reports per record in `output/biocheck/<record_id>/`.

## 4. Optional R route — EMODnetBiocheck / obistools

Once a dataset is in Darwin Core (from Task 3, or already published on an IPT), the official
EMODnet Biology checks can be run with [`R/run_biocheck.R`](R/run_biocheck.R):

```bash
Rscript R/run_biocheck.R --dwca output_dwca/12345_dwca.zip --out output/biocheck/12345
Rscript R/run_biocheck.R --dir  path/to/unzipped_dwca
Rscript R/run_biocheck.R --ipt  "https://ipt.vliz.be/eurobis/resource?r=ro_black_sea_fpk_2006_2015"
```

The script installs `EMODnetBiocheck` from GitHub (`remotes::install_github("EMODnet/EMODnetBiocheck")`)
and runs `checkdataset()`; if the installation fails it still runs the `obistools` checks
(`check_fields`, `check_onland`, `check_depth`, `check_eventdate`, `match_taxa`) and writes
`obistools_issues.csv`, `summary.txt` and `summary.json` (read by `zenodo_qc.py --dwca`).
In a DwC-A star schema the extensions link to the core through `id`; the script copies it into
`eventID` when missing, otherwise Biocheck reports every eMoF row as unlinked.

## 5. Latest run — stage 1 pre-screening (12 October 2026)

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
No dataset reaches *Pass* at this stage because none carries the full set of Darwin Core
identifiers (occurrenceID, scientificNameID, occurrenceStatus, basisOfRecord); these are
generated by Task 3, and the authoritative status comes from stage 2.

## 5b. Latest run — stage 2 final QC on the Darwin Core Archives (13 October 2026)

Input: the 153 DwC-A produced by Task 3 (32 *Converted*, 121 *Partial*), Python checks D01–D15
+ EMODnetBiocheck / obistools on every archive (`--dwca … --run-biocheck`).

| QC status | Archives | Main reason |
|---|---|---|
| Pass | 0 | — |
| Warnings | 16 | recommended terms only (depth, sampling instrument / effort), potential duplicates |
| Fail | 137 | see below |

| Main error | Archives | Fix |
|---|---|---|
| a required term absent from the whole dataset (eventDate, coordinates, scientificNameID) — all *Partial* archives | 121 | information is not in the Zenodo files (paper, README): add it by hand / mapping override |
| EMODnetBiocheck errors only (some rows with empty required values, `parentEventID` without a parent, depth below the bathymetry) | 9 | complete or drop the rows, review values |
| positions on land (D06) | 4 | review values |
| invalid dates, some empty core fields, taxa not in WoRMS | 3 | review values |

Average data score 84.5, average metadata score 91.2. The 16 *Warnings* archives: GeoEcoMar
meiobenthos BRIDGE-BS 2022–2024 and DOORS 2023, GeoEcoMar zooplankton BRIDGE-BS 2022–2024,
Svalbard benthic disturbance (3998049), harbour porpoise telemetry and acoustics (56708), moon
jellyfish observations, Zostera marina biomass and reflectance (19549677), Macaronesian islands
ecosystems, pelagic/benthic food pathways, Seine nursery 1995–2019, sandy beach meiofauna
(13201902), Pacific Arctic larval fish, North Atlantic deep-sea amphipods, seabird juvenile
survival, Fram Strait benthic fauna and Venice lagoon zooplankton (LTER-Italy).

`institutionCode` / `institutionID` are filled by the converter from an `institutionCode` column
of the data, from the EDMO code of a known institution among the creators (NIMRD "Grigore Antipa"
= 697, GeoEcoMar = 850, CNR-ISMAR Venice = 227), or from `institutions.json` set by the data
manager. Rows without a date take the date of the same sample/station; this moved the GeoEcoMar
zooplankton archive (15119469) to *Warnings*.

The converter fixes of 13 October 2026 (metadata / data-dictionary tables skipped, R-style and
hierarchy taxon names, relational tables joined, R `write.table` files realigned, date parsing)
took the archives from 167 to 152: 17 archives built from metadata sheets, trait tables or
non-taxon columns are no longer produced (they held no occurrences), 2 new archives were produced,
and 5 moved from *Partial* to *Converted* (MEDUSA jellyfish 20713997, deep-sea corals 18171047,
Baltic NIS 8357382, 3998049, 19549677). A re-run with `--ids` drops the earlier results of
archives that are no longer produced.

A second reconversion (hierarchy names, observation-table rule for the title taxon, derived
dates and coordinates) gave 153 archives: the false occurrences built from noise, vessel and
logger tables were removed, and meiofauna 13201902, harbour porpoise 56708, 20065629 and fin
whale 10785717 became *Converted* (the first two with *Warnings*).

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
