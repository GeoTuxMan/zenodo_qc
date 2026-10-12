# -----------------------------------------------------------------------------
# run_biocheck.R - optional R route for Task 2 (Zenodo QC)
#
# Runs the official EMODnet Biology checks on a dataset that is already in
# Darwin Core form (e.g. a DwC-A produced by Task 3 - zenodo_dwc_converter -
# or published on an IPT):
#
#   1. EMODnetBiocheck::checkdataset()  (EMODnet Biology / LifeWatch QC rules)
#   2. obistools checks                 (fields, on land, depth, dates, WoRMS)
#
# Usage (from the repository root):
#   Rscript R/run_biocheck.R --dwca path/to/dwca.zip          # zipped archive
#   Rscript R/run_biocheck.R --dir  path/to/unzipped_dwca     # folder with event/occurrence/emof
#   Rscript R/run_biocheck.R --ipt  "https://ipt.vliz.be/eurobis/resource?r=xxx"
#   Rscript R/run_biocheck.R --dir  path --out output/biocheck_xxx
#
# Outputs (in --out, default output/biocheck): biocheck_report.html (when
# EMODnetBiocheck is available), obistools_issues.csv, summary.txt
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
dwca <- get_arg("dwca")
dir  <- get_arg("dir")
ipt  <- get_arg("ipt")
out  <- get_arg("out", "output/biocheck")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

need <- function(pkg) suppressWarnings(requireNamespace(pkg, quietly = TRUE))

# --- install helpers (never block: the obistools route works without Biocheck)
if (!need("remotes")) try(install.packages("remotes", repos = "https://cloud.r-project.org"), silent = TRUE)
if (!need("obistools") && need("remotes")) try(remotes::install_github("iobis/obistools"), silent = TRUE)
if (!need("EMODnetBiocheck") && need("remotes")) {
  message("Installing EMODnetBiocheck from GitHub (EMODnet/EMODnetBiocheck)...")
  try(remotes::install_github("EMODnet/EMODnetBiocheck", upgrade = "never"), silent = TRUE)
}
has_biocheck <- need("EMODnetBiocheck")
has_obistools <- need("obistools")

# --- read the Darwin Core tables ---------------------------------------------
read_dwc <- function(folder) {
  f <- list.files(folder, full.names = TRUE, recursive = TRUE)
  pick <- function(pattern) {
    x <- f[grepl(pattern, basename(f), ignore.case = TRUE)]
    if (length(x)) read.delim(x[1], quote = "", stringsAsFactors = FALSE, encoding = "UTF-8",
                              check.names = FALSE, na.strings = c("", "NA")) else NULL
  }
  list(event = pick("^event"), occurrence = pick("^occurrence"),
       emof = pick("measurementorfact|emof"))
}

if (!is.null(dwca)) {
  dir <- file.path(tempdir(), "dwca")
  unzip(dwca, exdir = dir)
}
tables <- if (!is.null(dir)) read_dwc(dir) else list(event = NULL, occurrence = NULL, emof = NULL)

# In a DwC-A star schema the extensions link to the core through the "id" column.
# EMODnetBiocheck expects an explicit eventID in Occurrence / eMoF, so fill it from
# the core id when missing (otherwise every eMoF row is reported as unlinked).
if (!is.null(tables$event)) {
  for (nm in c("occurrence", "emof")) {
    t <- tables[[nm]]
    if (!is.null(t) && "id" %in% names(t) && (!"eventID" %in% names(t) || all(is.na(t$eventID)))) {
      t$eventID <- t$id
      tables[[nm]] <- t
    }
  }
}

summary_lines <- c(sprintf("run_biocheck.R - %s", format(Sys.time())))

# --- 1. EMODnetBiocheck ------------------------------------------------------
if (has_biocheck) {
  library(EMODnetBiocheck)
  res <- tryCatch({
    if (!is.null(ipt)) {
      checkdataset(IPTreport = ipt)
    } else {
      checkdataset(Event = tables$event, Occurrence = tables$occurrence, eMoF = tables$emof)
    }
  }, error = function(e) { message("checkdataset failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    saveRDS(res, file.path(out, "biocheck_result.rds"))
    if (!is.null(res$datacheck)) {
      writeLines(as.character(res$datacheck), file.path(out, "biocheck_datacheck.txt"))
    }
    for (nm in names(res)) {
      if (is.data.frame(res[[nm]])) write.csv(res[[nm]], file.path(out, paste0("biocheck_", nm, ".csv")),
                                              row.names = FALSE)
    }
    summary_lines <- c(summary_lines, "EMODnetBiocheck: OK (see biocheck_*.csv / .rds)")
  }
} else {
  summary_lines <- c(summary_lines, "EMODnetBiocheck: not installed - skipped (obistools checks only)")
}

# --- 2. obistools checks -----------------------------------------------------
if (has_obistools && !is.null(tables$occurrence)) {
  library(obistools)
  occ <- tables$occurrence
  if (!is.null(tables$event) && "id" %in% names(occ) && "id" %in% names(tables$event)) {
    # occurrence columns left empty are inherited from the event (DwC-A star schema)
    occ <- occ[, colSums(!is.na(occ)) > 0 | names(occ) == "id", drop = FALSE]
    occ <- merge(occ, tables$event[, setdiff(names(tables$event), setdiff(names(occ), "id"))],
                 by = "id", all.x = TRUE)
  }
  issues <- list()
  add <- function(name, x) {
    if (!is.null(x) && NROW(x)) issues[[name]] <<- data.frame(check = name, as.data.frame(x))
  }
  add("check_fields", tryCatch(check_fields(occ), error = function(e) NULL))
  if (all(c("decimalLongitude", "decimalLatitude") %in% names(occ))) {
    occ$decimalLongitude <- as.numeric(occ$decimalLongitude)
    occ$decimalLatitude  <- as.numeric(occ$decimalLatitude)
    land <- tryCatch(check_onland(occ, report = TRUE, buffer = 1000), error = function(e) NULL)
    add("check_onland", land)
    add("check_depth", tryCatch(check_depth(occ, report = TRUE), error = function(e) NULL))
  }
  if ("eventDate" %in% names(occ)) add("check_eventdate", tryCatch(check_eventdate(occ), error = function(e) NULL))
  if ("scientificName" %in% names(occ)) {
    names_u <- unique(na.omit(occ$scientificName))
    m <- tryCatch(match_taxa(head(names_u, 500), ask = FALSE), error = function(e) NULL)
    if (!is.null(m)) {
      unmatched <- names_u[seq_len(min(500, length(names_u)))][is.na(m$scientificNameID)]
      add("worms_unmatched", data.frame(scientificName = unmatched))
    }
  }
  keep <- c("check", "level", "field", "row", "message", "scientificName")
  all_issues <- if (length(issues)) do.call(rbind, lapply(issues, function(d) {
    d[] <- lapply(d, as.character)
    for (k in setdiff(keep, names(d))) d[[k]] <- NA_character_
    d[, keep, drop = FALSE]
  })) else data.frame()
  write.csv(all_issues, file.path(out, "obistools_issues.csv"), row.names = FALSE)
  summary_lines <- c(summary_lines, sprintf("obistools: %d issue rows (%s)", NROW(all_issues),
                                            paste(names(issues), collapse = ", ")))
} else {
  summary_lines <- c(summary_lines, "obistools: skipped (package or occurrence table missing)")
}

# --- machine-readable summary, read by zenodo_qc.py --dwca ----------------------
n_rows <- function(x) if (is.data.frame(x)) nrow(x) else 0L
bc <- if (exists("res") && !is.null(res)) res else NULL
general <- if (!is.null(bc)) bc$dtb$general_issues else NULL
summary_json <- list(
  biocheck = !is.null(bc),
  # general_issues = structural errors (missing / unlinked IDs, invalid fields)
  errors = n_rows(general),
  # eMoF vocabulary problems, non-exclusive unit IDs, marine taxa on land
  warnings = n_rows(bc$dtb$mof_issues) + sum(vapply(bc$exclusivity_fails, n_rows, 0L)) +
    n_rows(bc$MarTaxaonLand),
  messages = if (n_rows(general)) unique(paste0(general$table, ": ", general$field, " - ", general$message)) else list(),
  mof_messages = if (n_rows(bc$dtb$mof_issues)) unique(bc$dtb$mof_issues$message) else list(),
  marine_taxa_on_land = n_rows(bc$MarTaxaonLand),
  obistools_issue_rows = if (exists("all_issues")) NROW(all_issues) else NA
)
if (need("jsonlite")) {
  jsonlite::write_json(summary_json, file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
}

writeLines(summary_lines, file.path(out, "summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
