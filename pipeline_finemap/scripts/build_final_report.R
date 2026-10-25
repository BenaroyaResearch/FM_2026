library(data.table)

log_path <- snakemake@params[["log_path"]]
study_id <- snakemake@wildcards[["study"]]

associations <- fread(snakemake@input[["associations"]])
exclusions <- fread(snakemake@input[["exclusions"]])
# Snakemake hands back NULL rather than character(0) for a named input group that expanded
# to zero paths -- which happens when a study has zero eligible loci (every lead SNP excluded
# before summarizing, a real and already-supported outcome elsewhere in this pipeline, see the
# exclusions table below). dirname()/basename() error on NULL, so coerce here once.
summary_files <- snakemake@input[["summary_files"]]
if (is.null(summary_files)) summary_files <- character(0)

get_field <- function(lines, key) {
  pattern <- paste0("^", key, ":\\s*")
  match_line <- lines[grepl(pattern, lines)]
  if (length(match_line) == 0) {
    return(NA_character_)
  }
  trimws(sub(pattern, "", match_line[1]))
}

# Counts data rows in a gzipped TSV without loading it into R -- raw sumstats files can run
# into the millions of rows, so this shells out to zcat|wc -l and only reads back the count.
count_gz_rows <- function(path) {
  if (is.null(path) || length(path) == 0 || !isTRUE(file.exists(path)) || file.info(path)$size == 0) {
    return(NA_integer_)
  }
  raw_count <- tryCatch(
    system(sprintf("zcat %s | wc -l", shQuote(path)), intern = TRUE),
    error = function(e) NA_character_
  )
  n <- suppressWarnings(as.integer(trimws(raw_count[1])))
  if (is.na(n) || n < 1) NA_integer_ else n - 1L # subtract header row
}

# --- index summary files by locus ID, for fast lookup ---
summary_by_locus <- setNames(
  summary_files,
  basename(dirname(dirname(dirname(summary_files))))
)

# --- munging summary: SNP counts before/after munging ---
n_raw_snps <- count_gz_rows(snakemake@input[["raw_sumstats"]])
n_munged_snps <- count_gz_rows(snakemake@input[["munged_sumstats"]])

munging_summary_lines <- c(
  "MUNGING SUMMARY",
  strrep("-", 70),
  sprintf(
    "SNPS_BEFORE_MUNGING: %s",
    if (is.na(n_raw_snps)) {
      "NA (raw sumstats file missing or unreadable)"
    } else {
      format(n_raw_snps, big.mark = ",", scientific = FALSE)
    }
  ),
  sprintf(
    "SNPS_AFTER_MUNGING: %s",
    if (is.na(n_munged_snps)) {
      "NA (munged sumstats file missing or unreadable)"
    } else {
      format(n_munged_snps, big.mark = ",", scientific = FALSE)
    }
  ),
  ""
)

# --- top 20 QC_PASS SNPs by PIP, study-wide ---
susie_results_files <- snakemake@input[["susie_results"]]
if (is.null(susie_results_files)) susie_results_files <- character(0)
susie_summary_files <- snakemake@input[["susie_summary"]]
if (is.null(susie_summary_files)) susie_summary_files <- character(0)

susie_summary_by_locus <- setNames(
  susie_summary_files,
  basename(dirname(dirname(dirname(susie_summary_files))))
)
susie_results_by_locus <- setNames(
  susie_results_files,
  basename(dirname(dirname(dirname(susie_results_files))))
)

qc_pass_loci <- character(0)
qc_warn_fail_dt <- data.table(LOCUS = character(), QC_STATUS = character(), REASON = character())
for (loc in names(susie_summary_by_locus)) {
  loc_summary_lines <- readLines(susie_summary_by_locus[[loc]])
  qc_status_loc <- get_field(loc_summary_lines, "QC_STATUS")
  if (is.na(qc_status_loc)) next
  if (qc_status_loc == "QC_PASS") {
    qc_pass_loci <- c(qc_pass_loci, loc)
  } else if (qc_status_loc %in% c("QC_WARN", "QC_FAIL")) {
    qc_flags_loc <- get_field(loc_summary_lines, "qc_flags")
    qc_warn_fail_dt <- rbind(
      qc_warn_fail_dt,
      data.table(
        LOCUS = loc,
        QC_STATUS = qc_status_loc,
        REASON = if (is.na(qc_flags_loc)) "none" else qc_flags_loc
      )
    )
  }
}

top20_dt <- data.table(SNP = character(), LOCUS = character(), PIP = numeric())
if (length(qc_pass_loci) > 0) {
  all_pass_rows <- rbindlist(
    lapply(qc_pass_loci, function(loc) {
      f <- susie_results_by_locus[[loc]]
      if (is.null(f) || !file.exists(f) || file.info(f)$size == 0) {
        return(NULL)
      }
      tab <- fread(f, select = c("SNP", "PIP"))
      if (nrow(tab) == 0) {
        return(NULL)
      }
      tab[, LOCUS := loc]
      tab
    }),
    fill = TRUE
  )
  if (!is.null(all_pass_rows) && nrow(all_pass_rows) > 0) {
    setorder(all_pass_rows, -PIP)
    top20_dt <- head(all_pass_rows, 20)
  }
}

top20_lines <- c("TOP 20 QC_PASS SNPs BY PIP (study-wide)", strrep("-", 70))
if (nrow(top20_dt) == 0) {
  top20_lines <- c(
    top20_lines,
    "No QC_PASS loci with SNP-level PIP data available for this study."
  )
} else {
  top20_lines <- c(
    top20_lines,
    "SNP\tLOCUS\tPIP",
    sprintf("%s\t%s\t%s", top20_dt$SNP, top20_dt$LOCUS, top20_dt$PIP)
  )
}
top20_lines <- c(top20_lines, "")

# --- loci flagged QC_WARN or QC_FAIL, study-wide ---
if (nrow(qc_warn_fail_dt) > 0) {
  setorder(qc_warn_fail_dt, QC_STATUS) # "QC_FAIL" sorts before "QC_WARN" ascending -- more severe first
}
qc_warn_fail_lines <- c("LOCI WITH QC_WARN OR QC_FAIL", strrep("-", 70))
if (nrow(qc_warn_fail_dt) == 0) {
  qc_warn_fail_lines <- c(qc_warn_fail_lines, "No loci flagged QC_WARN or QC_FAIL for this study.")
} else {
  qc_warn_fail_lines <- c(
    qc_warn_fail_lines,
    "LOCUS\tQC_STATUS\tREASON",
    sprintf(
      "%s\t%s\t%s",
      qc_warn_fail_dt$LOCUS,
      qc_warn_fail_dt$QC_STATUS,
      qc_warn_fail_dt$REASON
    )
  )
}
qc_warn_fail_lines <- c(qc_warn_fail_lines, "")

report_lines <- c(
  sprintf("FINAL REPORT — STUDY %s", study_id),
  sprintf("Generated: %s", Sys.time()),
  sprintf("Total lead associations: %d", nrow(associations)),
  strrep("=", 70),
  "",
  munging_summary_lines,
  top20_lines,
  qc_warn_fail_lines
)

for (i in seq_len(nrow(associations))) {
  snp <- associations$SNP[i]
  chr <- associations$CHR[i]
  bp <- associations$BP[i]

  report_lines <- c(
    report_lines,
    sprintf("LOCUS: %s (chr%s:%s)", snp, chr, bp),
    strrep("-", 70)
  )

  excl_row <- exclusions[SNP == snp]

  # --- case 1: excluded before ever reaching summarization ---
  if (nrow(excl_row) > 0 && isTRUE(excl_row$EXCLUDE[1])) {
    report_lines <- c(
      report_lines,
      "FINE-MAPPED: NO",
      sprintf("REASON: %s", excl_row$EXCLUSION_REASON[1]),
      sprintf(
        "  LD_FLAG: %s | ALIGN_STATUS: %s",
        excl_row$LD_FLAG[1],
        excl_row$ALIGN_STATUS[1]
      ),
      ""
    )
    next
  }

  summary_file <- summary_by_locus[[snp]]

  # --- case 2: no summary file exists at all (shouldn't normally happen for a
  #     non-excluded locus, but handled explicitly rather than assumed impossible) ---
  if (is.null(summary_file) || !file.exists(summary_file)) {
    report_lines <- c(
      report_lines,
      "FINE-MAPPED: NO",
      "REASON: no summary record found for this locus (unexpected — investigate)",
      ""
    )
    next
  }

  summary_lines <- readLines(summary_file)
  status <- get_field(summary_lines, "STATUS")

  # --- case 3: SuSiE never ran, or the summary script itself errored ---
  if (is.na(status) || status %in% c("NOT_RUN", "ERROR")) {
    reason <- if (identical(status, "ERROR")) {
      get_field(summary_lines, "error")
    } else {
      get_field(summary_lines, "reason")
    }
    report_lines <- c(
      report_lines,
      "FINE-MAPPED: NO",
      sprintf("REASON: %s", if (is.na(reason)) "unknown" else reason),
      ""
    )
    next
  }

  # --- case 4: locus was fine-mapped and summarized — embed the summary report verbatim,
  #     dropping its own STATUS/LOCUS/Generated header lines since they duplicate the
  #     "LOCUS: ..." heading already written above ---
  body_lines <- summary_lines[!grepl("^(STATUS|LOCUS|Generated):", summary_lines)]

  report_lines <- c(
    report_lines,
    "FINE-MAPPED: YES",
    body_lines,
    ""
  )
}

writeLines(report_lines, snakemake@output[[1]])

cat(
  sprintf(
    "%s: final report written for %s (%d loci)\n",
    Sys.time(),
    study_id,
    nrow(associations)
  ),
  file = log_path,
  append = TRUE
)
