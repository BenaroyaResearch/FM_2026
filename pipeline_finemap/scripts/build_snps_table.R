library(data.table)

study_id <- snakemake@wildcards[["study"]]
trait <- snakemake@params[["trait"]]
window_bp <- as.integer(snakemake@params[["window_bp"]])

associations <- fread(snakemake@input[["associations"]])
results_files <- snakemake@input[["susie_results"]]

build_locus_rows <- function(f) {
  locus_id <- basename(dirname(dirname(dirname(f))))

  if (!file.exists(f) || file.info(f)$size == 0) {
    message(sprintf("build_snps_table: skipping locus %s — %s is missing or empty", locus_id, f))
    return(NULL)
  }
  tab <- fread(f)
  if (nrow(tab) == 0) {
    message(sprintf("build_snps_table: skipping locus %s — %s has a header but no data rows", locus_id, f))
    return(NULL)
  }

  lead_bp_row <- associations[SNP == locus_id]
  if (nrow(lead_bp_row) == 0) {
    message(sprintf("build_snps_table: skipping locus %s — no matching lead SNP in %s", locus_id, snakemake@input[["associations"]]))
    return(NULL)
  }
  lead_bp <- lead_bp_row$BP[1]

  cs_bps <- sprintf("%d, %d", lead_bp + window_bp, lead_bp - window_bp)

  data.table(
    TRAIT = trait,
    STUDY = study_id,
    LOCUS = locus_id,
    CS_BPS = cs_bps,
    SNP = tab$SNP,
    CHR = tab$CHR,
    BP = tab$BP,
    REF = tab$A1,
    ALT = tab$A2,
    MAF = tab$MAF,
    PVAL = tab$PVAL,
    Z = tab$Z,
    PIP = tab$PIP,
    CS = tab$CS,
    CS_SIZE = tab$CS_SIZE,
    LAMBDA = tab$LAMBDA
  )
}

all_rows <- rbindlist(lapply(results_files, build_locus_rows), fill = TRUE)

# Every per-locus file failing (missing/empty/no matching lead SNP) previously wrote a silent
# 1-byte blank output instead of erroring -- discovered when this was empty for 29 of 71
# finalized studies despite their real per-locus fine-mapping output being present and correct.
# write.table() ran unconditionally on a 0-row table with no check that aggregation actually
# produced anything. Fail loudly instead -- but only for the real bug case (>=1 locus file was
# provided and every one of them failed its own checks). A study with zero eligible loci at all
# (results_files itself empty, a real and already-confirmed state for e.g. GCST90475664/
# GCST90476372 -- both genuinely have 0 eligible loci, not a broken join) is not a failure and
# should still produce a valid, empty (header-only) table rather than error.
if (length(results_files) > 0 && (is.null(all_rows) || nrow(all_rows) == 0)) {
  stop(sprintf(
    "build_snps_table: no rows produced for study %s from %d input file(s) -- every locus was skipped (see messages above). Check the per-locus susie_results.tsv files and their join against %s.",
    study_id, length(results_files), snakemake@input[["associations"]]
  ))
}
if (is.null(all_rows) || ncol(all_rows) == 0) {
  all_rows <- data.table(
    TRAIT = character(), STUDY = character(), LOCUS = character(), CS_BPS = character(),
    SNP = character(), CHR = character(), BP = integer(), REF = character(), ALT = character(),
    MAF = numeric(), PVAL = numeric(), Z = numeric(), PIP = numeric(), CS = character(),
    CS_SIZE = integer(), LAMBDA = numeric()
  )
}

write.table(
  all_rows,
  snakemake@output[[1]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
