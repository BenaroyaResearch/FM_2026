# Purpose: for high-LD neighbors of the lead SNP that ld_neighbor_check found missing from the
# GWAS subset (dropped during munging QC — strand ambiguity, allele mismatch with reference,
# etc.), check whether that exact CHR:BP survived in one of the munge pipeline's cut-SNP side
# files, verify its alleles actually agree with the 1000G reference panel at that position (in
# either orientation), and splice it back into the subset only if so. Matches on CHR:BP; never
# invents data, and never rescues a SNP whose alleles disagree with the reference.

library(data.table)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]

gwas_subset <- fread(snakemake@input[["gwas_subset"]])
ld_lines <- readLines(snakemake@input[["ld_check"]])
cut_files <- snakemake@input[["cut_files"]]
bim <- fread(snakemake@input[["bim"]])

# CHR:BP -> reference alleles. A BP with more than one bim row (multi-allelic in the wider
# region window) can't be matched unambiguously, so it's dropped from the lookup entirely.
bim_lookup <- bim[, .(BP = V4, bim_A1 = V5, bim_A2 = V6)]
bim_dup_bps <- bim_lookup$BP[duplicated(bim_lookup$BP) | duplicated(bim_lookup$BP, fromLast = TRUE)]
bim_lookup <- bim_lookup[!BP %in% bim_dup_bps]

write_report <- function(dt) {
  fwrite(dt, snakemake@output[["rescued_snps_report"]], sep = "\t")
}

report_cols <- c("SNP", "CHR", "BP", "A1", "A2", "bim_A1", "bim_A2", "source_file")
empty_report <- data.table(
  SNP = character(),
  CHR = character(),
  BP = numeric(),
  A1 = character(),
  A2 = character(),
  bim_A1 = character(),
  bim_A2 = character(),
  source_file = character()
)

table_start <- which(grepl("^SNP_B\\tBP_B\\tR2", ld_lines))

# No missing-neighbor table at all (PASS / NOT_IN_REFERENCE / empty .ld) — nothing to rescue,
# pass the subset through unchanged.
if (length(table_start) == 0) {
  cat(
    sprintf(
      "%s: locus %s — no missing high-LD neighbors to rescue, passing GWAS subset through unchanged\n",
      Sys.time(),
      locus_id
    ),
    file = log_path,
    append = TRUE
  )
  fwrite(gwas_subset, snakemake@output[["gwas_subset_rescued"]], sep = "\t")
  write_report(empty_report)
  quit(status = 0)
}

data_lines <- ld_lines[(table_start + 1):length(ld_lines)]
missing <- rbindlist(lapply(strsplit(data_lines, "\t"), function(r) {
  data.table(SNP_B = r[1], BP_B = as.numeric(r[2]), R2 = as.numeric(r[3]))
}))

lead_chr <- gwas_subset$CHR[1]
required_cols <- names(gwas_subset)

# Read every cut-SNP file once, restrict immediately to this locus's chromosome, and skip any
# file that doesn't have the columns needed to rebuild a valid GWAS subset row. Keeps track of
# which file each row came from for the rescued-SNPs report.
cut_pool_parts <- lapply(cut_files, function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || !all(required_cols %in% names(dt))) {
    return(NULL)
  }
  dt <- dt[CHR == lead_chr, ..required_cols]
  dt[, source_file := basename(f)]
  dt
})
cut_pool_parts <- cut_pool_parts[!sapply(cut_pool_parts, is.null)]
cut_pool <- if (length(cut_pool_parts) > 0) {
  unique(rbindlist(cut_pool_parts), by = "BP")
} else {
  data.table()
}

candidates_dt <- if (nrow(cut_pool) > 0) {
  cut_pool[BP %in% missing$BP_B]
} else {
  cut_pool
}

n_missing <- nrow(missing)

# Verify every candidate's alleles actually agree with the 1000G reference at that position,
# in either orientation, before letting it back in. A candidate whose BP isn't in the reference
# at all (multi-allelic, or genuinely absent) is treated the same as an allele mismatch.
if (nrow(candidates_dt) > 0) {
  candidates_dt <- merge(candidates_dt, bim_lookup, by = "BP", all.x = TRUE)
  allele_match <- with(
    candidates_dt,
    !is.na(bim_A1) &
      ((A1 == bim_A1 & A2 == bim_A2) | (A1 == bim_A2 & A2 == bim_A1))
  )
  rescued_dt <- candidates_dt[allele_match]
  mismatched_dt <- candidates_dt[!allele_match]
} else {
  rescued_dt <- candidates_dt
  mismatched_dt <- candidates_dt
}

n_found <- nrow(rescued_dt)
n_mismatched <- nrow(mismatched_dt)

if (n_found > 0) {
  for (i in seq_len(n_found)) {
    cat(
      sprintf(
        "%s: locus %s — WARNING: grabbing cut SNP %s (%s:%d) from munge files\n",
        Sys.time(),
        locus_id,
        rescued_dt$SNP[i],
        rescued_dt$CHR[i],
        rescued_dt$BP[i]
      ),
      file = log_path,
      append = TRUE
    )
  }
}

if (n_mismatched > 0) {
  for (i in seq_len(n_mismatched)) {
    ref_alleles <- if (is.na(mismatched_dt$bim_A1[i])) {
      "not found in reference"
    } else {
      sprintf("%s/%s", mismatched_dt$bim_A1[i], mismatched_dt$bim_A2[i])
    }
    cat(
      sprintf(
        "%s: locus %s — WARNING: cut SNP %s (%s:%d) found in munge files but alleles (%s/%s) do not match 1000G reference (%s) — excluding from rescue\n",
        Sys.time(),
        locus_id,
        mismatched_dt$SNP[i],
        mismatched_dt$CHR[i],
        mismatched_dt$BP[i],
        mismatched_dt$A1[i],
        mismatched_dt$A2[i],
        ref_alleles
      ),
      file = log_path,
      append = TRUE
    )
  }
}

gwas_combined <- if (n_found > 0) {
  rbindlist(list(gwas_subset, rescued_dt[, ..required_cols]), use.names = TRUE)
} else {
  gwas_subset
}

cat(
  sprintf(
    "%s: locus %s — rescued %d of %d missing high-LD neighbor(s) from %d cut-SNP file(s) (%d not found in any cut file, %d excluded for allele mismatch with reference)\n",
    Sys.time(),
    locus_id,
    n_found,
    n_missing,
    length(cut_files),
    n_missing - n_found - n_mismatched,
    n_mismatched
  ),
  file = log_path,
  append = TRUE
)

fwrite(gwas_combined, snakemake@output[["gwas_subset_rescued"]], sep = "\t")

if (n_found > 0) {
  write_report(rescued_dt[, ..report_cols])
} else {
  write_report(empty_report)
}
