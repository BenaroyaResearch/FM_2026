# Purpose: produces the record of which loci got excluded from fine-mapping and why.
# Combines LD neighbor flag with alignment outcome

library(data.table)

log_path <- snakemake@params[["log_path"]]
associations <- fread(snakemake@input[["associations"]])
ld_check_files <- snakemake@input[["ld_checks"]]
alignment_files <- snakemake@input[["alignments"]]

# Extract LD flags from ld_check_files
ld_flags <- data.table(SNP = character(), LD_FLAG = character())
for (f in ld_check_files) {
  lines <- readLines(f)
  flag_line <- lines[grepl("^FLAG:", lines)]
  flag <- trimws(sub("^FLAG:\\s*", "", flag_line[1]))
  locus_id <- basename(dirname(dirname(dirname(f))))
  ld_flags <- rbind(ld_flags, data.table(SNP = locus_id, LD_FLAG = flag))
}

# Checks for alignment fail and records
align_status <- data.table(SNP = character(), ALIGN_STATUS = character())
for (f in alignment_files) {
  first_line <- readLines(f, n = 1)
  locus_id <- basename(dirname(dirname(f)))
  status <- if (grepl("^ALIGNMENT_STATUS: FAILED", first_line)) {
    "FAILED"
  } else {
    "OK"
  }
  align_status <- rbind(
    align_status,
    data.table(SNP = locus_id, ALIGN_STATUS = status)
  )
}

# Aligning LD status and align status with SNPs
exclusions <- merge(associations, ld_flags, by = "SNP", all.x = TRUE)
exclusions <- merge(exclusions, align_status, by = "SNP", all.x = TRUE)

# Loci absent from align_status never reached that step at all -- expected for the two flags
# that stop a locus at the filter_eligible_loci checkpoint (NOT_IN_REFERENCE, CHECK_FAILED);
# for any other flag, missing alignment output means something fell through a pipeline gap.
exclusions[is.na(LD_FLAG), LD_FLAG := "UNKNOWN"]
exclusions[
  is.na(ALIGN_STATUS) & !(LD_FLAG %in% c("NOT_IN_REFERENCE", "CHECK_FAILED")),
  ALIGN_STATUS := "NOT_REACHED"
]
exclusions[
  is.na(ALIGN_STATUS) & LD_FLAG %in% c("NOT_IN_REFERENCE", "CHECK_FAILED"),
  ALIGN_STATUS := "N/A"
]

exclusions[,
  EXCLUDE := (LD_FLAG %in% c("NOT_IN_REFERENCE", "CHECK_FAILED")) | (ALIGN_STATUS == "FAILED")
]

exclusions[,
  EXCLUSION_REASON := fifelse(
    LD_FLAG == "NOT_IN_REFERENCE",
    "lead SNP not in 1000G reference panel",
    fifelse(
      LD_FLAG == "CHECK_FAILED",
      "LD-neighbor check errored; missing-neighbor status could not be verified",
      fifelse(
        ALIGN_STATUS == "FAILED",
        "GWAS/LD-matrix dimension mismatch after alignment",
        fifelse(
          ALIGN_STATUS == "NOT_REACHED",
          "did not reach alignment step (check pipeline for gaps)",
          "none"
        )
      )
    )
  )
]

n_excluded <- sum(exclusions$EXCLUDE)
n_total <- nrow(exclusions)

cat(
  sprintf(
    "%s: exclusions table built — %d of %d loci excluded from fine-mapping\n",
    Sys.time(),
    n_excluded,
    n_total
  ),
  file = log_path,
  append = TRUE
)

if (any(exclusions$ALIGN_STATUS == "NOT_REACHED")) {
  not_reached <- exclusions[ALIGN_STATUS == "NOT_REACHED", SNP]
  cat(
    sprintf(
      "%s: *** WARNING: %d locus/loci passed LD check but never reached alignment: %s — investigate before trusting this table ***\n",
      Sys.time(),
      length(not_reached),
      paste(not_reached, collapse = ", ")
    ),
    file = log_path,
    append = TRUE
  )
  message(
    "!!! LOCI PASSED LD CHECK BUT MISSING ALIGNMENT OUTPUT: ",
    paste(not_reached, collapse = ", "),
    " !!!"
  )
}

write.table(
  exclusions,
  snakemake@output[[1]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
