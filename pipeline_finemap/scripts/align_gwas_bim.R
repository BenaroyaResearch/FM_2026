# Purpose: Makes sure GWAS data's row order exactly matches the LD matrix's row/column order.
# Makes sure alleles are oriented the same way in both the X/
# General workflow:
# 1. Split the GWAS in duplicate and non-duplicate positions
# 2. Recombines, forcing the final row order of the GWAS position to match that of the 1000g bim file
# 3. Checks row count matches LD matrix row count
# 4. Checks every SNP allele orientation agianst 1000g bim file. Flips A1/A2 and B wherever reference disagrees.
#    Flags anything that matches neither orientation as unresolved.
# 5. Writes final aligned table.

library(data.table)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]

gwas_subset <- fread(snakemake@input[["gwas_subset"]]) # GWAS positions
R <- as.matrix(fread(snakemake@input[["matrix"]], header = FALSE)) # LD matrix
bim <- fread(snakemake@input[["bim"]]) # Filtered 1000g bim

# Emits a WARNING when a vital SuSiE column (BETA, SE, or Z) is computed rather than read
# directly from the sumstats, so anyone reading the log knows this locus's effect-size
# columns aren't all first-hand data.
log_derived_warning <- function(column, from_columns) {
  cat(
    sprintf(
      "%s: locus %s — WARNING: %s derived from sumstat columns %s\n",
      Sys.time(),
      locus_id,
      column,
      paste(from_columns, collapse = ", ")
    ),
    file = log_path,
    append = TRUE
  )
}

# run_susie.R needs a Z-score column and nothing else -- it doesn't care how it got there.
# Only attempt to derive one (from BETA/SE, an OR-derived BETA, and/or an OR-CI-derived SE)
# when a native Z isn't already present -- otherwise a study with a perfectly good native Z
# still has BETA/SE derived from OR for no reason, which does nothing but add unused columns
# and log a misleading "derived from OR" warning for a locus whose real score was native.
if (!"Z" %in% names(gwas_subset)) {
  # Some studies' munged sumstats keep OR instead of BETA (never converted during munging).
  # log(OR) is the log-odds equivalent of BETA, and -log(OR) == log(1/OR), so deriving BETA
  # here up front makes the allele-flip logic below (which negates signed columns) correct
  # without any further change.
  if (!"BETA" %in% names(gwas_subset) && "OR" %in% names(gwas_subset)) {
    gwas_subset[, BETA := log(OR)]
    log_derived_warning("BETA", "OR")
  }

  # Some studies' munged sumstats drop SE entirely (e.g. it was 100% NA pre-munge and
  # MungeSumstats discarded it) but keep the OR's confidence interval. SE is recoverable from
  # that CI on the log-OR scale: for a 95% CI, upper/lower = exp(log(OR) +/- 1.96*SE), so
  # SE = (log(upper) - log(lower)) / (2*1.96). Scoped to the OR case specifically (not
  # attempted when only BETA is native) since there's no confirmed instance in this pipeline's
  # data of CI bounds on a beta/log scale, and assuming the wrong scale would silently
  # mis-compute SE rather than fail loudly.
  if (
    !"SE" %in% names(gwas_subset) &&
      "OR" %in% names(gwas_subset) &&
      all(c("CI_UPPER", "CI_LOWER") %in% names(gwas_subset))
  ) {
    gwas_subset[, SE := (log(CI_UPPER) - log(CI_LOWER)) / (2 * 1.96)]
    log_derived_warning("SE", c("CI_UPPER", "CI_LOWER"))
  }

  # Some studies' munged sumstats carry Z natively (handled above); others need it derived
  # from BETA/SE (including a BETA just derived from OR, and/or an SE just derived from the OR's
  # CI bounds, above). A study with neither a native SE nor usable CI bounds has no way to
  # recover a valid Z here, so that locus is excluded rather than silently mis-scored.
  if (all(c("BETA", "SE") %in% names(gwas_subset))) {
    gwas_subset[, Z := BETA / SE]
    log_derived_warning("Z", c("BETA", "SE"))
  } else {
    msg <- sprintf(
      "locus %s — ALIGNMENT_FAILED: no usable effect-size column (need Z, or BETA+SE, or OR+SE, or OR+CI_UPPER+CI_LOWER). EXCLUDED FROM FINE-MAPPING.",
      locus_id
    )
    cat(
      sprintf("%s: *** %s ***\n", Sys.time(), toupper(msg)),
      file = log_path,
      append = TRUE
    )
    message("!!! ", msg, " !!!")
    writeLines(
      c(
        "ALIGNMENT_STATUS: FAILED",
        "note: no usable effect-size column (Z, BETA+SE, OR+SE, or OR+CI_UPPER+CI_LOWER) found for this locus"
      ),
      snakemake@output[[1]]
    )
    quit(status = 0)
  }
}

gwas_subset$key <- paste(gwas_subset$CHR, gwas_subset$BP, sep = ":") # Build CHR:BP key for every row

# Catches every duplicated position in GWAS subset and separates into duplicated and non-duplicated.
# Have to check duplicated position against 1000g file.
dup_keys <- gwas_subset$key[
  duplicated(gwas_subset$key) | duplicated(gwas_subset$key, fromLast = TRUE)
] # Flags repeats after first occurrence and before last occurrence. Collects all rows sharing a position.
gwas_dups <- gwas_subset[gwas_subset$key %in% dup_keys, ]
gwas_nodups <- gwas_subset[!gwas_subset$key %in% dup_keys, ]

# Renaming bim rows to be able to match on column names with gwas subset
bim$key <- paste(bim$V1, bim$V4, sep = ":")
dedup_bim_alleles <- bim[, .(key, bim_A1 = V5, bim_A2 = V6)]

# Matches duplicate GWAS rows to the 1000g bim file on position, keeping whichever duplicate (if
# any) has alleles matching the reference panel's pair -- in EITHER orientation. A duplicate
# position usually means GWAS Catalog reported more than one alt allele at a multiallelic site as
# separate rows; only one of those rows corresponds to the specific biallelic pair the reference
# panel kept, and that row can just as validly be flip-oriented relative to it as an exact match
# (the flip itself is applied uniformly to every row, duplicate-origin or not, further down).
# Requiring an exact match here (as this used to) silently drops a duplicate row whose only
# "problem" is needing the same flip every other flipped SNP at this locus already gets -- leaving
# a position with a real GWAS entry looking like it has none once the row-order-preserving merge
# below fills the gap with a phantom, all-NA-except-BP row for it.
dup_candidates <- merge(gwas_dups, dedup_bim_alleles, by = "key")
matched_dups <- dup_candidates[
  (A1 == bim_A1 & A2 == bim_A2) | (A1 == bim_A2 & A2 == bim_A1)
]
matched_dups <- matched_dups[, names(gwas_subset), with = FALSE]

# Makes sure that column order between tables are the same for binding
gwas_nodups <- as.data.table(gwas_nodups)
matched_dups <- as.data.table(matched_dups)
setcolorder(matched_dups, names(gwas_nodups))

# Combines the cleaned duplicates and original non-duplicates
gwas_cleaned <- rbindlist(
  list(gwas_nodups, matched_dups),
  use.names = TRUE,
  fill = TRUE
)
gwas_cleaned$key <- NULL # Removes CHR:BP key

# Gathering duplicated SNPs in GWAS reference panel (~multi allelic SNPs)
bim_dup_bps <- bim$V4[duplicated(bim$V4) | duplicated(bim$V4, fromLast = TRUE)]

# Removing multi-allelic positions
if (length(bim_dup_bps) > 0) {
  cat(
    sprintf(
      "%s: locus %s — %d multi-allelic position(s) in reference panel, excluding from this locus: %s\n",
      Sys.time(),
      locus_id,
      length(bim_dup_bps),
      paste(unique(bim_dup_bps), collapse = ", ")
    ),
    file = log_path,
    append = TRUE
  )
  bim <- bim[!V4 %in% bim_dup_bps]
  gwas_cleaned <- gwas_cleaned[!BP %in% bim_dup_bps]
}

# Force gwas_cleaned into the exact same row order as bim/R. Allows or correct Z-score-SNP comparison downstream.
bim_order <- bim[, .(BP = V4)] # Captures original bim order
bim_order[, row_num := .I] # Assigns index
gwas_cleaned <- merge(bim_order, gwas_cleaned, by = "BP", all.x = TRUE) # Joins gwas_cleaned onto bim, matching rows where BP is equal in both
setorder(gwas_cleaned, row_num)
gwas_cleaned[, row_num := NULL]

# If the row counts of the LD matrix and the cleaned GWAS subset don't match, raise dimensional mismatch warning
# Proceed to fallback alignment: only keep GWAS rows with matching positions somewhere in .bim file
# Checks for same count but not order
if (nrow(data.frame(R)) != nrow(gwas_cleaned)) {
  cat(
    sprintf(
      "%s: locus %s — dimensional mismatch (R matrix: %d rows, GWAS cleaned: %d rows), attempting fallback alignment\n",
      Sys.time(),
      locus_id,
      nrow(data.frame(R)),
      nrow(gwas_cleaned)
    ),
    file = log_path,
    append = TRUE
  )
  gwas_subset_final <- subset(gwas_cleaned, BP %in% bim$V4)
} else {
  gwas_subset_final <- gwas_cleaned
}

# If alignment still doesn't match, exlcude from FM
if (nrow(gwas_subset_final) != nrow(R)) {
  msg <- sprintf(
    "locus %s — ALIGNMENT_FAILED: dimensions still mismatched after fallback (GWAS: %d rows, LD matrix: %d rows). EXCLUDED FROM FINE-MAPPING.",
    locus_id,
    nrow(gwas_subset_final),
    nrow(R)
  )
  cat(
    sprintf("%s: *** %s ***\n", Sys.time(), toupper(msg)),
    file = log_path,
    append = TRUE
  )
  message("!!! ", msg, " !!!")

  writeLines(
    c(
      "ALIGNMENT_STATUS: FAILED",
      sprintf("gwas_rows: %d", nrow(gwas_subset_final)),
      sprintf("ld_matrix_rows: %d", nrow(R)),
      "note: dimensions could not be reconciled, locus excluded from fine-mapping"
    ),
    snakemake@output[[1]]
  )
  quit(status = 0)
}

# Rebuild key and attach reference alleles
bim_alleles <- bim[, .(key, bim_A1 = V5, bim_A2 = V6)]
# Rebuild key
gwas_subset_final$key <- paste(
  gwas_subset_final$CHR,
  gwas_subset_final$BP,
  sep = ":"
)
# merge() sorts its output by the join key -- "key" is a CHR:BP string, so that sort is
# lexicographic (e.g. "10:100000235" < "10:99024977"), not numeric, and silently undoes the
# row order established above to match the LD matrix. Preserve it explicitly across the merge.
gwas_subset_final[, .orig_row := .I]
# Combine two tables on position so that both GWAS and bim alleles are present
gwas_subset_final <- merge(
  gwas_subset_final,
  bim_alleles,
  by = "key",
  all.x = TRUE
)
setorder(gwas_subset_final, .orig_row)
gwas_subset_final[, .orig_row := NULL]

# Marks rows where GWAS alleles are exactly swapped relative to reference. Flips Beta sign when reference alleles are flipped.
flip_idx <- with(
  gwas_subset_final,
  A1 == bim_A2 & A2 == bim_A1 & !is.na(bim_A1)
)
n_flipped <- sum(flip_idx, na.rm = TRUE) # Sum of SNPs with flipped alleles


if (n_flipped > 0) {
  # Flips whichever signed effect column(s) this locus actually has -- BETA (native or
  # derived from OR above) and/or Z (native or derived from BETA/SE above) -- rather than
  # assuming BETA specifically exists.
  for (signed_col in intersect(c("BETA", "Z"), names(gwas_subset_final))) {
    gwas_subset_final[[signed_col]][flip_idx] <- -gwas_subset_final[[signed_col]][flip_idx]
  }
  tmp_a1 <- gwas_subset_final$A1[flip_idx]
  gwas_subset_final$A1[flip_idx] <- gwas_subset_final$A2[flip_idx] # Assign A2 to A1
  gwas_subset_final$A2[flip_idx] <- tmp_a1 # Assign A1 to A2

  cat(
    sprintf(
      "%s: locus %s — flipped %d of %d SNPs (%.1f%%) to match reference panel allele orientation\n",
      Sys.time(),
      locus_id,
      n_flipped,
      nrow(gwas_subset_final),
      100 * n_flipped / nrow(gwas_subset_final)
    ),
    file = log_path,
    append = TRUE
  )
}

# Checks for SNPs that match neither the original nor flipped orientation. Doesn't do anything about them.
unresolved <- with(
  gwas_subset_final,
  !is.na(bim_A1) &
    !((A1 == bim_A1 & A2 == bim_A2) | (A1 == bim_A2 & A2 == bim_A1))
)
if (any(unresolved)) {
  cat(
    sprintf(
      "%s: locus %s — WARNING: %d SNP(s) match neither allele orientation in reference, check for strand issues\n",
      Sys.time(),
      locus_id,
      sum(unresolved)
    ),
    file = log_path,
    append = TRUE
  )
}

gwas_subset_final$key <- NULL
gwas_subset_final$bim_A1 <- NULL
gwas_subset_final$bim_A2 <- NULL

# run_susie.R pairs each row of this file with the corresponding row of R purely by position, so
# row order has to exactly match bim's (the order R was built in) -- not just row count. Checked
# directly rather than trusted, since the merges above are exactly the kind of operation (joins
# that silently re-sort by their key) that can undo the ordering established earlier without
# changing row count at all.
if (!identical(gwas_subset_final$BP, bim$V4)) {
  msg <- sprintf(
    "locus %s — ALIGNMENT_FAILED: final row order does not match the LD matrix's row order (bim). EXCLUDED FROM FINE-MAPPING.",
    locus_id
  )
  cat(sprintf("%s: *** %s ***\n", Sys.time(), toupper(msg)), file = log_path, append = TRUE)
  message("!!! ", msg, " !!!")
  writeLines(
    c(
      "ALIGNMENT_STATUS: FAILED",
      "note: final row order does not match the LD matrix row order, locus excluded from fine-mapping"
    ),
    snakemake@output[[1]]
  )
  quit(status = 0)
}

cat(
  sprintf(
    "%s: locus %s — aligned GWAS subset: %d SNPs, dimensions match LD matrix (%d x %d)\n",
    Sys.time(),
    locus_id,
    nrow(gwas_subset_final),
    nrow(R),
    ncol(R)
  ),
  file = log_path,
  append = TRUE
)

writeLines(
  c("ALIGNMENT_STATUS: OK", sprintf("n_snps: %d", nrow(gwas_subset_final))),
  paste0(snakemake@output[[1]], ".status")
)

write.table(
  gwas_subset_final,
  snakemake@output[[1]],
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)
