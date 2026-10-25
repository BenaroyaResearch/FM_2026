library(data.table)

log_path <- snakemake@params[["log_path"]]
study_id <- snakemake@wildcards[["study"]]

pip_full <- fread(snakemake@input[["pip_full"]])
cs_comp <- fread(snakemake@input[["cs_level"]])

summary_lines <- c(
  sprintf("OPENTARGETS ALIGNMENT SUMMARY — STUDY %s", study_id),
  sprintf("Generated: %s", Sys.time()),
  strrep("=", 70),
  ""
)

if (nrow(pip_full) == 0) {
  summary_lines <- c(
    summary_lines,
    "No OpenTargets or SuSiE PIP data found for this study."
  )
  writeLines(summary_lines, snakemake@output[["summary"]])
  quit(status = 0)
}

# --- PIP alignment stats, matched SNPs only ---
both <- pip_full[match_status == "both"]
n_both <- nrow(both)
n_mine_only <- sum(pip_full$match_status == "mine_only")
n_ot_only <- sum(pip_full$match_status == "OT_only")

summary_lines <- c(
  summary_lines,
  "PIP ALIGNMENT",
  sprintf("  SNPs present in both: %d", n_both),
  sprintf("  SNPs mine only: %d", n_mine_only),
  sprintf("  SNPs OpenTargets only: %d", n_ot_only),
  ""
)

if (n_both >= 3) {
  for (m in unique(both$ot_finemapping_method)) {
    sub <- both[ot_finemapping_method == m]
    if (nrow(sub) < 3) {
      summary_lines <- c(
        summary_lines,
        sprintf("  Method=%s: n=%d (too few for stats)", m, nrow(sub))
      )
      next
    }
    r <- cor(sub$my_pip, sub$ot_pip)
    r2 <- r^2
    mae <- mean(abs(sub$my_pip - sub$ot_pip))
    rmse <- sqrt(mean((sub$my_pip - sub$ot_pip)^2))
    summary_lines <- c(
      summary_lines,
      sprintf(
        "  Method=%s: n=%d, R=%.3f, R^2=%.3f, MAE=%.3f, RMSE=%.3f",
        m,
        nrow(sub),
        r,
        r2,
        mae,
        rmse
      )
    )
  }
  summary_lines <- c(summary_lines, "")

  both[, abs_diff := abs(my_pip - ot_pip)]
  worst <- both[order(-abs_diff)][1:min(10, .N)]
  summary_lines <- c(
    summary_lines,
    "Top disagreements (my_pip vs ot_pip):",
    sprintf(
      "  %s: mine=%.3f, OT=%.3f, diff=%.3f",
      worst$SNP,
      worst$my_pip,
      worst$ot_pip,
      worst$abs_diff
    )
  )
} else {
  summary_lines <- c(
    summary_lines,
    "Too few matched SNPs for alignment statistics."
  )
}
summary_lines <- c(summary_lines, "")

# --- per-locus R2 threshold summary ---
both_all <- pip_full[match_status == "both" & !is.na(locus)]

locus_r2_values <- sapply(unique(both_all$locus), function(loc) {
  sub <- both_all[locus == loc]
  if (nrow(sub) < 3) {
    return(NA_real_)
  }
  cor(sub$my_pip, sub$ot_pip)^2
})

n_loci_with_r2 <- sum(!is.na(locus_r2_values))
n_loci_above_08 <- sum(locus_r2_values > 0.8, na.rm = TRUE)
pct_above_08 <- if (n_loci_with_r2 > 0) {
  100 * n_loci_above_08 / n_loci_with_r2
} else {
  NA_real_
}

summary_lines <- c(
  summary_lines,
  "PER-LOCUS R2 SUMMARY",
  sprintf("  Loci with computable R2 (n>=3 matched SNPs): %d", n_loci_with_r2),
  sprintf(
    "  Loci with R2 > 0.8: %d of %d (%.1f%%)",
    n_loci_above_08,
    n_loci_with_r2,
    pct_above_08
  ),
  ""
)

# --- credible set alignment stats ---
summary_lines <- c(
  summary_lines,
  "CREDIBLE SET ALIGNMENT",
  sprintf("  My CS / OT CS pairs compared: %d", nrow(cs_comp)),
  sprintf(
    "  Distinct OT credible sets covered: %d",
    if (nrow(cs_comp) > 0) uniqueN(cs_comp$ot_cs_id) else 0
  )
)

if (nrow(cs_comp) > 0) {
  summary_lines <- c(
    summary_lines,
    sprintf(
      "  Mean Jaccard overlap: %.3f",
      mean(cs_comp$jaccard, na.rm = TRUE)
    ),
    sprintf(
      "  Median Jaccard overlap: %.3f",
      median(cs_comp$jaccard, na.rm = TRUE)
    ),
    sprintf(
      "  My mean CS size: %.1f",
      mean(cs_comp$my_cs_size, na.rm = TRUE)
    )
  )
}

writeLines(summary_lines, snakemake@output[["summary"]])
cat(
  sprintf("%s: alignment summary written\n", Sys.time()),
  file = log_path,
  append = TRUE
)
