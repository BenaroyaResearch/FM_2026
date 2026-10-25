library(data.table)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]
cfg <- snakemake@config[["locus_summary_thresholds"]]

# --- small helper: pull "key: value" out of a vector of summary lines ---
get_field <- function(lines, key) {
  pattern <- paste0("^", key, ":\\s*")
  match_line <- lines[grepl(pattern, lines)]
  if (length(match_line) == 0) {
    return(NA_character_)
  }
  trimws(sub(pattern, "", match_line[1]))
}

format_list <- function(x) {
  if (length(x) == 0) "none" else paste(x, collapse = ", ")
}

result <- tryCatch(
  {
    susie_lines <- readLines(snakemake@input[["susie_summary"]])
    qc_status <- get_field(susie_lines, "QC_STATUS")

    # --- hard stop: SuSiE never actually produced a real result for this locus ---
    if (is.na(qc_status) || qc_status == "NOT_RUN") {
      cat(
        sprintf(
          "%s: locus %s — NOT_RUN: SuSiE did not run (skipped or failed upstream)\n",
          Sys.time(),
          locus_id
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "STATUS: NOT_RUN",
          sprintf("LOCUS: %s", locus_id),
          "reason: SuSiE did not produce a result for this locus"
        ),
        snakemake@output[[1]]
      )
      quit(status = 0)
    }

    susie_results <- fread(snakemake@input[["susie_results"]])

    # === 1. convergence ===
    niter <- suppressWarnings(as.integer(get_field(susie_lines, "niter")))
    converged_lt50 <- if (is.na(niter)) NA else niter < 50

    # === 2/3. missing LD neighbors, bucketed by R2 ===
    ld_lines <- readLines(snakemake@input[["ld_check"]])
    table_start <- which(grepl("^SNP_B\\tBP_B\\tR2", ld_lines))

    n_missing_high <- 0
    n_missing_low <- 0

    if (length(table_start) > 0 && table_start < length(ld_lines)) {
      data_lines <- ld_lines[(table_start + 1):length(ld_lines)]
      r2_vals <- as.numeric(sapply(strsplit(data_lines, "\t"), `[`, 3))
      n_missing_high <- sum(r2_vals >= 0.8, na.rm = TRUE)
      n_missing_low <- sum(r2_vals >= 0.6 & r2_vals < 0.8, na.rm = TRUE)
    }

    # === 4. SNPs rescued ===
    rescued_input <- snakemake@input[["rescued_snps"]]
    n_rescued <- 0
    if (!is.null(rescued_input) && length(rescued_input) > 0 && file.exists(rescued_input) && file.info(rescued_input)$size > 0) {
      n_rescued <- nrow(fread(rescued_input))
    }

    # === 5/6. credible sets and their sizes ===
    n_cs <- as.integer(get_field(susie_lines, "n_credible_sets"))
    cs_table <- unique(susie_results[!is.na(CS), .(CS, CS_SIZE)])
    cs_sizes_str <- if (nrow(cs_table) > 0) {
      paste(sprintf("%s=%d", cs_table$CS, cs_table$CS_SIZE), collapse = ", ")
    } else {
      "none"
    }

    # === 7. SNPs with PIP ~= 1 ===
    pip_at_one_threshold <- if (!is.null(cfg$pip_approx_1_threshold)) {
      cfg$pip_approx_1_threshold
    } else {
      1.00
    }
    pip_one_snps <- susie_results[PIP >= pip_at_one_threshold, SNP]

    # === 8. lambda ===
    lambda_est <- as.numeric(get_field(susie_lines, "lambda"))

    # === 9. percent SNPs failing QC ===
    qc_failed_rate <- as.numeric(get_field(susie_lines, "qc_failed_rate"))

    # === 10. red SNPs (per-SNP kriging outliers) ===
    red_snps <- if ("QC_OUTLIER" %in% names(susie_results)) {
      susie_results[QC_OUTLIER == TRUE, SNP]
    } else {
      character()
    }

    # === 11/12. OpenTargets comparison ===
    pip_full <- fread(snakemake@input[["pip_full"]])
    locus_bps <- susie_results$BP
    span <- range(locus_bps, na.rm = TRUE)

    # R2 vs OT: same convention as build_locus_pip_plots.R (locus column is reliable for
    # match_status == "both", since it's populated from the "my" side of the join)
    both_rows <- pip_full[locus == locus_id & match_status == "both"]
    r2_vs_ot <- if (nrow(both_rows) >= 3) {
      sprintf("%.4f", cor(both_rows$my_pip, both_rows$ot_pip)^2)
    } else {
      "not available (fewer than 3 matched SNPs)"
    }

    ot_high_pip_threshold <- if (!is.null(cfg$opentargets_high_pip_threshold)) {
      cfg$opentargets_high_pip_threshold
    } else {
      0.9
    }

    # NOTE: pip_full is a genome-wide outer join on BP (see build_opentargets_comparison.R).
    # The `locus` column is only populated from the "my" side of that join, so it is NA for
    # every match_status == "OT_only" row -- filtering on `locus == locus_id` would silently
    # return zero OT-only SNPs for every locus. Instead, OT-only SNPs are scoped to this locus
    # by genomic position: this locus's own BP span (min/max BP of its own SNPs), not a true
    # locus join. This is an accepted approximation (BP-only, no CHR key, consistent with how
    # the rest of the pipeline already joins on BP) -- loci with wide or overlapping windows
    # could in principle double-count or misattribute an OT-only SNP.
    locus_rows <- pip_full[BP %in% locus_bps]
    ot_only_in_locus <- pip_full[
      match_status == "OT_only" & BP >= span[1] & BP <= span[2]
    ]

    ot_absent_from_locus <- ot_only_in_locus[ot_pip > ot_high_pip_threshold]
    ot_present_mine_low <- locus_rows[
      match_status == "both" & ot_pip > ot_high_pip_threshold & my_pip <= ot_high_pip_threshold
    ]
    mine_absent_from_ot <- locus_rows[
      match_status == "mine_only" & my_pip > ot_high_pip_threshold
    ]
    mine_present_ot_low <- locus_rows[
      match_status == "both" & my_pip > ot_high_pip_threshold & ot_pip <= ot_high_pip_threshold
    ]

    # === SNPs belonging to a credible set (final table) ===
    cs_snps <- susie_results[!is.na(CS)]

    # === write report ===
    out_lines <- c(
      "STATUS: OK",
      sprintf("LOCUS: %s", locus_id),
      sprintf("Generated: %s", Sys.time()),
      "",
      "--- Convergence ---",
      sprintf(
        "SUSIE_CONVERGED_UNDER_50_ITER: %s",
        if (is.na(converged_lt50)) "NA" else converged_lt50
      ),
      sprintf("niter: %s", if (is.na(niter)) "NA" else niter),
      "",
      "--- LD neighbors ---",
      sprintf("MISSING_LD_NEIGHBORS_R2_GE_0.8: %d", n_missing_high),
      sprintf("MISSING_LD_NEIGHBORS_R2_0.6_TO_0.8: %d", n_missing_low),
      sprintf("SNPS_RESCUED: %d", n_rescued),
      "",
      "--- Credible sets ---",
      sprintf("N_CREDIBLE_SETS: %s", if (is.na(n_cs)) "NA" else n_cs),
      sprintf("CS_SIZES: %s", cs_sizes_str),
      "",
      "--- PIP ---",
      sprintf("N_SNPS_PIP_APPROX_1: %d", length(pip_one_snps)),
      sprintf("SNPS_PIP_APPROX_1: %s", format_list(pip_one_snps)),
      "",
      "--- Model fit ---",
      sprintf("LAMBDA: %s", if (is.na(lambda_est)) "NA" else sprintf("%.4f", lambda_est)),
      "",
      "--- QC ---",
      sprintf(
        "PCT_SNPS_FAILING_QC: %s",
        if (is.na(qc_failed_rate)) "NA" else sprintf("%.2f%%", 100 * qc_failed_rate)
      ),
      sprintf("N_RED_SNPS: %d", length(red_snps)),
      sprintf("RED_SNPS: %s", format_list(red_snps)),
      "",
      "--- OpenTargets comparison ---",
      sprintf("R2_VS_OPENTARGETS: %s", r2_vs_ot),
      "",
      sprintf(
        "OT_HIGH_PIP_NOT_PICKED_UP_BY_ME (OT_pip>%.1f, absent from this locus entirely): %d",
        ot_high_pip_threshold,
        nrow(ot_absent_from_locus)
      ),
      sprintf(
        "OT_HIGH_PIP_PRESENT_BUT_MY_PIP_LOW (both present, ot_pip>%.1f, my_pip<=%.1f): %d",
        ot_high_pip_threshold,
        ot_high_pip_threshold,
        nrow(ot_present_mine_low)
      ),
      sprintf("  SNPS: %s", format_list(ot_present_mine_low$SNP)),
      sprintf(
        "MY_HIGH_PIP_NOT_IN_OPENTARGETS (my_pip>%.1f, absent from OT entirely): %d",
        ot_high_pip_threshold,
        nrow(mine_absent_from_ot)
      ),
      sprintf("  SNPS: %s", format_list(mine_absent_from_ot$SNP)),
      sprintf(
        "MY_HIGH_PIP_PRESENT_BUT_OT_PIP_LOW (both present, my_pip>%.1f, ot_pip<=%.1f): %d",
        ot_high_pip_threshold,
        ot_high_pip_threshold,
        nrow(mine_present_ot_low)
      ),
      sprintf("  SNPS: %s", format_list(mine_present_ot_low$SNP)),
      "",
      "--- SNPs, PIPs, and credible sets ---",
      "SNP\tCS\tPIP",
      sprintf("%s\t%s\t%s", cs_snps$SNP, cs_snps$CS, cs_snps$PIP)
    )

    writeLines(out_lines, snakemake@output[[1]])

    cat(
      sprintf(
        "%s: locus %s — summary report written (n_cs=%s, niter=%s, lambda=%s, %d red SNPs)\n",
        Sys.time(),
        locus_id,
        if (is.na(n_cs)) "NA" else n_cs,
        if (is.na(niter)) "NA" else niter,
        if (is.na(lambda_est)) "NA" else sprintf("%.4f", lambda_est),
        length(red_snps)
      ),
      file = log_path,
      append = TRUE
    )
  },
  error = function(e) {
    cat(
      sprintf(
        "%s: locus %s — ERROR during locus summary: %s\n",
        Sys.time(),
        locus_id,
        conditionMessage(e)
      ),
      file = log_path,
      append = TRUE
    )
    writeLines(
      c(
        "STATUS: ERROR",
        sprintf("LOCUS: %s", locus_id),
        sprintf("error: %s", conditionMessage(e))
      ),
      snakemake@output[[1]]
    )
  }
)
