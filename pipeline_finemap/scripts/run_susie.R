# Purpose: reads the aligned GWAS and LD matrix for one locus, checks for upstream failure markers,
# computes Z-scores and the lambda diagnostic, fits SuSiE, classifies QC status, writes results/plots,
# and handles every failure path (alignment failure, dimension mismatch, unexpected error) by producing
# placeholder outputs rather than crashing the rule.

library(data.table)
library(susieR)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]
n_samples <- as.numeric(snakemake@params[["n"]])
L <- as.integer(snakemake@params[["L"]])
#z_ld_weight <- as.numeric(snakemake@params[["z_ld_weight"]])
lambda_warn <- as.numeric(snakemake@params[["lambda_warn_threshold"]])
lambda_fail <- as.numeric(snakemake@params[["lambda_fail_threshold"]])
plot_path <- snakemake@output[["plot"]]
lambda_plot_path <- snakemake@output[["lambda_plot"]]
finemapping_method <- snakemake@params[["finemapping_method"]]

write_placeholder <- function(path, msg) {
  pdf(path)
  plot.new()
  text(0.5, 0.5, msg)
  dev.off()
}

cat(
  sprintf(
    "%s: locus %s — running fine-mapping with method: %s\n",
    Sys.time(),
    locus_id,
    finemapping_method
  ),
  file = log_path,
  append = TRUE
)

result <- tryCatch(
  {
    first_line <- readLines(snakemake@input[["aligned_gwas"]], n = 1)
    if (grepl("^ALIGNMENT_STATUS: FAILED", first_line)) {
      cat(
        sprintf(
          "%s: locus %s — skipping SuSiE, alignment previously failed\n",
          Sys.time(),
          locus_id
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "SUSIE_STATUS: SKIPPED",
          "QC_STATUS: NOT_RUN",
          "reason: alignment failed upstream"
        ),
        snakemake@output[["summary"]]
      )
      file.create(snakemake@output[["results"]])
      write_placeholder(plot_path, "No plot — alignment failed upstream")
      write_placeholder(lambda_plot_path, "No plot — alignment failed upstream")
      quit(status = 0)
    }

    gwas <- fread(snakemake@input[["aligned_gwas"]])
    R <- as.matrix(fread(snakemake@input[["matrix"]], header = FALSE))
    freq <- fread(snakemake@input[["freq"]])
    matrix_snplist <- readLines(snakemake@input[["snplist"]])

    if (nrow(gwas) != nrow(R)) {
      cat(
        sprintf(
          "%s: locus %s — ERROR: dimension mismatch at SuSiE step (GWAS: %d, LD: %d rows)\n",
          Sys.time(),
          locus_id,
          nrow(gwas),
          nrow(R)
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "SUSIE_STATUS: FAILED",
          "QC_STATUS: NOT_RUN",
          "reason: dimension mismatch"
        ),
        snakemake@output[["summary"]]
      )
      file.create(snakemake@output[["results"]])
      write_placeholder(plot_path, "No plot — dimension mismatch")
      write_placeholder(lambda_plot_path, "No plot — dimension mismatch")
      quit(status = 0)
    }

    # matrix_snplist records the exact row/column order the LD matrix R was built in (plink's
    # --write-snplist, one "CHR:BP" per line). z_scores below is paired with R purely by row
    # position, so if aligned_gwas's row order ever drifts from this -- matching row counts but
    # not row order, exactly what a stray order-scrambling merge upstream produces silently --
    # every SNP's Z gets checked against the wrong SNP's LD row. estimate_s_rss/lambda cannot be
    # trusted to catch this in general (a locus can drift only partially and still fit), so the
    # order is verified directly here rather than inferred from the diagnostic it would corrupt.
    gwas_key <- paste(gwas$CHR, gwas$BP, sep = ":")
    if (!identical(gwas_key, matrix_snplist)) {
      n_mismatch <- sum(gwas_key != matrix_snplist)
      cat(
        sprintf(
          "%s: locus %s — ERROR: aligned GWAS row order does not match the LD matrix's row order (%d of %d positions mismatched). Not safe to pair Z-scores with R.\n",
          Sys.time(),
          locus_id,
          n_mismatch,
          length(matrix_snplist)
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "SUSIE_STATUS: FAILED",
          "QC_STATUS: NOT_RUN",
          "reason: aligned GWAS row order does not match the LD matrix row order"
        ),
        snakemake@output[["summary"]]
      )
      file.create(snakemake@output[["results"]])
      write_placeholder(plot_path, "No plot — row order mismatch between GWAS and LD matrix")
      write_placeholder(lambda_plot_path, "No plot — row order mismatch between GWAS and LD matrix")
      quit(status = 0)
    }

    # A reference-panel SNP with a genuine allele frequency but zero genotype variance (e.g.
    # every sample heterozygous -- the classic signature of a segmental-duplication/paralogous-
    # mapping artifact, which easily clears the --maf 0.01 filter since allele frequency and
    # genotype variance are not the same thing) makes its Pearson correlation with every other
    # SNP mathematically undefined, so plink's --r square silently writes NaN for its entire row
    # and column. susie_rss/estimate_s_rss then fail with a generic "infinite or missing values
    # in 'x'", giving no indication of which SNP or why -- checked directly here instead so the
    # locus fails with a specific, actionable reason.
    # A degenerate SNP's own self-correlation (the diagonal) is 0/0 = NaN, same as its
    # correlation with everything else -- checking the diagonal identifies exactly the source
    # SNP(s), rather than every row that merely *contains* one of their NaN entries (which, in a
    # symmetric matrix, ends up being nearly the whole matrix once even one SNP is degenerate).
    bad_snp_idx <- which(!is.finite(diag(R)))
    if (length(bad_snp_idx) > 0) {
      bad_snps <- gwas$SNP[bad_snp_idx]
      cat(
        sprintf(
          "%s: locus %s — ERROR: LD matrix has non-finite (NaN/Inf) values for %d SNP(s), likely a reference-panel artifact (zero genotype variance despite nonzero allele frequency): %s\n",
          Sys.time(),
          locus_id,
          length(bad_snps),
          paste(bad_snps, collapse = ", ")
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "SUSIE_STATUS: FAILED",
          "QC_STATUS: NOT_RUN",
          sprintf(
            "reason: LD matrix has non-finite (NaN/Inf) values for %d SNP(s), likely a reference-panel artifact (zero genotype variance despite nonzero allele frequency): %s",
            length(bad_snps),
            paste(bad_snps, collapse = ", ")
          )
        ),
        snakemake@output[["summary"]]
      )
      file.create(snakemake@output[["results"]])
      write_placeholder(plot_path, "No plot — LD matrix contains non-finite (NaN/Inf) values")
      write_placeholder(lambda_plot_path, "No plot — LD matrix contains non-finite (NaN/Inf) values")
      quit(status = 0)
    }

    # align_gwas_bim.R already guarantees a correctly-signed Z column here -- computed from
    # whichever of native Z, BETA/SE, or OR/SE this study's sumstats actually had, and
    # already flipped for allele-swapped SNPs. No need to re-derive it from BETA/SE.
    z_scores <- gwas$Z
    lambda_est <- estimate_s_rss(z_scores, R, n = n_samples)

    finemapping_method <- snakemake@params[["finemapping_method"]]

    unmappable_effects <- switch(
      finemapping_method,
      "susieR" = "none",
      "susie-inf" = "inf",
      "susie-ash" = "ash",
      stop(sprintf("Unknown finemapping_method: %s", finemapping_method))
    )

    if (unmappable_effects == "none") {
      fitted <- susie_rss(
        z = z_scores,
        R = R,
        n = n_samples,
        L = L
      )
    } else {
      fitted <- susie_rss(
        z = z_scores,
        R = R,
        n = n_samples,
        L = L,
        unmappable_effects = unmappable_effects,
        estimate_residual_variance = as.logical(snakemake@params[[
          "estimate_residual_variance"
        ]])
      )
    }

    niter <- if (!is.null(fitted$niter)) as.integer(fitted$niter) else NA_integer_

    converged <- isTRUE(fitted$converged)

    pip <- fitted$pip
    cs_sets <- fitted$sets$cs

    cs_membership <- rep(NA_character_, nrow(gwas))
    cs_size <- rep(NA_integer_, nrow(gwas))
    if (!is.null(cs_sets) && length(cs_sets) > 0) {
      for (cs_name in names(cs_sets)) {
        idx <- cs_sets[[cs_name]]
        cs_membership[idx] <- cs_name
        cs_size[idx] <- length(idx)
      }
    }

    largest_cs_size <- if (length(cs_sets) > 0) {
      max(sapply(cs_sets, length))
    } else {
      NA_integer_
    }

    converged <- isTRUE(fitted$converged)
    convergence_method <- if (unmappable_effects == "none") "ELBO" else "PIP"

    freq[, BP := as.integer(sub(".*:", "", SNP))] # extract position from the chr:pos ID
    maf_lookup <- setNames(freq$MAF, freq$BP)
    maf_vals <- maf_lookup[as.character(gwas$BP)]

    out <- data.table(
      SNP = gwas$SNP,
      CHR = gwas$CHR,
      BP = gwas$BP,
      A1 = gwas$A1,
      A2 = gwas$A2,
      MAF = maf_vals,
      PVAL = gwas$P,
      Z = z_scores,
      PIP = pip,
      CS = cs_membership,
      CS_SIZE = cs_size,
      LAMBDA = lambda_est
    )
    out[, logLR := NA_real_]
    out[, z_std_diff := NA_real_]
    out[, KRIGING_Z := NA_real_]
    out[, QC_OUTLIER := NA]

    n_cs <- length(cs_sets)

    # --- QC status: judgment-call thresholds, adjustable via config ---
    qc_flags <- character()
    if (!converged) {
      qc_flags <- c(qc_flags, "did not converge")
    }
    if (lambda_est > lambda_fail) {
      qc_status <- "QC_FAIL"
      qc_flags <- c(
        qc_flags,
        sprintf(
          "lambda %.3f exceeds fail threshold %.3f",
          lambda_est,
          lambda_fail
        )
      )
    } else if (lambda_est > lambda_warn || !converged) {
      qc_status <- "QC_WARN"
      if (lambda_est > lambda_warn) {
        qc_flags <- c(
          qc_flags,
          sprintf(
            "lambda %.3f exceeds warn threshold %.3f",
            lambda_est,
            lambda_warn
          )
        )
      }
    } else {
      qc_status <- "QC_PASS"
    }

    # major_warnings starts as a copy of qc_flags (convergence/lambda issues) and picks up
    # additional warnings below (e.g. kriging diagnostic problems) that aren't QC-status
    # inputs themselves but are still worth surfacing to a reader of the summary file.
    major_warnings <- qc_flags

    # --- PIP plot (existing) ---
    pdf(plot_path, width = 8, height = 6)
    susie_plot(
      fitted,
      y = "PIP",
      main = sprintf(
        "%s (n_cs=%d, lambda=%.3f, %s)",
        locus_id,
        n_cs,
        lambda_est,
        qc_status
      )
    )
    dev.off()

    # --- Lambda / z-score consistency diagnostic plot ---
    kriging_result <- tryCatch(
      kriging_rss(z_scores, R, n = n_samples),
      error = function(e) NULL
    )

    # --- extract red-flagged (outlier) SNP count from kriging_rss diagnostic ---
    # Criterion verified against susieR's actual source (R/susie_rss_utils.R,
    # kriging_rss): it colors a point red on the diagnostic plot when
    # `logLR > 2 & abs(z) > 2`, where `z` is kriging_rss's own internally
    # PVE-adjusted observed z-score -- NOT `z_std_diff` (the LD-inconsistency
    # residual, (z - condmean)/sqrt(condvar)), which is a different quantity and
    # selects a different SNP set. `z_std_diff` is still captured below as an
    # informational LD-inconsistency diagnostic; it just no longer drives QC_OUTLIER.

    qc_failed_snps <- NA_integer_
    qc_failed_rate <- NA_real_

    if (!is.null(kriging_result)) {
      diag_data <- kriging_result$plot$data

      if (
        !is.null(diag_data) &&
          all(c("logLR", "z", "z_std_diff") %in% names(diag_data))
      ) {
        is_outlier <- diag_data$logLR > 2 & abs(diag_data$z) > 2
        qc_failed_snps <- sum(is_outlier, na.rm = TRUE)
        qc_failed_rate <- qc_failed_snps / nrow(diag_data)

        if (nrow(diag_data) == nrow(out)) {
          out[, logLR := diag_data$logLR]
          out[, z_std_diff := diag_data$z_std_diff]
          out[, KRIGING_Z := diag_data$z]
          out[, QC_OUTLIER := is_outlier]
        } else {
          cat(
            sprintf(
              "%s: locus %s — WARNING: kriging diagnostic row count (%d) does not match SuSiE results row count (%d); leaving logLR/z_std_diff/KRIGING_Z/QC_OUTLIER as NA for this locus\n",
              Sys.time(),
              locus_id,
              nrow(diag_data),
              nrow(out)
            ),
            file = log_path,
            append = TRUE
          )
          major_warnings <- c(
            major_warnings,
            sprintf(
              "kriging diagnostic row count (%d) does not match SuSiE results row count (%d)",
              nrow(diag_data),
              nrow(out)
            )
          )
        }
      } else {
        cat(
          sprintf(
            "%s: locus %s — WARNING: kriging_rss diagnostic data missing expected columns (logLR, z, z_std_diff); available columns: %s. Cannot compute QC-failed SNP count.\n",
            Sys.time(),
            locus_id,
            if (!is.null(diag_data)) {
              paste(names(diag_data), collapse = ", ")
            } else {
              "NULL"
            }
          ),
          file = log_path,
          append = TRUE
        )
        major_warnings <- c(
          major_warnings,
          sprintf(
            "kriging_rss diagnostic data missing expected columns (logLR, z, z_std_diff); available columns: %s",
            if (!is.null(diag_data)) paste(names(diag_data), collapse = ", ") else "NULL"
          )
        )
      }
    }

    if (!is.null(kriging_result) && !is.null(kriging_result$plot)) {
      ggplot2::ggsave(
        lambda_plot_path,
        kriging_result$plot,
        width = 7,
        height = 6
      )
    } else {
      write_placeholder(
        lambda_plot_path,
        "kriging_rss diagnostic unavailable for this locus"
      )
    }

    cat(
      sprintf(
        "%s: locus %s [%s] — SuSiE completed, %s, %d credible set(s), max PIP %.4f, lambda=%.4f, converged=%s (%s)\n",
        Sys.time(),
        locus_id,
        finemapping_method,
        qc_status,
        n_cs,
        max(pip, na.rm = TRUE),
        lambda_est,
        converged,
        convergence_method
      ),
      file = log_path,
      append = TRUE
    )

    if (length(qc_flags) > 0) {
      cat(
        sprintf(
          "%s: locus %s — QC flags: %s\n",
          Sys.time(),
          locus_id,
          paste(qc_flags, collapse = "; ")
        ),
        file = log_path,
        append = TRUE
      )
    }

    writeLines(
      c(
        "SUSIE_STATUS: OK",
        sprintf("QC_STATUS: %s", qc_status),
        sprintf("convergence_method: %s", convergence_method),
        sprintf("n_credible_sets: %d", n_cs),
        sprintf(
          "largest_cs_size: %s",
          if (is.na(largest_cs_size)) "NA" else largest_cs_size
        ),
        sprintf("n_snps: %d", nrow(gwas)),
        sprintf("lambda: %.4f", lambda_est),
        sprintf("converged: %s", converged),
        sprintf("niter: %s", if (is.na(niter)) "NA" else niter),
        sprintf(
          "qc_failed_snps: %s",
          if (is.na(qc_failed_snps)) "NA" else qc_failed_snps
        ),
        sprintf(
          "qc_failed_rate: %s",
          if (is.na(qc_failed_rate)) "NA" else sprintf("%.4f", qc_failed_rate)
        ),
        sprintf(
          "qc_flags: %s",
          if (length(qc_flags) > 0) paste(qc_flags, collapse = "; ") else "none"
        ),
        sprintf(
          "WARNINGS: %s",
          if (length(major_warnings) > 0) paste(major_warnings, collapse = "; ") else "none"
        )
      ),
      snakemake@output[["summary"]]
    )

    write.table(
      out,
      snakemake@output[["results"]],
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  },
  error = function(e) {
    cat(
      sprintf(
        "%s: locus %s — ERROR during SuSiE: %s\n",
        Sys.time(),
        locus_id,
        conditionMessage(e)
      ),
      file = log_path,
      append = TRUE
    )
    writeLines(
      c(
        "SUSIE_STATUS: FAILED",
        "QC_STATUS: NOT_RUN",
        sprintf("error: %s", conditionMessage(e))
      ),
      snakemake@output[["summary"]]
    )
    file.create(snakemake@output[["results"]])
    if (!file.exists(plot_path)) {
      write_placeholder(plot_path, paste("No plot —", conditionMessage(e)))
    }
    if (!file.exists(lambda_plot_path)) {
      write_placeholder(
        lambda_plot_path,
        paste("No plot —", conditionMessage(e))
      )
    }
  }
)
