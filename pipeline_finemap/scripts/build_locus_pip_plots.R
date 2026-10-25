library(data.table)
library(ggplot2)

log_path <- snakemake@params[["log_path"]]
finemapping_method <- snakemake@params[["finemapping_method"]]
window_label <- snakemake@params[["window_label"]]
study_id <- snakemake@wildcards[["study"]]

full_pip <- fread(snakemake@input[["pip_full"]])
both <- full_pip[match_status == "both" & !is.na(locus)]

loci <- unique(both$locus)

# --- compute per-locus R2 up front, before building any output text ---
locus_results <- list()
for (loc in loci) {
  sub <- both[locus == loc]
  n <- nrow(sub)
  r2 <- if (n >= 3) cor(sub$my_pip, sub$ot_pip)^2 else NA_real_
  locus_results[[loc]] <- list(n = n, r2 = r2)
}

r2_values <- sapply(locus_results, function(x) x$r2)
r2_values <- r2_values[!is.na(r2_values)]

n_loci_with_r2 <- length(r2_values)
n_loci_above_08 <- sum(r2_values > 0.8)
pct_above_08 <- if (n_loci_with_r2 > 0) {
  100 * n_loci_above_08 / n_loci_with_r2
} else {
  NA_real_
}

top_summary_line <- sprintf(
  "Loci with R2 > 0.8: %d of %d (%.1f%%)",
  n_loci_above_08,
  n_loci_with_r2,
  pct_above_08
)

# --- header, now including the same summary line ---
summary_lines <- c(
  sprintf("PER-LOCUS OPENTARGETS PIP ALIGNMENT — STUDY %s", study_id),
  sprintf("Generated: %s", Sys.time()),
  strrep("=", 70),
  sprintf("Loci with computable R2 (n>=3 matched SNPs): %d", n_loci_with_r2),
  top_summary_line,
  strrep("=", 70),
  sprintf("%-20s %8s %10s", "LOCUS", "N_SNPS", "R2"),
  strrep("-", 70)
)

# Derived from output$summary (already correctly rooted, real vs. a window-test's
# locus_output_root) rather than a hardcoded "results/..." literal, which silently wrote every
# locus's plot into the real results/ tree regardless of which root a run targeted -- discovered
# when a window-size test's plot overwrote the real GCST004131 result. Computed once, before the
# loop below shadows the base `sub()` function with a same-named data.table variable.
# Anchored on "/<study_id>/checking/" specifically (not a bare "/checking/") because a window
# test's own locus_output_root can itself contain the literal substring "checking"
# (e.g. results/checking/window_test_350kb), which a naive strip would cut at instead.
study_root <- base::sub(paste0("(/", study_id, ")/checking/.*$"), "\\1", snakemake@output[["summary"]])

# --- per-locus table + plots ---
for (loc in loci) {
  n <- locus_results[[loc]]$n
  r2 <- locus_results[[loc]]$r2

  if (is.na(r2)) {
    summary_lines <- c(
      summary_lines,
      sprintf("%-20s %8d %10s", loc, n, "NA (too few)")
    )
    next
  }

  summary_lines <- c(summary_lines, sprintf("%-20s %8d %10.3f", loc, n, r2))

  sub <- both[locus == loc]
  plot_path <- file.path(
    study_root,
    loc,
    window_label,
    finemapping_method,
    sprintf("%s_opentargets_pip_alignment.pdf", loc)
  )

  p <- ggplot(sub, aes(x = ot_pip, y = my_pip)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(alpha = 0.6, color = "steelblue") +
    annotate(
      "text",
      x = 0.05,
      y = 0.95,
      hjust = 0,
      size = 4,
      label = sprintf("R^2 = %.3f\nn = %d", r2, n)
    ) +
    labs(
      title = sprintf("%s — PIP alignment vs OpenTargets", loc),
      x = "OpenTargets PIP",
      y = "My PIP"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_minimal()

  ggsave(plot_path, p, width = 6, height = 5)
}

writeLines(summary_lines, snakemake@output[["summary"]])

cat(
  sprintf(
    "%s: per-locus PIP alignment — %d loci with matched SNPs, plots written to individual locus folders\n",
    Sys.time(),
    length(loci)
  ),
  file = log_path,
  append = TRUE
)
