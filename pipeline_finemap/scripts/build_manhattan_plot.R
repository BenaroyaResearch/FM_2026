# Purpose: LD-binned regional (Manhattan) plot for one locus -- each SNP colored by its LD r^2
# to the locus's own lead SNP (the standard LocusZoom-style bins), credible-set SNPs labeled,
# lead SNP marked with a diamond. Ported from the working-but-commented-out draft in
# presentation_figures_IBD/presentation_figs.R, which already solved the r2-to-lead computation
# and binning against the same LD_matrix files run_susie already produces.

library(data.table)
library(ggplot2)
library(ggrepel)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]
window_bp <- as.integer(snakemake@params[["window_bp"]])
plot_path <- snakemake@output[[1]]

write_placeholder <- function(path, msg) {
  pdf(path)
  plot.new()
  text(0.5, 0.5, msg)
  dev.off()
}

results_path <- snakemake@input[["susie_results"]]

if (!file.exists(results_path) || file.info(results_path)$size == 0) {
  cat(
    sprintf(
      "%s: locus %s — skipping manhattan plot, no susie results (upstream skipped or failed)\n",
      Sys.time(),
      locus_id
    ),
    file = log_path,
    append = TRUE
  )
  write_placeholder(plot_path, "No plot — no susie results (upstream skipped or failed)")
  quit(status = 0)
}

results <- tryCatch(fread(results_path), error = function(e) NULL)
if (is.null(results) || nrow(results) == 0) {
  cat(
    sprintf(
      "%s: locus %s — skipping manhattan plot, empty susie results\n",
      Sys.time(),
      locus_id
    ),
    file = log_path,
    append = TRUE
  )
  write_placeholder(plot_path, "No plot — empty susie results")
  quit(status = 0)
}

results[, neg_log10_p := -log10(PVAL)]

# --- LD to lead SNP, from the same matrix/snplist run_susie already consumes ---
# --r square (not --r2) produces signed Pearson r, so the extracted row must be squared to get r^2.
snplist <- readLines(snakemake@input[["snplist"]])
R_mat <- as.matrix(fread(snakemake@input[["matrix"]], header = FALSE))

# snplist entries are always plain "CHR:BP" (matching --write-snplist's plink output), never an
# rsID and never a "chrN:pos"-style locus label -- so locus_id itself never appears in it, and
# the lead SNP's own row can't be relied on to still exist in `results` either (it may have been
# dropped during alignment/QC even though the rest of the locus fine-mapped fine, and for a
# chr:BP-labeled locus the BP embedded in locus_id is often the pre-liftover position, not the
# resolved one actually used downstream). The one place a locus's true, resolved CHR:BP always
# lives is the associations table get_associations.R wrote -- the same source compute_ld_r2
# already uses for its own --ld-snp seed, for exactly this reason.
associations <- fread(snakemake@input[["associations"]])
lead_row <- associations[SNP == locus_id]
lead_chr <- if (nrow(lead_row) > 0) lead_row$CHR[1] else NA_character_
lead_bp <- if (nrow(lead_row) > 0) lead_row$BP[1] else NA_integer_

lead_idx <- if (!is.na(lead_bp)) which(snplist == paste0(lead_chr, ":", lead_bp)) else integer(0)

r2_available <- length(snplist) == nrow(R_mat) && length(lead_idx) == 1

if (r2_available) {
  r2_to_lead <- R_mat[lead_idx, ]^2
  bp_from_snplist <- suppressWarnings(as.integer(sub(".*:", "", snplist)))
  if (!all(is.na(bp_from_snplist))) {
    ld_dt <- data.table(BP = bp_from_snplist, r2_to_lead = r2_to_lead)
    results <- merge(results, ld_dt, by = "BP", all.x = TRUE)
  } else {
    ld_dt <- data.table(SNP = snplist, r2_to_lead = r2_to_lead)
    results <- merge(results, ld_dt, by = "SNP", all.x = TRUE)
  }
} else {
  results[, r2_to_lead := NA_real_]
  cat(
    sprintf(
      "%s: locus %s — WARNING: could not locate lead SNP in LD matrix/snplist, LD binning unavailable for this plot\n",
      Sys.time(),
      locus_id
    ),
    file = log_path,
    append = TRUE
  )
}

ld_bin_levels <- c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0.0-0.2", "no LD info")
results[, ld_bin := fifelse(
  is.na(r2_to_lead), "no LD info",
  fifelse(r2_to_lead >= 0.8, "0.8-1.0",
  fifelse(r2_to_lead >= 0.6, "0.6-0.8",
  fifelse(r2_to_lead >= 0.4, "0.4-0.6",
  fifelse(r2_to_lead >= 0.2, "0.2-0.4", "0.0-0.2")))))]
results[, ld_bin := factor(ld_bin, levels = ld_bin_levels)]

ld_colors <- setNames(
  c("#D62728", "#FF7F0E", "#2CA02C", "#1F77B4", "#9467BD", "grey70"),
  ld_bin_levels
)
legend_breaks <- setdiff(ld_bin_levels, "no LD info")

# --- window: the full locus window (lead SNP +/- window_bp), not just the significant span --
# matches the actual region that was fetched and fine-mapped, so every locus's plot shows the
# same, predictable extent rather than being cropped to whatever happened to reach r2 >= 0.2.
plot_xlim <- if (!is.na(lead_bp)) {
  c(lead_bp - window_bp, lead_bp + window_bp)
} else {
  range(results$BP)
}

# fread() infers an all-NA column (a locus with zero credible sets) as logical rather than
# character, and fifelse() requires its "yes"/"no" arguments to share a type -- cast explicitly
# so a locus with no credible sets doesn't crash the plot.
results[, CS_label := fifelse(is.na(CS), "Not in a credible set", as.character(CS))]
cs_snps <- results[CS_label != "Not in a credible set"]
# Same reasoning as the LD lookup above: the lead SNP's own row may not carry the locus_id label
# (or may be entirely absent from results), so fall back to matching on its resolved position.
lead_snp_row <- results[SNP == locus_id]
if (nrow(lead_snp_row) == 0 && !is.na(lead_bp)) {
  lead_snp_row <- results[BP == lead_bp]
}

p <- ggplot(results, aes(x = BP, y = neg_log10_p, color = ld_bin)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_point(
    data = lead_snp_row,
    shape = 23,
    size = 4,
    fill = "black",
    color = "black",
    inherit.aes = FALSE,
    mapping = aes(x = BP, y = neg_log10_p)
  )

# Above ~10 SNPs needing labels, ggrepel (even with max.overlaps = Inf) has to fling labels into
# whatever empty space it can find and connect them back with long, crossing leader lines -- the
# result is unreadable rather than merely crowded, so skip labeling entirely past that point
# instead of producing a plot that looks broken. A large, tightly-correlated credible set (dozens
# of near-identical SNPs) is a legitimate result, not a labeling problem to solve.
if (nrow(cs_snps) > 0 && nrow(cs_snps) <= 10) {
  p <- p +
    ggrepel::geom_text_repel(
      data = cs_snps,
      aes(label = SNP),
      size = 3,
      color = "black",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.5,
      segment.color = "grey50",
      min.segment.length = 0,
      show.legend = FALSE
    )
} else if (nrow(cs_snps) > 10) {
  cat(
    sprintf(
      "%s: locus %s — %d credible-set SNPs, skipping labels (readable only at <=10)\n",
      Sys.time(),
      locus_id,
      nrow(cs_snps)
    ),
    file = log_path,
    append = TRUE
  )
}

p <- p +
  scale_color_manual(
    values = ld_colors,
    name = expression(paste("LD (", r^2, ") to lead SNP")),
    breaks = legend_breaks
  ) +
  coord_cartesian(xlim = plot_xlim) +
  labs(
    title = sprintf("Regional Association Plot — %s", locus_id),
    x = "Position (BP)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(plot_path, p, width = 10, height = 7)

cat(
  sprintf("%s: locus %s — manhattan plot written\n", Sys.time(), locus_id),
  file = log_path,
  append = TRUE
)
