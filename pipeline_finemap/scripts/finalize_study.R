library(data.table)

log_path <- snakemake@params[["log_path"]]
study_id <- snakemake@wildcards[["study"]]

exclusions <- fread(snakemake@input[["exclusions"]])
summary_files <- snakemake@input[["summary_files"]]
if (is.null(summary_files)) summary_files <- character(0)

cat(
  sprintf("%s: ==== STUDY %s PIPELINE COMPLETE ====\n", Sys.time(), study_id),
  file = log_path,
  append = TRUE
)
cat(
  sprintf(
    "%s: total loci: %d | excluded: %d | fine-mapped: %d\n",
    Sys.time(),
    nrow(exclusions),
    sum(exclusions$EXCLUDE),
    length(summary_files)
  ),
  file = log_path,
  append = TRUE
)
cat(
  sprintf("%s: ======================================\n", Sys.time()),
  file = log_path,
  append = TRUE
)

finemapping_method <- snakemake@config[["finemapping_method"]]

# Per-locus summaries are many files, not one -- report as a pattern + count rather than
# enumerating every path.
per_locus_summary_pattern <- sprintf(
  "results/%s/<locus>/%s/<locus>_summary.txt",
  study_id,
  finemapping_method
)

manifest_lines <- c(
  "",
  "---- OUTPUT MANIFEST ----",
  sprintf("Exclusions table:               %s", snakemake@input[["exclusions"]]),
  sprintf(
    "Per-locus summaries:            %s  (%d loci)",
    per_locus_summary_pattern,
    length(summary_files)
  ),
  sprintf("SNPs table:                     %s", snakemake@input[["snps_table"]]),
  sprintf("Final report:                   %s", snakemake@input[["final_report"]]),
  sprintf("OpenTargets alignment summary:  %s", snakemake@input[["alignment_summary"]]),
  sprintf("Locus R2 summary:               %s", snakemake@input[["locus_r2_summary"]]),
  sprintf("Finalize sentinel:              %s", snakemake@output[[1]]),
  "--------------------------"
)

cat(paste0(manifest_lines, "\n"), file = log_path, append = TRUE, sep = "")

writeLines("DONE", snakemake@output[[1]])
