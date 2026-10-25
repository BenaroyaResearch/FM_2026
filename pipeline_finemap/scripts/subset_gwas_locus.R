# Subset the sum stats to contain only SNPs within {bp window} from lead association

library(data.table)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]
window_bp <- as.integer(snakemake@params[["window_bp"]]) # Established in config file

associations <- fread(snakemake@input[["associations"]])
lead_row <- associations[SNP == locus_id]

# Confirms that locus actually exists in associations table
# (In case associations table got regenerated or someone created/deleted .locus marker file
# so there is a mismatch between true associations table and Snakemake associations list)
if (nrow(lead_row) == 0) {
  cat(
    sprintf(
      "%s: ERROR - locus %s not found in associations table\n",
      Sys.time(),
      locus_id
    ),
    file = log_path,
    append = TRUE
  )
  quit(status = 1)
}

lead_chr <- lead_row$CHR[1]
lead_bp <- lead_row$BP[1]
window_lo <- lead_bp - window_bp
window_hi <- lead_bp + window_bp

# Filters via zcat+awk before the data ever reaches R -- munged files can run to tens of
# millions of rows (multi-GB once loaded into a data.table), and this rule runs once per locus
# with no reuse across a study's loci, which is what caused this workspace's OOM crashes. awk
# streams the gzip'd file line-by-line (O(1) memory) and hands R only the already-windowed
# subset. CHR/BP are located by header name rather than a fixed column position, since studies'
# munged files vary in which other columns they carry (BETA vs OR vs Z, etc.).
munged_path <- snakemake@input[["munged"]]
filter_cmd <- sprintf(
  "zcat %s | awk -F'\\t' -v chr=%s -v lo=%d -v hi=%d 'NR==1{for(i=1;i<=NF;i++){if($i==\"CHR\")c=i;if($i==\"BP\")b=i} print; next} $c==chr && $b>=lo && $b<=hi'",
  shQuote(munged_path),
  shQuote(lead_chr),
  window_lo,
  window_hi
)
subset_dt <- fread(cmd = filter_cmd)

cat(
  sprintf(
    "%s: locus %s (%s:%d) — subset window +/-%dbp, %d SNPs retained\n",
    Sys.time(),
    locus_id,
    lead_chr,
    lead_bp,
    window_bp,
    nrow(subset_dt)
  ),
  file = log_path,
  append = TRUE
)

# Checks if no SNPs were found within the window
if (nrow(subset_dt) == 0) {
  cat(
    sprintf("%s: WARNING - locus %s subset is empty\n", Sys.time(), locus_id),
    file = log_path,
    append = TRUE
  )
}

fwrite(subset_dt, snakemake@output[[1]], sep = "\t")
