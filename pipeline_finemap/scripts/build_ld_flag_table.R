# Purpose: Pulls LD flag value for each locus and puts into table {study}_ld_flags.tsv. Used to filter eligible FM loci.

library(data.table)

# Pulls in LD flags and full associations list
log_path <- snakemake@params[["log_path"]]
associations <- fread(snakemake@input[["associations"]])
ld_check_files <- snakemake@input[["ld_checks"]]

flags <- data.table(SNP = character(), LD_FLAG = character())

# For each locus file, read LD flag and add row to flags
for (f in ld_check_files) {
  lines <- readLines(f)
  flag_line <- lines[grepl("^FLAG:", lines)]
  flag <- trimws(sub("^FLAG:\\s*", "", flag_line[1]))

  locus_id <- basename(dirname(dirname(dirname(f))))
  flags <- rbind(flags, data.table(SNP = locus_id, LD_FLAG = flag))
}

# Merge rows into table
flag_table <- merge(associations, flags, by = "SNP", all.x = TRUE)

# Counts how many loci fell into each category
cat(
  sprintf(
    "%s: LD flag summary — PASS: %d, MOD_LOW: %d, MOD_HIGH: %d, HIGH: %d, NOT_IN_REFERENCE: %d, CHECK_FAILED: %d\n",
    Sys.time(),
    sum(flag_table$LD_FLAG == "PASS", na.rm = TRUE),
    sum(flag_table$LD_FLAG == "MOD_LOW", na.rm = TRUE),
    sum(flag_table$LD_FLAG == "MOD_HIGH", na.rm = TRUE),
    sum(flag_table$LD_FLAG == "HIGH", na.rm = TRUE),
    sum(flag_table$LD_FLAG == "NOT_IN_REFERENCE", na.rm = TRUE),
    sum(flag_table$LD_FLAG == "CHECK_FAILED", na.rm = TRUE)
  ),
  file = log_path,
  append = TRUE
)

write.table(
  flag_table,
  snakemake@output[[1]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
