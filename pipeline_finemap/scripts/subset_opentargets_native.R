# subset_opentargets_native.R
#
# Filters the shared, batch-wide OpenTargets cache down to just this one
# study's rows. Cheap operation — the cache itself is the expensive part
# (built once by build_opentargets_cache.R), this just does a fast filter.

library(data.table)

study_id <- snakemake@wildcards[["study"]]
cache <- fread(snakemake@input[["cache"]])

subset_rows <- cache[studyId == study_id]

fwrite(subset_rows, snakemake@output[[1]], sep = "\t")

cat(sprintf(
  "%s: found %d OpenTargets rows for study %s\n",
  Sys.time(),
  nrow(subset_rows),
  study_id
))
