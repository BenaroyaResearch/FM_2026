library(readr)
library(data.table)

raw_table <- snakemake@params[["raw_table"]]
study_table_path <- snakemake@input[["study_table"]]
output_path <- snakemake@output[[1]]

study_table <- fread(study_table_path)
my_gcst_ids <- unique(study_table$gcst_id)

cat(sprintf(
  "%s: filtering for %d specific study IDs\n",
  Sys.time(),
  length(my_gcst_ids)
))

keep_cols <- c(
  "studyId",
  "variantId",
  "chromosome",
  "position",
  "locusStart",
  "locusEnd",
  "credibleSetIndex",
  "finemappingMethod",
  "purityMeanR2",
  "purityMinR2",
  "qualityControls",
  "confidence",
  "ldSet",
  "locus"
)

# Columns known to contain multi-line embedded strings (Python-style list-of-dicts).
# Newlines/carriage returns here are purely cosmetic and safe to collapse to a
# single space -- doing so guarantees every row is exactly one physical line,
# so ANY downstream tool (fread, read_tsv, awk, wc -l) parses the file identically.
sanitize_newlines <- function(dt) {
  char_cols <- names(dt)[sapply(dt, is.character)]
  for (col in char_cols) {
    x <- dt[[col]]
    x <- gsub("[\r\n]+", " ", x)
    # Collapse a closing brace, whitespace, opening brace (with no comma
    # between them) into a properly comma-separated pair. This repairs
    # list-of-dict strings where the original line break WAS the only
    # separator between entries, which the newline->space substitution
    # above turns into invalid JSON (missing comma) rather than valid
    # JSON (comma already present, newline just cosmetic).
    x <- gsub("\\}\\s+\\{", "}, {", x)
    dt[[col]] <- x
  }
  dt
}

matched_rows <- list()

filter_chunk <- function(chunk, pos) {
  setDT(chunk)
  is_match <- chunk$studyId %in% my_gcst_ids
  if (any(is_match)) {
    sub <- chunk[is_match, ..keep_cols]
    sub <- sanitize_newlines(sub)
    matched_rows[[length(matched_rows) + 1]] <<- sub
  }
}

cat(sprintf("%s: starting scan of %s\n", Sys.time(), raw_table))

read_tsv_chunked(
  raw_table,
  callback = SideEffectChunkCallback$new(filter_chunk),
  chunk_size = 50000,
  col_types = cols(.default = "c")
)

cache <- if (length(matched_rows) > 0) {
  rbindlist(matched_rows, fill = TRUE)
} else {
  data.table()
}

fwrite(cache, output_path, sep = "\t")

n_found <- if (nrow(cache) > 0) uniqueN(cache$studyId) else 0
cat(sprintf(
  "%s: done. Cached %d rows across %d of your %d studies.\n",
  Sys.time(),
  nrow(cache),
  n_found,
  length(my_gcst_ids)
))
