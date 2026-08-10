# Script for standardizing GRCh38 format and munging via MungeSumstats
# Also applies three-tier fine-mapping decision logic based on final SNP count and sample size

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1] # Sumstats
output_path <- args[2] # Munged file
decision_path <- args[3] # FM decision
log_path <- args[4] # Log path
nthread <- as.integer(args[5]) # Threads
n_samples <- as.numeric(args[6]) # Sample size
n_high <- as.numeric(args[7]) # "High" sample size value
snp_high <- as.numeric(args[8]) # "High" SNP value
n_mid_low <- as.numeric(args[9]) # COME BACK
snp_mid_low <- as.numeric(args[10])
filename_string <- args[11]

library(MungeSumstats)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

# Append to log file
cat(
  sprintf("%s: starting munge (nThread=%d)\n", Sys.time(), nthread),
  file = log_path,
  append = TRUE
)

# Check the resolved filename for build hints before attempting inference.
# ".h." and "38" both indicate GRCh38; "37" indicates GRCh37. If the
# filename contains contradictory signals (e.g. both "37" and ".h."),
# do not guess — fall through to inference instead.
cat(
  sprintf(
    "%s: checking filename for build hints: %s\n",
    Sys.time(),
    filename_string
  ),
  file = log_path,
  append = TRUE
)

# Checking for hints on genome build in sum stats file name to avoid adding unnecessary inference step
has_38_token <- grepl(
  "\\.h\\.|GRCh38|hg38|[_.-]b38",
  filename_string,
  ignore.case = TRUE
)
has_37_token <- grepl(
  "GRCh37|hg19|[_.-]b37",
  filename_string,
  ignore.case = TRUE
)

says_38 <- has_38_token
says_37 <- has_37_token

if (says_38 && says_37) {
  cat(
    sprintf(
      "%s: filename contains conflicting build signals (37 and 38/.h.), cannot trust filename, proceeding to infer genome build\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
  inferred_build <- NULL # falls through to inference block below
} else if (says_38) {
  inferred_build <- "GRCh38"
  cat(
    sprintf(
      "%s: filename indicates GRCh38, skipping inference\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
} else if (says_37) {
  inferred_build <- "GRCh37"
  cat(
    sprintf(
      "%s: filename indicates GRCh37, skipping inference\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
} else {
  cat(
    sprintf(
      "%s: no build hint in filename, proceeding to infer genome build\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
  inferred_build <- NULL
}

# If filename gave no reliable signal, falls back to get_genome_build() MungeSumstats tool
if (is.null(inferred_build)) {
  inferred_build <- tryCatch(
    {
      result <- get_genome_builds(
        sumstats_list = list(study = input_path),
        nThread = nthread
      )
      result[["study"]]
    },
    error = function(e) {
      cat(
        sprintf(
          "%s: ERROR - build inference failed: %s\n",
          Sys.time(),
          conditionMessage(e)
        ),
        file = log_path,
        append = TRUE
      )
      quit(status = 4) # exit code 4: "couldn't determine build at all"
    }
  )

  # Checks if build inference succeeded but returned empty
  if (
    is.null(inferred_build) ||
      length(inferred_build) == 0 ||
      is.na(inferred_build)
  ) {
    cat(
      sprintf(
        "%s: ERROR - build inference returned no result, cannot proceed\n",
        Sys.time()
      ),
      file = log_path,
      append = TRUE
    )
    quit(status = 4) # exit code 4: "couldn't determine build at all"
  }
}

cat(
  sprintf("%s: genome build for this run: %s\n", Sys.time(), inferred_build),
  file = log_path,
  append = TRUE
)

# Attempt to munge with following flags. Uses inferred genome build and lifts to GRCh38 if needed.

munge_log_dir <- dirname(output_path)

reformatted <- tryCatch(
  {
    format_sumstats(
      path = input_path,
      ref_genome = inferred_build,
      convert_ref_genome = "GRCh38",
      indels = TRUE, # Indels exist in data
      drop_indels = FALSE, # Do not drop the indels ** Potential place for filtering change
      bi_allelic_filter = FALSE, # Keep multi-allelic SNPs ** Potential place for filtering change
      flip_frq_as_biallelic = TRUE, # Treat multi-allelic SNPs as biallelic and flip frq if needed ** Potential place for filtering change
      strand_ambig_filter = TRUE, # Removes AT and CG SNPs
      nThread = nthread, # Number of CPU threads available to munge
      return_data = TRUE, # Return as R object
      log_folder_ind = TRUE, # Write detailed logs
      log_folder = munge_log_dir,
      log_mungesumstats_msgs = TRUE # Capture MungeSumstat's messages
    )
  },
  # Captures any errors thrown during munging
  error = function(e) {
    cat(
      sprintf(
        "%s: ERROR - munge failed: %s\n",
        Sys.time(),
        conditionMessage(e)
      ),
      file = log_path,
      append = TRUE
    )
    quit(status = 2) # error 2: munge error
  }
)

# Save the raw return value immediately, no matter what it turns out to be structurally.
# This protects the actual munge output from any bug in the code that comes after
# (especially any failure that may occur while saving sum stats as tsv)
saveRDS(reformatted, sub("\\.tsv\\.gz$", "_raw.rds", output_path))
cat(
  sprintf(
    "%s: raw munge result saved (class: %s)\n",
    Sys.time(),
    paste(class(reformatted), collapse = ",")
  ),
  file = log_path,
  append = TRUE
)

# Count number of SNPs after munging
munged <- reformatted$sumstats
n_snps <- nrow(munged)

# In case SNP count fails (munge result isn't a proper table)
if (is.null(n_snps) || length(n_snps) == 0) {
  cat(
    sprintf(
      "%s: ERROR - could not determine row count (class of $sumstats: %s)\n",
      Sys.time(),
      paste(class(munged), collapse = ",")
    ),
    file = log_path,
    append = TRUE
  )
  quit(status = 3) # error 3: munge successful but SNP count failed
}

# Check to make sure output isn't empty
if (n_snps == 0) {
  cat(
    sprintf(
      "%s: ERROR - munge completed but produced 0 variants\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
  quit(status = 3) #  # error 3: munge successful, but returned no SNPs
}

# Munge result
fwrite(munged, output_path, sep = "\t")
cat(
  sprintf("%s: munge completed, %d SNPs retained\n", Sys.time(), n_snps),
  file = log_path,
  append = TRUE
)

# Munge log appended to preexisting study log
msg_log_path <- reformatted$log_files$MungeSumstats_log_msg
if (!is.null(msg_log_path) && file.exists(msg_log_path)) {
  msg_lines <- readLines(msg_log_path)
  cat(
    sprintf("%s: ---- MungeSumstats log messages ----\n", Sys.time()),
    file = log_path,
    append = TRUE
  )
  con <- file(log_path, open = "a")
  writeLines(msg_lines, con)
  close(con)
  cat(
    sprintf("%s: ---- end MungeSumstats log messages ----\n", Sys.time()),
    file = log_path,
    append = TRUE
  )
}

# Decision logic
# Is the SNP count below absolute minimum (designated in config (currently 1,000,000))? If yes, do not FM.
if (n_snps < snp_mid_low) {
  decision <- "none"
  reason <- sprintf(
    "SNP count (%d) below minimum threshold (%d) — not enough coverage to fine-map",
    n_snps,
    snp_mid_low
  )
  # If high SNP count and high sample size, eligible for functional and normal FM. Thresholds also designated in config.
  # ** Potential place for filtering change **
} else if (n_samples > n_high && n_snps > snp_high) {
  decision <- "both"
  reason <- sprintf(
    "N (%.0f) > %d and SNPs (%d) > %d — sufficient for both regular and functional FM",
    n_samples,
    n_high,
    n_snps,
    snp_high
  )
  # If doesn't satisfy either of the above conditions, proceed with normal FM only.
} else {
  decision <- "regular_only"
  reason <- sprintf(
    "N (%.0f), SNPs (%d) — adequate coverage but does not meet both-path thresholds; defaulting to regular FM only",
    n_samples,
    n_snps
  )
}

# Record FM decision. Will be read by pipeline_finemap to determine next steps.
cat(
  sprintf("%s: DECISION=%s. %s\n", Sys.time(), decision, reason),
  file = log_path,
  append = TRUE
)

writeLines(decision, decision_path)
