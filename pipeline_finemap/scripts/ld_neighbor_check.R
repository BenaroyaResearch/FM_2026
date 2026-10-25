# Determines if subset of GWAS summary stats surrounding lead SNP is missing any high-LD neighbors

library(data.table)

log_path <- snakemake@params[["log_path"]]
locus_id <- snakemake@wildcards[["locus"]]
output_path <- snakemake@output[[1]]

result <- tryCatch(
  {
    ld_file_first_line <- readLines(snakemake@input[["ld_file"]], n = 1)

    # Checks if the lead SNP isn't in reference panel (recorded in previous step, checked here)
    if (
      length(ld_file_first_line) > 0 &&
        ld_file_first_line == "LEAD_SNP_NOT_IN_REFERENCE"
    ) {
      cat(
        sprintf(
          "%s: locus %s — NOT_IN_REFERENCE: lead SNP absent from 1000G panel, cannot verify LD neighbors, excluding from fine-mapping\n",
          Sys.time(),
          locus_id
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(
        c(
          "FLAG: NOT_IN_REFERENCE",
          "missing_neighbors_count: NA",
          "note: lead SNP not present in reference panel, LD could not be assessed"
        ),
        output_path
      )
      quit(status = 0)
    }

    # Checks if previous rule found no SNPs in high-LD with lead SNP (therefore producing empty.ld file of size 0)
    # Passes LD neighbor check
    if (file.info(snakemake@input[["ld_file"]])$size == 0) {
      cat(
        sprintf(
          "%s: locus %s — PASS: no SNPs in high LD with lead SNP found (empty .ld file)\n",
          Sys.time(),
          locus_id
        ),
        file = log_path,
        append = TRUE
      )
      writeLines(c("FLAG: PASS", "missing_neighbors_count: 0"), output_path)
      quit(status = 0)
    }

    # Reading GWAS positions present at this locus, real PLINK LD output, and reference .bim for indel checking
    locus_snps <- fread(snakemake@input[["gwas_subset"]])
    ld_data <- fread(snakemake@input[["ld_file"]])
    bim_file <- fread(snakemake@input[["bim"]])

    missing_neighbors <- data.table(
      SNP_B = character(),
      BP_B = numeric(),
      R2 = numeric()
    )

    # Check every high-LD neighbor
    for (j in seq_len(nrow(ld_data))) {
      bp <- ld_data$BP_B[j]

      # Grab its position and check if the position is found anywhere in GWAS subset
      if (bp %in% locus_snps$BP) {
        next
      }

      # If neighbor is absent from GWAS subset, grab its ID and R2 for the record
      snp_b <- ld_data$SNP_B[j]
      r2_val <- ld_data$R2[j]
      bim_row <- bim_file[bim_file$V2 == snp_b, ]

      # Defensive check: if high LD neighbor isn't found in bim (extremely unlikely)
      if (nrow(bim_row) == 0) {
        missing_neighbors <- rbind(
          missing_neighbors,
          data.table(SNP_B = snp_b, BP_B = bp, R2 = r2_val)
        )
        next
      }

      # Located high LD neighbor in the .bim file and checks if it's an indel (allele count is greater than one).
      # Indels were deliberately kept in munging but can't be checked against the reference in the same way.
      # Absence is not counted as a missing high-LD neighbor
      allele1 <- bim_row$V5[1]
      allele2 <- bim_row$V6[1]
      is_indel <- nchar(allele1) > 1 || nchar(allele2) > 1

      if (is_indel) {
        next
      }

      missing_neighbors <- rbind(
        missing_neighbors,
        data.table(SNP_B = snp_b, BP_B = bp, R2 = r2_val)
      )
    }

    # No real missing neighbors --> Pass
    # Missing one neighbor with 0.6 < LD < 0.8 --> Mod low
    # Missing one neighbor with 0.8 <= LD --> Mod high
    # Missing multiple high LD neighbors --> High
    n_missing <- nrow(missing_neighbors)

    if (n_missing == 0) {
      flag <- "PASS"
    } else if (n_missing == 1) {
      if (missing_neighbors$R2[1] >= 0.8) {
        flag <- "MOD_HIGH"
      } else {
        flag <- "MOD_LOW"
      }
    } else {
      flag <- "HIGH"
    }

    cat(
      sprintf(
        "%s: locus %s — %s (%d missing neighbor(s))\n",
        Sys.time(),
        locus_id,
        flag,
        nrow(missing_neighbors)
      ),
      file = log_path,
      append = TRUE
    )

    # Always writes the flag and count of missing neighbors.
    # Appends list of missing neighbors.
    out_lines <- c(
      sprintf("FLAG: %s", flag),
      sprintf("missing_neighbors_count: %d", nrow(missing_neighbors))
    )

    if (nrow(missing_neighbors) > 0) {
      out_lines <- c(
        out_lines,
        "SNP_B\tBP_B\tR2",
        apply(missing_neighbors, 1, function(r) paste(r, collapse = "\t"))
      )
    }

    writeLines(out_lines, output_path)
  },
  # If any part of the LD neighbor check itself goes awry, we have no verified information
  # about whether this locus is missing high-LD neighbors -- that's worse than a known HIGH
  # flag (which at least reflects a completed check), so it's excluded rather than waved through.
  error = function(e) {
    cat(
      sprintf(
        "%s: locus %s — ERROR during LD-neighbor check: %s. Excluding locus (missing-neighbor status could not be verified).\n",
        Sys.time(),
        locus_id,
        conditionMessage(e)
      ),
      file = log_path,
      append = TRUE
    )
    writeLines(
      c(
        "FLAG: CHECK_FAILED",
        "missing_neighbors_count: NA",
        sprintf("error: %s", conditionMessage(e))
      ),
      output_path
    )
  }
)
