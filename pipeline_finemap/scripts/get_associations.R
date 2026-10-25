# Resolves each lead SNP's rsID to a GRCh38 CHR:BP by looking it up directly in dbSNP, rather
# than trusting study_table.tsv's lead_snp_positions column (which carries no build label).
# Checks the given position against both GRCh37 and GRCh38 dbSNP records for that rsID and
# converts to GRCh38 if it matched GRCh37. A lead SNP whose position matches neither build is
# excluded rather than silently mismatched downstream.

library(data.table)

log_path <- snakemake@params[["log_path"]]
decision <- trimws(readLines(snakemake@input[["decision"]]))
ref_lib <- snakemake@params[["ref_lib"]]
chain_file <- snakemake@params[["chain_file"]]

# Initializes associations table.
empty_table <- function() {
  data.table(SNP = character(), CHR = character(), BP = integer())
}

# Hypothetically, a study that failed pipeline_munge should not proceed to FM and
# therefore would be excluded from the config of pipeline_finemap. However, if this is
# not the case, pipeline_finemap will gracefully handle a study with fine-map decision = none.

# If decision to FM was none, produce an empty associations table and exit
if (decision == "none") {
  cat(
    sprintf(
      "%s: fm_decision=none, refusing to proceed. Writing empty associations table.\n",
      Sys.time()
    ),
    file = log_path,
    append = TRUE
  )
  write.table(
    empty_table(),
    snakemake@output[[1]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  quit(status = 0)
}

# If FM decision was not none, proceed.
cat(
  sprintf(
    "%s: fm_decision=%s, proceeding with associations table\n",
    Sys.time(),
    decision
  ),
  file = log_path,
  append = TRUE
)

# study_table.tsv stores each study's lead SNPs and positions as comma-separated strings
# in a single cell (column headers = "lead_snps", "lead_positions").
# splits IDs and positions on the commas
lead_snp_ids <- trimws(unlist(strsplit(snakemake@params[["lead_snps"]], ",")))
lead_positions <- trimws(unlist(strsplit(
  snakemake@params[["lead_positions"]],
  ","
)))
# Further splits into chr and bp on the colon
# sapply(strsplit(...), [, 1) common R idiom meaning "split every string, then take the 1st element
# of the resulting list"
lead_chr <- sapply(strsplit(lead_positions, ":"), `[`, 1)
lead_bp <- as.integer(sapply(strsplit(lead_positions, ":"), `[`, 2))

# Some loci have no known rsID and use a "chrN:pos" placeholder in lead_snps instead (see
# gather_loci/get_loci). snpsById() errors on anything that isn't a real rsID, so those can't be
# build-checked by ID -- but their position CAN be build-checked directly (see the
# !is_rsid[i] branch below, which queries dbSNP by position instead of by rsID).
is_rsid <- grepl("^rs[0-9]+$", lead_snp_ids) # Are lead SNPs listed by rsid? T/F
has_chrbp <- !is_rsid & !is.na(lead_bp) #

# SNPlocs.Hsapiens.dbSNP155.GRCh37/GRCh38 are too large to ship in the container image and live
# in the shared renv library instead (the same packages pipeline_munge uses to build-check the
# full sumstats file); loading them pulls in a good deal of Bioconductor machinery, so this is
# skipped entirely for a study with nothing to build-check either way (no rsIDs and no chr:BP
# entries with a usable position).

# Are there any lead positions to build check on?
if (any(is_rsid) || any(has_chrbp)) {
  required_pkgs <- c(
    "SNPlocs.Hsapiens.dbSNP155.GRCh37",
    "SNPlocs.Hsapiens.dbSNP155.GRCh38"
  )
  missing_pkgs <- required_pkgs[!dir.exists(file.path(ref_lib, required_pkgs))]
  if (length(missing_pkgs) > 0) {
    msg <- sprintf(
      "reference package(s) missing from %s: %s",
      ref_lib,
      paste(missing_pkgs, collapse = ", ")
    )
    cat(
      sprintf("%s: ERROR - %s\n", Sys.time(), msg),
      file = log_path,
      append = TRUE
    )
    stop(msg)
  }
  .libPaths(c(ref_lib, .libPaths()))
  suppressMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh37))
  suppressMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh38))

  # Looks up every lead SNP's canonical position in both builds in one batched call each (far
  # fewer queries than one snpsById() per rsID), then matches results back onto individual SNPs
  # by the returned RefSNP_id. ifnotfound="drop" can return fewer rows than requested, so
  # relying on row order/index to realign would silently mismatch SNPs to the wrong position.
  lookup_build <- function(snplocs, ids) {
    hits <- snpsById(snplocs, unique(ids), ifnotfound = "drop")
    data.table(
      SNP = mcols(hits)$RefSNP_id,
      CHR = as.character(seqnames(hits)),
      BP = pos(hits)
    )
  }

  if (any(is_rsid)) {
    hits38 <- lookup_build(
      SNPlocs.Hsapiens.dbSNP155.GRCh38,
      lead_snp_ids[is_rsid]
    )
    hits37 <- lookup_build(
      SNPlocs.Hsapiens.dbSNP155.GRCh37,
      lead_snp_ids[is_rsid]
    )
  }

  # chr:BP entries have no rsID to key a batched lookup on (unlike the rsID path above), so their
  # build is checked one position at a time below via snpsByOverlaps(); the liftover chain is
  # only needed if one of them actually turns out to be GRCh37.
  if (any(has_chrbp)) {
    chain <- rtracklayer::import.chain(chain_file)
  }
} else {
  cat(
    sprintf(
      "%s: no rsIDs or chr:BP position(s) with usable coordinates among %d lead SNP(s) for this study, skipping SNPlocs build-check entirely\n",
      Sys.time(),
      length(lead_snp_ids)
    ),
    file = log_path,
    append = TRUE
  )
}

results <- list()
n_converted <- 0
n_excluded <- 0

for (i in seq_along(lead_snp_ids)) {
  snp_id <- lead_snp_ids[i]

  if (!is_rsid[i]) {
    pos_valid <- !is.na(lead_bp[i])
    if (!pos_valid) {
      cat(
        sprintf(
          "%s: WARNING: lead SNP %s has no rsID and no valid position, excluding from associations\n",
          Sys.time(),
          snp_id
        ),
        file = log_path,
        append = TRUE
      )
      n_excluded <- n_excluded + 1
      next
    }

    # No rsID means the batched-by-ID lookup above can't be used, but the position itself can
    # still be build-checked directly: does a known dbSNP variant exist at this exact CHR:BP in
    # GRCh38, in GRCh37, in both, or in neither? A hit in exactly one build is the same kind of
    # confident evidence the rsID path above requires (an exact position match), so the same
    # exclude-when-ambiguous philosophy applies here too, rather than the old behavior of
    # trusting study_table.tsv's position as GRCh38 with no verification at all.
    pos_gr <- GenomicRanges::GRanges(
      lead_chr[i],
      IRanges::IRanges(lead_bp[i], width = 1)
    )
    hit38_pos <- tryCatch(
      snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, pos_gr),
      error = function(e) NULL
    )
    hit37_pos <- tryCatch(
      snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh37, pos_gr),
      error = function(e) NULL
    )
    found38 <- !is.null(hit38_pos) && length(hit38_pos) > 0
    found37 <- !is.null(hit37_pos) && length(hit37_pos) > 0

    if (found38 && !found37) {
      cat(
        sprintf(
          "%s: lead SNP %s (%s:%d) has no rsID; confirmed GRCh38 by position (matches known dbSNP variant %s)\n",
          Sys.time(),
          snp_id,
          lead_chr[i],
          lead_bp[i],
          mcols(hit38_pos)$RefSNP_id[1]
        ),
        file = log_path,
        append = TRUE
      )
      results[[length(results) + 1]] <- data.table(
        SNP = snp_id,
        CHR = lead_chr[i],
        BP = lead_bp[i]
      )
    } else if (found37 && !found38) {
      lifted <- tryCatch(
        unlist(rtracklayer::liftOver(
          GenomicRanges::GRanges(
            paste0("chr", lead_chr[i]),
            IRanges::IRanges(lead_bp[i], width = 1)
          ),
          chain
        )),
        error = function(e) NULL
      )
      if (!is.null(lifted) && length(lifted) == 1) {
        new_chr <- sub("^chr", "", as.character(seqnames(lifted)))
        new_bp <- start(lifted)
        cat(
          sprintf(
            "%s: converted lead SNP %s from GRCh37 (%s:%d) to GRCh38 (%s:%d) via liftOver -- no rsID, but position matched a known GRCh37 dbSNP variant (%s)\n",
            Sys.time(),
            snp_id,
            lead_chr[i],
            lead_bp[i],
            new_chr,
            new_bp,
            mcols(hit37_pos)$RefSNP_id[1]
          ),
          file = log_path,
          append = TRUE
        )
        results[[length(results) + 1]] <- data.table(
          SNP = snp_id,
          CHR = new_chr,
          BP = new_bp
        )
        n_converted <- n_converted + 1
      } else {
        cat(
          sprintf(
            "%s: WARNING: lead SNP %s (%s:%d) matched a known GRCh37 dbSNP variant but liftOver to GRCh38 failed or was ambiguous (%s), excluding from associations\n",
            Sys.time(),
            snp_id,
            lead_chr[i],
            lead_bp[i],
            if (is.null(lifted)) {
              "liftOver error"
            } else {
              sprintf("%d results, expected 1", length(lifted))
            }
          ),
          file = log_path,
          append = TRUE
        )
        n_excluded <- n_excluded + 1
      }
    } else {
      reason <- if (found38 && found37) {
        "matches a known variant in both GRCh37 and GRCh38 at this exact position"
      } else {
        "matches no known dbSNP variant in either build at this exact position"
      }
      cat(
        sprintf(
          "%s: WARNING: lead SNP %s (%s:%d) has no rsID and its genome build could not be confidently resolved (%s), excluding from associations\n",
          Sys.time(),
          snp_id,
          lead_chr[i],
          lead_bp[i],
          reason
        ),
        file = log_path,
        append = TRUE
      )
      n_excluded <- n_excluded + 1
    }
    next
  }

  hit38 <- hits38[SNP == snp_id]
  hit37 <- hits37[SNP == snp_id]
  pos_valid <- !is.na(lead_bp[i])

  is_grch38 <- pos_valid &&
    nrow(hit38) > 0 &&
    hit38$CHR[1] == lead_chr[i] &&
    hit38$BP[1] == lead_bp[i]
  is_grch37 <- !is_grch38 &&
    pos_valid &&
    nrow(hit37) > 0 &&
    nrow(hit38) > 0 &&
    hit37$CHR[1] == lead_chr[i] &&
    hit37$BP[1] == lead_bp[i]

  if (is_grch38) {
    results[[length(results) + 1]] <- data.table(
      SNP = snp_id,
      CHR = hit38$CHR[1],
      BP = hit38$BP[1]
    )
  } else if (is_grch37) {
    cat(
      sprintf(
        "%s: converted lead SNP %s from GRCh37 (%s:%d) to GRCh38 (%s:%d)\n",
        Sys.time(),
        snp_id,
        lead_chr[i],
        lead_bp[i],
        hit38$CHR[1],
        hit38$BP[1]
      ),
      file = log_path,
      append = TRUE
    )
    results[[length(results) + 1]] <- data.table(
      SNP = snp_id,
      CHR = hit38$CHR[1],
      BP = hit38$BP[1]
    )
    n_converted <- n_converted + 1
  } else {
    cat(
      sprintf(
        "%s: WARNING: lead SNP %s (%s:%s) does not match either GRCh37 or GRCh38 dbSNP position for this rsID, excluding from associations\n",
        Sys.time(),
        snp_id,
        lead_chr[i],
        ifelse(pos_valid, lead_bp[i], "NA")
      ),
      file = log_path,
      append = TRUE
    )
    n_excluded <- n_excluded + 1
  }
}

matched <- if (length(results) > 0) rbindlist(results) else empty_table()

cat(
  sprintf(
    "%s: %d of %d lead SNP(s) resolved to GRCh38 (%d converted from GRCh37, %d excluded - build unresolved)\n",
    Sys.time(),
    nrow(matched),
    length(lead_snp_ids),
    n_converted,
    n_excluded
  ),
  file = log_path,
  append = TRUE
)

write.table(
  matched,
  snakemake@output[[1]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
