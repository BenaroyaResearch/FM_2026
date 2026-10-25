library(data.table)
library(jsonlite)
library(ggplot2)

log_path <- snakemake@params[["log_path"]]
cs_size_threshold <- as.integer(snakemake@params[["cs_size_exclude_threshold"]])

associations <- fread(snakemake@input[["associations"]])
ot_native <- fread(snakemake@input[["ot_native"]])
susie_files <- snakemake@input[["susie_results"]]

# ============================================================
# SHARED PARSING (unchanged from before, confirmed working)
# ============================================================
parse_credible_set <- function(x) {
  if (is.na(x) || x == "" || x == "[]") {
    return(NULL)
  }
  cleaned <- gsub("[\r\n]+", " ", x)
  cleaned <- gsub("\\}\\s+\\{", "}, {", cleaned)
  cleaned <- gsub("'", "\"", cleaned)
  cleaned <- gsub("\\bTrue\\b", "true", cleaned)
  cleaned <- gsub("\\bFalse\\b", "false", cleaned)
  cleaned <- gsub("\\bNone\\b", "null", cleaned)
  tryCatch(fromJSON(cleaned, simplifyDataFrame = TRUE), error = function(e) {
    NULL
  })
}
get_members <- function(row) {
  parsed <- parse_credible_set(row$ldSet)
  if (
    !is.null(parsed) &&
      is.data.frame(parsed) &&
      "posteriorProbability" %in% names(parsed)
  ) {
    return(parsed)
  }
  parsed <- parse_credible_set(row$locus)
  if (
    !is.null(parsed) &&
      is.data.frame(parsed) &&
      "posteriorProbability" %in% names(parsed)
  ) {
    return(parsed)
  }
  NULL
}

load_my_locus <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(NULL)
  }
  tab <- fread(f)
  if (nrow(tab) == 0) {
    return(NULL)
  }
  tab[, locus := basename(dirname(dirname(dirname(f))))]
  tab
}
# susie_files is NULL rather than character(0) when a study has zero eligible loci (every
# lead SNP excluded before scoring) -- rbindlist on an all-NULL/empty list yields a data.table
# with no columns at all, so referencing CS below would error. Guarded rather than assumed
# non-empty since "every locus excluded" is a real, already-supported outcome elsewhere in
# this pipeline (see the exclusions table).
my_all <- if (!is.null(susie_files) && length(susie_files) > 0) {
  rbindlist(lapply(susie_files, load_my_locus), fill = TRUE)
} else {
  data.table()
}
if (nrow(my_all) > 0) {
  my_all[, cs_key := fifelse(is.na(CS), NA_character_, paste(locus, CS, sep = "_"))]
}

all_ot_members <- list()
n_ot_skipped <- 0
if (nrow(ot_native) > 0) {
  for (i in seq_len(nrow(ot_native))) {
    result <- tryCatch(
      {
        members <- get_members(ot_native[i])
        if (
          is.null(members) ||
            !is.data.frame(members) ||
            nrow(members) == 0 ||
            !("variantId" %in% names(members))
        ) {
          NULL
        } else {
          members_dt <- as.data.table(members)
          members_dt[,
            c("m_chr", "m_bp", "m_ref", "m_alt") := tstrsplit(
              variantId,
              "_",
              fixed = TRUE
            )
          ]
          members_dt[, m_bp := as.integer(m_bp)]
          members_dt[,
            ot_cs_id := paste0(
              ot_native$studyId[i],
              "_",
              ot_native$credibleSetIndex[i]
            )
          ]
          members_dt[, ot_credible_set_index := ot_native$credibleSetIndex[i]]
          members_dt[, ot_finemapping_method := ot_native$finemappingMethod[i]]
          members_dt[, ot_quality_controls := ot_native$qualityControls[i]]
          members_dt[, ot_confidence := ot_native$confidence[i]]
          members_dt[,
            ot_purity_mean_r2 := as.numeric(ot_native$purityMeanR2[i])
          ]
          members_dt[, ot_purity_min_r2 := as.numeric(ot_native$purityMinR2[i])]
          members_dt
        }
      },
      error = function(e) NULL
    )
    if (!is.null(result)) {
      all_ot_members[[length(all_ot_members) + 1]] <- result
    } else {
      n_ot_skipped <- n_ot_skipped + 1
    }
  }
}
ot_flat <- if (length(all_ot_members) > 0) {
  rbindlist(all_ot_members, fill = TRUE)
} else {
  data.table()
}
if (n_ot_skipped > 0) {
  cat(
    sprintf(
      "%s: skipped %d of %d OpenTargets rows (unparseable credible-set structure)\n",
      Sys.time(),
      n_ot_skipped,
      nrow(ot_native)
    ),
    file = log_path,
    append = TRUE
  )
}

# ============================================================
# PART 1 — PIP ALIGNMENT. No size filter anywhere in this section.
# Full outer join: every SNP present in EITHER side.
# ============================================================

my_pips <- if (nrow(my_all) > 0) {
  my_all[, .(
    SNP,
    CHR,
    BP,
    locus,
    my_pip = PIP,
    my_cs = CS,
    my_cs_size = CS_SIZE
  )]
} else {
  # Typed rather than a bare data.table() -- the merge below joins on BP/m_bp, and this side
  # being empty must not mean that column is missing when the other side has real rows.
  data.table(
    SNP = character(),
    CHR = character(),
    BP = integer(),
    locus = character(),
    my_pip = numeric(),
    my_cs = integer(),
    my_cs_size = integer()
  )
}
ot_pips <- if (nrow(ot_flat) > 0) {
  unique(ot_flat[, .(
    m_bp,
    ot_pip = posteriorProbability,
    ot_cs_id,
    ot_finemapping_method
  )])
} else {
  data.table(
    m_bp = integer(),
    ot_pip = numeric(),
    ot_cs_id = character(),
    ot_finemapping_method = character()
  )
}

full_pip <- merge(my_pips, ot_pips, by.x = "BP", by.y = "m_bp", all = TRUE)

full_pip[,
  match_status := fifelse(
    !is.na(my_pip) & !is.na(ot_pip),
    "both",
    fifelse(!is.na(my_pip) & is.na(ot_pip), "mine_only", "OT_only")
  )
]

cat(
  sprintf(
    "%s: PIP alignment — %d SNPs total (%d both, %d mine only, %d OT only)\n",
    Sys.time(),
    nrow(full_pip),
    sum(full_pip$match_status == "both"),
    sum(full_pip$match_status == "mine_only"),
    sum(full_pip$match_status == "OT_only")
  ),
  file = log_path,
  append = TRUE
)

fwrite(full_pip, snakemake@output[["pip_full"]], sep = "\t")

# --- Plot A: SNPs present in both only ---
both_only <- full_pip[match_status == "both"]
if (nrow(both_only) >= 3) {
  r2 <- cor(both_only$my_pip, both_only$ot_pip)^2
  p1 <- ggplot(both_only, aes(x = ot_pip, y = my_pip)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(alpha = 0.5, color = "steelblue") +
    annotate(
      "text",
      x = 0.05,
      y = 0.95,
      hjust = 0,
      size = 4,
      label = sprintf("R^2 = %.3f\nn = %d", r2, nrow(both_only))
    ) +
    labs(
      title = "PIP alignment — SNPs present in both",
      x = "OpenTargets PIP",
      y = "My PIP"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_minimal()
} else {
  p1 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Too few matched SNPs") +
    theme_void()
}
ggsave(snakemake@output[["pip_matched_plot"]], p1, width = 7, height = 6)

# --- Plot B: everything, three categories, distinct shapes ---
full_pip[,
  plot_pip := fifelse(
    match_status == "mine_only",
    my_pip,
    fifelse(match_status == "OT_only", ot_pip, my_pip)
  )
]
full_pip[, plot_pip_y := fifelse(match_status == "OT_only", ot_pip, my_pip)]
full_pip[, plot_pip_x := fifelse(match_status == "mine_only", my_pip, ot_pip)]
# for "both": x=ot_pip, y=my_pip (standard). for mine_only: plotted on y-axis at x=NA-safe placeholder handled via shape/legend instead of forcing an axis position
plot_data <- full_pip[!is.na(match_status)]
plot_data[, x_val := fifelse(match_status == "mine_only", my_pip, ot_pip)]
plot_data[, y_val := fifelse(match_status == "OT_only", ot_pip, my_pip)]

p2 <- ggplot(
  plot_data,
  aes(x = x_val, y = y_val, shape = match_status, color = match_status)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5) +
  scale_shape_manual(
    values = c(both = 16, mine_only = 17, OT_only = 15),
    labels = c(
      both = "In both",
      mine_only = "Mine only",
      OT_only = "OpenTargets only"
    )
  ) +
  scale_color_manual(
    values = c(
      both = "steelblue",
      mine_only = "forestgreen",
      OT_only = "firebrick"
    ),
    labels = c(
      both = "In both",
      mine_only = "Mine only",
      OT_only = "OpenTargets only"
    )
  ) +
  labs(
    title = "PIP alignment — all SNPs from either pipeline",
    x = "PIP (OT axis for matched/OT-only; my PIP for mine-only)",
    y = "PIP (my axis for matched/mine-only; OT PIP for OT-only)",
    shape = "SNP category",
    color = "SNP category"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal()
ggsave(snakemake@output[["pip_all_plot"]], p2, width = 8, height = 6)

# ============================================================
# PART 2 — CREDIBLE SET ALIGNMENT. Size filter applies ONLY here.
# Anchor: SNPs fine-mapped by me AND present in OT.
# ============================================================

anchor_snps <- full_pip[match_status == "both" & !is.na(my_cs)]

cs_rows <- list()
processed_pairs <- character()

if (nrow(anchor_snps) > 0 && nrow(ot_flat) > 0) {
  anchor_snps[, cs_key := my_cs] # my_cs here is bare CS label; need locus-qualified key
  anchor_bps <- unique(anchor_snps$BP)

  # find which OT credible set(s) each anchor SNP belongs to
  ot_membership <- ot_flat[m_bp %in% anchor_bps, .(m_bp, ot_cs_id)]

  # rebuild locus-qualified key for anchor SNPs from my_all
  anchor_with_locus <- merge(
    anchor_snps[, .(BP, my_cs)],
    my_all[, .(BP, locus, CS)],
    by.x = c("BP", "my_cs"),
    by.y = c("BP", "CS")
  )
  anchor_with_locus[, cs_key := paste(locus, my_cs, sep = "_")]

  pairs <- unique(merge(
    anchor_with_locus[, .(BP, cs_key)],
    ot_membership,
    by.x = "BP",
    by.y = "m_bp"
  )[, .(cs_key, ot_cs_id)])

  for (p in seq_len(nrow(pairs))) {
    cs_key_val <- pairs$cs_key[p]
    ot_cs_id_val <- pairs$ot_cs_id[p]
    pair_id <- paste(cs_key_val, ot_cs_id_val, sep = "||")
    if (pair_id %in% processed_pairs) {
      next
    }
    processed_pairs <- c(processed_pairs, pair_id)

    my_size_check <- nrow(my_all[cs_key == cs_key_val])
    ot_size_check <- nrow(ot_flat[ot_cs_id == ot_cs_id_val])

    if (
      my_size_check > cs_size_threshold || ot_size_check > cs_size_threshold
    ) {
      next # size filter applies ONLY in this section
    }

    my_bps <- my_all[cs_key == cs_key_val, BP]
    ot_bps <- ot_flat[ot_cs_id == ot_cs_id_val, m_bp]
    n_shared <- length(intersect(my_bps, ot_bps))
    union_n <- length(union(my_bps, ot_bps))
    jaccard <- if (union_n > 0) n_shared / union_n else 0

    locus_val <- unique(my_all[cs_key == cs_key_val, locus])[1]
    ot_row_ref <- ot_flat[ot_cs_id == ot_cs_id_val][1]

    cs_rows[[length(cs_rows) + 1]] <- data.table(
      locus = locus_val,
      my_cs = cs_key_val,
      my_cs_size = my_size_check,
      ot_cs_id = ot_cs_id_val,
      ot_credible_set_index = ot_row_ref$ot_credible_set_index,
      ot_finemapping_method = ot_row_ref$ot_finemapping_method,
      ot_cs_size = ot_size_check,
      n_shared_snps = n_shared,
      jaccard = jaccard
    )
  }
}

cs_out <- if (length(cs_rows) > 0) {
  rbindlist(cs_rows, fill = TRUE)
} else {
  data.table(
    locus = character(),
    my_cs = character(),
    my_cs_size = integer(),
    ot_cs_id = character(),
    ot_credible_set_index = character(),
    ot_finemapping_method = character(),
    ot_cs_size = integer(),
    n_shared_snps = integer(),
    jaccard = numeric()
  )
}

cat(
  sprintf(
    "%s: credible set alignment — %d (my CS, OT CS) pairs compared, covering %d distinct OT credible sets\n",
    Sys.time(),
    nrow(cs_out),
    uniqueN(cs_out$ot_cs_id)
  ),
  file = log_path,
  append = TRUE
)

fwrite(cs_out, snakemake@output[["cs_level"]], sep = "\t")
