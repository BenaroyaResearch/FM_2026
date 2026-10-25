library(data.table)
library(stringr)

log_path <- snakemake@params[["log_path"]]
gcst_id <- snakemake@wildcards[["study"]]

associations <- fread(snakemake@input[["associations"]])
susie_files <- snakemake@input[["susie_results"]]
susie_lookup <- setNames(susie_files, basename(dirname(dirname(dirname(susie_files)))))

ref_table <- fread(
  snakemake@params[["lead_variant_table"]],
  select = c("SNPS", "finemappingGWAS", "PIP")
)

extract_own_pip <- function(finemap_gwas_str, pip_str, target_gcst) {
  if (is.na(finemap_gwas_str) || finemap_gwas_str == "") {
    return(NA_real_)
  }
  gwas_list <- trimws(str_split(finemap_gwas_str, ",")[[1]])
  pip_list <- trimws(str_split(pip_str, ",")[[1]])
  idx <- which(gwas_list == target_gcst)
  if (length(idx) == 0) {
    return(NA_real_)
  }
  as.numeric(pip_list[idx[1]])
}

ref_table[,
  reference_pip := mapply(
    extract_own_pip,
    finemappingGWAS,
    PIP,
    MoreArgs = list(target_gcst = gcst_id)
  )
]
ref_lookup <- ref_table[!is.na(reference_pip), .(SNPS, reference_pip)]

get_my_lead_pip <- function(locus_id) {
  f <- susie_lookup[[locus_id]]
  if (is.null(f) || !file.exists(f)) {
    return(NA_real_)
  }
  tab <- fread(f)
  row <- tab[SNP == locus_id]
  if (nrow(row) == 0) {
    return(NA_real_)
  }
  row$PIP[1]
}

associations[, my_pip := sapply(SNP, get_my_lead_pip)]
comparison <- merge(
  associations,
  ref_lookup,
  by.x = "SNP",
  by.y = "SNPS",
  all.x = TRUE
)

n_matched <- sum(!is.na(comparison$reference_pip) & !is.na(comparison$my_pip))
n_close <- sum(
  abs(comparison$reference_pip - comparison$my_pip) < 0.2,
  na.rm = TRUE
)

cat(
  sprintf(
    "%s: reference comparison — %d loci matched to prior fine-mapping, %d within 0.2 PIP agreement\n",
    Sys.time(),
    n_matched,
    n_close
  ),
  file = log_path,
  append = TRUE
)

write.table(
  comparison,
  snakemake@output[[1]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
