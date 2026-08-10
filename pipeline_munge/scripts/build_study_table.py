import pandas as pd

raw_csv = "config/07272026_library_leadVariantTable.csv" # Original MPRA variant table
df = pd.read_csv(raw_csv, low_memory=False)

df = df[df["SUMSTATS_AVAILABLE"] == True] # Only take associations with sumstats available
df = df[~df["SUPERPOPULATION"].str.contains("/", na=False)] # For now, only take associations that are of one ancestry (EUR or EAS)

# Group associations in the same study and calculate which fraction of rows does not have fine-mapping
finemap_fraction = df.groupby("STUDY_KEY")["HAS_FINEMAPPING"].apply(
    lambda x: (x == False).mean()
)
eligible_study_keys = finemap_fraction[finemap_fraction >= 0.5].index # Eligible studies are those that are missing FM for at least 50% of their associations

print(f"{len(finemap_fraction)} studies before finemapping-fraction filter")
print(f"{len(eligible_study_keys)} studies retained (>=50% of associations not yet finemapped)")

# Filter down to eligible studies
df = df[df["STUDY_KEY"].isin(eligible_study_keys)]

# Checks that all associations for a given study has matching information
for col in ["SUPERPOPULATION", "TOTAL_SAMPLE_SIZE", "SUMMARY STATS LOCATION", "PUBMEDID"]:
    bad = df.groupby("STUDY_KEY")[col].nunique()
    bad = bad[bad > 1]
    if len(bad) > 0:
        print(f"WARNING: {len(bad)} studies have inconsistent {col} — check before proceeding:")
        print(bad.index.tolist())

# For a given study, drops any duplicate associations. 
# Also, builds CHR:BP string and lists in the same order as lead associations.
def build_position_string(study_key):
    rows = df[df["STUDY_KEY"] == study_key].drop_duplicates(subset="SNPS")
    return ",".join(f"{r['CHR_ID']}:{r['CHR_POS']}" for _, r in rows.iterrows())

# Collapses every study's rows into one row. 
grouped = df.groupby("STUDY_KEY").agg(
    Trait=("DISEASE/TRAIT", "first"),
    Ancestry=("SUPERPOPULATION", "first"),
    n=("TOTAL_SAMPLE_SIZE", "first"),
    GWAS_dir_url=("SUMMARY STATS LOCATION", "first"),
    gcst_id=("GCST_ACCESSION_ID", "first"),
    lead_snps=("SNPS", lambda x: ",".join(x.dropna().unique())),
    pubmed_id=("PUBMEDID", "first"),
).reset_index()

# Add fraction of SNPs fine-mapped
grouped["pct_not_finemapped"] = grouped["STUDY_KEY"].map(finemap_fraction)

# Add lead SNP positions
grouped["lead_snp_positions"] = grouped["STUDY_KEY"].apply(build_position_string)

# Checks that # of lead SNPs and # of lead SNP positions are the same. 
# Flags and drops studies with misalignment.
mismatch_mask = grouped["lead_snps"].str.count(",") != grouped["lead_snp_positions"].str.count(",")
n_mismatched = mismatch_mask.sum()

if n_mismatched > 0:
    print(f"DROPPED {n_mismatched} studies with misaligned SNP/position lists:")
    print(grouped.loc[mismatch_mask, ["STUDY_KEY", "gcst_id"]].to_string(index=False))
    grouped = grouped[~mismatch_mask].reset_index(drop=True)
    
# Sort studies by descending sample size
grouped = grouped.sort_values("n", ascending=False).reset_index(drop=True)

# Create separate study ID for organizational purposes
grouped["study_id"] = ["study" + str(i+1).zfill(3) for i in range(len(grouped))]

# Create log file path
grouped["log_file"] = grouped["gcst_id"].apply(lambda g: f"results/{g}/{g}.log")

grouped.to_csv("config/study_table.tsv", sep="\t", index=False)
print(f"{len(grouped)} studies written, sorted by sample size descending")
print(grouped[["study_id", "gcst_id", "Trait", "n", "pct_not_finemapped"]].head(10))