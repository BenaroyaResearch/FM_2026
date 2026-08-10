import pandas as pd


def load_study_info(config, study_table_path="config/study_table.tsv"):
    study_table = pd.read_csv(study_table_path, sep="\t")
    study_table = study_table[study_table["gcst_id"].isin(config["studies"])]
    return study_table.set_index("gcst_id").to_dict(orient="index")
