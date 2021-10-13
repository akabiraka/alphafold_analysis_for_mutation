import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd
import subprocess

# configs
# stabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_stabilizing.xlsx"
# destabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_destabilizing.xlsx"

# stabilizing_downsampled_dfs = pd.read_excel(stabilizing_downsampled_dfs_filepath)
# destabilizing_downsampled_dfs = pd.read_excel(destabilizing_downsampled_dfs_filepath)

def compute_data_distribution(input_file="data/ssym_downsampled_stabilizing.xlsx", output_file="data/stabilizing_downsampled_data_distribution.xlsx"):
    dfs = pd.read_excel(input_file)
    pdb_ids = dfs["pdb_id"].unique().tolist()
    x = pd.DataFrame()

    # print(stabilizing_downsampled_dfs[stabilizing_downsampled_dfs["pdb_id"]=="1bni"]["seq_len"].unique()[0])
    for i, pdb_id in enumerate(pdb_ids):
        # print(pdb_id)
        # print(dfs[dfs["pdb_id"]==pdb_id]["seq_len"].unique()[0])
        # print(' '.join(dfs[dfs["pdb_id"]==pdb_id]["inverse_pdb_id"].tolist()))
        x.loc[i, "pdb_id"] = pdb_id
        x.loc[i, "&1"] = "&"
        x.loc[i, "seq_len"] = dfs[dfs["pdb_id"]==pdb_id]["seq_len"].unique()[0]
        x.loc[i, "&2"] = "&"
        x.loc[i, "inverse_pdb_ids"] = ' '.join(dfs[dfs["pdb_id"]==pdb_id]["inverse_pdb_id"].tolist())
        x.loc[i, "newline"] = "\\\\\hline"
    # break

    x.sort_values(by=["seq_len"]).to_excel(output_file, index=False)

compute_data_distribution(input_file="data/ssym_downsampled_stabilizing.xlsx", output_file="outputs/stabilizing_downsampled_data_distribution.xlsx")
compute_data_distribution(input_file="data/ssym_downsampled_destabilizing.xlsx", output_file="outputs/destabilizing_downsampled_data_distribution.xlsx")
# destabilizing_downsampled_dfs.sort_values(by=["pdb_id", "seq_len"], inplace=True, ignore_index=True)
# print(destabilizing_downsampled_dfs.head())
# destabilizing_downsampled_dfs.to_excel(destabilizing_downsampled_dfs_filepath, index=False)
# print(destabilizing_downsampled_dfs[["pdb_id", "seq_len"]].drop_duplicates().sort_values(by=["seq_len"]))