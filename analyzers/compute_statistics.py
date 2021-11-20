import sys
sys.path.append("../alphafold_analysis_for_mutation")

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def get_statistics(score_df, column_name='tms'):
    # print(score_df.head())
    pdb_id = score_df["pdb_id"].unique()[0]
    seq_len = score_df["seq_len"].dropna().unique()[0]
    min = score_df[column_name].min(axis=0)
    max = score_df[column_name].max(axis=0)
    avg = score_df[column_name].mean(axis=0)
    std = score_df[column_name].std(axis=0)
    return pdb_id, seq_len, min, max, avg, std

inp_dir="outputs/scoring/"
out_dir="outputs/score_statistics/"
score_types = ["rms", "tms", "maxsub", "gdt_ts", "gdt_ha", "plddt", 'ptm_score']
    
for score_type in score_types:
    print("Computing statistics for: ", score_type)

    score_data = []
    for file in os.listdir(inp_dir):
        score_df = pd.read_csv(inp_dir+file)
        pdb_id, seq_len, min, max, avg, std = get_statistics(score_df, score_type)
        score_data.append([pdb_id, seq_len, min, max, avg, std])

    score_data.sort(key=lambda a_protein_scores: a_protein_scores[1]) # sorted on seq_len
    score_data_df = pd.DataFrame(score_data, columns=["pdb_id", "seq_len", score_type+"_min", score_type+"_max", score_type+"_avg", score_type+"_std"])
    
    score_data_df.to_csv(out_dir+score_type+".csv", index=False)