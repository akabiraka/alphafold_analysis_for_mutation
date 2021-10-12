import sys
sys.path.append("../alphafold_analysis_for_mutation")

import numpy as np
import pandas as pd
from os import listdir
import matplotlib.pyplot as plt

def get_statistics(score_df, column_name='tms'):
    pdb_id = score_df["pdb_id"].unique()[0]
    seq_len = int(score_df["seq_len"].unique()[0])
    min = score_df[column_name].min(axis=0)
    max = score_df[column_name].max(axis=0)
    avg = score_df[column_name].mean(axis=0)
    std = score_df[column_name].std(axis=0)
    return pdb_id, seq_len, min, max, avg, std
    
def plot_min_scores_vs_seq_len(wt_pdbs_score_df, mt_pdbs_score_df, column_name="rms", title="image", xlabel="Dataset", ylabel="Min RMS ($\AA$)"):
    # two plots vertically aligned
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)

    wt_x = wt_pdbs_score_df["seq_len"]
    wt_y = np.array(wt_pdbs_score_df[column_name+"_min"]).astype(float)
    axs[0].scatter(wt_x, wt_y, color='blue', label="WT")
    axs[0].legend()


    mt_x = mt_pdbs_score_df["seq_len"]
    mt_y = np.array(mt_pdbs_score_df[column_name+"_min"]).astype(float)
    axs[1].scatter(mt_x, mt_y, color='red', label="Variant")
    axs[1].legend()

    for ax in axs:
        ax.label_outer()

    plt.legend(loc='best')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.show()
    plt.savefig("output_images/min_statistics_vs_seq_len_scatter_plots/{}".format(title), dpi=300, format="pdf")
    plt.close()

# columns = ['tms', 'rms', 'maxsub', 'gdt_ts', 'gdt_ha', 'plddt', 'ptm_score']

configs = [
    {
        # the lower the better
        "column_name": "rms",
        "title": "min_rms_vs_seq_len_scatter_plot.pdf",
        "xlabel": "Dataset",
        "ylabel": "Min RMS ($\AA$)"
    },
    {
        # M-score has the value in (0,1], where 1 indicates a perfect match between two structures.
        "column_name":"tms",
        "title":"min_tms_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min TM-score"
    },
    {
        # MaxSub score has a range between 0 and 1, where 1 is an identical pair of structures.
        "column_name":"maxsub",
        "title":"min_maxsub_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min MaxSub-score"
    },
    {
        # a percentage, ranging from 0 to 100. In general, the higher the GDT_TS score, the more closely a model approximates a given reference structure.
        "column_name":"gdt_ts",
        "title":"min_gdt_ts_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min GDT_TS"
    },
    {
        # same as gdt-ts. only higher accuracy version of gdt-ts
        "column_name":"gdt_ha",
        "title":"min_gdt_ha_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min GDT_HA"
    },
    {
        # 100 is best
        "column_name":"plddt",
        "title":"min_plddt_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min PLDDT"
    },
    {
        # same as gdt-ts. only higher accuracy version of gdt-ts
        "column_name":"ptm_score",
        "title":"min_ptm_score_vs_seq_len_scatter_plot.pdf",
        "xlabel":"Dataset",
        "ylabel":"Min PTM-score"
    },
]

for config in configs:
    column_name=config["column_name"]
    title=config["title"]
    xlabel=config["xlabel"]
    ylabel=config["ylabel"]
    dir = "output_all_scores/"

    print("plotting ", title)

    score_data = []
    for file in listdir(dir):
        score_df = pd.read_excel(dir+file)
        pdb_id, seq_len, min, max, avg, std = get_statistics(score_df, column_name)
        score_data.append([pdb_id, seq_len, min, max, avg, std])

    score_data.sort(key=lambda a_protein_scores: a_protein_scores[1]) # sorted on seq_len
    score_data_df = pd.DataFrame(score_data, columns=["pdb_id", "seq_len", column_name+"_min", column_name+"_max", column_name+"_avg", column_name+"_std"])
    # print(score_data_df.head()) 

    ssym_df = pd.read_excel("data/ssym_classified_full.xlsx") 

    # computing for each MT protein
    unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]+ssym_df["chain_id"]).unique(), columns=["unique_wt_pdb_ids"])
    wt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_wt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_wt_pdb_ids")

    # computing for each WT protein
    unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inverse_pdb_id"]+ssym_df["inverse_chain_id"]).unique(), columns=["unique_mt_pdb_ids"])
    mt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_mt_pdb_ids")

    # computing for each stabilizing mutations
    # stabilizing_dfs = ssym_df[ssym_df["is_destabilizing"] == False]

    # computing for each destabilizing mutations
    # destabilizing_dfs = ssym_df[ssym_df["is_destabilizing"] == True]

    print(score_data_df.shape, wt_pdbs_score_df.shape, mt_pdbs_score_df.shape)
    plot_min_scores_vs_seq_len(wt_pdbs_score_df, mt_pdbs_score_df, column_name=column_name, title=title, xlabel=xlabel, ylabel=ylabel)
    
# plotting two plots over each other
# fig, ax = plt.subplots(2, sharex=True, sharey=True)
# wt_x = wt_pdbs_score_df["seq_len"]
# wt_y = np.array(wt_pdbs_score_df[column_name+"_min"]).astype(float)
# ax[0].scatter(wt_x, wt_y, color='blue', label="WT")

# mt_x = mt_pdbs_score_df["seq_len"]
# mt_y = np.array(mt_pdbs_score_df[column_name+"_min"]).astype(float)
# ax[1].scatter(mt_x, mt_y, color='red', label="Variant")

# plt.xticks(rotation="vertical")
# plt.legend(loc='best')
# plt.xlabel("Dataset")
# plt.ylabel("Min RMSD ($\AA$)")
# plt.title("WT-proteins vs. minimum TM-score")
# plt.show()

# print(wt_pdbs_score_df.head())
# print(wt_pdbs_score_df["seq_len"].unique())
# print(mt_pdbs_score_df.head())
# print(mt_pdbs_score_df["seq_len"].unique())

