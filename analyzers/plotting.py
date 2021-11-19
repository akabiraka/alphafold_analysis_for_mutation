import sys

from numpy.core.fromnumeric import mean
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
    axs[0].scatter(wt_x, wt_y, color='green', label="WT")
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
def plot_error_bar_statistics_vs_seq_len(wt_pdbs_score_df, mt_pdbs_score_df, column_name="rms", filename="image", xlabel="Dataset", ylabel="Min RMS ($\AA$)"):
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    
    # axs[0].scatter(x_wt_mutation_site, y_wt_plddt_stat, color='green', label="Wildtype")
    # axs[0].legend()
    # axs[1].scatter(x_mt_mutation_site, y_mt_plddt_stat, color='red', label="Variant")
    # axs[1].legend()
    
    seq_lens = wt_pdbs_score_df["seq_len"]
    pdb_ids = wt_pdbs_score_df["pdb_id"]
    mins = wt_pdbs_score_df[column_name+"_min"]
    maxs = wt_pdbs_score_df[column_name+"_max"]
    means = wt_pdbs_score_df[column_name+"_avg"]
    stds = wt_pdbs_score_df[column_name+"_std"]
    # print(pdb_ids.size, means.shape)
    indexes = range(pdb_ids.size)
    axs[0].errorbar(indexes, means, stds, color="green", fmt=".", ecolor="lightgreen", elinewidth=4, label="Wildtype")
    axs[0].legend()
    axs[0].errorbar(indexes, means, yerr=[means-mins, maxs-means], fmt='.', color="green", ecolor='lightgreen', lw=1, capsize=2)
    
    stabilizing_mt_pdbs_score_df = mt_pdbs_score_df[mt_pdbs_score_df["is_destabilizing"]==0]
    pdb_ids = stabilizing_mt_pdbs_score_df["pdb_id"]
    seq_lens = stabilizing_mt_pdbs_score_df["seq_len"]
    mins = stabilizing_mt_pdbs_score_df[column_name+"_min"]
    maxs = stabilizing_mt_pdbs_score_df[column_name+"_max"]
    means = stabilizing_mt_pdbs_score_df[column_name+"_avg"]
    stds = stabilizing_mt_pdbs_score_df[column_name+"_std"]
    indexes = range(pdb_ids.size)
    axs[1].errorbar(indexes, means, stds, color="orange", fmt=".", ecolor="lightgreen", elinewidth=4, label="Stabilizing variant")
    axs[1].legend()
    axs[1].errorbar(indexes, means, yerr=[means-mins, maxs-means], fmt='.', color="orange", ecolor='lightgreen', lw=1, capsize=2)
    
    
    destabilizing_mt_pdbs_score_df = mt_pdbs_score_df[mt_pdbs_score_df["is_destabilizing"]==1]
    pdb_ids = destabilizing_mt_pdbs_score_df["pdb_id"]
    seq_lens = destabilizing_mt_pdbs_score_df["seq_len"]
    mins = destabilizing_mt_pdbs_score_df[column_name+"_min"]
    maxs = destabilizing_mt_pdbs_score_df[column_name+"_max"]
    means = destabilizing_mt_pdbs_score_df[column_name+"_avg"]
    stds = destabilizing_mt_pdbs_score_df[column_name+"_std"]
    indexes = range(pdb_ids.size)
    axs[1].errorbar(indexes, means, stds, color="red", fmt=".", ecolor="lightgreen", elinewidth=4, label="Destabilizing variant")
    axs[1].legend()
    axs[1].errorbar(indexes, means, yerr=[means-mins, maxs-means], fmt='.', color="red", ecolor='lightgreen', lw=1, capsize=2)
    

    # # categories = np.array(mt_pdbs_score_df["is_destabilizing"])
    # # colormap = np.array(['yellow', 'red'])
    # pdb_ids = mt_pdbs_score_df["pdb_id"]
    # seq_lens = mt_pdbs_score_df["seq_len"]
    # mins = mt_pdbs_score_df[column_name+"_min"]
    # maxs = mt_pdbs_score_df[column_name+"_max"]
    # means = mt_pdbs_score_df[column_name+"_avg"]
    # stds = mt_pdbs_score_df[column_name+"_std"]
    # indexes = range(pdb_ids.size)
    # axs[1].errorbar(indexes, means, stds, color="red", fmt=".", ecolor="lightgreen", elinewidth=4, label="Variant")
    # axs[1].legend()
    # axs[1].errorbar(indexes, means, yerr=[means-mins, maxs-means], fmt='.', color="red", ecolor="lightgreen", lw=1, capsize=2)
    
    for ax in axs:
        ax.label_outer()
    
    plt.legend(loc="best")
    plt.xticks([], [])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.show()
    plt.savefig("output_images/global_structure_quality/{}.pdf".format(filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
            
configs = [
    {
        # the lower the better
        "column_name": "rms",
        "filename": "rmsd_vs_dataset.pdf",
        "xlabel": "Dataset",
        "ylabel": "RMSD ($\AA$)"
    },
    {
        # M-score has the value in (0,1], where 1 indicates a perfect match between two structures.
        "column_name":"tms",
        "filename":"tms_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"TM-Score"
    },
    {
        # MaxSub score has a range between 0 and 1, where 1 is an identical pair of structures.
        "column_name":"maxsub",
        "filename":"maxsub_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"MaxSub-Score"
    },
    {
        # a percentage, ranging from 0 to 100. In general, the higher the GDT_TS score, the more closely a model approximates a given reference structure.
        "column_name":"gdt_ts",
        "filename":"gdt_ts_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"GDT-TS"
    },
    {
        # same as gdt-ts. only higher accuracy version of gdt-ts
        "column_name":"gdt_ha",
        "filename":"gdt_ha_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"GDT-HA"
    },
    {
        # 100 is best
        "column_name":"plddt",
        "filename":"plddt_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"PLDDT"
    },
    {
        # same as gdt-ts. only higher accuracy version of gdt-ts
        "column_name":"ptm_score",
        "filename":"ptm_score_vs_dataset.pdf",
        "xlabel":"Dataset",
        "ylabel":"PTM-Score"
    },
]

for config in configs:
    column_name, filename, xlabel, ylabel, dir=config["column_name"], config["filename"], config["xlabel"], config["ylabel"], "output_all_scores/"
    print("plotting ", filename)

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
    # unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inverse_pdb_id"]+ssym_df["inverse_chain_id"]).unique(), columns=["unique_mt_pdb_ids"])
    ssym_df["unique_mt_pdb_ids"] = ssym_df["inverse_pdb_id"]+ssym_df["inverse_chain_id"]
    unique_mt_pdb_ids_df = ssym_df[["unique_mt_pdb_ids", "is_destabilizing"]]
    mt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_mt_pdb_ids")

    # computing for each stabilizing mutations
    # stabilizing_dfs = ssym_df[ssym_df["is_destabilizing"] == False]

    # computing for each destabilizing mutations
    # destabilizing_dfs = ssym_df[ssym_df["is_destabilizing"] == True]

    print(score_data_df.shape, wt_pdbs_score_df.shape, mt_pdbs_score_df.shape)
    # plot_min_scores_vs_seq_len(wt_pdbs_score_df, mt_pdbs_score_df, column_name=column_name, title=title, xlabel=xlabel, ylabel=ylabel)
    wt_pdbs_score_df.sort_values(by="seq_len", inplace=True, ignore_index=True)
    mt_pdbs_score_df.sort_values(by="seq_len", inplace=True, ignore_index=True)
    # print(wt_pdbs_score_df.head(), mt_pdbs_score_df.head())
    plot_error_bar_statistics_vs_seq_len(wt_pdbs_score_df, mt_pdbs_score_df, column_name, filename, xlabel, ylabel)
    # break


