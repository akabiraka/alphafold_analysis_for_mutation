import sys
sys.path.append("../alphafold_analysis_for_mutation")

import matplotlib.pyplot as plt
import pandas as pd

def plot_helper_1(wt_plddt_statistics_df, 
                                    mt_plddt_statistics_df, 
                                    ylabel="0-neighborhood"):
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    
    pdb_ids = wt_plddt_statistics_df["pdb_id"]
    mutation_sites = wt_plddt_statistics_df["mutation_site"].astype(str)
    # pdb_id_mutation_sites = wt_plddt_statistics_df["pdb_id"] + wt_plddt_statistics_df["mutation_site"].astype(str)
    mins = wt_plddt_statistics_df["min"]
    maxs = wt_plddt_statistics_df["max"]
    means = wt_plddt_statistics_df["avg"]
    stds = wt_plddt_statistics_df["std"]
    axs[0].errorbar(mutation_sites, means, stds, color="green", fmt=".", ecolor="lightgreen", elinewidth=4, label="Wildtype")
    axs[0].legend()
    axs[0].errorbar(mutation_sites, means, yerr=[means-mins, maxs-means], fmt='.', color="green", ecolor='lightgreen', lw=1, capsize=2)
    
    stabilizing_mt_plddt_statistics_df = mt_plddt_statistics_df[mt_plddt_statistics_df["mutation_type"]=="stabilizing"]
    pdb_ids = stabilizing_mt_plddt_statistics_df["pdb_id"]
    mutation_sites = stabilizing_mt_plddt_statistics_df["mutation_site"].astype(str)
    mins = stabilizing_mt_plddt_statistics_df["min"]
    maxs = stabilizing_mt_plddt_statistics_df["max"]
    means = stabilizing_mt_plddt_statistics_df["avg"]
    stds = stabilizing_mt_plddt_statistics_df["std"]
    axs[1].errorbar(mutation_sites, means, stds, color="orange", fmt=".", ecolor="lightgreen", elinewidth=4, label="Stabilizing variant")
    axs[1].legend()
    axs[1].errorbar(mutation_sites, means, yerr=[means-mins, maxs-means], fmt='.', color="orange", ecolor='lightgreen', lw=1, capsize=2)
    
    destabilizing_mt_plddt_statistics_df = mt_plddt_statistics_df[mt_plddt_statistics_df["mutation_type"]=="destabilizing"]
    pdb_ids = destabilizing_mt_plddt_statistics_df["pdb_id"]
    mutation_sites = destabilizing_mt_plddt_statistics_df["mutation_site"].astype(str)
    mins = destabilizing_mt_plddt_statistics_df["min"]
    maxs = destabilizing_mt_plddt_statistics_df["max"]
    means = destabilizing_mt_plddt_statistics_df["avg"]
    stds = destabilizing_mt_plddt_statistics_df["std"]
    axs[1].errorbar(mutation_sites, means, stds, color="red", fmt=".", ecolor="lightgreen", elinewidth=4, label="Destabilizing variant")
    axs[1].legend()
    axs[1].errorbar(mutation_sites, means, yerr=[means-mins, maxs-means], fmt='.', color="red", ecolor='lightgreen', lw=1, capsize=2)
    # pdb_ids = mt_plddt_statistics_df["pdb_id"]
    # mutation_sites = mt_plddt_statistics_df["mutation_site"].astype(str)
    # # pdb_id_mutation_sites = mt_plddt_statistics_df["pdb_id"] + mt_plddt_statistics_df["mutation_site"].astype(str)
    # mins = mt_plddt_statistics_df["min"]
    # maxs = mt_plddt_statistics_df["max"]
    # means = mt_plddt_statistics_df["avg"]
    # stds = mt_plddt_statistics_df["std"]
    # axs[1].errorbar(mutation_sites, means, stds, color="red", fmt=".", ecolor="lightgreen", elinewidth=4, label="Variant")
    # axs[1].legend()
    # axs[1].errorbar(mutation_sites, means, yerr=[means-mins, maxs-means], fmt='.', color="red", ecolor='lightgreen', lw=1, capsize=2)
    
    for ax in axs:
        ax.label_outer()
    
    plt.legend(loc="best")
    plt.xticks([], [])
    # plt.xticks(rotation='vertical')
    plt.xlabel("Mutation site")
    plt.ylabel("PLDDT ({})".format(ylabel))
    # plt.show()
    plt.savefig("outputs/images/local_plddt_conf/{}.pdf".format(ylabel), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()


def plot_helper_2(grouped_plddt_stat_by_wtpdbid_df, ylabel):
    pdb_ids = grouped_plddt_stat_by_wtpdbid_df.index.values
    mins = grouped_plddt_stat_by_wtpdbid_df["min"]
    maxs = grouped_plddt_stat_by_wtpdbid_df["max"]
    means = grouped_plddt_stat_by_wtpdbid_df["avg"]
    stds = grouped_plddt_stat_by_wtpdbid_df["std"]
    plt.errorbar(pdb_ids, means, stds, color="green", fmt=".", ecolor="lightgreen", elinewidth=4, label="Wildtype")
    plt.errorbar(pdb_ids, means, yerr=[means-mins, maxs-means], fmt='.', color="green", ecolor='lightgreen', lw=1, capsize=2)
    
    plt.legend(loc="best")
    plt.xlabel("Dataset (Wildtype)")
    plt.xticks([], [])
    plt.ylabel("Mutaion site PLDDT in variants ({})".format(ylabel))
    # plt.show()
    plt.savefig("outputs/images/local_plddt_conf_grouped_by_wt/{}.pdf".format(ylabel), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
    
    

def plot_mutation_site_vs_plddt_confidence(unique_wt_pdb_ids_df, unique_mt_pdb_ids_df):
    for neighbor in range(4):
        local_plddt_conf_stat_df = pd.read_csv("outputs/score_statistics/local_plddt_conf_{}_neighbor.csv".format(neighbor))
        wt_local_plddt_conf_stat_df = pd.merge(left=local_plddt_conf_stat_df, right=unique_wt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="wt_pdb_id")
        mt_local_plddt_conf_stat_df = pd.merge(left=local_plddt_conf_stat_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="mt_pdb_id")
        print(wt_local_plddt_conf_stat_df.shape, mt_local_plddt_conf_stat_df.shape)
        plot_helper_1(wt_local_plddt_conf_stat_df, mt_local_plddt_conf_stat_df, ylabel="{}-neighborhood".format(neighbor))
        # break
        
def plot_wt_vs_mutation_site_in_mt_plddt_confidence(ssym_df, unique_mt_pdb_ids_df):
    for neighbor in range(4):
        local_plddt_conf_stat_df = pd.read_csv("outputs/score_statistics/local_plddt_conf_{}_neighbor.csv".format(neighbor))
        mt_local_plddt_conf_stat_df = pd.merge(left=local_plddt_conf_stat_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="mt_pdb_id")
        mt_plddt_stat_with_wtpdbid_df = pd.merge(left=mt_local_plddt_conf_stat_df, right=ssym_df[["wt_pdb_id", "mt_pdb_id"]], how="left",
                                                   left_on="pdb_id", right_on="mt_pdb_id", left_index=False, right_index=False)
        grouped_plddt_stat_by_wtpdbid_df = mt_plddt_stat_with_wtpdbid_df.groupby(by="wt_pdb_id", sort=False).mean()
        plot_helper_2(grouped_plddt_stat_by_wtpdbid_df, ylabel="{}-neighborhood".format(neighbor))
        # break


if __name__=="__main__":
    ssym_df = pd.read_csv("data/ssym_684_classified.csv") 
    ssym_df["wt_pdb_id"] = ssym_df["pdb_id"]+ssym_df["chain_id"]
    unique_wt_pdb_ids_df = pd.DataFrame(ssym_df["wt_pdb_id"].unique(), columns=["wt_pdb_id"])
    ssym_df["mt_pdb_id"] = ssym_df["inv_pdb_id"]+ssym_df["inv_chain_id"]
    unique_mt_pdb_ids_df = ssym_df[["mt_pdb_id", "mutation_type"]] 
    
    plot_mutation_site_vs_plddt_confidence(unique_wt_pdb_ids_df, unique_mt_pdb_ids_df)
    plot_wt_vs_mutation_site_in_mt_plddt_confidence(ssym_df, unique_mt_pdb_ids_df)