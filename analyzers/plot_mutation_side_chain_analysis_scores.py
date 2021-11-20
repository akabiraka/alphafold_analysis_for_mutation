
import pandas as pd
import matplotlib.pyplot as plt

def plot(wt_df, mt_df):
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    
    x_wt_mutation_site = wt_df["mutation_site"].astype(int)
    y_wt_rmsd = wt_df["rmsd"].astype(float)
    axs[0].scatter(x_wt_mutation_site, y_wt_rmsd, color='green', label="Wildtype")
    axs[0].legend()
    
    stabilizing_mt_df = mt_df[mt_df["mutation_type"]=="stabilizing"]
    x_mt_mutation_site = stabilizing_mt_df["mutation_site"].astype(int)
    y_mt_rmsd = stabilizing_mt_df["rmsd"].astype(float)
    axs[1].scatter(x_mt_mutation_site, y_mt_rmsd, color='orange', label="Stabilizing variant")
    axs[1].legend()
    
    destabilizing_mt_df = mt_df[mt_df["mutation_type"]=="destabilizing"]
    x_mt_mutation_site = destabilizing_mt_df["mutation_site"].astype(int)
    y_mt_rmsd = destabilizing_mt_df["rmsd"].astype(float)
    axs[1].scatter(x_mt_mutation_site, y_mt_rmsd, color='red', label="Destabilizing variant")
    axs[1].legend()
    
    for ax in axs:
        ax.label_outer()
    plt.legend(loc="best")
    plt.xlabel("Mutation site")
    plt.ylabel("RMSD($\AA$)")
    # plt.show()
    plt.savefig("outputs/images/mutation_site_vs_rmsd_plots/{}-neighborhood.pdf".format(0), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
        
if __name__=="__main__":
    score_data_df = pd.read_csv("outputs/score_statistics/side_chain_analysis.csv")
    ssym_df = pd.read_csv("data/ssym_684_classified.csv") 

    unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]+ssym_df["chain_id"]).unique(), columns=["unique_wt_pdb_ids"])
    wt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_wt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_wt_pdb_ids")

    ssym_df["unique_mt_pdb_ids"] = ssym_df["inv_pdb_id"]+ssym_df["inv_chain_id"]
    unique_mt_pdb_ids_df = ssym_df[["unique_mt_pdb_ids", "mutation_type"]]
    mt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_mt_pdb_ids")
    
    print(wt_pdbs_score_df.shape, mt_pdbs_score_df.shape)
    plot(wt_pdbs_score_df, mt_pdbs_score_df)