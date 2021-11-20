import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

def plot_correlation_coefficient_wt(df, output_file_name):
    """Plots PCC and SCC values with scatter plot between 
    average PLDDT and RMSD values.

    Args:
        df (dataframe): input dataframe that has at least two column, named "avg" and "rmsd"
        output_file_name (str): the filename to save the figure.
    """
    plddt_avg = df["avg"].astype(float)
    rmsd = df["rmsd"].astype(float)

    p_r, p_p_value = pearsonr(plddt_avg, rmsd)
    s_r, s_p_value = spearmanr(rmsd, plddt_avg)
    # print(r, p_value)
    plt.scatter(plddt_avg, rmsd, color="green", label="Wildtype")
    plt.title("PCC={:.2f}(p-value={:.2f})\nSCC={:.2f}(p-value={:.2f})".format(p_r, p_p_value, s_r, s_p_value))
    plt.xlabel("PLDDT confidence")
    plt.ylabel("RMSD ($\AA$)")
    plt.legend()
    # plt.show()
    plt.savefig("outputs/images/correlation_between_plddt_and_rmsd/{}.pdf".format(output_file_name), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
    
def plot_correlation_coefficient_mt(df, output_file_name):
    """Plots PCC and SCC values with scatter plot between 
    average PLDDT and RMSD values. Also for stabilizing and destabilizing 
    mutations it puts two colors orange and red.

    Args:
        df (dataframe): input dataframe that has at least two column, named "avg" and "rmsd"
        output_file_name (str): the filename to save the figure.
    """
    plddt_avg = df["avg"].astype(float)
    rmsd = df["rmsd"].astype(float)

    p_r, p_p_value = pearsonr(plddt_avg, rmsd)
    s_r, s_p_value = spearmanr(rmsd, plddt_avg)
    # print(r, p_value)
    stabilizing_plddt_avg = df[df["mutation_type"]=="stabilizing"]["avg"]
    stabilizing_rmsd = df[df["mutation_type"]=="stabilizing"]["rmsd"]
    destabilizing_plddt_avg = df[df["mutation_type"]=="destabilizing"]["avg"]
    destabilizing_rmsd = df[df["mutation_type"]=="destabilizing"]["rmsd"]
    
    plt.scatter(stabilizing_plddt_avg, stabilizing_rmsd, color="orange", label="Stabilizing variant")
    plt.scatter(destabilizing_plddt_avg, destabilizing_rmsd, color="red", label="Destabilizing variant")
    plt.legend()
    plt.title("PCC={:.2f}(p-value={:.2f})\nSCC={:.2f}(p-value={:.2f})".format(p_r, p_p_value, s_r, s_p_value))
    plt.xlabel("PLDDT confidence")
    plt.ylabel("RMSD ($\AA$)")
    # plt.show()
    plt.savefig("outputs/images/correlation_between_plddt_and_rmsd/{}.pdf".format(output_file_name), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()


score_data_df = pd.read_csv("outputs/score_statistics/side_chain_analysis.csv")
ssym_df = pd.read_csv("data/ssym_684_classified.csv") 

# plotting correlation coefficient for wildtype proteins.
unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]+ssym_df["chain_id"]).unique(), columns=["unique_wt_pdb_ids"])
wt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_wt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_wt_pdb_ids")
plot_correlation_coefficient_wt(wt_pdbs_score_df, "wild-type")

# plotting correlation coefficient for mutant proteins.
ssym_df["unique_mt_pdb_ids"] = ssym_df["inv_pdb_id"]+ssym_df["inv_chain_id"]
unique_mt_pdb_ids_df = ssym_df[["unique_mt_pdb_ids", "mutation_type"]]
mt_pdbs_score_df = pd.merge(left=score_data_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="unique_mt_pdb_ids")
plot_correlation_coefficient_mt(mt_pdbs_score_df, "variant")
