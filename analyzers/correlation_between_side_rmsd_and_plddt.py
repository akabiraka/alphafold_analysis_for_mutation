import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

def plot_correlation_coefficient_wt(df, output_file_name, color):
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
    plt.savefig("output_images/correlation_between_plddt_and_rmsd/{}.pdf".format(output_file_name), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
    
def plot_correlation_coefficient_mt(df, output_file_name):
    plddt_avg = df["avg"].astype(float)
    rmsd = df["rmsd"].astype(float)

    p_r, p_p_value = pearsonr(plddt_avg, rmsd)
    s_r, s_p_value = spearmanr(rmsd, plddt_avg)
    # print(r, p_value)
    stabilizing_plddt_avg = df[df["is_destabilizing"]==0]["avg"]
    stabilizing_rmsd = df[df["is_destabilizing"]==0]["rmsd"]
    destabilizing_plddt_avg = df[df["is_destabilizing"]==1]["avg"]
    destabilizing_rmsd = df[df["is_destabilizing"]==1]["rmsd"]
    
    plt.scatter(stabilizing_plddt_avg, stabilizing_rmsd, color="orange", label="Stabilizing variant")
    plt.scatter(destabilizing_plddt_avg, destabilizing_rmsd, color="red", label="Destabilizing variant")
    plt.legend()
    plt.title("PCC={:.2f}(p-value={:.2f})\nSCC={:.2f}(p-value={:.2f})".format(p_r, p_p_value, s_r, s_p_value))
    plt.xlabel("PLDDT confidence")
    plt.ylabel("RMSD ($\AA$)")
    # plt.show()
    plt.savefig("output_images/correlation_between_plddt_and_rmsd/{}.pdf".format(output_file_name), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()

wt_df = pd.read_excel("outputs/wt_mutation_plddt_statistics_0_neighbor.xlsx")
plot_correlation_coefficient_wt(wt_df, "wild-type", "green")

ssym_df = pd.read_excel("data/ssym_classified_full.xlsx")[["inverse_pdb_id", "inverse_chain_id", "is_destabilizing"]]
mt_df = pd.read_excel("outputs/mt_mutation_plddt_statistics_0_neighbor.xlsx")
mt_df = pd.merge(left=mt_df, right=ssym_df, how="left", left_on="pdb_id", right_on="inverse_pdb_id")


plot_correlation_coefficient_mt(mt_df, "variant")

