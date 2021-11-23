import sys
sys.path.append("../alphafold_analysis_for_mutation")

import matplotlib.pyplot as plt
import pandas as pd


def plot_col_histogram(filepath, colums, labels, xlabel, ylabel, output_filename, 
                       bins=15, save=False):
    df = pd.read_csv(filepath)
    for i, col_name in enumerate(colums):
        df[col_name].hist(bins=bins, grid=False, label=labels[i], alpha=0.5)
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if save:
        plt.savefig("outputs/images/data_distribution/{}.pdf".format(output_filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
        plt.close()
    else: plt.show()
    
plot_col_histogram("data/ssym_684_classified.csv", bins=40, save=True, labels=[""], colums=["wild_seq_len"],
                   xlabel="Sequence length", ylabel="Number", output_filename="seq_len_vs_num_of_mutations")

plot_col_histogram("data/ssym_684_classified.csv", bins=30, save=True, labels=[""], colums=["mutant_reside_num"],
                   xlabel="Mutation site", ylabel="Number", output_filename="mutation_site_vs_num_of_mutations")

plot_col_histogram("data/ssym_684_classified.csv", bins=40, save=True, labels=[""], colums=["ddG"],
                   xlabel="$\Delta \Delta G (kcal/mol)$", ylabel="Number", output_filename="ddg_vs_num_of_mutations")

# plot_col_histogram("data/ssym_684_classified.csv", bins=20, save=False, labels=[""], colums=["mutation_type"],
#                    xlabel="Mutation type", ylabel="Number of mutations", output_filename="ddg_vs_num_of_mutations")

# plot_col_histogram("data/ssym_684_classified.csv", bins=20, save=False, labels=[""], colums=["wild_ss_type"],
#                    xlabel="Secondary strucutre types", ylabel="Number of mutations", output_filename="ddg_vs_num_of_mutations")
  

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

# compute_data_distribution(input_file="data/ssym_downsampled_stabilizing.xlsx", output_file="outputs/stabilizing_downsampled_data_distribution.xlsx")
# compute_data_distribution(input_file="data/ssym_downsampled_destabilizing.xlsx", output_file="outputs/destabilizing_downsampled_data_distribution.xlsx")
