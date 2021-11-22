import sys
sys.path.append("../alphafold_analysis_for_mutation")
import pandas as pd
import matplotlib.pyplot as plt

sec_struct_dict = {"H": "Alpha helix",
                   "B": "Beta-bridge",
                   "E": "Strand",
                   "G": "3-10 helix",
                   "I": "Pi helix",
                   "T": "Turn",
                   "S": "Bend",
                   "-": "None"}

score_data_df = pd.read_csv("outputs/score_statistics/ss_scores.csv")
grouped_by_ss = score_data_df.groupby(by="ground_truth_ss", sort=False)

indexes = score_data_df["ground_truth_ss"].unique()
means = grouped_by_ss.mean()["acc"]
stds = grouped_by_ss.std()["acc"]

ss_names, mins, maxs = [], [], []
for index in indexes:
    ss_names.append(sec_struct_dict.get(index))
    mins.append(score_data_df[score_data_df["ground_truth_ss"]==index]["acc"].min())
    maxs.append(score_data_df[score_data_df["ground_truth_ss"]==index]["acc"].max())

print(mins, maxs)
plt.errorbar(ss_names, means, stds, color="green", fmt=".", ecolor="lightgreen", elinewidth=4, label="Seconday structure accuracy")
plt.errorbar(ss_names, means, yerr=[means-mins, maxs-means], fmt='.', color="green", ecolor='salmon', lw=1, capsize=2)
plt.xticks(rotation=45)
plt.legend(loc="best")
# plt.show()
plt.savefig("outputs/images/secondary_structure_acc/acc_at_mutation_site.pdf", dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
plt.close()

