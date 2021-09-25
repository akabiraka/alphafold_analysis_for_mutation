import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd

# configs
input_filepath = "data/ssym_classified_full.xlsx"
stabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_stabilizing.xlsx"
destabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_destabilizing.xlsx"

def get_freq_by_seq_len(dfs):
    """Returns the frequency of rows by unique sequence lengths

    Args:
        dfs (DataFrame): a DataFrame object

    Returns:
        key-value pair: sequence lengths, frequencey
    """
    unique_seq_lens = dfs["seq_len"].unique().tolist()
    freq = []
    for seq_len in unique_seq_lens:
        freq.append(len(dfs[dfs["seq_len"]==seq_len]))
    return {"unique_seq_lens": unique_seq_lens,
            "freq": freq}
    
def downsample(dfs, n_samples):
    """Downsample the data by n_samples without replacement.

    Args:
        dfs (DataFrace): a DataFrame object
        n_samples (int): number of sample we want

    Returns:
        [type]: [description]
    """
    unique_mutations = set()
    entries = []    
    unique_seq_lens = dfs["seq_len"].unique().tolist()
    unique_seq_lens.sort(reverse=True)  
    while True:
        for seq_len in unique_seq_lens:
            if len(unique_mutations)==n_samples: 
                return pd.concat(entries, ignore_index=True)
            
            temp_df = dfs[dfs["wild_structure_seq_len"]==seq_len].sample(n=1, replace=False)
            mutation_key = temp_df["pdb_id"].values[0] + temp_df["inverse_pdb_id"].values[0]
            
            if mutation_key not in unique_mutations:
                entries.append(temp_df)
                unique_mutations.add(mutation_key)
                print(unique_mutations)
            
            # return 0
    
    

ssym_classified_dfs = pd.read_excel(input_filepath)
# print(ssym_classified_dfs["is_destabilizing"])

stabilizing_dfs = ssym_classified_dfs[ssym_classified_dfs["is_destabilizing"] == False ]
destabilizing_dfs = ssym_classified_dfs[ssym_classified_dfs["is_destabilizing"] == True]
print(get_freq_by_seq_len(stabilizing_dfs))
print(get_freq_by_seq_len(destabilizing_dfs))

# downsample
stabilizing_downsampled_dfs = downsample(stabilizing_dfs, n_samples=40)
destabilizing_downsampled_dfs = downsample(destabilizing_dfs, n_samples=40)

# save to excel
stabilizing_downsampled_dfs.to_excel(stabilizing_downsampled_dfs_filepath, index=False)
destabilizing_downsampled_dfs.to_excel(destabilizing_downsampled_dfs_filepath, index=False)

# print heads
print(stabilizing_downsampled_dfs.head())
print(destabilizing_downsampled_dfs.head())
