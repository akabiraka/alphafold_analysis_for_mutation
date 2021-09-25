import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd
import subprocess

# configs
stabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_stabilizing.xlsx"
destabilizing_downsampled_dfs_filepath = "data/ssym_downsampled_destabilizing.xlsx"

stabilizing_downsampled_dfs = pd.read_excel(stabilizing_downsampled_dfs_filepath)
destabilizing_downsampled_dfs = pd.read_excel(destabilizing_downsampled_dfs_filepath)

# listing the unique pdb_ids 
stabilizing_downsampled_pdb_list =  stabilizing_downsampled_dfs["pdb_id"].unique().tolist()
stabilizing_downsampled_inverse_pdb_list = stabilizing_downsampled_dfs["inverse_pdb_id"].unique().tolist()
destabilizing_downsampled_pdb_list = destabilizing_downsampled_dfs["pdb_id"].unique().tolist()
destabilizing_downsampled_inverse_pdb_list = destabilizing_downsampled_dfs["inverse_pdb_id"].unique().tolist()

# print the num of unique pdbs and pdb_ids to fold using AlphaFold
print(len(stabilizing_downsampled_pdb_list), stabilizing_downsampled_pdb_list)
print(len(stabilizing_downsampled_inverse_pdb_list), stabilizing_downsampled_inverse_pdb_list)
print(len(destabilizing_downsampled_pdb_list), destabilizing_downsampled_pdb_list)
print(len(destabilizing_downsampled_inverse_pdb_list), destabilizing_downsampled_inverse_pdb_list)
# 13 ['1amq', '2lzm', '1l63', '1ihb', '2rn2', '1iob', '1ey0', '1lz1', '4lyz', '1oh0', '1bni', '1vqb', '5pti']
# 40 ['5eaa', '1l11', '1l85', '1mx4', '1gob', '1hib', '2eyo', '1lhk', '1ir8', '1dmm', '1bse', '1vqj', '1bpt', '1l46', '235l', '1mx6', '3aa5', '2f0i', '2hea', '1lsn', '1bns', '1vqg', '1nag', '1qit', '3c8s', '240l', '2eyp', '1gb5', '1iot', '1bsc', '1d2w', '2eyl', '1yan', '1flu', '1fan', '1qir', '1dye', '2ey5', '1gf9', '1flq']
# 10 ['2lzm', '1l63', '1ihb', '2rn2', '1ey0', '1lz1', '4lyz', '1cey', '1oh0', '1rn1']
# 40 ['1l33', '255l', '1mx2', '3aa3', '1syc', '1gfu', '1hem', '1e6l', '1dmq', '1rgc', '1l57', '1l59', '1rda', '2ey6', '1inu', '1l19', '115l', '1lav', '1ey5', '1e6m', '1l20', '1l68', '1kvc', '1kab', '149l', '1rbu', '1ios', '110l', '1rbs', '1sye', '1lhh', '1rdc', '1gby', '109l', '1e6k', '1l42', '3aa2', '1ey7', '1ouc', '1dyg']


def run_command(command):
   return subprocess.getoutput(command)

def copy_pdb_fasta(pdb_list):
   for pdb_id in pdb_list:
      print(pdb_id)
      pdb_copy_command = "cp ../protein_data/data/pdbs_clean/{}* data/pdbs/".format(pdb_id)
      run_command(command=pdb_copy_command)
      fasta_copy_command = "cp ../protein_data/data/fastas/{}_from_structure.fasta* data/fastas/".format(pdb_id)
      run_command(command=fasta_copy_command)
    
copy_pdb_fasta(stabilizing_downsampled_pdb_list)
copy_pdb_fasta(stabilizing_downsampled_inverse_pdb_list)
copy_pdb_fasta(destabilizing_downsampled_pdb_list)
copy_pdb_fasta(destabilizing_downsampled_inverse_pdb_list)