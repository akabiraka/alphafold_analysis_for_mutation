import sys
sys.path.append("../alphafold_analysis_for_mutation")
import os
import pandas as pd
import glob
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
from Bio.PDB import *

input_file_path = "data/ssym_684_classified.csv"
out_file_path = "outputs/score_statistics/ss_scores.csv"
af_pred_dir = "data/alphafold2_predicted_pdbs/"
cln_pdb_dir = "data/pdbs_clean/"
n_proteins_to_skip = 0
n_proteins_to_evalutate = 30000

PDBData = PDBData(cln_pdb_dir)
dfs = pd.read_csv(input_file_path)


def get_row_items(row):
    pdb_id=row["pdb_id"]
    chain_id=row["chain_id"]
    mutation=row["mutation"]
    wild_residue=row["wild_residue"]
    mutant_residue=row["mutant_resiude"]
    mutation_site=row["mutant_reside_num"]
    ddg=row["ddG"]
    mutant_pdb_id=row["inv_pdb_id"]
    mutant_chain_id=row["inv_chain_id"]
    return pdb_id, chain_id, mutant_pdb_id, mutant_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg

def get_one_based_mutation_site(cln_pdb_file, chain_id, mutation_site):
    residue_ids_dict = PDBData.get_residue_ids_dict(pdb_file=cln_pdb_file, chain_id=chain_id)
    return residue_ids_dict.get(mutation_site)+1
    
def clean(inp_pdb_file, chain_id, out_pdb_file):
    chains=list(PDBParser(QUIET=True).get_structure("", inp_pdb_file)[0].get_chains())
    if len(chains)==1: return inp_pdb_file, chains[0].id
    
    structure = PDBParser(QUIET=True).get_structure("", inp_pdb_file)
    io=PDBIO()
    io.set_structure(structure)
    io.save(out_pdb_file, select=ChainAndAminoAcidSelect(chain_id))
    return out_pdb_file, chain_id



def generate_ss_prediction(pdb_id, chain_id, cln_pdb_file, af_pred_dir, mutation_site, one_based_mutation_site):
    tmr_working_dir="data/tmp/"
    acc=0
    wrong_ss_items=[]
    ground_truth_ss = PDBData.get_secondary_structure_at_residue(pdb_id, chain_id, cln_pdb_file, mutation_site)
    for pred_pdb_file in glob.iglob(af_pred_dir+"*.pdb"):
        pred_file_name = pred_pdb_file.split("/")[-1].split(".")[0]
        cln_pred_pdb_file, pred_chain_id = clean(pred_pdb_file, chain_id, tmr_working_dir+pred_file_name+".pdb")
        pred_ss = PDBData.get_secondary_structure_at_residue(pdb_id, pred_chain_id, cln_pred_pdb_file, one_based_mutation_site)
        
        if ground_truth_ss==pred_ss: acc=acc+1
        else: wrong_ss_items.append({pred_file_name: pred_ss})
            
    return ground_truth_ss, acc, wrong_ss_items

pdb_id_set = set()
ss_scores=[]
for i, row in dfs.iterrows():
    if i+1 <= n_proteins_to_skip: continue
    pdb_id, chain_id, mutant_pdb_id, mutant_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg = get_row_items(row)
    print("Row no: ", i+1, pdb_id, chain_id, mutant_pdb_id, mutant_chain_id)
    
    if pdb_id not in pdb_id_set:
        wild_cln_pdb_file=cln_pdb_dir+pdb_id+chain_id+".pdb"
        one_based_mutation_site = get_one_based_mutation_site(wild_cln_pdb_file, chain_id, mutation_site)
            
        wild_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+pdb_id+"_*/")
        if len(wild_af_pred_dir)!=0: # remove this line later
            wild_af_pred_dir = wild_af_pred_dir[0]
            ground_truth_ss, acc, wrong_ss_ids = generate_ss_prediction(pdb_id, chain_id, wild_cln_pdb_file, wild_af_pred_dir, mutation_site, one_based_mutation_site)
            ss_scores.append([pdb_id+chain_id, ground_truth_ss, mutation_site, acc, wrong_ss_ids])
            print(pdb_id+chain_id, ground_truth_ss, mutation_site, acc, wrong_ss_ids)
        pdb_id_set.add(pdb_id)
    
    if mutant_pdb_id not in pdb_id_set:
        mutant_cln_pdb_file=cln_pdb_dir+mutant_pdb_id+mutant_chain_id+".pdb"
        one_based_mutation_site = get_one_based_mutation_site(mutant_cln_pdb_file, mutant_chain_id, mutation_site)
        
        mutant_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+mutant_pdb_id+"_*/")
        if len(mutant_af_pred_dir)!=0: # remove this line later
            mutant_af_pred_dir = mutant_af_pred_dir[0]
            ground_truth_ss, acc, wrong_ss_ids = generate_ss_prediction(mutant_pdb_id, mutant_chain_id, mutant_cln_pdb_file, mutant_af_pred_dir, mutation_site, one_based_mutation_site)
            ss_scores.append([mutant_pdb_id+mutant_chain_id, ground_truth_ss, mutation_site, acc, wrong_ss_ids])
            print(mutant_pdb_id+mutant_chain_id, ground_truth_ss, mutation_site, acc, wrong_ss_ids)
        pdb_id_set.add(mutant_pdb_id)
        
    print()
    if i+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
        break 
    
ss_scores_df =  pd.DataFrame(ss_scores, columns=["pdb_id", "ground_truth_ss", "mutation_site", "acc", "wrong_ss_ids"])
ss_scores_df.to_csv(out_file_path, index=False)  