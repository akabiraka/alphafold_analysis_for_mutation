import sys
sys.path.append("../alphafold_analysis_for_mutation")

import matplotlib.pyplot as plt
from os import listdir
import pandas as pd
import numpy as np
import pickle
import glob
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
from Bio.PDB import *

af_pred_dir = "data/alphafold2_predicted_pdbs/"
n_proteins_to_skip = 0
n_proteins_to_evalutate = 10000
PDBData = PDBData()


def clean(inp_pdb_file, chain_id, out_pdb_file):
    chains=list(PDBParser(QUIET=True).get_structure("", inp_pdb_file)[0].get_chains())
    if len(chains)==1: return inp_pdb_file, chains[0].id
    
    structure = PDBParser(QUIET=True).get_structure("", inp_pdb_file)
    io=PDBIO()
    io.set_structure(structure)
    io.save(out_pdb_file, select=ChainAndAminoAcidSelect(chain_id))
    return out_pdb_file, chain_id
    
def parse_alphafold_pdb(pdb_filepath):
    """Given a AlphaFold predicted pdb filepath, it parses the file 
    and returns a dictionary of residue number and corresponding plddt score.

    Args:
        pdb_filepath (string): alphafold pdb filepath

    Returns:
        dictionary: residue num vs plddt score
    """
    pdb_filepath = open(pdb_filepath, 'r')
    lines = pdb_filepath.readlines()
    residue_plddt_dict = {}
    for line in lines:
        if line.startswith("ATOM"):
            an_atom_info = line.split()
            residue_num = int(an_atom_info[5])
            plddt = float(an_atom_info[10])
            residue_plddt_dict[residue_num] = plddt
            # print(residue_num, plddt)
    return residue_plddt_dict

def compute_plddts_of_a_protein(pdb_id, chain_id, af_pred_dir):

    """Given a alphafold prediction directory, 
    it calls parse_alphafold_pdb for each pdb file,
    and saves the plddt scores protein-wise.
    and returns a list of plddt scores for each protein
        plddt_scores (list): list of (alphafold_predicted_pdb_filename, plddt_dict)
    
    Args:
        dir (string): a alphafold prediction directory that contains n(=45) pdb files.
    """
    tmr_working_dir="data/tmp/"
    plddt_scores=[]
    for pred_pdb_file in glob.iglob(af_pred_dir+"*.pdb"):
        pred_file_name = pred_pdb_file.split("/")[-1].split(".")[0]
        cln_pred_pdb_file, pred_chain_id = clean(pred_pdb_file, chain_id, tmr_working_dir+pred_file_name+".pdb")
        residue_plddt_dict = parse_alphafold_pdb(pdb_filepath=cln_pred_pdb_file)
        plddt_scores.append([pdb_id, pred_chain_id, pred_file_name, residue_plddt_dict])
                        
    with open("outputs/plddt_scores/{}{}.pkl".format(pdb_id.lower(), chain_id), "wb") as f:
        pickle.dump(plddt_scores, f)

def get_row_items(row):
    pdb_id=PDBData.validate_pdb_id(row["pdb_id"])
    chain_id=row["chain_id"]
    mutation=row["mutation"]
    wild_residue=row["wild_residue"]
    mutant_residue=row["mutant_resiude"]
    mutation_site=row["mutant_reside_num"]
    ddg=row["ddG"]
    mutant_pdb_id=PDBData.validate_pdb_id(row["inv_pdb_id"])
    mutant_chain_id=row["inv_chain_id"]
    return pdb_id, chain_id, mutant_pdb_id, mutant_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg   
                 

input_file_path = "data/ssym_684_classified.csv"
dfs = pd.read_csv(input_file_path)

for i, row in dfs.iterrows():
    if i+1 <= n_proteins_to_skip: continue
    pdb_id, chain_id, mutant_pdb_id, mutant_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg = get_row_items(row)
    print("Row no: ", i+1, pdb_id, chain_id, mutant_pdb_id, mutant_chain_id)
    
    wild_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+pdb_id+"_*/")
    if len(wild_af_pred_dir)!=0: # remove this line later
        wild_af_pred_dir = wild_af_pred_dir[0]
        compute_plddts_of_a_protein(pdb_id, chain_id, wild_af_pred_dir)
    
    
    mutant_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+mutant_pdb_id+"_*/")
    if len(mutant_af_pred_dir)!=0: # remove this line later
        mutant_af_pred_dir = mutant_af_pred_dir[0]
        compute_plddts_of_a_protein(mutant_pdb_id, mutant_chain_id, mutant_af_pred_dir)
    
    print()
    if i+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
        break  

