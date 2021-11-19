import sys
sys.path.append("../alphafold_analysis_for_mutation")
import os
import pandas as pd
import glob
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
from utils.CleanSlate import CleanSlate
from Bio.PDB import *

from GDTScore import GDTScore
from RMSScore import RMSScore
from TMScore import TMScore

input_file_path = "data/ssym_684_classified.csv"
cln_pdb_dir = "data/pdbs_clean/"
af_pred_dir = "data/alphafold2_predicted_pdbs/"
out_dir = "outputs/scoring/"
n_proteins_to_skip = 233
n_proteins_to_evalutate = 1#0000
dfs = pd.read_csv(input_file_path)

# object initialization
gdt_score = GDTScore()
rms_score = RMSScore()
tm_score = TMScore()
PDBData = PDBData(cln_pdb_dir)
cln = CleanSlate()

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

def run_score_pipeline(ref_pdb_filepath, ref_chain_id, sample_pdb_filepath):
    # print(ref_pdb_filepath, sample_pdb_filepath)
    tms, maxsub, gdt_ts, gdt_ha = gdt_score.get_by_TMscore_software(ref_pdb_filepath, sample_pdb_filepath)
    rms = rms_score.get_by_SVDSuperImposer(ref_pdb_filepath, ref_chain_id, sample_pdb_filepath)
    # tms = tm_score.get(ref_pdb_filepath, sample_pdb_filepath)
    return tms, rms, maxsub, gdt_ts, gdt_ha

def get_af_scores(filepath):
    settings_handle = open(filepath, 'r')
    lines = settings_handle.readlines()
    seq = lines.__getitem__(1).rstrip('\n')[9:]
    data = []
    for line in lines:
        if line.startswith("rank_") and not line.startswith("rank_by"):
            line = line[:-1]
            precicted_pdb_file_name, plddt, ptm_score = line.split(" ")
            plddt = plddt.split(":")[1]
            ptm_score = ptm_score.split(":")[1]
            data.append([precicted_pdb_file_name+"_unrelaxed", len(seq), plddt, ptm_score])
    return data

def clean(inp_pdb_file, chain_id, out_pdb_file):
    structure = PDBParser(QUIET=True).get_structure("", inp_pdb_file)
    io=PDBIO()
    io.set_structure(structure)
    io.save(out_pdb_file, select=ChainAndAminoAcidSelect(chain_id))
    return out_pdb_file

def get_strucutre_comparison_scores(pdb_id, chain_id, ref_pdb_filepath, af_pred_dir):
    tmr_working_dir="data/tmp/"
    scores=[]
    for pred_pdb_file in glob.iglob(af_pred_dir+"*.pdb"):
        pred_file_name = pred_pdb_file.split("/")[-1].split(".")[0]
        cln_pred_pdb_file = clean(pred_pdb_file, chain_id, tmr_working_dir+pred_file_name+".pdb")
        # print(cln_pred_pdb_file)
        tms, maxsub, rms, gdt_ts, gdt_ha = run_score_pipeline(ref_pdb_filepath=ref_pdb_filepath, ref_chain_id=chain_id, sample_pdb_filepath=cln_pred_pdb_file)
        scores.append([pdb_id, chain_id, pred_file_name, tms, rms, maxsub, gdt_ts, gdt_ha])
    cln.clean_all_files(tmr_working_dir, ".pdb")
    return scores

def generate_scores(pdb_id, chain_id, cln_pdb_file, af_pred_dir):
    scores = get_strucutre_comparison_scores(pdb_id, chain_id, cln_pdb_file, af_pred_dir)
    scores_df =  pd.DataFrame(scores, columns=["pdb_id", "chain_id", "af_pdb_filename", "tms", "rms", "maxsub", "gdt_ts", "gdt_ha"])
    
    af_scores_data = get_af_scores(filepath=af_pred_dir+"/settings.txt")        
    af_scores_dfs = pd.DataFrame(af_scores_data, columns=["af_pdb_filename", "seq_len", "plddt", "ptm_score"])
    
    all_scores_dfs = pd.merge(left=scores_df, right=af_scores_dfs, how='left', left_on='af_pdb_filename', right_on='af_pdb_filename')
    all_scores_dfs.to_csv(out_dir+pdb_id+chain_id+".csv", index=False)    
     
    
for i, row in dfs.iterrows():
    if i+1 <= n_proteins_to_skip: continue
    pdb_id, chain_id, mutant_pdb_id, mutant_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg = get_row_items(row)
    print("Row no: ", i+1, pdb_id, chain_id, mutant_pdb_id, mutant_chain_id)
    
    wild_cln_pdb_file=cln_pdb_dir+pdb_id+chain_id+".pdb"
    wild_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+pdb_id+"_*/")
    if len(wild_af_pred_dir)==0: continue # remove this line later
    wild_af_pred_dir = wild_af_pred_dir[0]
    generate_scores(pdb_id, chain_id, wild_cln_pdb_file, wild_af_pred_dir)
    
    
    mutant_cln_pdb_file=cln_pdb_dir+mutant_pdb_id+mutant_chain_id+".pdb"
    mutant_af_pred_dir=glob.glob(af_pred_dir+"prediction_"+mutant_pdb_id+"_*/")
    if len(mutant_af_pred_dir)==0: continue # remove this line later
    mutant_af_pred_dir = mutant_af_pred_dir[0]
    generate_scores(mutant_pdb_id, mutant_chain_id, mutant_cln_pdb_file, mutant_af_pred_dir)
    
    print()
    if i+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
        break     