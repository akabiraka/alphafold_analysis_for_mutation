import sys
sys.path.append("../alphafold_analysis_for_mutation")
from os import listdir
import pandas as pd

from GDTScore import GDTScore
from RMSScore import RMSScore
from TMScore import TMScore

pdb_dir = "data/pdbs/"
alphafold_predicted_pdb_dir = "data/pdbs_alphafold_predicted/"
n_proteins_to_evaluate = 10

gdt_score = GDTScore()
rms_score = RMSScore()
tm_score = TMScore()

def run_score_pipeline(ref_pdb_filepath, sample_pdb_filepath):
    # print(ref_pdb_filepath, sample_pdb_filepath)
    _, maxsub, gdt_ts, gdt_ha = gdt_score.get_by_TMscore_software(ref_pdb_filepath, sample_pdb_filepath)
    rms = rms_score.get_by_SVDSuperImposer(ref_pdb_filepath, sample_pdb_filepath)
    tms = tm_score.get(ref_pdb_filepath, sample_pdb_filepath)
    return tms, rms, maxsub, gdt_ts, gdt_ha

def parse_settings_file(filepath):
    settings_handle = open(filepath, 'r')
    lines = settings_handle.readlines()
    seq = lines.__getitem__(1)[9:]
    data = []
    for line in lines:
        if line.startswith("rank_") and not line.startswith("rank_by"):
            line = line[:-1]
            precicted_pdb_file_name, plddt, ptm_score = line.split(" ")
            plddt = plddt.split(":")[1]
            ptm_score = ptm_score.split(":")[1]
            data.append([precicted_pdb_file_name+"_unrelaxed", len(seq), plddt, ptm_score])
    return data

# # sample usage 
# settings_filepath = "data/pdbs_alphafold_predicted/prediction_1AMQ_e9684/settings.txt"
# parse_settings_file(settings_filepath)


ith_protein=1
# for each pdb file
for pdb_file in listdir(pdb_dir):
    structure_comparison_data = []
    if pdb_file.endswith(".pdb"):
        pdb_id = pdb_file[:4].lower()
        chain_id = pdb_file[4:5]
        print("Start processing: ",ith_protein, pdb_id+chain_id)
        # for each alphafold predicted collection pdbs directory
        for pdb_alphafold_dir in listdir(alphafold_predicted_pdb_dir):
            if "_"+pdb_id+"_" in pdb_alphafold_dir.lower():
                # for each pdb predicted by alphafold corresponding to the pdb file
                for pdb_alphafold_file in listdir(alphafold_predicted_pdb_dir+pdb_alphafold_dir+"/"):
                    if pdb_alphafold_file.endswith(".pdb"):
                        # print(pdb_id, pdb_alphafold_dir, pdb_alphafold_file)
                        # generate the file paths
                        ref_pdb_filepath = pdb_dir+pdb_id+chain_id+".pdb" #"data/pdbs/1amqA.pdb"
                        sample_pdb_filepath = alphafold_predicted_pdb_dir+pdb_alphafold_dir+"/"+pdb_alphafold_file #"data/pdbs_alphafold_predicted/prediction_1AMQ_e9684/rank_1_model_3_ptm_seed_0_unrelaxed.pdb"
                        # compute the scores
                        tms, maxsub, rms, gdt_ts, gdt_ha = run_score_pipeline(ref_pdb_filepath, sample_pdb_filepath)
                        structure_comparison_data.append([pdb_file[:-4], pdb_alphafold_file[:-4], tms, rms, maxsub, gdt_ts, gdt_ha])
                
                # parse data from settings, i.e plddt and ptmScore, for each alphafold predicted pdbs directory 
                settings_filepath = alphafold_predicted_pdb_dir+pdb_alphafold_dir+"/settings.txt"
                settings_data = parse_settings_file(filepath=settings_filepath)
                
                # creating dfs
                structure_comparison_dfs =  pd.DataFrame(structure_comparison_data, columns=["pdb_id", "pdb_alphafold_filename", "tms", "rms", "maxsub", "gdt_ts", "gdt_ha"])
                settings_dfs = pd.DataFrame(settings_data, columns=["pdb_alphafold_filename", "seq_len", "plddt", "ptm_score"])
                # print(structure_comparison_dfs.shape)
                # print(settings_dfs.shape)
                # print(structure_comparison_dfs.head())
                # print(settings_dfs.head())
                # do left join to merge using pdb_alphafold_filename
                
                print("Finished processing: ", pdb_id+chain_id)
                all_scores_dfs = pd.merge(left=structure_comparison_dfs, right=settings_dfs, how='left', left_on='pdb_alphafold_filename', right_on='pdb_alphafold_filename')
                all_scores_dfs.to_excel("outputs/"+pdb_id+chain_id+".xlsx", index=False)                
        #         break    
        # break 
    if ith_protein == n_proteins_to_evaluate:
        break    
    ith_protein=ith_protein+1
         