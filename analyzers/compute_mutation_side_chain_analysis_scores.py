import sys
sys.path.append("../alphafold_analysis_for_mutation")

from Bio.PDB import *
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib.pyplot import cla
import numpy as np
import pandas as pd
from os.path import exists
import glob
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from objects.Selector import ChainAndAminoAcidSelect
from utils.CleanSlate import CleanSlate
cln = CleanSlate()


def clean(inp_pdb_file, chain_id, out_pdb_file):
    chains=list(PDBParser(QUIET=True).get_structure("", inp_pdb_file)[0].get_chains())
    if len(chains)==1: return inp_pdb_file, chains[0].id
    
    structure = PDBParser(QUIET=True).get_structure("", inp_pdb_file)
    io=PDBIO()
    io.set_structure(structure)
    io.save(out_pdb_file, select=ChainAndAminoAcidSelect(chain_id))
    return out_pdb_file, chain_id

def get_atom_coords(target_residue, predicted_residue):
    """Map atoms and return 3d coordinates of mapped atoms.

    Args:
        target_residue (Bio.PDB.Residue): a ref Bio.PDB.Residue
        predicted_residue (Bio.PDB.Residue): a predicted residue

    Returns:
        two np array: 2 nx3 np array, n indicates number of mapped atoms.
    """
    target_side_chain_atom_coords = []
    predicted_side_chain_atom_coords = []
    
    target_main_chain_atom_coords = []
    predicted_main_chain_atom_coords = []
    for atom in target_residue.get_atoms():
        # print(atom.name, predicted_residue[atom.name].name)
        # print(atom.get_coord(), predicted_residue[atom.name].get_coord())
        # print(atom.name)
        if atom.name in ["N", "CA", "C"]:
            # print("Main chain atom: ", atom.name)
            target_main_chain_atom_coords.append(atom.get_coord())
            predicted_main_chain_atom_coords.append(predicted_residue[atom.name].get_coord())
            
        if atom.name not in ["N", "CA", "C"] and predicted_residue.has_id(atom.name):
            # print("Side chain atom: ", atom.name)
            target_side_chain_atom_coords.append(atom.get_coord())
            predicted_side_chain_atom_coords.append(predicted_residue[atom.name].get_coord())
            
    print(np.array(target_main_chain_atom_coords).shape, np.array(predicted_main_chain_atom_coords).shape, 
            np.array(target_side_chain_atom_coords).shape, np.array(predicted_side_chain_atom_coords).shape)
    return np.array(target_main_chain_atom_coords), np.array(predicted_main_chain_atom_coords), np.array(target_side_chain_atom_coords), np.array(predicted_side_chain_atom_coords)


def compute_side_chain_rmsd_for_mutation_site_(target_pdb_file, target_chain_id, predicted_pdb_file, pred_chain_id, mutation_site):
    target_residue = PDBParser(QUIET=True).get_structure("protein", target_pdb_file)[0][target_chain_id][mutation_site]
    predicted_residue = PDBParser(QUIET=True).get_structure("protein", predicted_pdb_file)[0][pred_chain_id][mutation_site]
    target_main_chain_atom_coords, predicted_main_chain_atom_coords, target_side_chain_atom_coords, predicted_side_chain_atom_coords=get_atom_coords(target_residue, predicted_residue)
    svd_super_imposer = SVDSuperimposer()
    svd_super_imposer.set(reference_coords=target_main_chain_atom_coords, coords=predicted_main_chain_atom_coords)
    svd_super_imposer.run()
    rot, tran = svd_super_imposer.get_rotran()
    predicted_side_chain_atom_coords_transformed = np.dot(predicted_side_chain_atom_coords, rot) + tran
    mse = mean_squared_error(y_true=target_side_chain_atom_coords, y_pred=predicted_side_chain_atom_coords_transformed)
    rmse = np.sqrt(mse)
    return rmse
    
def compute_all_mutation_side_chain_rmsd(inp_file_path, out_file_path):
    mutation_site_df = pd.read_csv(inp_file_path)
    mutation_site_df["rmsd"] = 0.0
    # print(mutation_site_df.head())
    for i, row in mutation_site_df.iterrows():
        pdb_id = row["pdb_id"]
        target_chain_id = pdb_id[-1]
        mutation_site = int(row["mutation_site"])
        target_pdb_file = "data/pdbs_clean/{}.pdb".format(pdb_id)
        predicted_pdb_file_path = glob.glob("data/alphafold2_predicted_pdbs/prediction_{}_*/rank_1_*_unrelaxed.pdb".format(pdb_id[0:4]))
        print(pdb_id[0:4], mutation_site, target_pdb_file, predicted_pdb_file_path)
        tmr_working_dir="data/tmp/"
        if len(predicted_pdb_file_path) > 0:
            pred_pdb_file = predicted_pdb_file_path[0]
            pred_file_name = pred_pdb_file.split("/")[-1].split(".")[0]
            cln_pred_pdb_file, pred_chain_id = clean(pred_pdb_file, target_chain_id, tmr_working_dir+pred_file_name+".pdb")
            rmsd = compute_side_chain_rmsd_for_mutation_site_(target_pdb_file, target_chain_id, cln_pred_pdb_file, pred_chain_id, mutation_site)
            mutation_site_df.loc[i, "rmsd"] = rmsd
            # print(rmsd)
            # break
        cln.clean_all_files(tmr_working_dir, ".pdb") 
    mutation_site_df.to_csv(out_file_path, index=False)    
    # print(mutation_site_df.tail())

    
if __name__ == "__main__":
    compute_all_mutation_side_chain_rmsd(inp_file_path="outputs/score_statistics/local_plddt_conf_0_neighbor.csv",
                                         out_file_path="outputs/score_statistics/side_chain_analysis.csv")