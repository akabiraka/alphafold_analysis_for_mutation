import Bio.PDB
import numpy as np
import pandas as pd
from os.path import exists
import glob

pdb_parser = Bio.PDB.PDBParser(QUIET=True)
super_imposer = Bio.PDB.Superimposer()
svd_super_imposer = Bio.SVDSuperimposer.SVDSuperimposer()


def __get_atom_coords(target_residue, predicted_residue):
    """Map atoms and return 3d coordinates of mapped atoms.

    Args:
        target_residue (Bio.PDB.Residue): a ref Bio.PDB.Residue
        predicted_residue (Bio.PDB.Residue): a predicted residue

    Returns:
        two np array: 2 nx3 np array, n indicates number of mapped atoms.
    """
    target_side_chain_atom_coords = []
    predicted_side_chain_atom_coords = []
    for atom in target_residue.get_atoms():
        # print(atom.name, predicted_residue[atom.name].name)
        # print(atom.get_coord(), predicted_residue[atom.name].get_coord())
        if predicted_residue.has_id(atom.name):
            target_side_chain_atom_coords.append(atom.get_coord())
            predicted_side_chain_atom_coords.append(predicted_residue[atom.name].get_coord())
    return np.array(target_side_chain_atom_coords), np.array(predicted_side_chain_atom_coords)

def compute_mutation_side_chain_rmsd(target_pdb_file, predicted_pdb_file, mutation_site, chain_id = "A"):
    """Computes RMS using Bio.SVDSuperimposer of two side chains.

    Args:
        target_pdb_file (str): a pdb file path
        predicted_pdb_file (str): a pdb file path
        mutation_site (int): a mutated residue number
        chain_id (str, optional): chain id. Defaults to "A".

    Returns:
        float: SVD super imposed rms value
    """
    # extracting the mutation site specific residue
    target_residue = pdb_parser.get_structure("protein", target_pdb_file)[0][chain_id][mutation_site]
    predicted_residue = pdb_parser.get_structure("protein", predicted_pdb_file)[0]["A"][mutation_site]
    
    target_coords, predicted_coords = __get_atom_coords(target_residue, predicted_residue)
    
    # apply the imposer
    svd_super_imposer.set(reference_coords=target_coords, coords=predicted_coords)
    svd_super_imposer.run()
    # print(svd_super_imposer.get_init_rms(), svd_super_imposer.get_rms())
    return svd_super_imposer.get_rms()

# target_pdb_file = "data/pdbs/3aa5X.pdb"
# predicted_pdb_file = glob.glob("data/pdbs_alphafold_predicted/prediction_3AA5_*/rank_1_*_unrelaxed.pdb")[0]
# mutation_site = 44
# chain_id = "X"
# rms = compute_mutation_side_chain_rmsd(target_pdb_file, predicted_pdb_file, mutation_site, chain_id)
# print(rms)


def compute_all_mutation_side_chain_rmsd(file_path="outputs/mt_mutation_site_plddt_statistics.xlsx"):
    mutation_site_df = pd.read_excel(file_path)
    # print(mutation_site_df.head())
    for i, row in mutation_site_df.iterrows():
        pdb_id = row["pdb_id"]
        chain_id = row["chain_id"]
        mutation_site = int(row["mutation_site"])
        target_pdb_file = "data/pdbs/{}{}.pdb".format(pdb_id, chain_id)
        predicted_pdb_file = glob.glob("data/pdbs_alphafold_predicted/prediction_{}_*/rank_1_*_unrelaxed.pdb".format(pdb_id.upper()))[0]
        # print(pdb_id, mutation_site, target_pdb_file, predicted_pdb_file)
        # print(exists(target_pdb_file), exists(predicted_pdb_file))
        rmsd = compute_mutation_side_chain_rmsd(target_pdb_file, predicted_pdb_file, mutation_site, chain_id)
        # print(pdb_id, chain_id, mutation_site, rmsd)
        mutation_site_df["rmsd"][i] = rmsd
        # break

    mutation_site_df.to_excel(file_path, index=False)    
    # print(mutation_site_df.head())

compute_all_mutation_side_chain_rmsd(file_path="outputs/mt_mutation_site_plddt_statistics.xlsx")
compute_all_mutation_side_chain_rmsd(file_path="outputs/wt_mutation_site_plddt_statistics.xlsx")