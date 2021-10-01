import sys
sys.path.append("../alphafold_analysis_for_mutation")

import Bio.PDB
import Bio.SVDSuperimposer
import numpy as np

class RMSScore(object):
    def __init__(self) -> None:
        super().__init__()
        self.pdb_parser = Bio.PDB.PDBParser(QUIET = True)
        self.super_imposer = Bio.PDB.Superimposer()
        self.svd_super_imposer = Bio.SVDSuperimposer.SVDSuperimposer()

    def __get_CA_atoms(self, chain):
        atoms = []
        for residue in chain:
            atoms.append(residue['CA'])
        return atoms

    def __get_CA_atom_coordinates(self, chain):
        atom_coords = []
        for residue in chain:
            atom_coords.append(residue['CA'].get_coord())
        return np.array(atom_coords)
                
    def get_by_SVDSuperImposer(self, ref_pdb_filepath, sample_pdb_filepath):
        # reading the pdb structures
        ref_structure = self.pdb_parser.get_structure("reference", ref_pdb_filepath)
        sample_structure = self.pdb_parser.get_structure("sample", sample_pdb_filepath)
        
        # get the CA atom coordinates
        ref_atom_coords = self.__get_CA_atom_coordinates(chain=ref_structure[0]["A"])
        sample_atom_coords = self.__get_CA_atom_coordinates(chain=sample_structure[0]["A"])
        
        # apply the imposer
        self.svd_super_imposer.set(reference_coords=ref_atom_coords, coords=sample_atom_coords)
        self.svd_super_imposer.run()
        return self.svd_super_imposer.get_rms()
        
    def get_by_SuperImposer(self, ref_pdb_filepath, sample_pdb_filepath):
        # reading the pdb structures
        ref_structure = self.pdb_parser.get_structure("reference", ref_pdb_filepath)
        sample_structure = self.pdb_parser.get_structure("sample", sample_pdb_filepath)
        
        # get the CA atoms
        ref_atoms = self.__get_CA_atoms(chain=ref_structure[0]["A"])
        sample_atoms = self.__get_CA_atoms(chain=sample_structure[0]["A"])
        
        # apply the imposer
        self.super_imposer.set_atoms(fixed=ref_atoms, moving=sample_atoms)
        self.super_imposer.apply(sample_atoms)
        return self.super_imposer.rms


# sample usage
# rms_score = RMSScore()

# 1mx4 prediction_1MX4_21f84 rank_10_model_4_ptm_seed_4_unrelaxed.pdb
# ref_pdb_filepath = "data/pdbs/1mx4A.pdb"
# sample_pdb_filepath = "data/pdbs_alphafold_predicted/prediction_1MX4_21f84/rank_10_model_4_ptm_seed_4_unrelaxed.pdb"

# # call the imposers
# rms = rms_score.get_by_SuperImposer(ref_pdb_filepath, sample_pdb_filepath)
# print(rms)

# rms = rms_score.get_by_SVDSuperImposer(ref_pdb_filepath, sample_pdb_filepath)
# print(rms)