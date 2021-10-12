import Bio.PDB
from matplotlib.pyplot import cla
import numpy as np
import pandas as pd
from os.path import exists
import glob
import matplotlib.pyplot as plt

class MutationSideChainAnalysis(object):
    def __init__(self) -> None:
        super().__init__()
        self.pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        self.super_imposer = Bio.PDB.Superimposer()
        self.svd_super_imposer = Bio.SVDSuperimposer.SVDSuperimposer()

    def plot(self):
        color = ["green", "blue", "red"]
        for neighbor in range(3):
            wt_df = pd.read_excel("outputs/wt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            mt_df = pd.read_excel("outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            
            x_wt_mutation_site = wt_df["mutation_site"].astype(int)
            y_wt_rmsd = wt_df["rmsd"].astype(float)
            x_mt_mutation_site = mt_df["mutation_site"].astype(int)
            y_mt_rmsd = mt_df["rmsd"].astype(float)
            
            
            fig = plt.figure()
            gs = fig.add_gridspec(2, hspace=0)
            axs = gs.subplots(sharex=True, sharey=True)
            axs[0].scatter(x_wt_mutation_site, y_wt_rmsd, color='green', label="Wild-type")
            axs[0].legend()
            axs[1].scatter(x_mt_mutation_site, y_mt_rmsd, color='red', label="Variant")
            axs[1].legend()
            for ax in axs:
                ax.label_outer()
            plt.legend(loc="best")
            plt.xlabel("Mutation site")
            plt.ylabel("RMSD($\AA$) ({}-neighborhood)".format(neighbor))
            # plt.show()
            plt.savefig("output_images/mutation_site_rmsd_plots/mutation_site_rmsd_{}_neighbor.pdf".format(neighbor), dpi=300, format="pdf")
            plt.close()


    def __get_atom_coords(self, target_residue, predicted_residue):
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

    def compute_mutation_side_chain_rmsd(self, target_pdb_file, predicted_pdb_file, mutation_site, neighbor=0, chain_id = "A"):
        """Computes RMS using Bio.SVDSuperimposer of two side chains.

        Args:
            target_pdb_file (str): a pdb file path
            predicted_pdb_file (str): a pdb file path
            mutation_site (int): a mutated residue number
            chain_id (str, optional): chain id. Defaults to "A".

        Returns:
            float: SVD super imposed rms value
        """
        all_rmsds = []
        for i in range(-neighbor, neighbor+1):
            new_mutation_site = mutation_site+i  
            chain = self.pdb_parser.get_structure("protein", target_pdb_file)[0][chain_id]
            if chain.has_id(new_mutation_site)!=1: continue
                
            # extracting the mutation site specific residue
            target_residue = self.pdb_parser.get_structure("protein", target_pdb_file)[0][chain_id][new_mutation_site]
            predicted_residue = self.pdb_parser.get_structure("protein", predicted_pdb_file)[0]["A"][new_mutation_site]
        
            target_coords, predicted_coords = self.__get_atom_coords(target_residue, predicted_residue)
        
            # apply the imposer
            self.svd_super_imposer.set(reference_coords=target_coords, coords=predicted_coords)
            self.svd_super_imposer.run()
            # print(svd_super_imposer.get_init_rms(), svd_super_imposer.get_rms())
            all_rmsds.append(self.svd_super_imposer.get_rms())
        
        return all_rmsds

    def compute_all_mutation_side_chain_rmsd(self, neighbor=0, file_path="outputs/mt_mutation_plddt_statistics_0_neighbor.xlsx"):
        mutation_site_df = pd.read_excel(file_path)
        mutation_site_df["rmsd"] = 0.0
        # print(mutation_site_df.head())
        for i, row in mutation_site_df.iterrows():
            pdb_id = row["pdb_id"]
            chain_id = row["chain_id"]
            mutation_site = int(row["mutation_site"])
            target_pdb_file = "data/pdbs/{}{}.pdb".format(pdb_id, chain_id)
            predicted_pdb_file = glob.glob("data/pdbs_alphafold_predicted/prediction_{}_*/rank_1_*_unrelaxed.pdb".format(pdb_id.upper()))[0]
            print(pdb_id, chain_id, mutation_site, target_pdb_file, predicted_pdb_file)
            rmsd = np.mean(self.compute_mutation_side_chain_rmsd(target_pdb_file, predicted_pdb_file, mutation_site, neighbor, chain_id))
            mutation_site_df.loc[i, "rmsd"] = rmsd
            print(rmsd)
            # break

        mutation_site_df.to_excel(file_path, index=False)    
        # print(mutation_site_df.tail())

import unittest
class TestMutationSideChainAnalysis(unittest.TestCase):
    def setUp(self) -> None:
        self.mutation_side_chain_analysis =  MutationSideChainAnalysis()
    
    @unittest.skipIf(True, "Takes time.")
    def test_compute_mutation_side_chain_rmsd(self):
        #, target_pdb_file, predicted_pdb_file, mutation_site, chain_id = "A"):
        target_pdb_file = "data/pdbs/3aa5A.pdb"
        predicted_pdb_file = glob.glob("data/pdbs_alphafold_predicted/prediction_3AA5_*/rank_1_*_unrelaxed.pdb")[0]
        mutation_site = 44
        chain_id = "A"
        neighbor = 3
        rmds = self.mutation_side_chain_analysis.compute_mutation_side_chain_rmsd(target_pdb_file, 
                                                                                 predicted_pdb_file, 
                                                                                 mutation_site, 
                                                                                 neighbor,
                                                                                 chain_id)
        print(rmds)

    @unittest.skipIf(True, "Takes time.")        
    def test_compute_all_mutation_side_chain_rmsd(self):
        for neighbor in range(3):
            self.mutation_side_chain_analysis.compute_all_mutation_side_chain_rmsd(neighbor=neighbor, file_path="outputs/wt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            self.mutation_side_chain_analysis.compute_all_mutation_side_chain_rmsd(neighbor=neighbor, file_path="outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            # break
            
    def test_plot(self):
        self.mutation_side_chain_analysis.plot()
        
if __name__ == "__main__":
    unittest.main()        