import Bio.PDB
from matplotlib.pyplot import cla
import numpy as np
import pandas as pd
from os.path import exists
import glob
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

class MutationSideChainAnalysis(object):
    def __init__(self) -> None:
        super().__init__()
        self.pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        self.super_imposer = Bio.PDB.Superimposer()
        self.svd_super_imposer = Bio.SVDSuperimposer.SVDSuperimposer()

    def plot(self):
        for neighbor in range(1):
            ssym_df = pd.read_excel("data/ssym_classified_full.xlsx")[["inverse_pdb_id", "inverse_chain_id", "is_destabilizing"]]
            wt_df = pd.read_excel("outputs/wt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            mt_df = pd.read_excel("outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(neighbor))
            
            fig = plt.figure()
            gs = fig.add_gridspec(2, hspace=0)
            axs = gs.subplots(sharex=True, sharey=True)
            
            
            
            x_wt_mutation_site = wt_df["mutation_site"].astype(int)
            y_wt_rmsd = wt_df["rmsd"].astype(float)
            axs[0].scatter(x_wt_mutation_site, y_wt_rmsd, color='green', label="Wildtype")
            axs[0].legend()
            
            mt_df = pd.merge(left=mt_df, right=ssym_df, how="left", left_on="pdb_id", right_on="inverse_pdb_id")
            
            stabilizing_mt_df = mt_df[mt_df["is_destabilizing"]==0]
            x_mt_mutation_site = stabilizing_mt_df["mutation_site"].astype(int)
            y_mt_rmsd = stabilizing_mt_df["rmsd"].astype(float)
            axs[1].scatter(x_mt_mutation_site, y_mt_rmsd, color='orange', label="Stabilizing variant")
            axs[1].legend()
            
            destabilizing_mt_df = mt_df[mt_df["is_destabilizing"]==1]
            x_mt_mutation_site = destabilizing_mt_df["mutation_site"].astype(int)
            y_mt_rmsd = destabilizing_mt_df["rmsd"].astype(float)
            axs[1].scatter(x_mt_mutation_site, y_mt_rmsd, color='red', label="Destabilizing variant")
            axs[1].legend()
            
            for ax in axs:
                ax.label_outer()
            plt.legend(loc="best")
            plt.xlabel("Mutation site")
            plt.ylabel("RMSD($\AA$)")
            # plt.show()
            plt.savefig("output_images/mutation_site_vs_rmsd_plots/{}-neighborhood.pdf".format(neighbor), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
            plt.close()


    def get_atom_coords(self, target_residue, predicted_residue):
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

    
    def compute_side_chain_rmsd_for_mutation_site_(self, target_pdb_file, predicted_pdb_file, mutation_site):
        target_residue = self.pdb_parser.get_structure("protein", target_pdb_file)[0]["A"][mutation_site]
        predicted_residue = self.pdb_parser.get_structure("protein", predicted_pdb_file)[0]["A"][mutation_site]
        target_main_chain_atom_coords, predicted_main_chain_atom_coords, target_side_chain_atom_coords, predicted_side_chain_atom_coords = self.get_atom_coords(target_residue, predicted_residue)
        self.svd_super_imposer.set(reference_coords=target_main_chain_atom_coords, coords=predicted_main_chain_atom_coords)
        self.svd_super_imposer.run()
        rot, tran = self.svd_super_imposer.get_rotran()
        predicted_side_chain_atom_coords_transformed = np.dot(predicted_side_chain_atom_coords, rot) + tran
        mse = mean_squared_error(y_true=target_side_chain_atom_coords, y_pred=predicted_side_chain_atom_coords_transformed)
        rmse = np.sqrt(mse)
        return rmse
        
    def compute_all_mutation_side_chain_rmsd(self, inp_file_path, out_file_path):
        mutation_site_df = pd.read_csv(inp_file_path)
        mutation_site_df["rmsd"] = 0.0
        # print(mutation_site_df.head())
        for i, row in mutation_site_df.iterrows():
            pdb_id = row["pdb_id"]
            mutation_site = int(row["mutation_site"])
            target_pdb_file = "data/pdbs_clean/{}.pdb".format(pdb_id)
            predicted_pdb_file_path = glob.glob("data/alphafold2_predicted_pdbs/prediction_{}_*/rank_1_*_unrelaxed.pdb".format(pdb_id[0:4]))
            print(pdb_id[0:4], mutation_site, target_pdb_file, predicted_pdb_file_path)
            if len(predicted_pdb_file_path) > 0:
                predicted_pdb_file = predicted_pdb_file_path[0]
                rmsd = self.compute_side_chain_rmsd_for_mutation_site_(target_pdb_file, predicted_pdb_file, mutation_site)
                mutation_site_df.loc[i, "rmsd"] = rmsd
                # print(rmsd)
                # break

        mutation_site_df.to_csv(out_file_path, index=False)    
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
        neighbor = 1
        rmsd = self.mutation_side_chain_analysis.compute_mutation_side_chain_rmsd(target_pdb_file, 
                                                                                 predicted_pdb_file, 
                                                                                 mutation_site, 
                                                                                 neighbor,
                                                                                 chain_id)
        print(rmsd)

    @unittest.skipIf(True, "Takes time.")        
    def test_compute_all_mutation_side_chain_rmsd(self):
        self.mutation_side_chain_analysis.compute_all_mutation_side_chain_rmsd(file_path="outputs/wt_mutation_plddt_statistics_{}_neighbor.xlsx".format(0))
        self.mutation_side_chain_analysis.compute_all_mutation_side_chain_rmsd(file_path="outputs/mt_mutation_plddt_statistics_{}_neighbor.xlsx".format(0))
        # break
    
    @unittest.skipIf(True, "Takes time.")            
    def test_plot(self):
        self.mutation_side_chain_analysis.plot()
        
if __name__ == "__main__":
    # unittest.main()   
    # ssym_df = pd.read_csv("data/ssym_684_classified.csv") 
    # local_plddt_conf_stat_df = pd.read_csv("outputs/score_statistics/local_plddt_conf_{}_neighbor.csv".format(0))
    
    # ssym_df["wt_pdb_id"] = ssym_df["pdb_id"]+ssym_df["chain_id"]
    # unique_wt_pdb_ids_df = pd.DataFrame(ssym_df["wt_pdb_id"].unique(), columns=["wt_pdb_id"])
    # ssym_df["mt_pdb_id"] = ssym_df["inv_pdb_id"]+ssym_df["inv_chain_id"]
    # unique_mt_pdb_ids_df = ssym_df[["mt_pdb_id", "mutation_type"]]
    
    # wt_local_plddt_conf_stat_df = pd.merge(left=local_plddt_conf_stat_df, right=unique_wt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="wt_pdb_id")
    # mt_local_plddt_conf_stat_df = pd.merge(left=local_plddt_conf_stat_df, right=unique_mt_pdb_ids_df, how="inner", left_on="pdb_id", right_on="mt_pdb_id")
    # print(wt_local_plddt_conf_stat_df.shape, mt_local_plddt_conf_stat_df.shape)
    
    mutation_side_chain_analysis =  MutationSideChainAnalysis()
    mutation_side_chain_analysis.compute_all_mutation_side_chain_rmsd(inp_file_path="outputs/score_statistics/local_plddt_conf_0_neighbor.csv",
                                                                      out_file_path="outputs/score_statistics/side_chain_analysis.csv")
    # mutation_side_chain_analysis.plot()     
    
