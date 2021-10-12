import sys
sys.path.append("../alphafold_analysis_for_mutation")

import matplotlib.pyplot as plt
from os import listdir
import pandas as pd
import numpy as np
import pickle


class PLDDTScore(object):
    def __init__(self) -> None:
        super().__init__()
    
    def __read_pickel(self, filepath):
        with open(filepath, "rb") as f:
            return pickle.load(f)
    
    def get_wt_mt_unique_pdb_ids(self, with_chain_id=True):
        """Returns wild-type or mutant-type unique pdb_ids 
        with chain_id attached to it or not.

        Args:
            with_chain_id (bool, optional): [description]. Defaults to True.

        Returns:
            dataframe, dataframe: unique_wt_pdb_ids_df, unique_mt_pdb_ids_df
        """
        ssym_df = pd.read_excel("data/ssym_classified_full.xlsx") 
        if with_chain_id:
            # computing for each MT protein
            unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]+ssym_df["chain_id"]).unique(), columns=["unique_wt_pdb_ids"])
            # computing for each WT protein
            unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inverse_pdb_id"]+ssym_df["inverse_chain_id"]).unique(), columns=["unique_mt_pdb_ids"])
            return unique_wt_pdb_ids_df, unique_mt_pdb_ids_df
        
        else:
            # computing for each MT protein
            unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]).unique(), columns=["unique_wt_pdb_ids"])
            # computing for each WT protein
            unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inverse_pdb_id"]).unique(), columns=["unique_mt_pdb_ids"])
            return unique_wt_pdb_ids_df, unique_mt_pdb_ids_df

    def is_wild_type(self, pdb_id):
        """Check if the given pdb_id is wild-type or not

        Args:
            pdb_id (string): pdb id

        Returns:
            bool: is wild-type or not
        """
        unique_wt_pdb_ids_df, _ = self.get_wt_mt_unique_pdb_ids(with_chain_id=False)
        return unique_wt_pdb_ids_df["unique_wt_pdb_ids"].tolist().__contains__(pdb_id)

    def is_mutant_type(self, pdb_id):
        """Check if the given pdb_id is mutant-type or not

        Args:
            pdb_id (string): pdb id

        Returns:
            bool: is mutant-type or not
        """
        _, unique_mt_pdb_ids_df = self.get_wt_mt_unique_pdb_ids(with_chain_id=False)
        return unique_mt_pdb_ids_df["unique_mt_pdb_ids"].tolist().__contains__(pdb_id)

    def get_mutation_site(self, pdb_id):
        """Returns the mutation site(s) for a pdb_id from dataset.
        For a WT it should be a list of multiple residue numbers , b/c for each WT there are multiple MTs.
        For a MT it would be a list of single  residue number.

        Args:
            pdb_id (string): pdb id

        Returns:
            list: list of residue numbers (or mutation sites)
        """
        ssym_df = pd.read_excel("data/ssym_classified_full.xlsx") 
        if self.is_wild_type(pdb_id):
            mutation_site_df = ssym_df[ssym_df["pdb_id"]==pdb_id]["mutant_reside_num"]
            return mutation_site_df.unique().tolist()
        else:
            mutation_site_df = ssym_df[ssym_df["inverse_pdb_id"]==pdb_id]["mutant_reside_num"]
            return mutation_site_df.tolist()
        
    def plot_residue_plddt_dict(self, dict):
        """Given a dictionary it plots the keys vs. values

        Args:
            dict (dictionary): Dictionary
        """
        print(dict.values())
        print(dict.keys())
        # plt.scatter(x.keys(), x.values())
        plt.plot(dict.keys(), dict.values(), '-o')
        plt.show()
        
    def parse_alphafold_pdb(self, pdb_filepath):
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

    def compute_plddts_of_a_protein(self, pdb_id, dir):

        """Given a alphafold prediction directory, 
        it calls parse_alphafold_pdb for each pdb file
        and returns a list of plddt scores for each protein

        Args:
            dir (string): a alphafold prediction directory that contains n(=45) pdb files.
        Returns:
            plddt_scores (list): list of (alphafold_predicted_pdb, plddt_dict)
        """
        plddt_scores = []
        for pdb_alphafold_file in listdir(dir):
            if pdb_alphafold_file.endswith(".pdb"):
                pdb_filepath = dir+pdb_alphafold_file
                residue_plddt_dict = self.parse_alphafold_pdb(pdb_filepath=pdb_filepath)
                # plot_residue_plddt_dict(residue_plddt_dict)
                plddt_scores.append([pdb_id, pdb_alphafold_file, residue_plddt_dict])

        # saving plddt scores
        # plddt_scores_df = pd.DataFrame(plddt_scores, columns=["pdb_alphafold_file", "residue_plddt_dict"])
        # plddt_scores_df.to_excel("output_plddt_scores/{}.xlsx".format(pdb_id), index=False)                       
        with open("output_plddt_scores/{}.pkl".format(pdb_id), "wb") as f:
            pickle.dump(plddt_scores, f)
                    
    def compute_plddts_of_all_proteins(self):
        alphafold_predicted_pdb_dir = "data/pdbs_alphafold_predicted/"
        for pdb_alphafold_dir in listdir(alphafold_predicted_pdb_dir):
            pdb_id = pdb_alphafold_dir.split("_")[1]
            print(pdb_alphafold_dir.split("_")[1])
            self.compute_plddts_of_a_protein(pdb_id, dir=alphafold_predicted_pdb_dir+pdb_alphafold_dir+"/")
            # break

    def compute_mutation_site_statistics(self, mutation_site, plddt_scores):
        all_plddts_for_a_mutation_site = []
        for plddt_score in plddt_scores:
            residue_plddt_dict = plddt_score[2]
            plddt_for_a_mutation_site = residue_plddt_dict.get(int(mutation_site))
            all_plddts_for_a_mutation_site.append(plddt_for_a_mutation_site)
            
        # print(all_plddts_for_a_mutation_site)
        min = np.min(all_plddts_for_a_mutation_site)
        max = np.max(all_plddts_for_a_mutation_site)
        avg = np.mean(all_plddts_for_a_mutation_site)
        median = np.median(all_plddts_for_a_mutation_site)
        std = np.std(all_plddts_for_a_mutation_site)
        return min, max, avg, median, std
    
    def compute_mutation_site_specific_plddt_statistics(self):
        wt_mutation_site_plddt_statistics = []
        mt_mutation_site_plddt_statistics = []
        plddt_scores_dir="output_plddt_scores/"
        
        for plddt_score_file in listdir(plddt_scores_dir):
            pdb_id = plddt_score_file.split(".")[0].lower()
            print(plddt_score_file.split(".")[0].lower())
            plddt_scores = self.__read_pickel(plddt_scores_dir+plddt_score_file)
            mutation_sites = self.get_mutation_site(pdb_id)
            
            for mutation_site in mutation_sites:
                min, max, avg, median, std = self.compute_mutation_site_statistics(mutation_site, plddt_scores)
                if self.is_wild_type(pdb_id):
                    wt_mutation_site_plddt_statistics.append([pdb_id, mutation_site, min, max, avg, median, std])
                else:
                    mt_mutation_site_plddt_statistics.append([pdb_id, mutation_site, min, max, avg, median, std]) 
            # break

        wt_mutation_site_plddt_statistics_df = pd.DataFrame(wt_mutation_site_plddt_statistics, 
                                                                    columns=["pdb_id", "mutation_site", "min", "max", "avg", "median", "std"])    
        wt_mutation_site_plddt_statistics_df.to_excel("outputs/wt_mutation_site_plddt_statistics.xlsx", index=False)

        mt_mutation_site_plddt_statistics_df = pd.DataFrame(mt_mutation_site_plddt_statistics, 
                                                                    columns=["pdb_id", "mutation_site", "min", "max", "avg", "median", "std"])    
        mt_mutation_site_plddt_statistics_df.to_excel("outputs/mt_mutation_site_plddt_statistics.xlsx", index=False)
            
plddt_score = PLDDTScore()
# plddt_score.compute_plddts_of_all_proteins()
plddt_score.compute_mutation_site_specific_plddt_statistics()


