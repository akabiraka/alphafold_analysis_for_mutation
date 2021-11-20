
import os
import numpy as np
import pickle
import pandas as pd


def get_wt_mt_unique_pdb_ids(with_chain_id=True):
    """Returns wild-type or mutant-type unique pdb_ids 
    with chain_id attached to it or not.

    Args:
        with_chain_id (bool, optional): [description]. Defaults to True.

    Returns:
        dataframe, dataframe: unique_wt_pdb_ids_df, unique_mt_pdb_ids_df
    """
    ssym_df = pd.read_csv("data/ssym_684_classified.csv") 
    if with_chain_id:
        # computing for each MT protein
        unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]+ssym_df["chain_id"]).unique(), columns=["unique_wt_pdb_ids"])
        # computing for each WT protein
        unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inv_pdb_id"]+ssym_df["inv_chain_id"]).unique(), columns=["unique_mt_pdb_ids"])
        return unique_wt_pdb_ids_df, unique_mt_pdb_ids_df
    
    else:
        # computing for each MT protein
        unique_wt_pdb_ids_df = pd.DataFrame((ssym_df["pdb_id"]).unique(), columns=["unique_wt_pdb_ids"])
        # computing for each WT protein
        unique_mt_pdb_ids_df = pd.DataFrame((ssym_df["inv_pdb_id"]).unique(), columns=["unique_mt_pdb_ids"])
        return unique_wt_pdb_ids_df, unique_mt_pdb_ids_df


def is_wild_type(pdb_id):
    """Check if the given pdb_id is wild-type or not

    Args:
        pdb_id (string): pdb id

    Returns:
        bool: is wild-type or not
    """
    unique_wt_pdb_ids_df, _ = get_wt_mt_unique_pdb_ids(with_chain_id=False)
    return unique_wt_pdb_ids_df["unique_wt_pdb_ids"].tolist().__contains__(pdb_id)


def is_mutant_type(pdb_id):
    """Check if the given pdb_id is mutant-type or not

    Args:
        pdb_id (string): pdb id

    Returns:
        bool: is mutant-type or not
    """
    _, unique_mt_pdb_ids_df = get_wt_mt_unique_pdb_ids(with_chain_id=False)
    return unique_mt_pdb_ids_df["unique_mt_pdb_ids"].tolist().__contains__(pdb_id)


def __read_pickel(filepath):
    with open(filepath, "rb") as f:
        return pickle.load(f)
    

def get_mutation_site(pdb_id):
    """Returns the mutation site(s) for a pdb_id from dataset.
    For a WT it should be a list of multiple residue numbers , b/c for each WT there are multiple MTs.
    For a MT it would be a list of single  residue number.

    Args:
        pdb_id (string): pdb id

    Returns:
        DF: dataframe of (chain_id, mutant_reside_num)
    """
    ssym_df = pd.read_csv("data/ssym_684_classified.csv") 
    if is_wild_type(pdb_id):
        mutation_site_df = ssym_df[ssym_df["pdb_id"]==pdb_id]["mutant_reside_num"]
        return mutation_site_df.unique().tolist()
    else:
        mutation_site_df = ssym_df[ssym_df["inv_pdb_id"]==pdb_id]["mutant_reside_num"]
        return mutation_site_df.tolist()
    

def get_neighbors_plddt(mutation_site, residue_plddt_dict, neighbor=0):
    """Returns list of plddts of neighbors.

    Args:
        mutation_site (int): residue num, key of residue_plddt_dict
        residue_plddt_dict (dict): [description]
        neighbor (int, optional): Defaults to 0.

    Returns:
        list: list of plddt scores
    """
    neighbors_plddt_scores = []
    for i in range(-neighbor, neighbor+1):
        new_mutation_site = mutation_site+i  
        if new_mutation_site in residue_plddt_dict:
            neighbor_plddt_score = residue_plddt_dict.get(new_mutation_site)
            neighbors_plddt_scores.append(neighbor_plddt_score)
    return neighbors_plddt_scores


def compute_mutation_site_statistics(mutation_site, plddt_scores, neighbor=0):
    """It computes statistics of plddt scores of the mutation site.

    Args:
        mutation_site (int): residue number
        plddt_scores (dictionary): residue num vs plddt score

    Returns:
        float, float, float, float, float: statistics
    """
    all_plddts_for_a_mutation_site = []
    for plddt_score in plddt_scores:
        residue_plddt_dict = plddt_score[3]
        plddt_for_a_mutation_site = np.mean(get_neighbors_plddt(mutation_site, residue_plddt_dict, neighbor))
        all_plddts_for_a_mutation_site.append(plddt_for_a_mutation_site)
        
    # print(all_plddts_for_a_mutation_site)
    min = np.min(all_plddts_for_a_mutation_site)
    max = np.max(all_plddts_for_a_mutation_site)
    avg = np.mean(all_plddts_for_a_mutation_site)
    median = np.median(all_plddts_for_a_mutation_site)
    std = np.std(all_plddts_for_a_mutation_site)
    return min, max, avg, median, std


def compute_mutation_site_specific_plddt_statistics(neighbor=0):
    """This is a caller function to compute plddt statistics on mutation 
    site for all proteins and save them.
    """
    plddt_scores_dir="outputs/plddt_scores/"
    local_plddt_conf_stat_file = "outputs/score_statistics/local_plddt_conf_{}_neighbor.csv".format(neighbor)
    local_plddt_conf_stats = []
    
    for plddt_score_file in os.listdir(plddt_scores_dir):
        pdb_id, chain_id = plddt_score_file[0:4].lower(), plddt_score_file[4]
        print(pdb_id, chain_id)
        plddt_scores = __read_pickel(plddt_scores_dir+plddt_score_file)
        mutation_sites = get_mutation_site(pdb_id)
        
        for mutation_site in mutation_sites:
            min, max, avg, median, std = compute_mutation_site_statistics(int(mutation_site), plddt_scores, neighbor)
            local_plddt_conf_stats.append([pdb_id+chain_id, mutation_site, min, max, avg, median, std])

    local_plddt_conf_stats_df = pd.DataFrame(local_plddt_conf_stats, columns=["pdb_id", "mutation_site", "min", "max", "avg", "median", "std"])    
    local_plddt_conf_stats_df.to_csv(local_plddt_conf_stat_file, index=False)
        
        
if __name__=="__main__":
    compute_mutation_site_specific_plddt_statistics(neighbor=0)
    compute_mutation_site_specific_plddt_statistics(neighbor=1)
    compute_mutation_site_specific_plddt_statistics(neighbor=2)
    compute_mutation_site_specific_plddt_statistics(neighbor=3)