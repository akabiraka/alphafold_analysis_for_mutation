import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd
import numpy as np

from objects.Selector import ChainAndAminoAcidSelect
from objects.PDBData import PDBData
from objects.Mutation import Mutation
from objects.MutationCollection import MutationCollection
from utils.CleanSlate import CleanSlate

# configurations
pdb_dir = "data/pdbs/"
clean_pdb_dir = "data/pdbs_clean/"
fasta_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/ssym_684.xlsx"
ssym_classified_one_by_one_filepath = "data/ssym_classified_one_by_one.csv"
ssym_classified_full_filepath = "data/ssym_classified_full.xlsx"
n_proteins_to_skip = 0
n_proteins_to_evalutate = 1000

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
cln = CleanSlate()
mutation_collection = MutationCollection()

# helper functions
def get_mutation(row):
    mutation = Mutation(
        pdb_id=PDBData.validate_pdb_id(row.pdb_id),
        chain_id=row.chain_id,
        event=row.mutation,
        wild_residue=row.wild_residue,
        mutant_resiude=row.mutant_resiude,
        mutant_reside_num=row.mutant_reside_num,
        ddG=row.ddG,
        inverse_pdb_id=PDBData.validate_pdb_id(row.inv_pdb_id),
        inverse_chain_id=row.inv_chain_id
    )
    # print(mutation)
    return mutation


# data generation
dfs = pd.read_excel(input_file_path, skiprows=[1])
# print(dfs.head())

stabilizing_mutations = []
destabilizing_mutations = []
for row in dfs.itertuples():
    if row.Index+1 <= n_proteins_to_skip: continue
    print("Protein no: {}".format(row.Index+1))
    
    mutation = get_mutation(row)
    
    PDBData.download_structure(pdb_id=mutation.pdb_id)
    PDBData.download_structure(pdb_id=mutation.inverse_pdb_id)
    
    PDBData.download_fasta(pdb_id=mutation.pdb_id, is_save_file=True, fasta_dir=fasta_dir)
    PDBData.download_fasta(pdb_id=mutation.inverse_pdb_id, is_save_file=True, fasta_dir=fasta_dir)
    
    clean_wild_protein_structure = PDBData.clean(pdb_id=mutation.pdb_id, 
                          chain_id=mutation.chain_id,
                          selector=ChainAndAminoAcidSelect(mutation.chain_id))
    clean_mutant_protein_structure = PDBData.clean(pdb_id=mutation.inverse_pdb_id, 
                          chain_id=mutation.inverse_chain_id,
                          selector=ChainAndAminoAcidSelect(mutation.inverse_chain_id))
    
    wild_fasta_downloaded = PDBData.get_fasta_seq_record(pdb_id=mutation.pdb_id, fasta_dir=fasta_dir)
    mutant_fasta_downloaded = PDBData.get_fasta_seq_record(pdb_id=mutation.inverse_pdb_id, fasta_dir=fasta_dir)
    wild_fasta_from_structure, wild_seq_len_from_structure = PDBData.generate_fasta_from_pdb(pdb_id=mutation.pdb_id,
                                           chain_id=mutation.chain_id,
                                           input_pdb_filepath=clean_pdb_dir+mutation.pdb_id+mutation.chain_id+".pdb",
                                           save_as_fasta=True,
                                           output_fasta_dir=fasta_dir)
    mutant_fasta_from_structure, mutant_seq_len_from_structure = PDBData.generate_fasta_from_pdb(pdb_id=mutation.inverse_pdb_id,
                                           chain_id=mutation.inverse_chain_id,
                                           input_pdb_filepath=clean_pdb_dir+mutation.inverse_pdb_id+mutation.inverse_chain_id+".pdb",
                                           save_as_fasta=True,
                                           output_fasta_dir=fasta_dir)
    
    print(len(wild_fasta_downloaded.seq), len(mutant_fasta_downloaded.seq), wild_seq_len_from_structure, mutant_seq_len_from_structure)
    
    # setting sequence length from structures and fastas
    mutation.set_wild_structure_seq_len(wild_seq_len_from_structure)
    mutation.set_mutant_structure_seq_len(mutant_seq_len_from_structure)
    mutation.set_wild_downloaded_fasta_seq_len(len(wild_fasta_downloaded.seq))
    mutation.set_mutant_downloaded_fasta_seq_len(len(mutant_fasta_downloaded.seq))
    
    mutation_collection.save_one_by_one(mutation, filepath=ssym_classified_one_by_one_filepath)
    mutation_collection.append(mutation)
    
    print()
    if row.Index+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
        break
    

mutation_collection.sort()
mutation_collection.save(filepath=ssym_classified_full_filepath)
# print(mutation_collection)
# cln.clean_all()