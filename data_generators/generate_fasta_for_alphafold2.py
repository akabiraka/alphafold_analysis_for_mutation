import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd

from objects.Selector import AllChainAndStandardAminoAcidSelect
from objects.PDBData import PDBData
from objects.Mutation import Mutation
from objects.MutationCollection import MutationCollection

# configurations
pdb_dir = "data/for_alphafold2_pdbs/"
clean_pdb_dir = "data/for_alphafold2_pdbs_clean/"
fastas_dir = "data/for_alphafold2_fastas/"
CIF = "mmCif"
input_file_path = "data/ssym_684_copy.xlsx"
n_proteins_to_skip = 0
n_proteins_to_evalutate = 10000

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
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

# this downloads and generate fasta of all chains for AlphaFold experiments.
stabilizing_mutations = []
destabilizing_mutations = []
for i, row in dfs.iterrows():
    if i+1 <= n_proteins_to_skip: continue
    print("Row no: {}".format(i+1))
    
    mutation = get_mutation(row)
    
    PDBData.download_structure(pdb_id=mutation.pdb_id)
    PDBData.download_structure(pdb_id=mutation.inverse_pdb_id)
    
    clean_wild_protein_structure = PDBData.clean(pdb_id=mutation.pdb_id, 
                          chain_id="",
                          selector=AllChainAndStandardAminoAcidSelect())
    clean_mutant_protein_structure = PDBData.clean(pdb_id=mutation.inverse_pdb_id, 
                          chain_id="",
                          selector=AllChainAndStandardAminoAcidSelect())
    
    PDBData.generate_full_fasta_from_pdb(pdb_id=mutation.pdb_id, input_pdb_filepath=clean_pdb_dir+mutation.pdb_id+".pdb", output_fasta_dir=fastas_dir)
    PDBData.generate_full_fasta_from_pdb(pdb_id=mutation.inverse_pdb_id, input_pdb_filepath=clean_pdb_dir+mutation.inverse_pdb_id+".pdb", output_fasta_dir=fastas_dir)
    
    print()
    if i+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
        break