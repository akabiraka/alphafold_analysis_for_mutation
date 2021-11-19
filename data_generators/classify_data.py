import sys
sys.path.append("../alphafold_analysis_for_mutation")
import pandas as pd
from objects.Selector import ChainAndAminoAcidSelect
from objects.PDBData import PDBData

def get_row_items(row):
    pdb_id=PDBData.validate_pdb_id(row["pdb_id"])
    chain_id=row["chain_id"]
    mutation=row["mutation"]
    wild_residue=row["wild_residue"]
    mutant_residue=row["mutant_resiude"]
    mutation_site=row["mutant_reside_num"]
    ddg=row["ddG"]
    inv_pdb_id=PDBData.validate_pdb_id(row["inv_pdb_id"])
    inv_chain_id=row["inv_chain_id"]
    return pdb_id, chain_id, inv_pdb_id, inv_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg

# configurations
pdb_dir = "data/pdbs/"
clean_pdb_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/ssym_684_copy.xlsx"
out_file_path = "data/ssym_684_classified.csv"
n_proteins_to_skip = 233
n_proteins_to_evalutate = 10000

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)

# data generation
dfs = pd.read_excel(input_file_path, skiprows=[1])

if __name__=="__main__":
    """tasks: 
        1. download proteins
        2. clean pdbs using selected chain id and standard amino acids
        3. generate fasta
        4. compute analysis perspective: seq len, mutation type, sec structure type at mutation site
    """
    for i, row in dfs.iterrows():
        if i+1 <= n_proteins_to_skip: continue
        print("Row no: {}".format(i+1))
    
        pdb_id, chain_id, inv_pdb_id, inv_chain_id, wild_residue, mutation, mutant_residue, mutation_site, ddg = get_row_items(row)
        
        cln_wild_pdb_file=clean_pdb_dir+pdb_id+chain_id+".pdb"
        cln_mutant_pdb_file=clean_pdb_dir+inv_pdb_id+inv_chain_id+".pdb"
        wild_fasta_file=fastas_dir+pdb_id+chain_id+".fasta"
        mutant_fasta_file=fastas_dir+inv_pdb_id+inv_chain_id+".fasta"
        
        PDBData.download_structure(pdb_id=pdb_id)
        PDBData.download_structure(pdb_id=inv_pdb_id)
        
        clean_wild_protein_structure = PDBData.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))
        clean_mutant_protein_structure = PDBData.clean(pdb_id=inv_pdb_id, chain_id=inv_chain_id, selector=ChainAndAminoAcidSelect(inv_chain_id))
        
        wild_fasta, wild_seq_len = PDBData.generate_fasta_from_pdb(pdb_id=pdb_id, chain_id=chain_id, input_pdb_filepath=cln_wild_pdb_file, 
                                                                save_as_fasta=True, output_fasta_file=wild_fasta_file)
        mutant_fasta, mutant_seq_len = PDBData.generate_fasta_from_pdb(pdb_id=inv_pdb_id, chain_id=inv_chain_id, input_pdb_filepath=cln_mutant_pdb_file,
                                                                    save_as_fasta=True, output_fasta_file=mutant_fasta_file)
        
        wild_ss_type = PDBData.get_secondary_structure_at_residue(pdb_id, chain_id, cln_wild_pdb_file, mutation_site)
        mutant_ss_type = PDBData.get_secondary_structure_at_residue(inv_pdb_id, inv_chain_id, cln_mutant_pdb_file, mutation_site)
        
        dfs.loc[i, "wild_seq_len"] = wild_seq_len 
        dfs.loc[i, "mutant_seq_len"] = mutant_seq_len
        dfs.loc[i, "wild_ss_type"] = wild_ss_type
        dfs.loc[i, "mutant_ss_type"] = mutant_ss_type
        dfs.loc[i, "mutation_type"] = "stabilizing" if ddg>=0 else "destabilizing"
        # print(dfs.head())
        print()
        if i+1 == n_proteins_to_skip+n_proteins_to_evalutate: 
            break
        
    dfs.to_csv(out_file_path, index=False)
    
    
        

        