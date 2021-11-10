import sys
sys.path.append("../alphafold_analysis_for_mutation")

import pandas as pd

# from objects.Mutation import Mutation

class MutationCollection(object):
    def __init__(self) -> None:
        super(MutationCollection).__init__()
        self.mutations = []
        
    def append(self, mutation):
        """Append a Mutation object to the collection.

        Args:
            mutation (Mutation): Mutation type object
        """
        self.mutations.append(mutation)
        
    def sort(self):
        self.mutations.sort(key=lambda mutation: (mutation.is_destabilizing, mutation.wild_structure_seq_len, mutation.ddG))
    
    def save(self, filepath="data/ssym_classified_full.xlsx"):
        y = []
        for mutation in self.mutations:
            y.append(mutation.get_all_attributes())
        df = pd.DataFrame(y, columns = ["pdb_id", "chain_id", 
            "inverse_pdb_id", "inverse_chain_id",
            "wild_structure_seq_len", "mutant_structure_seq_len",
            "wild_downloaded_fasta_seq_len", "mutant_downloaded_fasta_seq_len",
            "event", "wild_residue", "mutant_resiude", "mutant_reside_num", 
            "is_destabilizing", "ddG"]) 
        df.to_excel(filepath, index=False)
        
    def save_one_by_one(self, mutation, filepath="data/ssym_classified.csv"):
        """Append a mutation line one.

        Args:
            mutation (Mutation): A Mutation object
            filepath (str, optional): Save the mutation in filepath. Defaults to "data/ssym_classified.csv".
        """
        y = []
        y.append(mutation.get_all_attributes())
        df = pd.DataFrame(y)
        df.to_csv(filepath, index=False, header=False, mode='a')
            
    def __str__(self):
        """Print the MutationCollection object or multiple mutations

        Returns:
            [type]: [description]
        """
        y = []
        for mutation in self.mutations:
            y.append(mutation.get_all_attributes())
        return y.__str__()