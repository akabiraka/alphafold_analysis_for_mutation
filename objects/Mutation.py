import sys
sys.path.append("../alphafold_analysis_for_mutation")

class Mutation(object):
    def __init__(self, pdb_id, chain_id, event, wild_residue, 
                 mutant_resiude, mutant_reside_num, ddG, 
                 inverse_pdb_id=None, inverse_chain_id=None):
        
        super(Mutation).__init__()
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.event = event
        self.wild_residue = wild_residue
        self.mutant_resiude = mutant_resiude
        self.mutant_reside_num = mutant_reside_num
        self.is_destabilizing = True if ddG<0 else False
        
        self.inverse_pdb_id = inverse_pdb_id
        self.inverse_chain_id = inverse_chain_id
        
        self.wild_structure_seq_len = 0
        self.mutant_structure_seq_len = 0
        
        self.wild_downloaded_fasta_seq_len = 0
        self.mutant_downloaded_fasta_seq_len = 0
        
        if type(ddG) == int or type(ddG) == float:
            self.ddG = ddG
        else:
            raise Exception("Reported ddG ", ddG, "is not a number.")
        
    def set_wild_downloaded_fasta_seq_len(self, seq_len):
        """Set the wild structure protein sequence length.

        Args:
            seq_len (int): sequence length
        """
        self.wild_downloaded_fasta_seq_len = seq_len
        
    def set_mutant_downloaded_fasta_seq_len(self, seq_len):
        """Set the wild structure protein sequence length.

        Args:
            seq_len (int): sequence length
        """
        self.mutant_downloaded_fasta_seq_len = seq_len
        

    def set_wild_structure_seq_len(self, seq_len):
        """Set the wild structure protein sequence length.

        Args:
            seq_len (int): sequence length
        """
        self.wild_structure_seq_len = seq_len
        
    def set_mutant_structure_seq_len(self, seq_len):
        """Set the mutant structure protein sequence length.

        Args:
            seq_len (int): sequence length
        """
        self.mutant_structure_seq_len = seq_len
            
    def get_all_attributes(self):
        """Returns all necessary attributes of a Mutation.

        Returns:
            array: 1D array
        """
        return [self.pdb_id, self.chain_id, 
                self.inverse_pdb_id, self.inverse_chain_id,
                self.wild_structure_seq_len, self.mutant_structure_seq_len, 
                self.wild_downloaded_fasta_seq_len, self.mutant_downloaded_fasta_seq_len,
                self.event, self.wild_residue, self.mutant_resiude, self.mutant_reside_num, 
                self.is_destabilizing, self.ddG]

    def __str__(self):
        """Overwritten the __str__ method to facilitate the print function.

        Returns:
            string: all attributes in string format.
        """
        return [self.pdb_id, self.chain_id, self.seq_len, 
                self.event, self.wild_residue, self.mutant_resiude, 
                self.mutant_reside_num, self.is_destabilizing,
                self.ddG, self.inverse_pdb_id, self.inverse_chain_id].__str__()
    
    # def get_iterable(self):
    #     return self.all_attributes
        