__author__ = "Anowarul Kabir"
__updated__ = "2020-07-25 14:50:23"

import sys
from numpy import mod
sys.path.append("../protein_data")
import requests
import os

from Bio.PDB import *
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Polypeptide import PPBuilder 
from Bio import SeqIO

class PDBData(object):
    def __init__(self, pdb_dir="data/pdbs/"):
        super(PDBData, self).__init__()
        self.pdb_dir = pdb_dir
        self.pdb_format = "mmCif"
        self.pdb_ext = ".cif"
        self.parser = MMCIFParser(QUIET=True)
        self.pdbl = PDBList()
        self.pdbio = PDBIO()
        self.ppb=PPBuilder()
        
    def validate_pdb_id(self, pdb_id):
        """
        Validate the given pdb_id which must be a 4 letter pdb_id. 
        Else raise an Exception.
        """
        if len(pdb_id) == 4:
            return pdb_id
        elif len(pdb_id) > 4:
            pdb_id = pdb_id[0:4].lower()
        else:
            raise Exception("Reported ", pdb_id, "is not a pdb id in 4 letter format.")
        return pdb_id
        
    def download_structure(self, pdb_id):
        """Download the pdb structure into mmCif format.

        Args:
            pdb_id (string): pdb_id
        """
        self.pdbl.retrieve_pdb_file(pdb_id, pdir=self.pdb_dir, file_format=self.pdb_format)
    
    def clean(self, pdb_id, chain_id, selector=None, clean_pdb_dir="data/pdbs_clean/"):
        """
            Given a Select instance, this function reads a mmCif protein data and saves into
            pdb format.
        Args:
            pdb_id ([string]): [description]
            chain_id ([character]): [description]
            selector: Subclass of Bio.PDB.Select class

        Returns:
            Structure: Return the cleaned structure
        """
        print("Cleaning {}:{} ... ..".format(pdb_id, chain_id))
        pdb_filename = self.pdb_dir + pdb_id + self.pdb_ext
        structure = self.parser.get_structure(pdb_id, pdb_filename)
        self.pdbio.set_structure(structure)
        pdb_filename = clean_pdb_dir + pdb_id+chain_id + ".pdb"
        if selector is None:
            self.pdbio.save(pdb_filename)
        else:
            self.pdbio.save(pdb_filename, select=selector)
        return PDBParser(QUIET=True).get_structure(pdb_id, pdb_filename)
            
    def get_chain_from_structure(self, structure, chain_id):
        """Given a structure and chain_id, it returns 
        corresponding chain of the structure

        Args:
            structure ([Bio.PDB.Structure]): [An instance of Bio.PDB.Structure]
            chain_id ([character]): [chain_id]

        Returns:
            chain: Bio.PDB.Chain
        """
        models = list(structure.get_models())
        chains = list(models[0].get_chains())
        for chain in chains:
            if chain.id == chain_id:
                return chain
            
    def get_chain_from_clean_pdb(self, pdb_id, chain_id, pdb_filename):
        """Return chain given a pdb file.

        Args:
            pdb_id (string)
            chain_id (character)
            pdb_filename (path): A pdb file path. i.e: "data/pdbs_clean/1c8cA.pdb"
        Returns:
            Chain: A chain of a Protein data
        """
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_filename)
        models = list(structure.get_models())
        chains = list(models[0].get_chains())
        # for each chain
        for chain in chains:
            if chain.id == chain_id:
                return chain
            
    def generate_fasta_from_pdb(self, pdb_id, chain_id, input_pdb_filepath, save_as_fasta=False, output_fasta_file=None):
        """Return sequence and length, and create fasta file from pdb data.

        Args:
            pdb_id (string)
            chain_id (string)
            input_pdb_filepath (string): the path must have a pdb file. i.e data/pdbs_clean/1amqA.pdb
        Returns:
            seq: primary sequence of the pdb file
        """
        print("Generating fasta {}:{} ... ..".format(pdb_id, chain_id))
        structure = PDBParser(QUIET=True).get_structure(pdb_id, input_pdb_filepath)
        residues = structure[0][chain_id].get_residues()
        seq = ""
        for residue in residues:
            seq += Polypeptide.three_to_one(residue.get_resname())
        
        if save_as_fasta and output_fasta_file is not None:
            with open(output_fasta_file, "w") as fasta_file_handle:
                fasta_file_handle.write(">{}:{}\n".format(pdb_id.upper(), chain_id))
                fasta_file_handle.write(seq)
        
        return seq, len(seq)
    
    def generate_full_fasta_from_pdb(self, pdb_id, input_pdb_filepath, output_fasta_dir=None, force=False):
        output_fasta_file = "{}{}.fasta".format(output_fasta_dir, pdb_id)
        if os.path.exists(output_fasta_file) and force==False: 
            print("Fasta is already set up for {}. To set-up again, set force=True.".format(pdb_id))
            return
        structure = PDBParser(QUIET=True).get_structure(pdb_id, input_pdb_filepath)
        seqs = []
        for model in structure:
            model_id = model.id
            for chain in model:
                chain_id = chain.id
                seq = ""
                for residue in chain:
                    seq += Polypeptide.three_to_one(residue.get_resname())
                seqs.append(seq)
                with open(output_fasta_file, "a") as fasta_file_handle:
                    fasta_file_handle.write(">{}_{}_{}\n".format(pdb_id, model_id, chain_id))
                    fasta_file_handle.write(seq+"\n")
                    
        seqs_dict = dict.fromkeys(seqs, 0) 
        for seq in seqs:
            seqs_dict[seq] +=1                               
        
        with open(output_fasta_file, "a") as fasta_file_handle:
            fasta_file_handle.write("\nAlphaFold2: sequence n-homooligomers\n")
            for i, seq_key in enumerate(seqs_dict.keys()):
                if i==0: fasta_file_handle.write(seq_key)
                else: fasta_file_handle.write(":"+seq_key)
            fasta_file_handle.write("\n")
            for i, seq_key in enumerate(seqs_dict.keys()):
                if i==0: fasta_file_handle.write(str(seqs_dict[seq]))
                else: fasta_file_handle.write(":"+str(seqs_dict[seq]))
        
    def __save_fasta(self, pdb_id, fasta_text, is_save_file=True, fasta_dir="data/fastas/"):
        """Private method for download_fasta method

        Args:
            pdb_id (string): pdb_id
            fasta_text (string): text in fasta format
            is_save_file (bool, optional): [description]. Defaults to True.
            fasta_dir (str, optional): Directory path. Defaults to "data/fastas/".
        """
        if is_save_file:
            filepath = "{}{}_downloaded.fasta".format(fasta_dir, pdb_id)
            with open(filepath, "w") as fasta_file_handle:
                fasta_file_handle.write(fasta_text)
                fasta_file_handle.close()
        
    def download_fasta(self, pdb_id, is_save_file=True, fasta_dir="data/fastas/"):
        """Download fasta from pdb using pdb_id

        Args:
            pdb_id (string): pdb_id
            is_save_file (bool, optional): Self-explained. Defaults to True.
            fasta_dir (str, optional): Directory path. Defaults to "data/fastas/".
        """
        print("Downloading fasta {} ... ..".format(pdb_id))
        http = requests.Session()
        save_fasta_hook = lambda response, *args, **kwargs: self.__save_fasta(pdb_id, response.text, is_save_file, fasta_dir)
        http.hooks["response"] = save_fasta_hook
        http.get('https://www.rcsb.org/fasta/entry/'+pdb_id)
        
    def get_fasta_seq_record(self, pdb_id, fasta_dir="data/fastas/"):
        """Returns a Bio.SeqRecord.

        Args:
            pdb_id (string): pdb_id
            fasta_dir (string): directory path

        Returns:
            [Bio.SeqRecord]: returns the 1st Bio.SeqRecord from fasta file.
        """
        filepath = "{}{}_downloaded.fasta".format(fasta_dir, pdb_id)
        return next(SeqIO.parse(filepath, "fasta"))
        
    
    def get_seq_from_pdb(self, input_pdb_filepath):   
        """Given a input pdb filepath, it returns the sequence

        Args:
            input_pdb_filepath (string): example: "data/pdbs_clean/1amqA.pdb"

        Returns:
            Bio.Seq.Seq
        """
        structure = PDBParser(QUIET=True).get_structure("id", input_pdb_filepath)
        polypeptide = self.ppb.build_peptides(structure).__getitem__(0)
        return polypeptide.get_sequence()
    
    def get_secondary_structure_at_residue(self, pdb_id, chain_id, cln_pdb_file, residue_pos):
        model = PDBParser(QUIET=True).get_structure(pdb_id, cln_pdb_file)[0]
        residue_id = model[chain_id][residue_pos].id
        dssp_dict = DSSP(model, cln_pdb_file, dssp="mkdssp")
        ss = dssp_dict[chain_id,residue_id][2]
        return ss
    
    def get_residue_ids_dict(self, pdb_file, chain_id):
        residues = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues()
        residue_ids_dict = {residue.id[1]:i for i, residue in enumerate(residues)}
        return residue_ids_dict
