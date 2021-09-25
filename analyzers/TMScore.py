import tmscoring

class TMScore(object):
    def __init__(self) -> None:
        super().__init__()
        
    def get(self, ref_pdb_filepath, sample_pdb_filepath):
        return tmscoring.get_tm(ref_pdb_filepath, sample_pdb_filepath)

# # sample usage
# tm_score = TMScore()
        
# ref_pdb_filepath = "data/pdbs/1amqA.pdb"
# sample_pdb_filepath = "data/pdbs_alphafold_predicted/prediction_1AMQ_e9684/rank_1_model_3_ptm_seed_0_unrelaxed.pdb"
# print(tm_score.get(ref_pdb_filepath, sample_pdb_filepath))