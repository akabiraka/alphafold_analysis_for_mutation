import sys
sys.path.append("../alphafold_analysis_for_mutation")
import subprocess
import os

class GDTScore(object):
    def __init__(self) -> None:
        super().__init__()
        
    def get_by_TMscore_software(self, ref_pdb_filepath, sample_pdb_filepath):
        result = subprocess.run(["./3rd_party_softwares/TMscore", sample_pdb_filepath, ref_pdb_filepath], stdout=subprocess.PIPE)
        output = result.stdout.decode('utf-8')
        # print(output)
        tm_score = output[913:919]
        maxsub = output[946:952]
        gdt_ts = output[979:985]
        gdt_ha = output[1056:1062]
        return tm_score, maxsub, gdt_ts, gdt_ha

# # sample usage
# gdt_score = GDTScore()
# ref_pdb_filepath = "data/pdbs/1amqA.pdb"
# sample_pdb_filepath = "data/pdbs_alphafold_predicted/prediction_1AMQ_e9684/rank_1_model_3_ptm_seed_0_unrelaxed.pdb"

# tm_score, maxsub, gdt_ts, gdt_ha = gdt_score.get_by_TMscore_software(ref_pdb_filepath, sample_pdb_filepath)
# print(tm_score, maxsub, gdt_ts, gdt_ha)