# Project Title

Analysis of AlphaFold2 for Modeling Structures of Wildtype and Variant Protein Sequences. Lets check this.

### Datafile

* ssym_684.xlsx
* ssym_684_copy.xlsx
* ssym_684_classified.csv

### Data generation

* `python data_generators/generate_fasta_for_alphafold2.py`:
* `python data_generators/classify_data.py`:

### Analysis

* `python analyzers/dataset_distribution.py`
* `python analyzers/compute_global_quality_scores.py`
* `python analyzers/compute_global_quality_stats.py`
* `python analyzers/plot_global_quality.py`
* `python analyzers/compute_local_plddt_scores.py`
* `python analyzers/compute_local_plddt_stats.py`
* `python analyzers/plot_local_plddt_conf.py`
* `python analyzers/compute_mutation_side_chain_analysis_scores.py`
* `python analyzers/plot_mutation_side_chain_analysis_scores.py`
* `python analyzers/compute_correlation_rmsd_and_plddt.py`

### Data issues

* Issue 1: 1briA's 1st residue is VAL which does not report CA atom, in that case consider C atom.

### Other software usage

* DSSP for secondary structure comparison.
  * Homepage: https://swift.cmbi.umcn.nl/gv/dssp/
  * Github link and installation guide: https://github.com/PDB-REDO/dssp
* TMscore (Zhang lab) to compute the globabl structure similarities: TM-Score, GDT-TS, GDT-HA and MaxSub-Score.
  * Homepage: https://zhanggroup.org/TM-score/
