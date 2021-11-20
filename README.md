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

* `python analyzers/compute_global_quality_scores.py`
* `python analyzers/compute_global_quality_stats.py`
* `python analyzers/plot_global_quality.py`
* `python analyzers/compute_local_plddt_score.py`
* `python analyzers/compute_local_plddt_stats.py`
* `python analyzers/plot_local_plddt_conf.py`
