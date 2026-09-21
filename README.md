# hcirv
hcirv analysis script
stat_mutate_type.py is not suitable for continuous base mutations in coding region. If continuous mutations are required, manual verification is required.

### `pairwise_identity.py`
Computes pairwise sequence identity matrix between two sets of sequences using global alignment (Needleman-Wunsch, via Biopython pairwise2).
Input: two FASTA files (target sequences and reference sequences).
Output: CSV matrix (rows = data sequences, columns = target sequences, values = identity fraction).
Dependency: Python ≥3.10, Biopython=1.78, pandas=1.5.3, numpy=1.19.1, seaborn=0.13.2.
Usage: `python scripts/pairwise_identity.py <aim.fasta> <data.fasta>`


### `cpg_island.sh`
Runs EMBOSS `cpgplot` to detect CpG islands in HCirV sequences.  
Parameters: window=50 bp, minlen=80 bp, CpG O/E ≥ 0.6, GC% ≥ 60%.  
Input: `*.fasta` files in the working directory.  
Output (per sequence): `.cpgplot` (raw scores), `.gff` (island coordinates), `.png` (profile plot).  
Dependency: EMBOSS cpgplot (version EMBOSS:6.6.0.0), GNU parallel (version 20160622).  
Usage: place FASTA files in working directory, then `bash scripts/cpg_island.sh`.

### `plot_orfs_and_CpGs.R`
Generates linear genome maps showing ORF positions and strand orientation for each HCirV sequence.
Input: tab-delimited file (`molecule`, `start`, `end`, `strand`, `genome_length`, etc.).
Output: PDF file (`<input>.pdf`) with one panel per sequence.
Dependency: R (≥4.0), base graphics only (no external packages required).
Usage: `Rscript scripts/plot_orfs_and_CpGs.R <input.tsv>`
