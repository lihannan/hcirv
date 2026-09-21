# hcirv
hcirv analysis script
stat_mutate_type.py is not suitable for continuous base mutations in coding region. If continuous mutations are required, manual verification is required.


### `pairwise_identity.py`
Computes pairwise sequence identity matrix between two sets of sequences using global alignment (Needleman-Wunsch, via Biopython pairwise2).
Input: two FASTA files (target sequences and reference sequences).
Output: CSV matrix (rows = data sequences, columns = target sequences, values = identity fraction).
Dependency: Python ≥3.10, Biopython=1.78, pandas=1.5.3, numpy=1.19.1, seaborn=0.13.2.
Usage: `python scripts/pairwise_identity.py <aim.fasta> <data.fasta>`


### `plot_pairwise_heatmaps.R`
Reads pairwise identity CSV matrices (output from `pairwise_identity.py`) and generates heatmap panels for each matrix using pheatmap. Also extracts self-identity and maximum off-diagonal identity per sequence.
Input: one or more `.csv` files in the working directory (rows = data sequences, columns = target sequences, values = identity fraction).
Output: `test.pdf` (grid of heatmap panels, 4 columns), `test.csv` (per-sequence self-identity and max off-diagonal identity).
Dependency: R (≥4.0), pheatmap=1.0.13, RColorBrewer=1.1-3, gridExtra=2.3.1.
Usage: place CSV files in working directory, then `Rscript plot_pairwise_heatmaps.R`


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


### `fasta_to_tsv_shannong.py`
Converts a multiple sequence alignment in FASTA format into a tab-delimited site-by-sequence table.
Each row represents one alignment position; each column represents one sequence.
Input: aligned FASTA file (all sequences same length).
Output: `<input>.seq` (tab-separated, one position per line, nucleotides separated by tabs).
Dependency: Python =3
Usage: `python fasta_to_seq.py <aligned.fasta>`


### `shannon_entropy.R`
Computes per-site Shannon entropy from a multiple sequence alignment and generates a diversity profile plot.
Input: tab-delimited file (`all.align.fasta.seq`) where each column is a sequence and each row is a site; values are A/T/G/C.
Output: on-screen plot (PDF export can be enabled).
Parameters: sliding window = 20 sites.
Dependency: R (≥4.0), posterior (1.7.0).
Usage: update `setwd()` path, then `Rscript scripts/shannon_entropy.R`


### `all.align.fasta.seq`
input file of `shannon_entropy.R`


### `ORFs_and_CpGs_data.tsv`
input file of `plot_orfs_and_CpGs.R`


