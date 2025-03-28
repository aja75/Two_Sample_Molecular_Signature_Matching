# Two_Sample_Molecular_Signature_Matching
Compare 2 Differential Expression Profiles in terms of their DEG Overlap % and Concordance validated by permutation
---
## Before Using (Windows):
1. If you want to compile, use auto-py-to-exe or an eqivalent (frontend.py is main)


2. If you don't want to compile, make sure python and dependencies are installed (double click on frontend.py or run from python terminal):
### Dependencies
1. numpy
2. pandas
3. PySimpleGUI-4-foss
4. matplotlib
5. scipy

## How To Use:
1. From Differential Expression (DE) files (use DEseq2, EdgeR, CuffDiff2, etc) select 2 different profiles
2. Specify if DE is formatted as log2FC or Z-Score 'column spelling must match these exactly'
3. Specify thresholds for the binary mapping default = {log2FC = 1, Z-Score = 2}
4. Enter the gene column name, ex. ENSEMBL_ID, GENE_SYMBOL, etc 'must be spelled exactly'
5. Enter permutation number, the more permutations, the longer it takes to run default{100}
6. Select p-value significance threshold, default = {a<=0.05}
7. Submit.
8. When Export button appears, select to write data to disk. you will be prompted to select output directory
