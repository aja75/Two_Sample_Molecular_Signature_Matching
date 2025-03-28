# Two_Sample_Molecular_Signature_Matching
Compare 2 Differential Expression Profiles in terms of their DEG Overlap % and Concordance validated by permutation
---

## How To Use:
1. From Differential Expression (DE) files (use DEseq2, EdgeR, CuffDiff2, etc) select 2 different profiles
2. Specify if DE is formatted as log2FC or Z-Score 'column spelling must match these exactly'
3. Specify thresholds for the binary mapping default = {log2FC = 1, Z-Score = 2}
4. Enter the gene column name, ex. ENSEMBL_ID, GENE_SYMBOL, etc 'must be spelled exactly'
5. Enter permutation number, the more permutations, the longer it takes to run default{100}
6. Select p-value significance threshold, default = {a<=0.05}
7. Submit.
8. When Export button appears, select to write data to disk. you will be prompted to select output directory
