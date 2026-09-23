# data/ — inputs (raw data are NOT redistributed)

## `data/raw/` (git-ignored)

Place the unpacked 10x-format count matrices here:

```
data/raw/
├── GSE140228/          # discovery: CD45+ immune cells from HCC patients (5 sites)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── GSE149614/          # validation: 10 HCC patients, 4 sites
│   └── ...
└── TF_names_v_1.01.txt # human TF list (Lambert et al. 2018)
```

### Download links

| Item | Source |
|---|---|
| GSE140228 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140228 |
| GSE149614 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614 |
| Human TF list | https://humantfs.ccbr.utoronto.ca/ |

Citations for the datasets:
* Zhang Q., He Y., Luo N., *et al.* Landscape and dynamics of single immune cells in hepatocellular carcinoma. *Cell* **2019**, 179, 829–845.
* Lu Y., Yang A., Quan C., *et al.* A single-cell atlas of the multicellular ecosystem of primary and metastatic hepatocellular carcinoma. *Nat. Commun.* **2022**, 13, 4077.

## `data/processed/` (git-ignored)

Intermediate `.rds` objects written by the pipeline (`*_myeloid.rds`,
`grn_edges.rds`, `grn_adjacency.rds`, `rwr_stationary.rds`,
`permutation_scores.rds`, `pseudotime.rds`, `candidate_tfs.txt`,
`grn_target_genes.txt`).
