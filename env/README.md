# env/ — computational environment records

`run_all.R` and each numbered script write a `sessionInfo()` dump here
(`sessionInfo.txt`, `sessionInfo_01.txt`, …) so that the exact package
versions behind every result table can be recovered.

Manuscript-reported versions: **R 4.5.3**, Seurat 4.3.0, GENIE3 1.20.0,
monocle3 1.3.1, pySCENIC 0.12.1.

Hardware note: the analysis was run on an institutional compute node
(8 cores, 64 GB RAM). GENIE3 with 1,000 trees is the dominant cost.
