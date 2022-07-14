# GSE109564_replication_Seurat4
This R script replicates the data processing done by Wu et al (2018).


Download GSE109564 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109564) into the working directory.
Run the R script to replicate the data processing done by Wu et al (https://pubmed.ncbi.nlm.nih.gov/29980650/).


Because the unsupervised clustering and cell type identification were originally done with Seurat V2, this R script delivers the same unsupervised clustering and cell type identification using Seurat V4 (https://satijalab.org/seurat/).


R version: 4.1.3
Bioconductor version: 3.14
All library packages are up-to-date for said versions.
