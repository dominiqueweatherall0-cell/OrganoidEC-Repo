This project uses:
- Week-6 fetal heart reference (GSE106118; GEO)
- Organoid endothelial cells (EMM2; 10x Genomics outputs)

Raw data are not included due to size, but can be downloaded from:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106118
- [Lab-provided EMM2 dataset]

Place downloaded files into `data/raw/` with the following structure:
- `data/raw/reference/` → GSE106118_barcode_information.txt.gz, GSE106118_UMI_count_merge.txt.gz
- `data/raw/EMM2/` → matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
