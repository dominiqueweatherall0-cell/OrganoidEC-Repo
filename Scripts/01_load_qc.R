## ---- setup ----
set.seed(1234)
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(dplyr); library(scales)
})

# Relative repo paths (adjust if needed)
ref_meta   <- "data/raw/reference/GSE106118_barcode_information.txt.gz"
ref_counts <- "data/raw/reference/GSE106118_UMI_count_merge.txt.gz"
emm2_dir   <- "data/raw/EMM2"  # contains matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz

fig_dir <- "results/figures"; dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
obj_dir <- "data/processed"; dir.create(obj_dir,  recursive = TRUE, showWarnings = FALSE)

## ---- helpers ----
dashify <- function(x){x <- gsub("_","-",x, fixed=TRUE); x <- gsub("\\s+","",x); make.unique(x)}
CreateSeuratFromCounts <- function(counts_df, project="Sample", meta=NULL, min.cells=3, min.features=200){
  stopifnot(ncol(counts_df) >= 2)
  genes_orig <- counts_df[[1]]
  mat <- as.matrix(data.frame(row.names=genes_orig, counts_df[,-1, drop=FALSE]))
  rownames(mat) <- dashify(rownames(mat))
  obj <- CreateSeuratObject(Matrix(mat, sparse=TRUE), project=project, meta.data=meta,
                            min.cells=min.cells, min.features=min.features)
  obj
}
CreateSeuratFrom10x <- function(dir, project=basename(dir), min.cells=3, min.features=200){
  raw <- ReadMtx(mtx=file.path(dir,"matrix.mtx.gz"),
                 features=file.path(dir,"features.tsv.gz"),
                 cells=file.path(dir,"barcodes.tsv.gz"))
  rownames(raw) <- dashify(rownames(raw))
  CreateSeuratObject(raw, project=project, min.cells=min.cells, min.features=min.features)
}

## ---- Week6Ref (GSE metadata filter HE6W) ----
metadata <- fread(ref_meta); counts <- fread(ref_counts)
get_col <- function(dt, nm){ idx <- which(tolower(colnames(dt)) == tolower(nm)); if(!length(idx)) return(NULL); dt[[idx]] }
processed_col <- get_col(metadata, "processed_file")
cell_col      <- get_col(metadata, "cell")
wk6_meta <- metadata[FALSE]
if(!is.null(processed_col)) wk6_meta <- rbind(wk6_meta, metadata[grepl("HE6W", processed_col, ignore.case=TRUE)])
if(!is.null(cell_col))      wk6_meta <- rbind(wk6_meta, metadata[grepl("HE6W", cell_col,      ignore.case=TRUE)])
wk6_meta <- unique(wk6_meta); if(nrow(wk6_meta) == 0) stop("No Week-6 rows found (HE6W).")

gene_col <- intersect(names(counts), c("gene","Gene","GENE","GENEID","GeneID","symbol","SYMBOL","V1"))
if(!length(gene_col)) gene_col <- names(counts)[1]
setnames(counts, gene_col[1], "gene")

wk6_cells <- unique(as.character(get_col(wk6_meta, "cell"))); wk6_cells <- wk6_cells[!is.na(wk6_cells)]
keep_cols <- intersect(names(counts), wk6_cells); if(!length(keep_cols)) stop("No matching Week-6 cells in counts.")
week6_counts <- counts[, c("gene", keep_cols), with = FALSE]

# quick count-based QC to drop extreme low-depth cells prior to Seurat
mat_dt  <- as.matrix(week6_counts[, -1, with = FALSE])
libsizes <- colSums(mat_dt); meanZero <- colMeans(mat_dt == 0)
keep <- (libsizes > 2000) & (meanZero < 0.98)
week6_counts_qc <- week6_counts[, c("gene", names(keep)[keep]), with = FALSE]

ref6 <- CreateSeuratFromCounts(week6_counts_qc, project="Week6Ref")
emm2 <- CreateSeuratFrom10x(emm2_dir, project="EMM2")

## ---- QC metrics + filtering (only ref6 & emm2) ----
compute_qc_metrics <- function(obj){
  DefaultAssay(obj) <- "RNA"
  feats <- rownames(obj)
  mt    <- grep("^MT-",  feats, value = TRUE)
  ribo  <- grep("^RP[SL]\\d", feats, value = TRUE)
  hb    <- grep("^HB[ABDEZ]\\d?", feats, value = TRUE)
  obj[["percent.mt"]]   <- if(length(mt))   PercentageFeatureSet(obj, mt)   else rep(NA_real_, ncol(obj))
  obj[["percent.ribo"]] <- if(length(ribo)) PercentageFeatureSet(obj, ribo) else rep(0, ncol(obj))
  obj[["percent.hb"]]   <- if(length(hb))   PercentageFeatureSet(obj, hb)   else rep(0, ncol(obj))
  obj
}
filter_qc <- function(obj, min_features=200, min_counts=500, max_mt=20, max_ribo=60, max_hb=5){
  keep <- obj$nFeature_RNA >= min_features & obj$nCount_RNA >= min_counts &
          (is.na(obj$percent.mt) | obj$percent.mt < max_mt) &
          obj$percent.ribo < max_ribo & obj$percent.hb < max_hb
  subset(obj, cells = colnames(obj)[keep])
}
ref6 <- compute_qc_metrics(ref6); emm2 <- compute_qc_metrics(emm2)
ref6 <- filter_qc(ref6);        emm2 <- filter_qc(emm2)

# QC violins (one png per dataset)
plot_qc <- function(o, nm){
  feats <- c("nFeature_RNA","nCount_RNA","percent.ribo","percent.hb")
  if(!all(is.na(o$percent.mt))) feats <- c(feats, "percent.mt")
  p <- VlnPlot(o, features = feats, ncol = length(feats), pt.size = 0.2) + ggtitle(nm)
  ggsave(file.path(fig_dir, paste0("qc_violin_", nm, ".png")), p, width = 7.5, height = 4.5, dpi = 300)
}
plot_qc(ref6, "Week6Ref"); plot_qc(emm2, "EMM2")

# save for next steps
saveRDS(ref6, file.path(obj_dir,"ref6_qc.rds"))
saveRDS(emm2, file.path(obj_dir,"emm2_qc.rds"))
