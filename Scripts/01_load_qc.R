# 01_load_qc.R — relative paths + auto-discovery; saves QC_organoid.png and EMM2_QC.png

set.seed(1234)
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(dplyr); library(cowplot)
})

# -------------------------
# Repo-relative paths
# -------------------------
raw_dir       <- file.path("data","raw")
ref_dir       <- file.path(raw_dir,"reference")
processed_dir <- file.path("data","processed")
fig_dir       <- "results"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Helpers
# -------------------------
find_one <- function(pattern, start = "."){
  hits <- list.files(start, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(hits) < 1) stop("Could not find file matching: ", pattern)
  hits[[1]]
}
dashify <- function(x){ x <- gsub("_","-", x, fixed=TRUE); x <- gsub("\\s+","", x); make.unique(x) }

CreateSeuratFromCounts <- function(counts_df, project="Sample", meta=NULL, min.cells=3, min.features=200){
  stopifnot(ncol(counts_df) >= 2)
  genes_orig <- counts_df[[1]]
  mat <- as.matrix(data.frame(row.names=genes_orig, counts_df[,-1, drop=FALSE]))
  genes_new <- dashify(rownames(mat)); rownames(mat) <- genes_new
  gene_map <- data.frame(original=genes_orig, seurat=genes_new, stringsAsFactors=FALSE)
  mat <- Matrix(mat, sparse=TRUE)
  obj <- CreateSeuratObject(mat, project=project, meta.data=meta, min.cells=min.cells, min.features=min.features)
  obj@misc$gene_symbol_map <- gene_map
  obj
}

CreateSeuratFrom10x <- function(dir, project=basename(dir), min.cells=3, min.features=200){
  mtx   <- file.path(dir, "matrix.mtx.gz")
  feats <- file.path(dir, "features.tsv.gz")
  bc    <- file.path(dir, "barcodes.tsv.gz")
  if (!all(file.exists(c(mtx,feats,bc))))
    stop("10x directory missing required files: ", dir)
  feat_df <- tryCatch({ data.table::fread(feats, header = FALSE) }, error=function(e) NULL)
  raw <- ReadMtx(mtx=mtx, features=feats, cells=bc)
  genes_orig <- rownames(raw); genes_new <- dashify(genes_orig); rownames(raw) <- genes_new
  obj <- CreateSeuratObject(raw, project=project, min.cells=min.cells, min.features=min.features)
  obj$condition <- project
  obj@misc$gene_symbol_map <- data.frame(original=genes_orig, seurat=genes_new, stringsAsFactors=FALSE)
  if (!is.null(feat_df) && ncol(feat_df) >= 2){
    colnames(feat_df)[1:2] <- c("feature_id","symbol"); obj@misc$feature_map <- feat_df
  } else obj@misc$feature_map <- NULL
  obj
}

# Env overrides allowed; otherwise auto-discover under data/raw/reference/
meta_path   <- Sys.getenv("ORG_REF_META",
                  if (dir.exists(ref_dir)) find_one("^GSE106118_barcode_information\\.txt\\.gz$", ref_dir) else find_one("^GSE106118_barcode_information\\.txt\\.gz$", raw_dir))
counts_path <- Sys.getenv("ORG_REF_COUNTS",
                  if (dir.exists(ref_dir)) find_one("^GSE106118_UMI_count_merge\\.txt\\.gz$", ref_dir) else find_one("^GSE106118_UMI_count_merge\\.txt\\.gz$", raw_dir))

# 10x dirs (env override first, else default to data/raw/<name>)
get_10x_dir <- function(name){
  envv <- Sys.getenv(paste0("ORG_10X_", toupper(name)))
  if (nzchar(envv)) return(envv)
  file.path(raw_dir, name)
}
dir_control <- get_10x_dir("Control")
dir_emm1    <- get_10x_dir("EMM1")
dir_emm2    <- get_10x_dir("EMM2")
dir_mm      <- get_10x_dir("MM")

# -------------------------
# Load reference metadata + counts
# -------------------------
metadata <- fread(meta_path)
counts   <- fread(counts_path)

get_col <- function(dt, nm){ idx <- which(tolower(colnames(dt)) == tolower(nm)); if (!length(idx)) return(NULL); dt[[idx]] }
processed_col <- get_col(metadata, "processed_file")
cell_col      <- get_col(metadata, "cell")

wk6_meta <- metadata[FALSE]
if (!is.null(processed_col)) wk6_meta <- rbind(wk6_meta, metadata[grepl("HE6W", processed_col, ignore.case=TRUE)])
if (!is.null(cell_col))      wk6_meta <- rbind(wk6_meta, metadata[grepl("HE6W", cell_col,      ignore.case=TRUE)])
wk6_meta <- unique(wk6_meta)
if (nrow(wk6_meta) == 0) stop("No Week-6 rows found (looked for 'HE6W').")

gene_col <- intersect(names(counts), c("gene","Gene","GENE","GENEID","GeneID","symbol","SYMBOL","V1"))
if (!length(gene_col)) gene_col <- names(counts)[1]
data.table::setnames(counts, gene_col[1], "gene")

wk6_cells <- unique(as.character(get_col(wk6_meta, "cell")))
wk6_cells <- wk6_cells[!is.na(wk6_cells)]
candidate_cols <- intersect(names(counts), wk6_cells)
if (!length(candidate_cols)) stop("No matching Week-6 cells among count columns.")
week6_counts <- counts[, c("gene", candidate_cols), with = FALSE]

# -------------------------
# Quick QC + objects
# -------------------------
mat_dt  <- as.matrix(week6_counts[, -1, with = FALSE])
libsizes <- colSums(mat_dt); meanZero <- colMeans(mat_dt == 0)
keep <- (libsizes > 2000) & (meanZero < 0.98)
week6_counts_qc <- week6_counts[, c("gene", names(keep)[keep]), with = FALSE]

ref6 <- CreateSeuratFromCounts(week6_counts_qc, project="Week6Ref")
ctrl <- CreateSeuratFrom10x(dir_control, "Control")
emm1 <- CreateSeuratFrom10x(dir_emm1,   "EMM1")
emm2 <- CreateSeuratFrom10x(dir_emm2,   "EMM2")
mm   <- CreateSeuratFrom10x(dir_mm,     "MM")

compute_qc_metrics <- function(obj){
  DefaultAssay(obj) <- "RNA"; obj <- JoinLayers(obj, assay="RNA")
  feats <- rownames(obj); ncell <- ncol(obj)
  mt_genes   <- grep("^MT-",  feats, value=TRUE)
  ribo_genes <- grep("^RP[SL]\\d", feats, value=TRUE)
  hb_genes   <- grep("^HB[ABDEZ]\\d?", feats, value=TRUE)
  obj[["percent.mt"]]   <- if (length(mt_genes))   PercentageFeatureSet(obj, mt_genes) else rep(NA_real_, ncell)
  obj[["percent.ribo"]] <- if (length(ribo_genes)) PercentageFeatureSet(obj, ribo_genes) else rep(0, ncell)
  obj[["percent.hb"]]   <- if (length(hb_genes))   PercentageFeatureSet(obj, hb_genes)   else rep(0, ncell)
  obj
}
filter_qc <- function(obj, min_features=200, min_counts=500, max_mt=20, max_ribo=60, max_hb=5){
  keep <- obj$nFeature_RNA >= min_features & obj$nCount_RNA >= min_counts &
    (is.na(obj$percent.mt) | obj$percent.mt < max_mt) & obj$percent.ribo < max_ribo & obj$percent.hb < max_hb
  subset(obj, cells = colnames(obj)[keep])
}

qc_list <- list(ref6=ref6, ctrl=ctrl, emm1=emm1, emm2=emm2, mm=mm)
qc_list <- lapply(qc_list, compute_qc_metrics)

# Save QC violins for Week6Ref and EMM2
plot_qc <- function(obj, nm, file_out){
  feats <- c("nFeature_RNA","nCount_RNA")
  if (!all(is.na(obj$percent.mt))) feats <- c(feats,"percent.mt")
  feats <- c(feats,"percent.ribo")
  p <- VlnPlot(obj, features = feats, ncol = length(feats), pt.size = 0.2) + ggtitle(nm)
  ggsave(file.path(fig_dir, file_out), p, width = 8.5, height = 4.8, dpi = 300)
}
plot_qc(qc_list$ref6, "Week6Ref", "QC_organoid.png")
plot_qc(qc_list$emm2, "EMM2",     "EMM2_QC.png")

qc_list <- lapply(qc_list, filter_qc)
ref6 <- qc_list$ref6; emm2 <- qc_list$emm2
saveRDS(ref6, file.path(processed_dir, "Week6Ref_qc.rds"))
saveRDS(emm2, file.path(processed_dir, "EMM2_qc.rds"))
