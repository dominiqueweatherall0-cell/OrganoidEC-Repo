# 02_gate_integrate_subtypes.R — EC gating, robust SCT integration, subtype scoring
# Saves: Endothelial_clusters.png, umap_condition.png, subtype_composition.png

set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(scales) })

fig_dir <- "results"; obj_dir <- file.path("data","processed")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ref6 <- readRDS(file.path(obj_dir, "Week6Ref_qc.rds"))
emm2 <- readRDS(file.path(obj_dir, "EMM2_qc.rds"))

dashify <- function(x){ x <- gsub("_","-", x, fixed=TRUE); x <- gsub("\\s+","", x); make.unique(x) }
gate_ec <- function(o){
  DefaultAssay(o) <- "RNA"
  marker_syms <- unique(gsub("-", "", toupper(c("PECAM1","CD31","CDH5","VE-CADHERIN","KDR","VEGFR2",
                                                "FLT1","ESAM","VWF","PLVAP","ENG","CD34"))))
  normalize_syms <- function(x) gsub("-", "", toupper(x))
  rn <- rownames(o)
  looks_id <- grepl("^(ENSG\\d+|AC\\d+|AL\\d+|AP\\d+|RP\\d+|LINC\\d+\\.?\\d*)$", rn)
  if (any(looks_id) && !is.null(o@misc$feature_map)){
    fmap <- o@misc$feature_map
    id_to_sym <- setNames(fmap$symbol, fmap$feature_id)
    new_syms  <- id_to_sym[o@misc$gene_symbol_map$original]
    ok <- !is.na(new_syms) & nzchar(new_syms)
    rn2 <- rn; rn2[ok] <- dashify(new_syms[ok])
    if (any(duplicated(rn2))){
      mat <- GetAssayData(o, assay="RNA", layer="counts")
      dup_groups <- split(seq_along(rn2), rn2)
      summed <- do.call(rbind, lapply(dup_groups, function(ix) Matrix::colSums(mat[ix, , drop=FALSE])))
      o <- CreateSeuratObject(as(summed,"dgCMatrix"), meta.data=o@meta.data, project=o@project.name)
    } else {
      mat <- GetAssayData(o, assay="RNA", layer="counts"); rownames(mat) <- rn2
      o[["RNA"]]@counts <- mat; o[["RNA"]]@data <- new("dgCMatrix")
    }
  }
  rn_sym_norm <- normalize_syms(rownames(o))
  present <- rn_sym_norm %in% marker_syms
  if (!any(present)) stop("No endothelial markers found after mapping.")
  ec_genes <- rownames(o)[present]
  o <- NormalizeData(o, verbose=FALSE)
  mat <- GetAssayData(o, assay="RNA", layer="counts")
  o$ECscore <- Matrix::colSums(mat[ec_genes, , drop=FALSE] > 0)
  subset(o, subset = ECscore >= 3)
}

ref6_ec <- gate_ec(ref6)
emm2_ec <- gate_ec(emm2)

ref6 <- RenameCells(ref6_ec, add.cell.id = "W6")
emm2 <- RenameCells(emm2_ec, add.cell.id = "E2")
ref6[["percent.mt"]] <- PercentageFeatureSet(ref6, pattern="^MT-")
emm2[["percent.mt"]] <- PercentageFeatureSet(emm2, pattern="^MT-")
objs <- list(ref6, emm2)
objs <- lapply(objs, function(x) SCTransform(x, variable.features.n=3000, vars.to.regress="percent.mt", verbose=FALSE))

safe_integrate_SCT <- function(objs){
  feats <- SelectIntegrationFeatures(objs, nfeatures = 3000)
  objs2 <- PrepSCTIntegration(objs, anchor.features = feats, verbose = FALSE)
  try1 <- try({ anchors <- FindIntegrationAnchors(objs2, normalization.method="SCT", anchor.features=feats, dims=1:30, verbose=FALSE)
                IntegrateData(anchors, normalization.method="SCT", dims=1:30, verbose=FALSE) }, silent=TRUE)
  if (!inherits(try1,"try-error")) return(try1)
  try2 <- try({ anchors <- FindIntegrationAnchors(objs2, normalization.method="SCT", anchor.features=feats, dims=1:30, k.filter=20, k.anchor=10, verbose=FALSE)
                IntegrateData(anchors, normalization.method="SCT", dims=1:30, k.weight=10, verbose=FALSE) }, silent=TRUE)
  if (!inherits(try2,"try-error")) return(try2)
  try3 <- try({ anchors <- FindIntegrationAnchors(objs2, normalization.method="SCT", anchor.features=feats, dims=1:20, k.filter=10, k.anchor=5, verbose=FALSE)
                IntegrateData(anchors, normalization.method="SCT", dims=1:20, k.weight=5, verbose=FALSE) }, silent=TRUE)
  if (!inherits(try3,"try-error")) return(try3)
  NULL
}
integrated <- safe_integrate_SCT(objs)
if (is.null(integrated)){
  merged <- merge(ref6, y = emm2, add.cell.ids = c("W6","E2"))
  DefaultAssay(merged) <- "RNA"; merged <- NormalizeData(merged, verbose=FALSE)
  merged <- FindVariableFeatures(merged, nfeatures=2000, verbose=FALSE)
  merged <- ScaleData(merged, verbose=FALSE); merged <- RunPCA(merged, npcs=30, verbose=FALSE)
  merged <- FindNeighbors(merged, dims=1:30, verbose=FALSE); merged <- FindClusters(merged, resolution=0.4, verbose=FALSE)
  merged <- RunUMAP(merged, dims=1:30, verbose=FALSE); merged$batch_mode <- "merged-noBatch"
  integrated <- merged
} else {
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose=FALSE); integrated <- RunPCA(integrated, npcs=30, verbose=FALSE)
  integrated <- FindNeighbors(integrated, dims=1:30, verbose=FALSE); integrated <- FindClusters(integrated, resolution=0.4, verbose=FALSE)
  integrated <- RunUMAP(integrated, dims=1:30, verbose=FALSE); integrated$batch_mode <- "SCT-integrated"
}
pref <- sub("^([A-Za-z0-9]+)_.*", "\\1", colnames(integrated))
integrated$condition <- factor(ifelse(pref=="E2","EMM2","Week6Ref"), levels=c("Week6Ref","EMM2"))

annotate_subtypes_cellwise <- function(obj){
  DefaultAssay(obj) <- "RNA"
  sig <- list(
    "Tip EC"                  = c("ESM1","DLL4","APLN","KDR","CXCR4"),
    "Arterial (stalk/Notch)"  = c("EFNB2","HEY1","JAG1","NOTCH1","GJA5"),
    "Venous EC"               = c("NR2F2","EPHB4","ACKR1"),
    "Activated Venule EC"     = c("SELE","ICAM1","VCAM1","ACKR1"),
    "Fenestrated EC"          = c("PLVAP","CA4","AQP1","EHD3"),
    "Continuous Capillary EC" = c("CLDN5","CDH5","TJP1","ICAM2"),
    "Lymphatic EC"            = c("PROX1","PDPN","LYVE1","FLT4"),
    "Endocardial EC"          = c("NFATC1","NPR3","KDR","HEY2")
  )
  sig <- lapply(sig, function(g) intersect(g, rownames(obj))); sig <- sig[vapply(sig, length, 1L) > 0]
  obj <- AddModuleScore(obj, features = sig, name = "SIG", assay="RNA", nbin=24)
  sig_cols <- grep("^SIG", colnames(obj@meta.data), value=TRUE)
  colnames(obj@meta.data)[match(sig_cols, colnames(obj@meta.data))] <- paste0(names(sig), "_score")
  md <- obj@meta.data; score_cols <- grep("_score$", colnames(md), value=TRUE)
  scmat <- as.matrix(md[, score_cols, drop=FALSE]); best_idx <- max.col(scmat, ties.method="first")
  best_lab <- sub("_score$","", score_cols)[best_idx]; maxval <- scmat[cbind(seq_len(nrow(scmat)), best_idx)]
  best_lab[maxval < quantile(maxval, 0.05, na.rm=TRUE)] <- "Low-signal"
  obj$subtype_cell <- best_lab; Idents(obj) <- "subtype_cell"; obj
}
ec_w6_emm2 <- annotate_subtypes_cellwise(integrated)
saveRDS(ec_w6_emm2, file = file.path(obj_dir, "ec_w6_emm2_final.rds"))

theme_pub <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face="bold"))

p_umap_subtype <- DimPlot(ec_w6_emm2, group.by="subtype_cell", label=FALSE) +
  labs(title=paste0("Endothelial subtypes (", unique(ec_w6_emm2$batch_mode), ")"), x="umap_1", y="umap_2") +
  theme_pub + theme(legend.title = element_blank())
ggsave(file.path(fig_dir, "Endothelial_clusters.png"), p_umap_subtype, width=11.5, height=4.8, dpi=300)

p_umap_cond <- DimPlot(ec_w6_emm2, group.by="condition", label=FALSE) +
  labs(title="UMAP: Batch correction check (condition)", x="umap_1", y="umap_2") + theme_pub
ggsave(file.path(fig_dir, "umap_condition.png"), p_umap_cond, width=8.4, height=4.6, dpi=300)

meta_tbl <- tibble::as_tibble(ec_w6_emm2@meta.data[, c("condition","subtype_cell")]) |>
  dplyr::rename(subtype = subtype_cell)
prop_df <- meta_tbl |>
  dplyr::count(condition, subtype, name="n") |>
  dplyr::group_by(condition) |>
  dplyr::mutate(prop = n/sum(n)) |> dplyr::ungroup()
p_comp <- ggplot(prop_df, aes(condition, prop, fill=subtype)) +
  geom_col(position="fill", width=0.8, color="white", linewidth=0.2) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  labs(title=paste0("EC subtype composition (", unique(ec_w6_emm2$batch_mode), "): Week6Ref vs EMM2"),
       x=NULL, y="Proportion of ECs", fill="EC subtype") + theme_pub
ggsave(file.path(fig_dir, "subtype_composition.png"), p_comp, width=12, height=6, dpi=300)
