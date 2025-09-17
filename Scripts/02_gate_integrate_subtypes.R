set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(scales) })

fig_dir <- "results/figures"; dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
obj_dir <- "data/processed"; dir.create(obj_dir,  recursive = TRUE, showWarnings = FALSE)

ref6 <- readRDS(file.path(obj_dir,"ref6_qc.rds"))
emm2 <- readRDS(file.path(obj_dir,"emm2_qc.rds"))

## ---- EC gating (symbol-only) ----
gate_ec <- function(o){
  DefaultAssay(o) <- "RNA"; o <- NormalizeData(o, verbose=FALSE)
  ec_genes <- intersect(c("PECAM1","CDH5","KDR","FLT1","ESAM","VWF","PLVAP","ENG","CD34"), rownames(o))
  mat <- GetAssayData(o, assay="RNA", slot="counts")
  o$ECscore <- Matrix::colSums(mat[ec_genes, , drop=FALSE] > 0)
  subset(o, subset = ECscore >= 3)
}
ref6_ec <- gate_ec(ref6); emm2_ec <- gate_ec(emm2)

## ---- SCT integration (ref6 + emm2) ----
ref6_ec <- RenameCells(ref6_ec, add.cell.id = "W6")
emm2_ec <- RenameCells(emm2_ec, add.cell.id = "E2")
objs <- list(ref6_ec, emm2_ec)
objs <- lapply(objs, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern="^MT-")
  SCTransform(x, variable.features.n=3000, vars.to.regress="percent.mt", verbose=FALSE)
})
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = objs, normalization.method="SCT",
                                  anchor.features=features, dims=1:30, verbose=FALSE)
integrated <- IntegrateData(anchorset = anchors, normalization.method="SCT", dims=1:30, verbose=FALSE)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose=FALSE)
integrated <- RunPCA(integrated, npcs=30, verbose=FALSE)
integrated <- FindNeighbors(integrated, dims=1:30, verbose=FALSE)
integrated <- FindClusters(integrated, resolution=0.4, verbose=FALSE)
integrated <- RunUMAP(integrated, dims=1:30, verbose=FALSE)
integrated$batch_mode <- "SCT-integrated"

pref <- sub("^([A-Za-z0-9]+)_.*","\\1", colnames(integrated))
integrated$condition <- factor(ifelse(pref=="E2","EMM2","Week6Ref"), levels=c("Week6Ref","EMM2"))

## ---- Subtype annotation (cellwise winner-takes-all) ----
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
  sig <- lapply(sig, function(g) intersect(g, rownames(obj)))
  obj <- AddModuleScore(obj, features = sig, name = "SIG", assay = "RNA", nbin = 24)
  score_cols <- grep("^SIG", colnames(obj@meta.data), value = TRUE)
  colnames(obj@meta.data)[match(score_cols, colnames(obj@meta.data))] <- paste0(names(sig), "_score")
  scmat <- as.matrix(obj@meta.data[, grep("_score$", colnames(obj@meta.data)), drop=FALSE])
  best_idx <- max.col(scmat, ties.method = "first")
  labs <- sub("_score$","", colnames(scmat))[best_idx]
  maxval <- scmat[cbind(seq_len(nrow(scmat)), best_idx)]
  labs[maxval < quantile(maxval, 0.05, na.rm=TRUE)] <- "Low-signal"
  obj$subtype_cell <- labs
  obj
}
ec_w6_emm2 <- annotate_subtypes_cellwise(integrated)

## ---- Plots: UMAP by condition + composition bar ----
theme_pub <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.major.y = element_line(linewidth = 0.2),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p_cond <- DimPlot(ec_w6_emm2, group.by="condition") +
  labs(title = "UMAP: Batch correction check (condition)", x="umap_1", y="umap_2") +
  theme_pub + theme(legend.title = element_blank())
ggsave(file.path(fig_dir,"umap_mix_condition.png"), p_cond, width=7.5, height=6.0, dpi=300)

meta_tbl <- ec_w6_emm2@meta.data[, c("condition","subtype_cell"), drop=FALSE] %>% as_tibble()
prop_df <- meta_tbl %>% count(condition, subtype_cell, name="n") %>% group_by(condition) %>%
  mutate(prop = n/sum(n)) %>% ungroup() %>% rename(subtype = subtype_cell)

p_comp <- ggplot(prop_df, aes(condition, prop, fill=subtype)) +
  geom_col(position="fill", width=0.8, color="white", linewidth=0.2) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits=c(0,1)) +
  labs(title = "EC subtype composition (SCT-integrated): Week6Ref vs EMM2",
       x=NULL, y="Proportion of ECs", fill="EC subtype") +
  theme_pub + guides(fill = guide_legend(nrow=2, byrow=TRUE))
ggsave(file.path(fig_dir,"ec_subtype_composition.png"), p_comp, width=10, height=6, dpi=300)

saveRDS(ec_w6_emm2, file.path(obj_dir,"ec_w6_emm2_final.rds"))
