# 04_mech_panel_de.R — mechanoreception panel DE within Tip ECs
# Saves: mechno_trnasuction_DE.png

set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

fig_dir <- "results"; obj_dir <- file.path("data","processed")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ec_w6_emm2 <- readRDS(file.path(obj_dir, "ec_w6_emm2_final.rds"))
DefaultAssay(ec_w6_emm2) <- "RNA"; ec_w6_emm2 <- JoinLayers(ec_w6_emm2, assay="RNA")

tip_cells <- WhichCells(ec_w6_emm2, expression = subtype_cell == "Tip EC")
tip_obj   <- subset(ec_w6_emm2, cells = tip_cells)
Idents(tip_obj) <- "condition"

mech_genes <- intersect(c("PIEZO1","KLF2","KLF4","MAPK7","MAP3K3","CAMK2D","PLVAP","ESM1"),
                        rownames(tip_obj))

if (length(mech_genes) > 0 && length(tip_cells) > 0){
  de <- FindMarkers(
    tip_obj, ident.1="Week6Ref", ident.2="EMM2",
    features=mech_genes, logfc.threshold=0, min.pct=0, assay="RNA", test.use="wilcox"
  ) %>% dplyr::mutate(gene = rownames(.), adj_p = p.adjust(p_val, method="BH"))
  de$gene <- factor(de$gene, levels = rev(mech_genes))

  p_bar <- ggplot(de, aes(x = avg_log2FC, y = gene)) +
    geom_col() +
    geom_text(aes(label = paste0("adj p=", formatC(adj_p, format="e", digits=2))),
              hjust = ifelse(de$avg_log2FC > 0, -0.02, 1.02), size = 3.2) +
    labs(title = "Tip ECs: mechanoreceptor-panel DE (Week6Ref – EMM2)",
         x = "avg_log2FC (positive = higher in Week6Ref)", y = NULL)
  ggsave(file.path(fig_dir, "mechno_trnasuction_DE.png"), p_bar, width=8.4, height=4.9, dpi=300)
} else {
  message("Mechanoreceptor genes or Tip cells not found; skipping DE plot.")
}
