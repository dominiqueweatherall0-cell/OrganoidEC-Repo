set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

fig_dir <- "results/figures"; obj_dir <- "data/processed"
ec_w6_emm2 <- readRDS(file.path(obj_dir,"ec_w6_emm2_final.rds"))
DefaultAssay(ec_w6_emm2) <- "RNA"

tip_cells <- WhichCells(ec_w6_emm2, expression = subtype_cell == "Tip EC")
tip_obj <- subset(ec_w6_emm2, cells = tip_cells)
Idents(tip_obj) <- "condition"

mech_genes <- intersect(c("KLF4","CAMK2D","PIEZO1","ESM1"), rownames(tip_obj))
de <- FindMarkers(tip_obj, ident.1="Week6Ref", ident.2="EMM2",
                  features=mech_genes, logfc.threshold=0, min.pct=0, assay="RNA")
de <- de %>% mutate(gene = rownames(.), padj = p.adjust(p_val, method="BH")) %>%
  arrange(desc(avg_log2FC))

theme_pub <- theme_minimal(base_size=14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face="bold"))
p_de <- ggplot(de, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC)) +
  geom_hline(yintercept = 0, linewidth=0.5, color="grey50") +
  geom_col(width=0.7, fill="#5a7bd8") +
  geom_text(aes(label=paste0("adj p=", formatC(padj, format="e", digits=2)),
                y=avg_log2FC + 0.06*sign(avg_log2FC))),
            vjust=ifelse(de$avg_log2FC>=0,-0.3,1.3), size=4) +
  coord_flip() +
  labs(title="Tip ECs: mechanoreceptor-panel DE (Week6Ref – EMM2)",
       x=NULL, y="avg_log2FC (positive = higher in Week6Ref)") +
  theme_pub
ggsave(file.path(fig_dir,"tip_mechanoreceptor_DE.png"), p_de, width=8.5, height=5.5, dpi=300)
