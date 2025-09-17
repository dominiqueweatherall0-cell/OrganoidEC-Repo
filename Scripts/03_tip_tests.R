set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(scales) })

fig_dir <- "results/figures"; obj_dir <- "data/processed"
ec_w6_emm2 <- readRDS(file.path(obj_dir,"ec_w6_emm2_final.rds"))
DefaultAssay(ec_w6_emm2) <- "RNA"; ec_w6_emm2$condition <- factor(ec_w6_emm2$condition, levels=c("Week6Ref","EMM2"))

# Tip module score (ESM1,DLL4,APLN,KDR,CXCR4)
tip_genes <- intersect(c("ESM1","DLL4","APLN","KDR","CXCR4"), rownames(ec_w6_emm2))
ec_w6_emm2@meta.data <- ec_w6_emm2@meta.data[, !grepl("^TipSig", colnames(ec_w6_emm2@meta.data)), drop=FALSE]
ec_w6_emm2 <- AddModuleScore(ec_w6_emm2, features=list(TipSig=tip_genes),
                             name="TipSig", nbin=24, assay="RNA", search=FALSE)
score_col <- grep("^TipSig", colnames(ec_w6_emm2@meta.data), value=TRUE)[1]
tip_cells <- WhichCells(ec_w6_emm2, expression = subtype_cell == "Tip EC")
md_tip <- ec_w6_emm2@meta.data[tip_cells, c("condition", score_col), drop=FALSE]

pval <- tryCatch(wilcox.test(md_tip[[score_col]] ~ md_tip$condition, exact=FALSE)$p.value, error=function(e) NA_real_)
lab  <- paste0("Wilcoxon p = ", ifelse(is.na(pval), "NA", formatC(pval, format="e", digits=2)))

p_tip <- ggplot(md_tip, aes(condition, .data[[score_col]], fill=condition)) +
  geom_violin(trim=TRUE, width=0.9, alpha=0.9) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.95) +
  labs(title="Tip EC module score by condition", x=NULL, y="Tip EC module score") +
  annotate("text", x=1.5, y=max(md_tip[[score_col]])*1.05, label=lab, size=5, fontface="bold") +
  scale_fill_manual(values=c("Week6Ref"="#f47f7f","EMM2"="#39b7b2")) +
  theme_minimal(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir,"tip_module_violin.png"), p_tip, width=7.5, height=5.5, dpi=300)

# Tip proportion (Fisher exact)
is_tip <- ec_w6_emm2$subtype_cell == "Tip EC"
tab <- table(condition = ec_w6_emm2$condition, tip = is_tip)
pf  <- tryCatch(fisher.test(tab)$p.value, error=function(e) NA_real_)
labf <- paste0("Fisher exact p = ", ifelse(is.na(pf), "NA", formatC(pf, format="e", digits=2)))
prop_df <- as.data.frame(tab) %>% mutate(tip = as.logical(tip)) %>%
  group_by(condition) %>% mutate(prop = Freq/sum(Freq)) %>% ungroup() %>% filter(tip)

p_prop <- ggplot(prop_df, aes(condition, prop, fill=condition)) +
  geom_col(width=0.7, color="white") +
  geom_text(aes(label = percent(prop, accuracy=0.1)), vjust=-0.5, size=4.5) +
  labs(title="Tip EC proportion by condition", x=NULL, y="Proportion of ECs") +
  annotate("text", x=1.5, y=max(prop_df$prop)*1.12, label=labf, size=5, fontface="bold") +
  scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0, max(prop_df$prop)*1.2)) +
  scale_fill_manual(values=c("Week6Ref"="#f47f7f","EMM2"="#39b7b2")) +
  theme_minimal(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir,"tip_proportion_bar.png"), p_prop, width=7.0, height=5.2, dpi=300)
