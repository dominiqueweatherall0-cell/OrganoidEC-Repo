# 03_tip_tests.R — Tip module violin (Wilcoxon) and Tip proportion (Fisher)
# Saves: tip_module_score.png, tip_proportion.png

set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(scales) })

fig_dir <- "results"; obj_dir <- file.path("data","processed")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ec_w6_emm2 <- readRDS(file.path(obj_dir, "ec_w6_emm2_final.rds"))
DefaultAssay(ec_w6_emm2) <- "RNA"; ec_w6_emm2 <- JoinLayers(ec_w6_emm2, assay="RNA")

tip_sig_genes <- intersect(c("ESM1","DLL4","APLN","KDR","CXCR4"), rownames(ec_w6_emm2))
ec_w6_emm2@meta.data <- ec_w6_emm2@meta.data[, !grepl("^TipSig", colnames(ec_w6_emm2@meta.data)), drop=FALSE]
ec_w6_emm2 <- AddModuleScore(ec_w6_emm2, features=list(TipSig=tip_sig_genes),
                             name="TipSig", nbin=24, assay="RNA", search=FALSE)
score_col <- grep("^TipSig", colnames(ec_w6_emm2@meta.data), value=TRUE)[1]

tip_cells <- WhichCells(ec_w6_emm2, expression = subtype_cell == "Tip EC")
md_tip <- ec_w6_emm2@meta.data[tip_cells, c("condition", score_col), drop=FALSE]

wil_p <- tryCatch(wilcox.test(
  x = md_tip[md_tip$condition=="Week6Ref", score_col, drop=TRUE],
  y = md_tip[md_tip$condition=="EMM2",    score_col, drop=TRUE],
  alternative="two.sided", exact=FALSE
)$p.value, error=function(e) NA_real_)
lab_wil <- paste0("Wilcoxon p = ", ifelse(is.na(wil_p),"NA", formatC(wil_p, format="e", digits=2)))

p_violin <- ggplot(md_tip, aes(condition, .data[[score_col]], fill=condition)) +
  geom_violin(trim=TRUE, width=0.9, alpha=0.9) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.95, position=position_dodge(0.9)) +
  labs(title="Tip EC module score by condition", x=NULL, y="Tip EC module score") +
  annotate("text", x=1.5, y=max(md_tip[[score_col]])*1.05, label=lab_wil, size=5.2, fontface="bold") +
  theme_minimal(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir, "tip_module_score.png"), p_violin, width=8.4, height=4.6, dpi=300)

is_tip <- ec_w6_emm2$subtype_cell == "Tip EC"
tab <- table(condition = ec_w6_emm2$condition, tip = is_tip)
fisher_p <- tryCatch(fisher.test(tab)$p.value, error=function(e) NA_real_)
lab_fisher <- paste0("Fisher exact p = ", ifelse(is.na(fisher_p),"NA", formatC(fisher_p, format="e", digits=2)))

prop_df <- as.data.frame(tab) |>
  dplyr::mutate(tip = as.logical(tip)) |>
  dplyr::group_by(condition) |>
  dplyr::mutate(prop = Freq / sum(Freq)) |>
  dplyr::ungroup() |>
  dplyr::filter(tip)

p_prop <- ggplot(prop_df, aes(x=condition, y=prop)) +
  geom_col(width=0.7) +
  geom_text(aes(label=scales::percent(prop, accuracy=0.1)), vjust=-0.5, size=4.5) +
  labs(title="Tip EC proportion by condition", x=NULL, y="Proportion of ECs") +
  annotate("text", x=1.5, y=max(prop_df$prop)*1.12, label=lab_fisher, size=5.2, fontface="bold") +
  scale_y_continuous(labels=scales::percent_format(accuracy=1), limits=c(0, max(prop_df$prop)*1.2)) +
  theme_minimal(base_size=14) + theme(legend.position="none")
ggsave(file.path(fig_dir, "tip_proportion.png"), p_prop, width=8.3, height=4.6, dpi=300)
