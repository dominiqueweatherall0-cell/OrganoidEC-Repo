# OrganoidEC-Repo
This repository presents an independent single-cell RNA-seq analysis of cardiac organoid endothelial cells (Aguirre lab, dataset GSE218582) benchmarked against a fetal heart reference (GSE106118). The analysis was motivated by Aguirre’s observation of reduced vascularization in EMM2/1 cardiac organoids, prompting an independent investigation into endothelial cell subtypes and mechanotransduction pathways.

Objectives:

- Independently process and analyze large-scale scRNA-seq data.

- Compare endothelial subtypes between fetal reference and cardiac organoids.

- Identify transcriptional pathways linked to tip endothelial cell (EC) maturation, with a focus on mechanotransduction.

Workflow:

- Quality Control – filtering and visualization of single-cell counts.

- Integration & Batch Correction – aligning fetal and organoid datasets.

- Subtype Annotation – classification of endothelial subtypes (tip, venous, fenestrated, endocardial, etc.).

- Tip EC Analysis – module scoring and statistical testing of tip ECs across conditions.

- Mechanotransduction Panel – differential expression analysis revealing a candidate shear stress–responsive pathway that may contribute to impaired tip EC maturation.

Repository Structure

Scripts/ – modular R scripts:

01_load_qc.R – load datasets and perform QC

02_gate_integrate_subtypes.R – endothelial gating, integration, subtype annotation

03_tip_tests.R – tip EC module scoring and proportion tests

04_mech_panel_de.R – mechanotransduction panel differential expression

Results/ – figures and summary outputs, organized by analysis stage.
