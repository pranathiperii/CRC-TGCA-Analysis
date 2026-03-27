# =============================================================================
# 01_bigquery_integration.R
# TCGA-COAD Integrative Analysis: Methylation + Gene Expression
#
# This script covers:
#   1. BigQuery SQL queries to pull TCGA-COAD data
#   2. Data loading and integration in R
#   3. Correlation analysis: methylation vs expression
#   4. Visualisation of epigenetic silencing
#
# Data source: TCGA-COAD via ISB-CGC BigQuery (isb-cgc-bq)
# Patients: 296 with complete RNA-seq + methylation + clinical data
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# SECTION 1: BigQuery SQL Queries
# =============================================================================
# Run these queries in the BigQuery console at bigquery.cloud.google.com
# Pin the project isb-cgc-bq to access TCGA data
# Save each result as a CSV file

# --- Query 1: RNA-seq expression for CRC-relevant genes ---
# Save as COAD_RNAseq.csv
#
# SELECT
#   case_barcode,
#   sample_barcode,
#   gene_name,
#   fpkm_uq_unstranded as expression
# FROM `isb-cgc-bq.TCGA_versioned.RNAseq_hg38_gdc_r42`
# WHERE project_short_name = 'TCGA-COAD'
#   AND sample_type_name = 'Primary Tumor'
#   AND gene_name IN (
#     'H19','PIGR','GNG4','CYP2B6','PRAP1','MAB21L2',
#     'APC','KRAS','TP53','SMAD4','PIK3CA','BRAF','MLH1',
#     'MSH2','MSH6','PMS2','CTNNB1','CDH1','VIM','SEPT9',
#     'SDC2','NDRG4','BMP3','TFPI2','SFRP1','SFRP2','WIF1',
#     'DKK1','DKK2','DKK3','RUNX3','CDKN2A','MGMT',
#     'MYC','EGFR','ERBB2','VEGFA','CDX2','CEACAM5','EPCAM'
#   )

# --- Query 2: DNA methylation beta values ---
# Save as COAD_methylation.csv
#
# SELECT
#   m.case_barcode,
#   m.probe_id,
#   m.beta_value
# FROM `isb-cgc-bq.TCGA_versioned.DNA_methylation_hg38_gdc_2017_01` m
# WHERE m.case_barcode IN (
#     SELECT DISTINCT case_barcode
#     FROM `isb-cgc-bq.TCGA_versioned.RNAseq_hg38_gdc_r42`
#     WHERE project_short_name = 'TCGA-COAD'
#     AND sample_type_name = 'Primary Tumor'
#   )
# AND m.probe_id IN (
#     'cg23479922','cg02839557','cg10673833','cg17501210',
#     'cg23606718','cg04420002','cg00486163','cg25305428',
#     'cg01837089','cg08697437','cg16366515','cg10673833',
#     'cg17272132','cg05022650','cg23606718','cg04398124',
#     'cg07376895','cg14584218','cg08972016','cg14190895'
# )

# --- Query 3: Clinical metadata ---
# Save as COAD_clinical.csv
#
# SELECT
#   submitter_id as case_barcode,
#   demo__gender as gender,
#   demo__vital_status as vital_status,
#   demo__days_to_death as days_to_death,
#   demo__age_at_index as age_at_diagnosis,
#   proj__project_id as project
# FROM `isb-cgc-bq.TCGA_versioned.clinical_gdc_r27`
# WHERE proj__project_id = 'TCGA-COAD'

# --- Query 4: Tumour vs Normal (for ML classifier) ---
# Save as COAD_tumour_normal.csv
#
# SELECT
#   case_barcode,
#   sample_barcode,
#   gene_name,
#   fpkm_uq_unstranded as expression,
#   sample_type_name
# FROM `isb-cgc-bq.TCGA_versioned.RNAseq_hg38_gdc_r42`
# WHERE project_short_name = 'TCGA-COAD'
#   AND sample_type_name IN ('Primary Tumor', 'Solid Tissue Normal')
#   AND gene_name IN (
#     'H19','PIGR','GNG4','CYP2B6','PRAP1','MAB21L2',
#     'APC','KRAS','TP53','SMAD4','PIK3CA','BRAF','MLH1',
#     'MSH2','MSH6','PMS2','CTNNB1','CDH1','VIM','SEPT9',
#     'SDC2','NDRG4','BMP3','TFPI2','SFRP1','SFRP2','WIF1',
#     'RUNX3','CDKN2A','MGMT','MYC','EGFR','ERBB2','VEGFA',
#     'CDX2','KRT20','MUC2','CEACAM5','EPCAM'
#   )

# =============================================================================
# SECTION 2: Load and Integrate Data
# =============================================================================
# Update these paths to wherever you saved your CSV files

rnaseq      <- read.csv("COAD_RNAseq.csv.csv")
methylation <- read.csv("COAD_methylation.csv.csv")
clinical    <- read.csv("COAD_clinical.csv..csv")

message(sprintf("RNA-seq: %d rows, %d columns", nrow(rnaseq), ncol(rnaseq)))
message(sprintf("Methylation: %d rows, %d columns", nrow(methylation), ncol(methylation)))
message(sprintf("Clinical: %d patients", nrow(clinical)))

# Reshape RNA-seq to wide format: one row per patient, one column per gene
rnaseq_wide <- rnaseq %>%
  select(case_barcode, gene_name, expression) %>%
  distinct(case_barcode, gene_name, .keep_all = TRUE) %>%
  pivot_wider(names_from = gene_name, values_from = expression)

# Reshape methylation to wide format: one row per patient, one column per probe
methyl_wide <- methylation %>%
  distinct(case_barcode, probe_id, .keep_all = TRUE) %>%
  pivot_wider(names_from = probe_id, values_from = beta_value)

# Merge all three by patient ID
integrated <- rnaseq_wide %>%
  inner_join(methyl_wide, by = "case_barcode") %>%
  inner_join(clinical,   by = "case_barcode")

message(sprintf("Integrated dataset: %d patients, %d features", nrow(integrated), ncol(integrated)))

# Log2 transform expression values
# Expression from RNA-seq FPKM values needs log2 transformation
# Beta values (methylation) are already on 0-1 scale, no transformation needed
expr_cols <- colnames(rnaseq_wide)[-1]
integrated[expr_cols] <- log2(integrated[expr_cols] + 1)

# =============================================================================
# SECTION 3: Correlation Analysis - Methylation vs Expression
# =============================================================================
# For each probe-gene pair, test whether higher promoter methylation
# is associated with lower gene expression (epigenetic silencing)

probe_gene_map <- data.frame(
  probe_id  = c("cg23479922", "cg02839557", "cg10673833",
                "cg17501210", "cg23606718"),
  gene_name = c("MLH1",       "APC",        "MLH1",
                "CDKN2A",     "CDH1")
)

results <- data.frame()
for (i in 1:nrow(probe_gene_map)) {
  probe <- probe_gene_map$probe_id[i]
  gene  <- probe_gene_map$gene_name[i]
  if (probe %in% colnames(integrated) & gene %in% colnames(integrated)) {
    ct <- cor.test(integrated[[probe]], integrated[[gene]], method = "pearson")
    results <- rbind(results, data.frame(
      probe     = probe,
      gene      = gene,
      r         = round(ct$estimate, 3),
      p_value   = round(ct$p.value, 4),
      n         = sum(!is.na(integrated[[probe]]) & !is.na(integrated[[gene]]))
    ))
  }
}

print(results)
write.csv(results, "methylation_expression_correlations.csv", row.names = FALSE)

# =============================================================================
# SECTION 4: Visualise CDH1 Epigenetic Silencing
# =============================================================================
# CDH1 (E-cadherin) shows the strongest negative correlation
# r = -0.24, p < 0.0001 (n = 296)
# Higher promoter methylation -> lower CDH1 expression
# This is direct evidence of epigenetic silencing

cdh1_plot <- ggplot(integrated, aes(x = cg23606718, y = CDH1)) +
  geom_point(alpha = 0.5, colour = "#2c7bb6", size = 1.5) +
  geom_smooth(method = "lm", colour = "#d7191c", se = TRUE, linewidth = 1) +
  labs(
    title    = "Epigenetic Silencing of CDH1 in TCGA-COAD",
    subtitle = "Promoter methylation (cg23606718) vs expression | n = 296 patients",
    x        = "Methylation beta value (cg23606718)",
    y        = "CDH1 expression (log2 FPKM-UQ)",
    caption  = "r = -0.24, p < 0.0001 | Pearson correlation"
  ) +
  theme_classic(base_size = 13)

ggsave("CDH1_epigenetic_silencing.png", cdh1_plot, width = 7, height = 5, dpi = 300)
message("Saved: CDH1_epigenetic_silencing.png")

# =============================================================================
# SECTION 5: Summary
# =============================================================================

message("=== Analysis Complete ===")
message(sprintf("Patients with complete data: %d", nrow(integrated)))
message("Correlation results saved to: methylation_expression_correlations.csv")
message("CDH1 plot saved to: CDH1_epigenetic_silencing.png")
message("For ML classifier, open CRC_TCGA_ML_Classifier.ipynb in Google Colab")
