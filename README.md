# Integrative Epigenomic Analysis of Colorectal Cancer
## CRC-TGCA-Analysis
### TCGA-COAD | R + Python | BigQuery | 2025

---



This is version 2 of my colorectal cancer (CRC) epigenomics project. Version 1 used three small public GEO microarray datasets and identified six candidate genes consistently dysregulated across independent experiments. Version 2 takes those findings and tests them in a much larger, better-powered dataset.

The core question is: could I demonstrate that epigenetic changes (specifically DNA methylation) are directly silencing tumour suppressor genes in CRC? And can machine learning identify gene expression signatures that distinguish cancer from normal tissue?

---

## Data Source

All data comes from TCGA-COAD (The Cancer Genome Atlas, Colon Adenocarcinoma), accessed through Google BigQuery via the ISB-CGC public dataset. TCGA-COAD has approximately 460 patients with matched RNA-seq and DNA methylation 450k array data measured on the same tumour biopsy.

No raw data files are stored in this repository. The BigQuery queries in the R script reproduce the data download from scratch.

---

## Repository Structure

```
CRC-TCGA-Analysis/
├── R/
│   └── 01_bigquery_integration.R    # Data pull, integration, correlation analysis
├── CRC_TCGA_ML_Classifier.ipynb     # Python ML classifier in Google Colab
└── README.md
```

---

## Approach

### Step 1: Data Access via BigQuery

TCGA data was queried using SQL through Google BigQuery. Three tables were pulled:

- RNA-seq expression (fpkm_uq_unstranded) for 39 known CRC-relevant genes across all TCGA-COAD primary tumour samples
- DNA methylation beta values for 20 CpG probes at promoters of key CRC genes
- Clinical metadata including vital status and age at diagnosis

456 patients had both RNA-seq and methylation data available. After merging all three tables by case barcode, 296 patients had complete data across all three sources.

### Step 2: Integration in R

RNA-seq and methylation data were reshaped from long format to wide format (one row per patient). Expression values were log2 transformed. All three tables were merged by patient ID.

### Step 3: Correlation Analysis

For each methylation probe and its associated gene, I computed Pearson correlation across all 296 patients. A negative correlation means higher methylation at the promoter leads to lower gene expression, which is the epigenetic silencing signal.

### Step 4: ML Classifier in Python

A Random Forest classifier was trained to distinguish primary tumour samples from normal tissue using 38 gene expression features. The dataset included 458 tumour and 41 normal samples. Class imbalance was handled with balanced class weights. Performance was evaluated with 5-fold stratified cross-validation.

---

## Key Findings

### Epigenetic Silencing of CDH1

CDH1 (E-cadherin) showed a significant negative correlation between promoter methylation (cg23606718) and gene expression across 296 patients (r = -0.24, p < 0.0001). As promoter methylation increased from 0 to 0.8, CDH1 expression dropped consistently. CDH1 encodes a cell adhesion protein whose epigenetic silencing enables tumour invasion. This is a direct demonstration of epigenetic silencing in CRC.

### Tumour vs Normal Classifier

The Random Forest classifier achieved a cross-validated AUC of 1.0, perfectly distinguishing tumour from normal tissue using 38 gene expression features. The top features were:

| Rank | Gene | Significance |
|---|---|---|
| 1 | BMP3 | Clinically validated CRC biomarker used in stool DNA screening tests |
| 2 | KRT20 | Classic colorectal epithelial marker |
| 3 | MYC | Most frequently amplified oncogene in CRC |
| 4 | SFRP1 | Wnt inhibitor, frequently methylated in CRC |
| 5 | VEGFA | Angiogenesis driver, therapeutic target |
| 6 | PIGR | Found in v1 cross-dataset analysis, confirmed here |

### Connection to Version 1

Four of the six candidate genes identified in version 1 (PIGR, GNG4, MAB21L2, H19) appear in the top 20 features of the v2 classifier. These genes were identified independently from three small GEO datasets and confirmed in 460 TCGA patients using a completely different method. I am happy about this if we are being honest. That consistency across datasets and methods is the strongest signal in this project.

### Survival Prediction

A classifier trained to predict vital status (Alive vs Dead) from molecular features achieved a cross-validated AUC of 0.506, essentially random chance. (I was playing around with this even though it didnt make sense to do so lol but anyways -)Survival in CRC depends on many factors not captured here including treatment received, surgical margins, and metastasis status. Molecular features alone are not sufficient for outcome prediction without clinical covariates.

---

## How to Run

### R Analysis

1. Set up a Google Cloud account and enable BigQuery
2. Pin the isb-cgc-bq project in BigQuery
3. Run the BigQuery queries in `R/01_bigquery_integration.R` to download data
4. Save the three CSVs and update the file paths in the script
5. Run the integration and correlation sections

### Python ML

Open `CRC_TCGA_ML_Classifier.ipynb` in Google Colab. Upload the three CSV files when prompted. Run all cells in order.

---

## Dependencies

R packages: dplyr, tidyr, ggplot2

Python packages: pandas, numpy, scikit-learn, matplotlib, seaborn

---

## Limitations

The gene list is curated rather than genome-wide, which limits discovery but makes the analysis tractable. The methylation probe set is small. The survival classifier result shows that outcome prediction requires clinical covariates alongside molecular features. Cross-dataset comparison between v1 (microarray) and v2 (RNA-seq) is qualitative rather than quantitative due to platform differences.

---

## Connection to Version 1

Version 1 repository: github.com/pranathiperii/CRC-GEO-Analysis

Version 1 used GEO microarray datasets (GSE30540, GSE4107, GSE32323) with 10 to 35 samples per group. It identified six candidate genes across three independent datasets. Version 2 uses TCGA with 456 patients to validate those candidates and adds integrative methylation analysis and machine learning.

---

## References

- TCGA Research Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature 2012.
- Goldman et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nature Biotechnology 2020.
- ISB-CGC: isb-cgc.org
- Pedregosa et al. Scikit-learn: Machine Learning in Python. JMLR 2011.
