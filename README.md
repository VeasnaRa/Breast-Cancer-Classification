# Breast-Cancer-Classification
A machine learning project for predicting breast cancer patient vital status using clinical variables and high-dimensional mRNA gene expression data through regularized regression methods.

## 📋 Project Overview

This project implements and compares multiple penalized regression models to classify breast cancer patients' vital status (Alive/Dead) using a dataset of 1,230 patients with 24 clinical features and 5,000 mRNA gene expression measurements. The analysis addresses key challenges in biomedical prediction: high dimensionality (p >> n), class imbalance, and feature correlation.

**Key Achievements:**
- Best model: **UniLasso + SMOTE** with Clinical + TOP20 genes
- **F1-Score:** 0.737 | **AUC:** 0.893 | **Recall:** 70%
- Effective handling of 5:1 class imbalance using SMOTE
- Sparse, interpretable solutions with biological relevance

## 🗂️ Repository Structure
```
.
├── final_implementation.Rmd      # Main analysis notebook (R Markdown)
├── final_implementation.html     # Rendered analysis notebook in html
├── final_implementation.pdf      # Rendered analysis notebook in pdf
├── function_model_fit.R          # Custom model fitting functions
├── mrr_bio.Rdata                 # Preprocessed dataset
├── lasso_extension/              # Extended Lasso implementations
├── model_metrics/                # Model evaluation results
├── slide_presentation.pdf        # Project presentation slides
└── summary_report.pdf            # Detailed technical report
```

## 📊 Dataset

**Clinical Data:** 1,231 patients × 24 features
- Demographics: age, weight, race, ethnicity
- Tumor characteristics: AJCC stage, morphology, diagnosis
- Treatment information: prior treatment, disease response

**Genomic Data:** 123 patients × 5,000 high-variance mRNA transcripts
- Log-transformed and variance-filtered gene expression
- Captures biological pathways: immune infiltration, stromal activation, metabolism

**Target Variable:** Vital Status (Alive/Dead)
- Class distribution: ~5:1 imbalance (Alive:Dead)

## 🔬 Methodology

### 1. Data Preprocessing
- **Missing value imputation:** Median for numerical, mode for categorical variables
- **Statistical testing:** t-tests and Wilcoxon rank-sum  tests for feature significance
- **Differential expression analysis:** Identified key genes (APOB, LYVE1, LINC01497, AC104211.1)

### 2. Feature Engineering
- Gene screening using univariate differential expression (limma)
- Nested feature subsets: TOP-20, TOP-50, TOP-100, TOP-500, TOP-1000, TOP-5000
- PCA analysis revealing expression-based clusters

### 3. Models Implemented

| Model | Description |
|-------|-------------|
| **Logistic Regression** | Baseline model with clinical predictors |
| **Ridge Regression** | L2 penalty for handling multicollinearity |
| **Lasso Regression** | L1 penalty for sparse feature selection |
| **Elastic Net** | Combined L1/L2 penalties |
| **Adaptive Lasso** | Weighted Lasso (Zou, 2006) for better variable selection |
| **UniLasso** | Two-stage sparse method (Chatterjee et al., 2025) |

### 4. Class Imbalance Handling
- **SMOTE** (Synthetic Minority Over-sampling Technique)
- Balanced training data to ~1:1 ratio (141 → 705 minority samples)
- k=5 nearest neighbors, test set kept original

## 📈 Results

### Best Model Performance (UniLasso + SMOTE with Clinical + TOP20)

| Metric | Value |
|--------|-------|
| F1-Score | 0.737 |
| AUC | 0.893 |
| Recall | 0.700 |
| Precision | 0.778 |
| Specificity | High |

### Key Findings

1. **Feature Set Optimization:** Clinical + TOP20 genes consistently outperformed larger gene sets
2. **SMOTE Impact:** Improved recall from 45% → 80% while maintaining high AUC
3. **Model Comparison:** UniLasso provided best precision-recall balance
4. **Overfitting Control:** Larger gene sets (>100) degraded test performance

### Model Comparison
```
Model              | F1    | Recall | Precision | AUC
-------------------|-------|--------|-----------|------
UniLasso + SMOTE   | 0.737 | 0.700  | 0.778     | 0.893
Adaptive Lasso     | 0.717 | 0.825  | 0.635     | 0.908
ElasticNet         | 0.717 | 0.825  | 0.635     | 0.905
Lasso              | 0.696 | 0.800  | 0.615     | 0.905
Logistic           | 0.711 | 0.800  | 0.640     | 0.902
Ridge              | 0.653 | 0.825  | 0.541     | 0.896
```

## 🚀 Getting Started

### Prerequisites
```r
# Required R packages
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyverse",
  "corrplot",
  "ggcorrplot",
  "naniar",
  "visdat",
  "rstatix",
  "DescTools",
  "car",
  "limma",
  "survival",
  "pheatmap",
  "diptest",
  "forcats",
  "glmnet",
  "caTools",
  "pROC",
  "gridExtra",
  "smotefamily",
  "reshape2"
))
```

### Running the Analysis

1. **Clone the repository:**
```bash
git clone https://github.com/VeasnaRa/Breast-Cancer-Classification.git
```
2. **Run the main analysis:**
```r
# Open and knit the R Markdown file
rmarkdown::render("final_implementation.Rmd")
```

4. **View results:**
Open `final_implementation.html` in your browser

## 📚 Key References

- Zou, H. (2006). "The Adaptive Lasso and Its Oracle Properties". *Journal of the American Statistical Association*
- Chatterjee, A., Hastie, T., & Tibshirani, R. (2025). "UniLasso: Two-stage sparse regression for high-dimensional genomic data"
- TCGA-BRCA cohort preprocessing methodology

## 👥 Authors

- **RA Veasna**
- **DIN Sokheng**

**Supervisor:** Juhyun PARK

**Institution:** École Nationale Supérieure d'Informatique pour l'Industrie et l'Entreprise (ENSIIE)

**Date:** December 15, 2025
