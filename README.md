# 🧬 BSArray

BSArray is an automated pipeline for Bulk Segregant Analysis (BSA) using Illumina Infinium SNP genotyping arrays.  
It enables rapid and cost-effective identification of genomic regions (QTL) associated with traits of interest.

---

## 🚀 Quick Start

### 1. Download the repository

Clone the repository:

```bash

git clone https://github.com/habibwdy/BSArray.git

```

or download it as a ZIP file from GitHub.

### 2. Open the project

Open **BSArray.Rproj** (recommended) or open the project folder in RStudio.

### 3. Install required packages

Run:

```r

install.packages(c("shiny"))

```

or install any additional packages if prompted.

### 4. Launch BSArray

```r

shiny::runApp()

```

Alternatively, open `app.R` in RStudio and click **Run App**.

### 5. Run your first analysis

BSArray includes example datasets in the repository.

Simply:

1. Load the example bulk genotype file

2. Load the example parental genotype file

3. Click **Run Analysis**

4. Explore the genomic landscape plots (the Yosemite plot ;) ) and candidate QTL

---

## 🚀 Key Features

- SNP array-based BSA (no sequencing required)
- Fast and computationally efficient
- Quality-weighted allele frequency analysis (W statistic)
- Sliding window smoothing for robust QTL detection
- Integrated visualization (genomic landscape plots)
- Supports multiple analysis modes (Default & Classic)

---

## 🧠 Conceptual Workflow

Theta → Allele Frequency → ΔAF → W statistic → Sliding window → QTL peaks 

---

## 📥 Input Requirements

### Bulk file (Standard Report from GenomeStudio)
- SNP Name  
- Sample ID  
- Allele1 – Top  
- Allele2 – Top  
- GC Score  
- GT Score  
- Theta  
- R  

### Parental file (Matrix Report)
- Top strand alleles  
- GenCall column excluded  

---

## ⚙️ Analysis Modes

### 🔹 Default Mode (W Statistic)
- Computes ΔAF between bulks  
- Applies genotype quality weighting  
- Uses sliding window smoothing  
- Detects QTL via empirical thresholds  

### 🔹 Classic Mode
- Direct genotype matching with parents  
- Supports:
  - Strict model  
  - Dominant model  
  - Recessive model  

---

## 📊 Output

- 📈 Genomic landscape plots  
- 📍 Candidate QTL regions  
- 🧬 Significant SNP markers  
- 📄 Excel output with full results  
- 🧪 Gene annotation (±10 kb)  

---

## 🧪 Step-by-Step Tutorial

BSArray automatically loads the required SNP and gene databases from the `data/` folder. Users should keep the original folder structure after downloading or cloning the repository.

1. Prepare input files from GenomeStudio  
2. Upload bulk and parent files  
3. Assign samples to bulks and parents  
4. Set parameters (default recommended)  
5. Run analysis  
6. Explore genomic landscape plots  
7. Identify peaks above thresholds  
8. Extract significant markers  
9. Interpret results  
10. Export results  

---

## 🎯 Summary

BSArray transforms SNP array data into actionable QTL insights by:
- Quantifying allele frequency differences  
- Adjusting for genotype quality  
- Aggregating signals across genomic regions  
- Detecting peaks exceeding empirical thresholds  

---

## 👤 Author
Habib Widyawan, Zenglu Li (UGA Soybean Breeding and Genetics Lab)
