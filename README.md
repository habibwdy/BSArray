# 🧬 BSArray

BSArray is an automated pipeline for Bulk Segregant Analysis (BSA) using Illumina Infinium SNP genotyping arrays.  
It enables rapid and cost-effective identification of genomic regions (QTL) associated with traits of interest.

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
