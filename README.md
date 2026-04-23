# 🧬 WGCNA + DESeq2 RNA-seq Analysis Pipeline

A reproducible, HPC-friendly pipeline for RNA-seq analysis combining **DESeq2 normalization** and **WGCNA (Weighted Gene Co-expression Network Analysis)** to identify gene modules and their association with traits.

---

## 📌 Overview

This pipeline performs:

- RNA-seq count data processing  
- Low-expression gene filtering  
- Variance Stabilizing Transformation (VST) using DESeq2  
- Co-expression network construction (WGCNA)  
- Module detection and merging  
- Module–trait correlation analysis  
- High-quality visualization outputs  

Designed for **large datasets** and optimized for **HPC (SLURM)** environments.

