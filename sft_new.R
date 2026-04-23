#!/usr/bin/env Rscript

### ==============================
###   WGCNA + DESeq2 Pipeline
###   (HPC-friendly / PNG plots)
### ==============================

library(DESeq2)
library(WGCNA)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
options(bitmapType = "cairo")  # Important for PNG output on HPC


outdir <- "results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)			
# ------------------------------
# 1. Load Count Data
# ------------------------------

counts <- read.csv("mega_Matrix_all_samples_copper.tsv",
                   header = TRUE,
                   sep = "\t",
                   row.names = 1)

counts <- na.omit(counts)

sample_names <- colnames(counts)
condition <- c(rep("Control", 6), rep("Patient", 6))   # adjust if sample order changes

colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Normalization
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)


# Save VST-normalized counts (TSV)
write.table(vst_counts,
            file = file.path(outdir, "vst_normalized_counts.tsv"),
            sep = "\t", quote = FALSE)


# ------------------------------
# 2. WGCNA Expression Preparation
# ------------------------------
datExpr <- t(vst_counts)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Save the matrix used for WGCNA
write.table(datExpr,
            file = file.path(outdir, "wgcna_input_matrix.tsv"),
            sep = "\t", quote = FALSE)
# ------------------------------
# 3. Sample Clustering Plot
# ------------------------------
sampleTree <- hclust(dist(datExpr), method = "average")
png("01_sample_clustering.png", width = 1800, height = 1200, res = 200)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5)
dev.off()

# ------------------------------
# 4. Soft-threshold Power Selection
# ------------------------------
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

png("02_soft_threshold.png", width = 1800, height = 1200, res = 200)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Fit (signed R^2)",
     type = "n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "blue")
dev.off()



softPower <- 12  # choose based on above plot

# ------------------------------
# 5. Network Construction
# ------------------------------
adjacency <- adjacency(datExpr, power = softPower,type = "signed")

TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
png("03_gene_clustering.png", width = 1800, height = 1200, res = 200)
plot(geneTree, main = "Gene clustering", xlab = "", sub = "")
dev.off()

# ------------------------------
# 8. Module identification (dynamic tree cut)
# ------------------------------
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)
moduleColors <- labels2colors(dynamicMods)

png(file.path(outdir, "04_module_dendrogram.png"), width = 1800, height = 1200, res = 200)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03)
dev.off()

# Save initial module assignment (use colnames(datExpr) as genes)
module_assignment_before <- data.frame(
  Gene = colnames(datExpr),
  ModuleColor = moduleColors,
  stringsAsFactors = FALSE
)
write.table(module_assignment_before,
            file = file.path(outdir, "module_assignments_before_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------
# 9. Module eigengenes and merging
# ------------------------------
MElist <- moduleEigengenes(datExpr, colors = moduleColors)
MEs <- MElist$eigengenes

# Save module eigengenes before merging
write.table(MEs,
            file = file.path(outdir, "module_eigengenes_before_merging.tsv"),
            sep = "\t", quote = FALSE)

ME.dissimilarity <- 1 - cor(MEs, use="complete")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")

png(file.path(outdir, "05_module_eigengene_clustering.png"), width = 1200, height = 800, res = 200)
plot(METree, main = "Clustering of module eigengenes")
abline(h = 0.25, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, moduleColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

png("06_Merged_Module_Dendrogram.png", 
    width = 1800, height = 1200, res = 150)

plotDendroAndColors(
  geneTree,
  cbind(moduleColors, mergedColors),
  c("Original Module", "Merged Module"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors (Original vs Merged)"
)

dev.off()

# Save merged module assignment (uses colnames(datExpr) as genes)
module_assignment_after <- data.frame(
  Gene = colnames(datExpr),
  MergedColor = mergedColors,
  stringsAsFactors = FALSE
)
write.table(module_assignment_after,
            file = file.path(outdir, "module_assignments_after_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save module eigengenes after merging
write.table(mergedMEs,
            file = file.path(outdir, "module_eigengenes_after_merging.tsv"),
            sep = "\t", quote = FALSE)

# ------------------------------
# 10. Save additional results and summaries
# ------------------------------
# gene dendrogram order (gene indices) and names
write.table(geneTree$order,
            file = file.path(outdir, "gene_dendrogram_order_indices.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE)
write.table(colnames(datExpr)[geneTree$order],
            file = file.path(outdir, "gene_dendrogram_order_genes.tsv"),
            sep = "\t", quote = FALSE, col.names = FALSE)

# TOM dissimilarity (warning: can be large)
# write.table(as.data.frame(as.matrix(dissTOM)),
#            file = file.path(outdir, "tom_dissimilarity.tsv"),
#            sep = "\t", quote = FALSE)

# Module sizes before/after
module_sizes_before <- as.data.frame(table(moduleColors))
colnames(module_sizes_before) <- c("Module", "NumGenes")
write.table(module_sizes_before,
            file = file.path(outdir, "module_sizes_before_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

module_sizes_after <- as.data.frame(table(mergedColors))
colnames(module_sizes_after) <- c("Module", "NumGenes")
write.table(module_sizes_after,
            file = file.path(outdir, "module_sizes_after_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Module gene lists before/after
module_info_before <- data.frame(Gene = colnames(datExpr), Module = moduleColors, stringsAsFactors = FALSE)
write.table(module_info_before,
            file = file.path(outdir, "module_gene_list_before_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

module_info_after <- data.frame(Gene = colnames(datExpr), Module = mergedColors, stringsAsFactors = FALSE)
write.table(module_info_after,
            file = file.path(outdir, "module_gene_list_after_merging.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save vst_counts and datExpr too
write.table(vst_counts, file = file.path(outdir, "vst_normalized_counts.tsv"), sep = "\t", quote = FALSE)
write.table(datExpr, file = file.path(outdir, "wgcna_input_matrix.tsv"), sep = "\t", quote = FALSE)



# 8. MODULE–TRAIT HEATMAP
# ------------------------------
# 8. MODULE–TRAIT HEATMAP
# ------------------------------

# Convert trait to numeric
#condition_num <- ifelse(condition == "Treated", 1, 0)
# FIX TRAIT ENCODING
condition_num <- ifelse(condition == "Patient", 1, 0)

datTraits <- data.frame(condition_num = condition_num,
                        row.names = sample_names)

nSamples <- nrow(datExpr)

# Module–trait correlation
moduleTraitCor <- cor(mergedMEs, datTraits$condition_num, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Text to display (correlation + p-value)
textMatrix <- paste0(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

# SAVE LARGE HIGH-QUALITY PNG
# SAVE LARGE HIGH-QUALITY PNG
png(
  "07_Module_Trait_Heatmap_HIGH_QUALITY.png",
  width  = 3800,   # increased width → wider boxes
  height = 3000,   # increased height → taller boxes
  res    = 350     # high resolution
)

# Increase margins so labels don't get crowded
par(mar = c(8, 10, 5, 5))  # bottom, left, top, right

# Draw heatmap
labeledHeatmap(
  Matrix = moduleTraitCor,
  
  xLabels = "Condition",
  yLabels = names(mergedMEs),
  ySymbols = names(mergedMEs),
  
  colorLabels = FALSE,
  colors = blueWhiteRed(100),
  textMatrix = textMatrix,
  
  # 🔽 Smaller numbers inside heatmap cells
  cex.text = 0.5,
  
  # 🔽 Reduce Y-axis (module name) text size
  cex.lab.y = 0.9,
  
  # 🔼 Keep X-axis readable
  cex.lab.x = 1.4,
  
  # Title size
  cex.main = 2.0,
  
  zlim = c(-1, 1),
  main = "Module–Trait Relationships"
)

dev.off()

cat("\n✔ WGCNA Pipeline Completed Successfully. Heatmap Saved.\n")

# ------------------------------
# 9. Export Module–Trait Results (TSV)
# ------------------------------

module_trait_results <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = as.numeric(moduleTraitCor[, 1]),
  P_value = as.numeric(moduleTraitPvalue[, 1]),
  Direction = ifelse(moduleTraitCor[, 1] > 0, "UP in Condition", "DOWN in Condition"),
  stringsAsFactors = FALSE
)

# Sort by absolute correlation (strongest first)
module_trait_results <- module_trait_results[
  order(abs(module_trait_results$Correlation), decreasing = TRUE),
]

# Save to TSV
write.table(
  module_trait_results,
  file = "07_Module_Trait_Correlation_Results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n✔ WGCNA Pipeline Completed Successfully. Results written to:", normalizePath(outdir), "\n")

