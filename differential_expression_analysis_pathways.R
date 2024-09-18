countsData<-round(kegg_pathway_table)

sampleInfo <- data.frame(
  row.names = sort(colnames(countsData)),
  treatment = c(rep("control", 21), rep("treatment", 21)),
  replicate = rep(c("C1", "C2", "C3", "T1", "T2", "T3"), each = 7),
  time = rep(c(0, 12, 24, 48, 96, 168, 336), 6)
)
sampleInfo <- sampleInfo[match(colnames(countsData), rownames(sampleInfo)),]

dds <- DESeqDataSetFromMatrix(countData=countsData, colData=sampleInfo, design=~treatment + time + treatment:time)

# Pre-filtering
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
results_treatment <- results(dds, name="treatment_treatment_vs_control")
results_time <- results(dds, name="time")
results_interaction <- results(dds, name="treatmenttreatment.time")

# Adjust for multiple testing
results_treatment$padj <- p.adjust(results_treatment$pvalue, method="BH")
results_time$padj <- p.adjust(results_time$pvalue, method="BH")
results_interaction$padj <- p.adjust(results_interaction$pvalue, method="BH")

# Sort by adjusted p-value and handle NAs
sorted_treatment <- results_treatment[order(results_treatment$padj, na.last=NA), ]
sorted_time <- results_time[order(results_time$padj, na.last=NA), ]
sorted_interaction <- results_interaction[order(results_interaction$padj, na.last=NA), ]

# Filter for significant genes
sig_treatment <- sorted_treatment[!is.na(sorted_treatment$padj) & sorted_treatment$padj < 0.05, ]
sig_time <- sorted_time[!is.na(sorted_time$padj) & sorted_time$padj < 0.05, ]
sig_interaction <- sorted_interaction[!is.na(sorted_interaction$padj) & sorted_interaction$padj < 0.05, ]

# Create a data frame with the desired columns
p_values_deseq_pathways <- data.frame(
    gene = rownames(results_treatment),
    treatment_fold_change = results_treatment$log2FoldChange,
    treatment_padj = results_treatment$padj,
    time_fold_change = results_time$log2FoldChange,
    time_padj = results_time$padj,
    interaction_fold_change = results_interaction$log2FoldChange,
    interaction_padj = results_interaction$padj
)
write.table(p_values_deseq_pathways, file=paste0(folder_path,"KEGG_DESEQ_P_vals_PROTISTS.tsv"), sep="\t",row.names=FALSE)
#results_treatment <- results(dds, name="treatment_treatment_vs_control")
#
## Initialize data frame with the treatment results
#combined_results <- data.frame(
#    gene = rownames(results_treatment),
#    treatment_fold_change = results_treatment$log2FoldChange,
#    treatment_padj = results_treatment$padj
#)
#
## Extract and combine time effects
#for (i in 2:7) {
#    time_contrast <- paste0("time_", i, "_vs_1")
#    res_time <- results(dds, name=time_contrast)
#    combined_results[[paste0("time", i, "_fold_change")]] = res_time$log2FoldChange
#    combined_results[[paste0("time", i, "_padj")]] = res_time$padj
#}
#
## Extract and combine interaction effects
#for (i in 2:7) {
#    interaction_contrast <- paste0("treatmenttreatment.time", i)
#    res_interaction <- results(dds, name=interaction_contrast)
#    combined_results[[paste0("interaction_time", i, "_fold_change")]] = res_interaction$log2FoldChange
#    combined_results[[paste0("interaction_time", i, "_padj")]] = res_interaction$padj
#}
#
#
#
#do_expr_analysis_edgeR<-function(counts){
#  library(edgeR)
#
#  # Define group information
#  group <- factor(c("control", "control", "control", "treatment", "treatment", "treatment"))
#
#  # Create DGEList object
#  y <- DGEList(counts=counts, group=group)
#
#  # Filter lowly expressed genes
#  keep <- filterByExpr(y)
#  y <- y[keep, , keep.lib.sizes=FALSE]
#
#  # Normalize the data
#  y <- calcNormFactors(y)
#
#  # Estimate dispersion
#  y <- estimateDisp(y, design=model.matrix(~group))
#
#  # Fit the model
#  fit <- glmQLFit(y, design=model.matrix(~group))
#
#  # Perform differential expression test
#  qlf <- glmQLFTest(fit, contrast=c(-1, 1))
#
#  # Extract and adjust results
#  results <- topTags(qlf, n=Inf)
#  results$table$FDR <- p.adjust(results$table$PValue, method="fdr")
#  return(results)
#}
#do_expr_analysis_deseq<-function(counts){
#  library(DESeq2)
#  # Create column data frame
#  colData <- data.frame(
#      row.names=colnames(counts),
#      condition=factor(c("control", "control", "control", "treatment", "treatment", "treatment"))
#  )
#
#  # Create DESeqDataSet object
#  dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
#
#  # Pre-filtering
#  dds <- dds[rowSums(counts(dds)) > 1, ]
#
#  # Run the DESeq2 pipeline
#  dds <- DESeq(dds)
#
#  # Extract results
#  results <- results(dds)
#
#  # Adjust for multiple testing
#  results$padj <- p.adjust(results$pvalue, method="BH")
#  return(results)
#}