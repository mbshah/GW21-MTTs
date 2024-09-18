# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)


countsData <- read.table("/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/KO_counts_HMM_version_trs_pred_0.001.tsv", header = TRUE, row.names = 1, sep = "\t")
countsData <- round(countsData)
sampleInfo <- data.frame(
  row.names = sort(colnames(countsData)),
  condition = c(rep("control", 21), rep("treatment", 21)),
  replicate = rep(c("C1", "C2", "C3", "T1", "T2", "T3"), each = 7),
  timePoint = rep(c("S1", "S2", "S3", "S4", "S5", "S6", "S7"), 6)
)
sampleInfo <- sampleInfo[match(colnames(countsData), rownames(sampleInfo)),]
#normalizing the counts first
dds <- DESeqDataSetFromMatrix(countData = countsData, colData = sampleInfo, design = ~condition + timePoint)
dds <- DESeq(dds)
normCounts <- counts(dds, normalized = TRUE)
normCounts <- round(normCounts)


# to do differential expression only on first time point to identify the genes brought in by wastewater:
sampleInfoS1 <- sampleInfo[grep("_S1", rownames(sampleInfo)),]
countsDataS1 <- normCounts[, rownames(sampleInfoS1)]
ddsS1 <- DESeqDataSetFromMatrix(countData = countsDataS1, colData = sampleInfoS1, design = ~condition)
ddsS1 <- DESeq(ddsS1)
res <- results(ddsS1, contrast = c("condition", "control", "treatment"))
# Order the results by the p-value (smallest first)
resOrdered <- res[order(res$log2FoldChange),] # sort by log fold chanfge and take significant ones

# Subset to get only significantly downregulated genes (here means the one significantly more in treatment meaning they were introduced by wastewater)
sigUp <- subset(resOrdered, padj < 0.05 & log2FoldChange < 0)
all_sigup<-rownames(sigUp)
pathways_info_sigup<- ko_pathway_final %>%
  filter(KEGG_Ortholog %in% all_sigup)
topKOs <- head(sigUp, 15) #taking top 15 most in wastewater treated sample
topKOsNames <- rownames(topKOs)
topCountsData <- countsData[topKOsNames,]
df <- topCountsData
new_columns <- list()

# Loop over the unique time points
for (i in 1:7) {
  # Calculate the mean for the 'C' columns
  new_columns[[paste0("C", i)]] <- rowMeans(df[, grepl(paste0("^C.*_S", i, "$"), names(df))], na.rm = TRUE)
  # Calculate the mean for the 'T' columns
  new_columns[[paste0("T", i)]] <- rowMeans(df[, grepl(paste0("^T.*_S", i, "$"), names(df))], na.rm = TRUE)
}
# Combine the new columns with the original data frame
df_means <- cbind(df, new_columns)
final_df_top_brought_in <- df_means %>% select(C1:T7)


# Convert row names to a column for ggplot2 compatibility
final_df_top_brought_in$KEGG_ortholog <- rownames(final_df_top_brought_in)

# Split the data frame into a list of data frames for each KEGG ortholog
list_of_dfs <- split(final_df_top_brought_in, final_df_top_brought_in$KEGG_ortholog)

# Initialize an empty list to store ggplot objects
plot_list <- list()
plot_list2<-list()


# Loop through each mini data frame and create a plot
for (i in 1:length(list_of_dfs)) {
  # Get the current data frame
  df <- list_of_dfs[[i]]

  # Melt the data frame from wide to long format for ggplot2
  long_df <- reshape2::melt(df, id.vars = 'KEGG_ortholog')

  # Separate the 'variable' into 'group' and 'time_point'
  long_df <- long_df %>%
    separate(variable, into = c("group", "time_point"), sep = 1)

  diff_df <- aggregate(value ~ time_point, data = long_df, FUN = function(x) diff(x))
  diff_df$time_point<-as.numeric(diff_df$time_point)
  diff_df$group<-ifelse(diff_df$value<0,"Less in Test", "More in Test")
  # Convert 'time_point' to a numeric value for plotting
  long_df$time_point <- as.numeric(gsub("C|T", "", long_df$time_point))
  kid <- long_df$KEGG_ortholog[1]
  koname <- get_kegg_name(kid)
  ktitle <- paste0(kid, " ", koname)
  kanno_vector<-ko_pathway_final[ko_pathway_final$KEGG_Ortholog == kid, ]$Pathway
  kanno_str<-""
  if(length(kanno_vector) >= 1){
    print(kanno_vector)
    #print(kanno_str)
    kanno_str<-paste(kanno_vector[[1]], collapse = "\n")
  }

  # Create a multiline plot for the current data frame without a legend
  p <- ggplot(diff_df, aes(x = time_point, y = value, color = group)) +
    geom_point() +
    #geom_bar(stat = "identity")+
    theme_minimal() +
    labs(title = ktitle, x = "Time Point", y = "Reads abundance difference")+xlim(0,12)+
    theme(legend.position = "none")+ annotate("text", label = kanno_str, x=9, y=300)

  # Add the plot to the list
  plot_list[[i]] <- p

  q <- ggplot(long_df, aes(x = time_point, y = value, color = group, fill= group)) +
    geom_line() +
    #geom_bar(stat = "identity")+
    theme_minimal() +
    labs(title = ktitle, x = "Time Point", y = "Reads abundance difference")+xlim(0,12)+
    theme(legend.position = "none")+ annotate("text", label = kanno_str, x=9, y=300)

  plot_list2[[i]] <- q
}


# Create a separate legend using the last plot
legend_plot <- get_legend(p + theme(legend.position = "bottom"))
legend_plot2<-get_legend(q + theme(legend.position = "bottom"))

# Combine all plots and the separate legend into one image
combined_plot <- do.call(gridExtra::grid.arrange, c(grobs = c(plot_list, list(legend_plot)), ncol = 3))
combined_plot <- do.call(gridExtra::grid.arrange, c(grobs = c(plot_list2, list(legend_plot2)), ncol = 3))

