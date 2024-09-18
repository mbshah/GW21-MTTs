###LOAD DATA USING 005-4 script
source("/home/hb0358/PycharmProjects/MTT_mbs/Statistics_GW21/scripts/load_required_fns_and_data.R")
#BiocManager::install("gage")
library(gage)
library(pathview)




p_values_t_no_na <- p_values_deseq_kegg[!is.na(p_values_deseq_kegg$interaction_padj), ]
significant_kos <- rownames(p_values_t_no_na[p_values_t_no_na$interaction_padj < 0.05, ])
kegg_table_filt<-kegg_table[significant_kos,]

countsData<-round(kegg_table_filt)

sampleInfo <- data.frame(
  row.names = sort(colnames(countsData)),
  treatment = c(rep("control", 21), rep("treatment", 21)),
  replicate = rep(c("C1", "C2", "C3", "T1", "T2", "T3"), each = 7),
  time = factor(rep(1:7, 6))
)

koProfile<-t(countsData)
Kegg_enrichment_per_time_output<-DataFrame()

current_time<-1
ko_treatment<-koProfile[row.names(koProfile) %in% row.names(sampleInfo[sampleInfo$treatment=='treatment' & sampleInfo$time==current_time,]), ]
ko_control<-koProfile[row.names(koProfile) %in% row.names(sampleInfo[sampleInfo$treatment=='control' & sampleInfo$time==current_time,]), ]
#ko_Medium<-koProfile[row.names(koProfile) %in% MAGs_table[MAGs_table$size_lev=='Medium',]$cluster_ID, ]
ko_differential_abundance<-data.frame(Treatment=double(),Control=double())#, Intermediate=double())
for (KID in colnames(koProfile)){
  y<-mean(ko_treatment[,KID])
  x<-mean(ko_control[,KID])
  #z<-mean(ko_Medium[,KID])
  ko_differential_abundance[KID,]<-c(x,y)
}
KOs<-ko_differential_abundance
KOs$KO<-row.names(KOs)


#KOs <- read.table("/home/hb0358/sciebo/mag+gen_spec/MAGs/MAGs_compiled_Paper/tables/gene_differential_presence.tsv", header=T, sep="\t", as.is = T)
# KO should be present in at least 10% of the genomes
KOs2 <- KOs[apply(KOs[,2:3], 1, function(x) any(x>0.1)),] # 2839
# presence in genomes at least twice at high
treatment <- KOs2[log2((KOs2$Treatment+0.0001)/(KOs2$Control+0.0001))>=1,'KO'] # 99
control <- KOs2[log2((KOs2$Treatment+0.0001)/(KOs2$Control+0.0001))<=(-1),'KO'] # 1617
both <- setdiff(setdiff(KOs2$KO, treatment), control) # 1123
all <- KOs$KO # 5777


background.ko <-all
kegg.ko  <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
threshold <- 0.05

sign.ko <- control
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.control<- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- treatment
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.treatment <- signKEGGEnrichment(KEGG.enrichment, threshold)

sign.ko <- both
KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

write.table(sign.KEGG.enrichment.control, file = paste0(folder_path,"KEGG_Enriched_control_sign.tsv"), col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.treatment, file = paste0(folder_path,"KEGG_Enriched_treatment_sign.tsv"), col.names = NA, sep = "\t")
write.table(sign.KEGG.enrichment.both, file = paste0(folder_path,"KEGG_Enriched_both_sign.tsv"), col.names = NA, sep = "\t")

write.table(control, file=paste0(folder_path,"KEGG_Enriched_control_list.tsv"))
write.table(treatment, file=paste0(folder_path,"KEGG_Enriched_treatment_list.tsv"))
write.table(both, file= paste0(folder_path,"KEGG_Enriched_both_list.tsv"))

###END KEGG ENRICHMENT ANALYSIS####

results <- list()
all_pathways <- character()
for (current_time in 1:7) {
  sprintf("Now Running %f",current_time)
  ko_treatment <- koProfile[row.names(koProfile) %in% row.names(sampleInfo[sampleInfo$treatment == 'treatment' & sampleInfo$time == current_time, ]), ]
  ko_control <- koProfile[row.names(koProfile) %in% row.names(sampleInfo[sampleInfo$treatment == 'control' & sampleInfo$time == current_time, ]), ]

  ko_differential_abundance <- data.frame(Treatment = double(), Control = double())
  for (KID in colnames(koProfile)) {
    y <- mean(ko_treatment[, KID])
    x <- mean(ko_control[, KID])
    ko_differential_abundance[KID, ] <- c(x, y)
  }
  KOs <- ko_differential_abundance
  KOs$KO <- row.names(KOs)

  KOs2 <- KOs[apply(KOs[, 2:3], 1, function(x) any(x > 0.1)), ]
  treatment <- KOs2[log2((KOs2$Treatment + 0.0001) / (KOs2$Control + 0.0001)) >= 1, 'KO']
  control <- KOs2[log2((KOs2$Treatment + 0.0001) / (KOs2$Control + 0.0001)) <= (-1), 'KO']
  both <- setdiff(setdiff(KOs2$KO, treatment), control)
  all <- KOs$KO

  background.ko <- all
  kegg.ko <- kegg.gsets(species = "ko", id.type = "kegg")$kg.sets
  threshold <- 0.05

  sign.ko <- control
  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.control <- signKEGGEnrichment(KEGG.enrichment, threshold)

  sign.ko <- treatment
  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.treatment <- signKEGGEnrichment(KEGG.enrichment, threshold)

  sign.ko <- both
  KEGG.enrichment <- KEGGEnrichment(sign.ko, background.ko, kegg.ko)
  sign.KEGG.enrichment.both <- signKEGGEnrichment(KEGG.enrichment, threshold)

  results[[current_time]] <- list(
    control = sign.KEGG.enrichment.control,
    treatment = sign.KEGG.enrichment.treatment,
    both = sign.KEGG.enrichment.both
  )

  if (is.matrix(sign.KEGG.enrichment.control)) {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.control[, "Pathway", drop = FALSE])[[1]])
  } else  {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.control["Pathway", drop = FALSE])[[1]])
  }
  if (is.matrix(sign.KEGG.enrichment.treatment)) {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.treatment[, "Pathway", drop = FALSE])[[1]])
  } else  {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.treatment["Pathway", drop = FALSE])[[1]])
  }
  if (is.matrix(sign.KEGG.enrichment.both)) {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.both[, "Pathway", drop = FALSE])[[1]])
  } else  {
    all_pathways <- append(all_pathways, list(sign.KEGG.enrichment.both["Pathway", drop = FALSE])[[1]])
  }
  all_pathways <- unique(all_pathways)
}

pathway_presence <- data.frame(Pathway = all_pathways)

for (time_point in 1:7) {
  pathway_presence[[paste0("Time", time_point, "_control")]] <- ""
  pathway_presence[[paste0("Time", time_point, "_treatment")]] <- ""
  pathway_presence[[paste0("Time", time_point, "_both")]] <- ""
}

for (time_point in 1:7) {
  if (is.matrix(results[[time_point]]$control)) {
    pathways_control <- results[[time_point]]$control[,"Pathway"]
  }else {
    pathways_control <- results[[time_point]]$control["Pathway"]
  }
  if (is.matrix(results[[time_point]]$treatment)) {
    pathways_treatment <- results[[time_point]]$treatment[,"Pathway"]
  }else {
    pathways_treatment <- results[[time_point]]$treatment["Pathway"]
  }
  if (is.matrix(results[[time_point]]$both)) {
    pathways_both <- results[[time_point]]$both[,"Pathway"]
  }else {
    pathways_both <- results[[time_point]]$both["Pathway"]
  }

  pathway_presence[[paste0("Time", time_point, "_control")]] <- ifelse(pathway_presence$Pathway %in% pathways_control, "x", "")
  pathway_presence[[paste0("Time", time_point, "_treatment")]] <- ifelse(pathway_presence$Pathway %in% pathways_treatment, "x", "")
  pathway_presence[[paste0("Time", time_point, "_both")]] <- ifelse(pathway_presence$Pathway %in% pathways_both, "x", "")
}


# Define a function to extract pathways and KOs from a sub-data structure
extract_pathways_and_kos <- function(sub_data) {
  if (is.null(sub_data)) return(data.frame(Pathway=character(), KOs=character(), stringsAsFactors=FALSE))

  if (is.matrix(sub_data)) {
    return(data.frame(Pathway=sub_data[,1], KOs=sub_data[,2], stringsAsFactors=FALSE))
  }else if (is.character(sub_data)) {
    return(data.frame(Pathway=sub_data["Pathway"], KOs=sub_data["KOs"], stringsAsFactors=FALSE))
  }

  return(data.frame(Pathway=character(), KOs=character(), stringsAsFactors=FALSE))
}

# Initialize an empty data frame to store all results
all_results <- data.frame(Pathway=character(), KOs=character(), stringsAsFactors=FALSE)

# Loop through each time point in the results list
for (time_point in results) {
  # Extract pathways and KOs from 'control', 'treatment', and 'both' if they exist
  if (!is.null(time_point$control)) {
    control_results <- extract_pathways_and_kos(time_point$control)
    all_results <- rbind(all_results, control_results)
  }

  if (!is.null(time_point$treatment)) {
    treatment_results <- extract_pathways_and_kos(time_point$treatment)
    all_results <- rbind(all_results, treatment_results)
  }

  if (!is.null(time_point$both)) {
    both_results <- extract_pathways_and_kos(time_point$both)
    all_results <- rbind(all_results, both_results)
  }
}

# Print the combined table of pathways and KOs
all_results_unique <- all_results %>%
  group_by(Pathway) %>%
  summarise(KOs_in_Data = paste(unique(unlist(strsplit(KOs, ", "))), collapse = ", ")) %>%
  ungroup()
merged_results <- merge(pathway_presence, all_results_unique, by = "Pathway", all.x = TRUE)
merged_results$Pathway_ID <- sapply(strsplit(merged_results$Pathway, " "), function(x) gsub("ko", "map", x[1]))
merged_results_filtered<-merged_results[merged_results$Pathway_ID %in% row.names(pathway_pedia),]
temp_pathway_pedia<-pathway_pedia
temp_pathway_pedia$Pathway_ID<-row.names(temp_pathway_pedia)
row.names(temp_pathway_pedia)<-NULL
merged_results_filtered<-merge(merged_results_filtered,temp_pathway_pedia,by="Pathway_ID", all.x = TRUE)
required_cols<-c("Pathway_ID","Subsystems","Pathway_Names",
                "Time1_control","Time2_control","Time3_control","Time4_control","Time5_control","Time6_control","Time7_control",
                "Time1_treatment","Time2_treatment","Time3_treatment","Time4_treatment","Time5_treatment","Time6_treatment","Time7_treatment",
                "Time1_both","Time2_both","Time3_both","Time4_both","Time5_both","Time6_both","Time7_both",
                "KOs_in_Data","KEGG_Orthologs")
merged_results_filtered<-merged_results_filtered[,required_cols]
write.csv(merged_results_filtered, file = paste0(folder_path, "KEGG_Enriched_Pathway_Presence_Protist_ver.csv"), row.names = FALSE)
