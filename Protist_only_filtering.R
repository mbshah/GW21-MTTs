library(dplyr)
#Load the kegg ortholog to gene table
ko_genes <- read.table("/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/ko_genes.list", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(ko_genes) <- c("KEGG_Ortholog", "Gene")
ko_genes$Org <- sapply(strsplit(ko_genes$Gene, ":"), `[`, 1)

#Filter genes list to have just protist genes using kos organisms list list manually made from using the KEGG website
protist_genomes_codes_from_kegg <- c(
  "mbr", "sre", "ddi", "dpp", "dfa", "ehi", "edi", "eiv", "acan", "pfa", "pfd",
  "pfh", "prei", "pgab", "pyo", "pcb", "pbe", "pvv", "pkn", "pvx", "pcy", "pmal",
  "prel", "pcot", "tan", "tpv", "tot", "beq", "bbo", "bmic", "bbig", "cpv", "cho",
  "tgo", "bbes", "tet", "ptm", "smin", "pti", "fcy", "tps", "ngd", "aaf", "pif",
  "psoj", "spar", "ehx", "gtt", "tbr", "tbg", "tcr", "lma", "lif", "ldo", "lmi",
  "lbz", "lpan", "ngr", "tva", "gla"
)
ko_genes_filter_by_org<-ko_genes[ko_genes$Org %in% protist_genomes_codes_from_kegg, ]
ko_genes_filter_by_org$KEGG_Ortholog <- gsub("ko:", "", ko_genes_filter_by_org$KEGG_Ortholog)
required_kos_list_by_org<-unique(sort(ko_genes_filter_by_org$KEGG_Ortholog))
writeLines(required_kos_list_by_org, "/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/KOs_just_protists_by_org")


#filter genes list to have just protist kos USING UNIPROT
protist_genes1<-readLines("/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/uniprotkb_taxonomy_id_2759_kegg_IDS_sort_uniq") #load list of protist genes
ko_genes_filter_by_gene<-ko_genes[ko_genes$Gene %in% protist_genes1, ]
ko_genes_filter_by_gene$KEGG_Ortholog <- gsub("ko:", "", ko_genes_filter_by_gene$KEGG_Ortholog)
required_kos_list_by_genes<-unique(sort(ko_genes_filter_by_gene$KEGG_Ortholog))
writeLines(required_kos_list_by_genes, "/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/KOs_just_protists_by_uniprot_genes")

length(unique(sort(ko_genes_filter_by_org$Org)))
length(protist_genomes_codes_from_kegg)
length(unique(sort(ko_genes_filter_by_gene$Org)))
setdiff(unique(sort(ko_genes_filter_by_gene$Org)),unique(sort(ko_genes_filter_by_org$Org))) #in gene based but not in org based: "ag"  "ccp" "cme" "gsl"
setdiff(unique(sort(ko_genes_filter_by_org$Org)),unique(sort(ko_genes_filter_by_gene$Org))) # in org based but not in gene based: "bbes" "mbr"  "pbe"  "pcot" "pfd"  "pfh"  "smin" "sre"


#OPen and import Pathway info file
ko_pathway <- read.table("/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/ko_pathway.list", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(ko_pathway) <- c("KEGG_Ortholog", "Pathway")
ko_pathway$KEGG_Ortholog <- gsub("ko:", "", ko_pathway$KEGG_Ortholog)
ko_pathway$Pathway <- gsub("path:[A-Za-z]+", "", ko_pathway$Pathway)
# Split the Pathway column into a list of pathways for each KO
ko_pathway_aggregated <- ko_pathway %>%
  group_by(KEGG_Ortholog) %>%
  summarise(Pathways = list(unique(Pathway))) %>%
  ungroup()

# Convert the list to a comma-separated string
ko_pathway_aggregated$Pathways <- sapply(ko_pathway_aggregated$Pathways, function(x) paste(x, collapse = ", "))
ko_pathway_aggregated_filtered<-ko_pathway_aggregated[ko_pathway_aggregated$KEGG_Ortholog %in% required_kos_list_by_org,]
ko_pathway_aggregated_filtered$Pathways <- strsplit(ko_pathway_aggregated_filtered$Pathways, ", ")
ko_pathway_aggregated_filtered$Pathways <- sapply(ko_pathway_aggregated_filtered$Pathways, function(x) {
  paste0("map", x)
})


pathway_data <- read.table("/home/hb0358/Disk2/Disk2/Data/MTT_mbs/ko/pathway_substem_map.tsv", sep = "\t", header = TRUE)
pathway_data$pathway <- gsub("[A-Za-z]+", "", pathway_data$pathway)
pathway_data$pathway <- paste0("map", pathway_data$pathway)


#combine_both_tables

# Expand the list into rows
library(tidyr)
ko_pathway_expanded <- ko_pathway_aggregated_filtered %>%
  unnest(Pathways)# %>%
  # Convert to numeric for joining
  #mutate(Pathways = as.numeric(Pathways))
detach("package:tidyr", unload=TRUE, character.only=TRUE)
# Join with 'pathway_data' to add subsystem and pathway name information
ko_pathway_full <- ko_pathway_expanded %>%
  left_join(pathway_data, by = c("Pathways" = "pathway"))

# Group by KEGG Ortholog and collapse the joined data back into lists ##CREATE A TABLE WITH KEGG ORTHOLOGS AS KEY
ko_pathway_filtered_final <- ko_pathway_full %>%
  group_by(KEGG_Ortholog) %>%
  summarise(
    Pathways = list(unique(Pathways)),
    Subsystems = list(unique(subsystem)),
    Pathway_Names = list(unique(pathway_name))
  ) %>%
  ungroup()

#CREATE TABLE WITH PATHWAYS AS KEY
pathway_ko_filtered_final <- ko_pathway_full %>%
  group_by(Pathways) %>%
  summarise(
    KEGG_Orthologs = list(unique(KEGG_Ortholog)),
    Subsystems = list(unique(subsystem)),
    Pathway_Names = list(unique(pathway_name))
  ) %>%
  ungroup()
patterns <- c("Human Diseases", "Organismal")
pattern <- paste(patterns, collapse = "|")
pathway_ko_filtered_final_NH<-pathway_ko_filtered_final[!grepl(pattern, pathway_ko_filtered_final$Subsystems),] #filtering non protist pathways for sure
pathway_ko_filtered_final_NH$KEGG_Orthologs<-sapply(pathway_ko_filtered_final_NH$KEGG_Orthologs, function(x) paste(x, collapse = "; "))
pathway_ko_filtered_final_NH$Subsystems<-sapply(pathway_ko_filtered_final_NH$Subsystems, function(x) paste(x, collapse = "; "))
pathway_ko_filtered_final_NH$Pathway_Names<-sapply(pathway_ko_filtered_final_NH$Pathway_Names, function(x) paste(x, collapse = "; "))
pathway_ko_filtered_final_NH_df<-as.data.frame(pathway_ko_filtered_final_NH)
write.table(pathway_ko_filtered_final_NH_df,file = paste0(folder_path, "pathway_pedia_only_protist.tsv"), row.names = F, sep="\t")
write.table(ko_pathway_filtered_final,file = paste0(folder_path, "all_KOs_protist.tsv"), row.names = F, sep="\t")



#USING THE ABOVE DATA TO CREATE PROTIST TABLES
#Create Filtered KEGGS Profile
folder_path<-"/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/"
ko_profile_file<-paste0(folder_path,"KO_counts_HMM_version_trs_pred_0.001.tsv")
kegg_table_full<-as.data.frame(read.csv(ko_profile_file,sep="\t",row.names = 1))
kegg_table_filtered<-kegg_table_full[row.names(kegg_table_full) %in% required_kos_list_by_org,]


#Create New PAthways PRofile
samples <- colnames(kegg_table_filtered)
pathway_counts <- data.frame(matrix(ncol = length(samples), nrow = length(pathway_ko_filtered_final_NH$Pathways)))
colnames(pathway_counts) <- samples
rownames(pathway_counts) <- pathway_ko_filtered_final_NH$Pathways
for (i in 1:nrow(pathway_ko_filtered_final_NH)) {
  pathway <- pathway_ko_filtered_final_NH$Pathways[i]
  kegg_orthologs <- unlist(pathway_ko_filtered_final_NH$KEGG_Orthologs[i])

  # Sum the counts for each sample
  counts <- colSums(kegg_table_filtered[rownames(kegg_table_filtered) %in% kegg_orthologs, , drop = FALSE])

  # Add the counts to the results data frame
  pathway_counts[pathway, ] <- counts
}


#WRITING BOTH TABLES
kegg_table_filtered$KO_ID <- row.names(kegg_table_filtered) #adding coloum to have header in final file
kegg_table_filtered <- kegg_table_filtered[, c("KO_ID", setdiff(names(kegg_table_filtered), "KO_ID"))] #moving the rownames col to the begining
write.table(kegg_table_filtered, file = paste0(folder_path, "KO_counts_HMM_version_trs_pred_0.001_PROTIST.tsv"), sep = "\t", row.names = FALSE)
kegg_table_filtered$KO_ID <- NULL #removing the coloumn added earlier to be able to use this table again


pathway_counts$Pathway_ID<-row.names(pathway_counts)#adding coloum to have header in final file
pathway_counts<-pathway_counts[,c("Pathway_ID",setdiff(names(pathway_counts),"Pathway_ID"))] #moving the rownames col to the begining
write.table(pathway_counts, file = paste0(folder_path, "Pathway_counts_HMM_version_trs_pred_0.001_PROTIST.tsv"), sep = "\t", row.names = FALSE)
pathway_counts$Pathway_ID<-NULL #removing the coloumn added earlier to be able to use this table again