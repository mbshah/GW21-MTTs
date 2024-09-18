
#if you wanto use KO table

p_values_t_no_na <- p_values_deseq_kegg[!is.na(p_values_deseq_kegg$interaction_padj), ]
significant_kos_pval <- rownames(p_values_t_no_na[p_values_t_no_na$treatment_padj < 0.05, ])
kegg_table_filtered<-kegg_table[significant_kos_pval,]
#kegg_table_filtered$C1_S2<-NULL

#if you want to use pathway table
p_values_t_no_na <- p_values_deseq_pathways[!is.na(p_values_deseq_pathways$treatment_padj), ]
significant_pathways_pval <- rownames(p_values_t_no_na[p_values_t_no_na$treatment_padj < 0.05, ])
kegg_pathway_table_filtered<-kegg_pathway_table[significant_pathways_pval,]
#
#cluster_table<-as.data.frame(read.csv(pathway_cluster_file,sep="\t",row.names = 1))
#cluster_table_no_na <- cluster_table[!is.na(cluster_table$Sig_cluster), ]
#significant_pathways_clus <- rownames(cluster_table_no_na[cluster_table_no_na$Sig_cluster == 3, ])


counts<-round(kegg_pathway_table)





#CHOOSE NORMALIZATION METHOD HERE
data.hel.norm<-deseq_norm(counts)
data.hel.norm.t<-t(data.hel.norm)
data.hel.norm.t<-data.hel.norm.t[order(row.names(data.hel.norm.t)),]

#calcualting and visualizing PCA
library(factoextra)
library(ggbiplot)
library(plotly)
data.hel.norm.t.pca<-prcomp(data.hel.norm.t,center = T)

#assigning groups to sorted datafirst 21 rows control last 21 treatment
groups.biplot <- c(rep("Control", 21), rep("Wastewater",21))

##Viasualization of PCA bloc, uncomment below to view
fviz_eig(data.hel.norm.t.pca) #plot to visualize and decide which pca axis to use
##visualize PCA
ggbiplot(data.hel.norm.t.pca, labels=rownames(data.hel.norm.t),
         ellipse.prob = 0.9,choices=c(1,3),var.axes=FALSE,
         ellipse=TRUE, groups=groups.biplot)
eig_val<-get_eig(data.hel.norm.t.pca)
row.names(eig_val)<-colnames(data.hel.norm.t.pca[["x"]])
hd_pca_components<-cbind(data.frame(data.hel.norm.t.pca[["x"]]),groups.biplot)


ax12<-ggplot(hd_pca_components, aes(x=PC1,y=PC2, color=groups.biplot,label=row.names(hd_pca_components)))+
  geom_text()+#geom_smooth(method = "gam", se = FALSE)+#stat_ellipse()+
  theme(legend.position = "none")+labs()
ax14<-ggplot(hd_pca_components, aes(x=PC1,y=PC4, color=groups.biplot,label=row.names(hd_pca_components)))+
  geom_text()+#geom_smooth(method = "gam", se = FALSE)+#stat_ellipse()+
  theme(legend.position = "bottom")
ax42<-ggplot(hd_pca_components, aes(x=PC4,y=PC2, color=groups.biplot,label=row.names(hd_pca_components)))+
  geom_text()+#geom_smooth(method = "loess", se = FALSE)+#stat_ellipse()+
  theme(legend.position = "none")
pca_3_axis<-ggarrange(ax12,ax42,ax14,heights = c(2,0.7), widths=c(2,0.7),ncol=2,nrow=2)
pca_3_axis

pca_out_file_png<-"/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/graphics/PCA/3xis_pca2.png"
pca_out_file_svg<-"/home/hb0358/Disk2/Disk2/Data/MTT_mbs/Statistics_GW21/graphics/PCA/3xis_pca2.svg"
#ggsave(file=pca_out_file_svg,device="svg",plot=pca_3_axis,width=1892, height=901, units="px",dpi=-150)
#ggsave(file=pca_out_file_png,device="png",plot=pca_3_axis,width=1892, height=901, units="px",dpi=-150)


plot_ly(hd_pca_components,x=~PC4, y=~PC2, z=~PC3, color=~groups.biplot,
        colors = c('#636EFA','#EF553B')) %>%
  add_markers(size = 45,text=row.names(hd_pca_components))



#calculate distance matrix and analyze decay


library(gplots)
library(reshape2)

data_dist <- vegdist(data.hel.norm.t, method = "bray")
data_dist.df<-as.data.frame(as.matrix(data_dist))
heatmap.2(as.matrix(data_dist.df), Rowv=F, Colv="Rowv", scale='none', symm =T)

data_dist.df_prep<-cbind(row.names(data_dist.df),data_dist.df)
colnames(data_dist.df_prep)[1]<-"sample_names"
data_dist_melt<-melt(data_dist.df_prep, id= 1)
data_dist_melt$time1<-""
data_dist_melt$time2<-""
data_dist_melt$couple<-""
data_dist_melt$time_dif<-""

time_dict<-c("S1"=1,"S2"=12, "S3"=24,"S4"=48,"S5"=96,"S6"=168,"S7"=336)
#Adding the time difference to the frame
for (i in seq_len(NROW(data_dist_melt))){
  sample_time<-strsplit(data_dist_melt$sample_names[i],"_")[[1]][2]
  data_dist_melt$time1[i]<-time_dict[sample_time]
  variable_time<-strsplit(toString(data_dist_melt$variable[i]),"_")[[1]][2]
  data_dist_melt$time2[i]<-time_dict[variable_time]
  data_dist_melt$couple[i]<-paste0(substr(data_dist_melt$sample_names[i],1,1),substr(toString(data_dist_melt$variable[i]),1,1))
  data_dist_melt$time_dif[i]<-data_dist_melt$time1[i]
}
data_dist_melt$time_dif<-abs(as.numeric(data_dist_melt$time_dif))
data_dist_melt<-data_dist_melt %>% filter(sample_names != variable) #remove all self comparision wich will be 0
data_dist_melt<-data_dist_melt %>% filter(couple != "TC") #remove all rows with TC because TC=CT
data_dist_melt<-data_dist_melt %>% filter(time1==time2) #keeping only comparisions of same time point
data_dist_melt_distinct<-distinct(data_dist_melt, value,time_dif, .keep_all= TRUE) # removing non disticnt comparisions
table(data_dist_melt$couple)

ggplot(data_dist_melt_distinct, aes(x = time_dif, y = value, color=couple))+
  geom_point()+
  geom_smooth(method = "lm", se = F, formula=y ~ log(x))+
  labs(title="time_decay_MTT")+
  xlab("time in hours") +
  ylab("bray dis")

##Calculating SNR
full_name<-c("CC_1","CT_1","CC_12","CT_12","CC_24","CT_24","CC_48","CT_48",
             "CC_96","CT_96","CC_168","CT_168","CC_336","CT_336")
time_scale<-c(1,1,12,12,24,24,48,48,96,96,168,168,336,336)
pair<-rep(c("CC","CT"),7)
MTT_final<-data.frame(full_name,time_scale,pair)
MTT_final$mean<-0
MTT_final$sd<-0
MTT_final$dataset<-"MTT"
MTT_final$mean_difference<-0
#added next two to confor to guido's tables so i can use his table
#MTT_final$factor<-0

#add missing things
for (i in seq_len(NROW(MTT_final))){
  temp_df<-data_dist_melt_distinct %>% filter(couple == MTT_final$pair[i] & time_dif==MTT_final$time_scale[i])
  MTT_final$mean[i]<-as.numeric(mean(temp_df$value))
  MTT_final$sd[i]<-sd(temp_df$value)
}

for (i in seq_len(NROW(MTT_final))){
  if (MTT_final$pair[i]=="CT"){
    MTT_final$mean_difference[i]<-MTT_final[MTT_final$time_scale== MTT_final$time_scale[i] & MTT_final$pair == "CT", "mean"] - MTT_final[MTT_final$time_scale== MTT_final$time_scale[i] & MTT_final$pair == "CC", "mean"]
    MTT_final$factor[i]<-round(MTT_final[MTT_final$time_scale== MTT_final$time_scale[i] & MTT_final$pair == "CT", "mean"] / MTT_final[MTT_final$time_scale== MTT_final$time_scale[i] & MTT_final$pair == "CC", "mean"], digits=2)
  }else{
    MTT_final$mean_difference[i]<-MTT_final[MTT_final$time_scale== MTT_final$time_scale[i] & MTT_final$pair == "CC", "mean"]
  }
}

ggplot(MTT_final, aes(x = as.numeric(interaction(dataset,time_scale)), y = mean_difference, fill = pair,group=dataset)) +
  geom_bar(stat = "identity",color="white") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),labels=c("1","12","24","48","96","168","336"))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, group=pair),position=position_dodge(),
                stat="identity",width=0.5,size=0.5)+
  #ylim(0,0.1)+
  geom_text(aes(label=factor),  position = position_stack(vjust=0.5))+
  labs(title="MTT", x ="Time (Hours)", y = "Bray-Curtis dissimilarity")


#testing normality
  normality.table<- shapiro.test(data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CT" & data_dist_melt_distinct$time_dif==1])
leveneTest(c(data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CT" & data_dist_melt_distinct$time_dif==1],
                   data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CC" & data_dist_melt_distinct$time_dif==1]),
           g=as.factor(c(rep("CT",9),rep("CC",3))))
bartlett.test(c(data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CT" & data_dist_melt_distinct$time_dif==1],
                   data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CC" & data_dist_melt_distinct$time_dif==1]),
              g=as.factor(c(rep("CT",9),rep("CC",3))))
t.test(data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CC" & data_dist_melt_distinct$time_dif==1], data_dist_melt_distinct$value[data_dist_melt_distinct$couple=="CT" & data_dist_melt_distinct$time_dif==1], var.equal=T)

##CALCULAING T-TEST PVALUES FORALL TIME POINTS
time_diffs <- c(1, 12, 24, 48,96, 168, 336)

# Initialize a data frame to store the p-values
MTT_significance_p_vals <- data.frame(
  time_diff = time_diffs,
  p_value = numeric(length(time_diffs))
)

# Loop through each time difference, calculate the p-value, and store it in the data frame
for (i in 1:length(time_diffs)) {
  time_diff <- time_diffs[i]
  p_val <- t.test(
    data_dist_melt_distinct$value[data_dist_melt_distinct$couple == "CC" & data_dist_melt_distinct$time_dif == time_diff],
    data_dist_melt_distinct$value[data_dist_melt_distinct$couple == "CT" & data_dist_melt_distinct$time_dif == time_diff],
    var.equal = TRUE
  )$p.value
  MTT_significance_p_vals$p_value[i] <- p_val
}