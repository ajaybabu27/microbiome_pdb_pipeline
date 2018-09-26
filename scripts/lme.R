#Linear Mixed Effect Modeling
library(lme4)
library(MuMIn)
library(nlme)
library(ggplot2)
library(htmlwidgets)
library(qiime2R)
library(tidyverse)
library(reshape2)

setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00007')

set.seed(100)

mapping_file<-read.table(file='mapping_final_combined.tsv',sep='\t',header = T,comment.char = '')

df_R2_full_raw<-NULL
df_R2_reduced_raw<-NULL

#Shannon metrics (raw)
file_ls <- as.character(unzip("2_core-metrics-results/shannon_vector.qza", list = TRUE)$Name)
shannon_stats=unz(filename = file_ls[grep("alpha-diversity",file_ls)], description = "2_core-metrics-results/shannon_vector.qza")
shannon_stats_df_raw=read.table(shannon_stats,sep='\t',header = T)

mapping_file_shannon<-merge(mapping_file,shannon_stats_df_raw,by=1)
mapping_file_shannon[is.na(mapping_file_shannon$eRAP_ID),'eRAP_ID']<-0
#mapping_file_shannon<-mapping_file_shannon[!is.na(mapping_file_shannon$caseERAPID),]

cdi_biom.model<-lmer(shannon ~ bacteria_reads + CDItestResult + CaseControlAnnot + ABX_admin_24hrprior_sample_collection_bool + (1|eRAP_ID), data=mapping_file_shannon, REML=F)
summary(cdi_biom.model)
full_R2f<-r.squaredGLMM(cdi_biom.model)[1]
full_R2m<-r.squaredGLMM(cdi_biom.model)[2]-r.squaredGLMM(cdi_biom.model)[1]
full_R2<-r.squaredGLMM(cdi_biom.model)[2]

df_R2_full_raw<-rbind(df_R2_full_raw,c(full_R2,full_R2f,full_R2m))
colnames(df_R2_full_raw)<-c('R2','R2-fixed','R2-mixed')

dependent="shannon ~ "
independent=c('bacteria_reads','CDItestResult','CaseControlAnnot','ABX_admin_24hrprior_sample_collection_bool')

for (y in independent){
  formula_right=paste(independent[! independent %in% y],collapse = '+')  
  formula_right=paste(formula_right,'(1|eRAP_ID)',sep='+')
  formula=as.formula(paste(dependent,formula_right,sep=''))
  
  cdi_biom.model.reduced<-lmer(formula, data=mapping_file_shannon, REML=F)
  
  reduced_R2<-r.squaredGLMM(cdi_biom.model.reduced)[2]
  cont=full_R2-reduced_R2
  pval<-anova(cdi_biom.model.reduced,cdi_biom.model)$`Pr(>Chisq)`[2]
  
  df_R2_reduced_raw<-rbind(df_R2_reduced_raw,c(y,reduced_R2,cont,pval))
  
}

colnames(df_R2_reduced_raw)<-c('cov_removed','reduced_R2','contribution','pval')

#hist(mapping_file_shannon$shannon)
#qqnorm(mapping_file_shannon$shannon)

#Shannon metric (rarefied)
file_ls <- as.character(unzip("3_alpha_diversity/alpha-rarefaction_1000.qzv", list = TRUE)$Name)
shannon_stats=unz(filename = file_ls[grep("shannon.csv",file_ls)], description = "3_alpha_diversity/alpha-rarefaction_1000.qzv")
shannon_stats_df=read.table(shannon_stats,sep=',',header = T)

#######################################################
#raw Vs Rarefaction shannon

shannon_stats_df_raw_norm<-merge(shannon_stats_df_raw,shannon_stats_df,by=1)
shannon_stats_df_raw_norm<-shannon_stats_df_raw_norm[,c('X',"eRAP_ID","caseERAPID","CaseControlAnnot","combo","SamplingDate",colnames(shannon_stats_df_raw_norm)[grepl("4000",colnames(shannon_stats_df_raw_norm))])]
shannon_stats_df_raw_norm<-shannon_stats_df_raw_norm[order(shannon_stats_df_raw_norm$caseERAPID,shannon_stats_df_raw_norm$eRAP_ID,shannon_stats_df_raw_norm$SamplingDate),]
order_rows<-order(shannon_stats_df_raw_norm$caseERAPID,shannon_stats_df_raw_norm$CaseControlAnnot,shannon_stats_df_raw_norm$eRAP_ID,shannon_stats_df_raw_norm$SamplingDate)
shannon_stats_df_raw_norm<-shannon_stats_df_raw_norm[order_rows,]
shannon_stats_df_raw_norm$X<-factor(as.character(shannon_stats_df_raw_norm$X),levels=unique(shannon_stats_df_raw_norm$X[order_rows]))
shannon_stats_df_raw_norm$caseERAPID<-as.character(shannon_stats_df_raw_norm$caseERAPID)
shannon_stats_df_raw_norm$eRAP_ID<-as.character(shannon_stats_df_raw_norm$eRAP_ID)

shannon_stats_df_raw<-merge(shannon_stats_df_raw,shannon_stats_df_raw_norm[,c('X','eRAP_ID','caseERAPID')],by=1)

shannon_stats_df_raw_norm_melt<-melt(shannon_stats_df_raw_norm)

png('lme_results/rawVSrarefiedShannonvalues.png',height = 2000,width = 2000,res=200)
ggplot()+geom_boxplot(data=shannon_stats_df_raw_norm_melt,aes(x=X, y=value))+geom_point(data=shannon_stats_df_raw,aes(x=X,y=shannon),color='red')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+facet_wrap(~caseERAPID,scales='free_x')
dev.off()


#######################################################


mapping_file_shannon_filt<-mapping_file_shannon[mapping_file_shannon$bacteria_reads<100000,]
mapping_file_shannon_filt$bins<-NA
ggplot(data = mapping_file_shannon_filt,aes(x=bacteria_reads))+geom_histogram(colour='red')
mapping_file_shannon_filt[mapping_file_shannon_filt$bacteria_reads<100000,'bins']<-'high'
mapping_file_shannon_filt[mapping_file_shannon_filt$bacteria_reads<40000,'bins']<-'low'
ggplot()+geom_boxplot(data=mapping_file_shannon_filt,aes(x=bins,y=shannon))
ggplot(data=mapping_file_shannon_filt,aes(x=bacteria_reads,y=shannon))+geom_point(aes(colour=bins))
cor(mapping_file_shannon_filt$bacteria_reads,mapping_file_shannon_filt$shannon)
table(mapping_file_shannon_filt$bins)
mapping_file_shannon_filt$bins<-as.factor(mapping_file_shannon_filt$bins)
wilcox.test(shannon ~ bins, data = mapping_file_shannon_filt) 


#uterate through rarefied shannon values
df_R2_full<-NULL
df_R2_reduced<-NULL
for (x in 1:1000){
  shannon_stats_df_subset=shannon_stats_df[,c('sample.id',paste('depth.4000_iter.',x,sep=''))]
  colnames(shannon_stats_df_subset)<-c('sample.id','shannon')
  
  mapping_file_shannon<-merge(mapping_file,shannon_stats_df_subset,by=1)
  mapping_file_shannon[is.na(mapping_file_shannon$eRAP_ID),'eRAP_ID']<-0
  mapping_file_shannon<-mapping_file_shannon[!is.na(mapping_file_shannon$caseERAPID),]
  
  
  cdi_biom.model<-lmer(shannon ~ bacteria_reads + CDItestResult + CaseControlAnnot + ABX_admin_24hrprior_sample_collection_bool + (1|eRAP_ID), data=mapping_file_shannon, REML=F)
  summary(cdi_biom.model)
  full_R2f<-r.squaredGLMM(cdi_biom.model)[1]
  full_R2m<-r.squaredGLMM(cdi_biom.model)[2]-r.squaredGLMM(cdi_biom.model)[1]
  full_R2<-r.squaredGLMM(cdi_biom.model)[2]
  
  df_R2_full<-rbind(df_R2_full,c(x,full_R2,full_R2f,full_R2m))
  
  
  dependent="shannon ~ "
  independent=c('bacteria_reads','CDItestResult','CaseControlAnnot','ABX_admin_24hrprior_sample_collection_bool')
  
  for (y in independent){
    formula_right=paste(independent[! independent %in% y],collapse = '+')  
    formula_right=paste(formula_right,'(1|eRAP_ID)',sep='+')
    formula=as.formula(paste(dependent,formula_right,sep=''))
    
    cdi_biom.model.reduced<-lmer(formula, data=mapping_file_shannon, REML=F)
    
    reduced_R2<-r.squaredGLMM(cdi_biom.model.reduced)[2]
    cont=full_R2-reduced_R2
    pval<-anova(cdi_biom.model.reduced,cdi_biom.model)$`Pr(>Chisq)`[2]
    
    df_R2_reduced<-rbind(df_R2_reduced,c(x,y,reduced_R2,cont,pval))
    
  }
}

colnames(df_R2_full)<-c('iter','R2','R2-fixed','R2-mixed')
colnames(df_R2_reduced)<-c('iter','cov_removed','reduced_R2','contribution','pval')
df_R2_reduced<-as.data.frame(df_R2_reduced)
df_R2_reduced$reduced_R2<-as.numeric(as.character(df_R2_reduced$reduced_R2))
df_R2_reduced$contribution<-as.numeric(as.character(df_R2_reduced$contribution))
df_R2_reduced$pval<-as.numeric(as.character(df_R2_reduced$pval))

df_R2_full_melt<-melt(df_R2_full[,c(2:4)])
df_R2_full_melt$value<-round(100*df_R2_full_melt$value,2)
df_R2_full_raw_melt<-melt(df_R2_full_raw)
df_R2_full_raw_melt$value<-round(100*df_R2_full_raw_melt$value,2)

#Plot overall R2 values
png('lme_results/R2overall_stats.png',width = 1000,height = 1500,res = 200)
ggplot()+geom_boxplot(data=df_R2_full_melt,aes(x=Var2,y=value))+geom_point(data=df_R2_full_raw_melt,aes(x=Var2,y=value),color='red')+
  geom_text(data=df_R2_full_raw_melt,aes(x=Var2,y=value,label=value), hjust=0, vjust=-1)+labs(x='R2 Breakdown',y='%')
dev.off()

df_R2_reduced_melt<-melt(df_R2_reduced[,c(2:4)])
df_R2_reduced_melt$value<-round(100*df_R2_reduced_melt$value,2)
df_R2_reduced_raw<-as.data.frame(df_R2_reduced_raw)
df_R2_reduced_raw$reduced_R2<-as.numeric(as.character(df_R2_reduced_raw$reduced_R2))
df_R2_reduced_raw$contribution<-as.numeric(as.character(df_R2_reduced_raw$contribution))
df_R2_reduced_raw$pval<-as.numeric(as.character(df_R2_reduced_raw$pval))
df_R2_reduced_raw_melt<-melt(df_R2_reduced_raw[,c(1:3)])
df_R2_reduced_raw_melt$value<-round(100*df_R2_reduced_raw_melt$value,2)
df_R2_reduced_raw$pval<-round(df_R2_reduced_raw$pval,2)

#Plot overall reduced R2 distribution
png('lme_results/reduced_R2_stats.png',width = 1000,height = 1500,res = 200)
ggplot()+geom_boxplot(data=df_R2_reduced_melt,aes(x=variable,y=value))+ facet_wrap(~ cov_removed, ncol=1)+
  geom_point(data=df_R2_reduced_raw_melt,aes(x=variable,y=value),color='red')+
  geom_text(data=df_R2_reduced_raw_melt,aes(x=variable,y=value,label=value),hjust=0, vjust=1)+labs(y='%')
dev.off()

#Plot ANOVA p-values
png('lme_results/reduced_modelVSfull_model_pvals.png',width = 2000,height = 1500,res = 200)
ggplot()+geom_boxplot(data=df_R2_reduced[,c(2,5)],aes(x=cov_removed,y=pval))+geom_hline(yintercept = 0.05,color='red')+
  geom_point(data=df_R2_reduced_raw[,c(1,4)],aes(x=cov_removed,y=pval),color='red')+
  geom_text(data=df_R2_reduced_raw[,c(1,4)],aes(x=cov_removed,y=pval,label=pval),hjust=0, vjust=2)
dev.off()

########################################################################
###The following code Will be moved to seperate script once finalized###
########################################################################

#phyloloseq Analysis

setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00007')

SVs<-read_qza("0_merged_OTU/table.qza")
otu_table<-as.data.frame(SVs$data)


sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])
feature_count<-as.data.frame(sort(rowSums(otu_table)))
feature_count2<-feature_count[feature_count$`sort(rowSums(otu_table))`<10,]
hist(feature_count2)
p<-ggplot(data=feature_count,aes(x=feature_count$`sort(rowSums(otu_table))`))+geom_histogram(bins = 200)
ggplotly(p)
sample_count<-as.data.frame(sort(colSums(otu_table)))
ggplot(sample_count,aes(x=sample_count$`sort(colSums(otu_table))`))+geom_histogram(bins = 100)
quantile(sample_count$`sort(colSums(otu_table))`)
mean(sample_count$`sort(colSums(otu_table))`)

file_ls <- as.character(unzip("/home/ajay/Desktop/minerva/sc/orga/scratch/kumara22/test_microbiome_analysis/0_merged_OTU/table.qza", list = TRUE)$Name)
biom_file=unz(filename = paste(file_ls[1],'/data/feature-table.biom',sep=''),description = "table.qza")


library("phyloseq")

SVs<-read_qza("0_merged_OTU/table.qza")
biom_tsv<-as.data.frame(SVs$data)

taxon<-read_qza("4_taxons/taxonomy.qza")
taxonomy_tsv<-taxon$data
taxonomy_tsv$Confidence<-NULL
colnames(taxonomy_tsv)<-c('Feature.ID','taxonomy')

biom_tsv=merge(biom_tsv,taxonomy_tsv,by.x = 0, by.y = 1)

write.table(biom_tsv,file = '0_merged_OTU/merged_taxon_biome.tsv',sep='\t',row.names = F)

#bash command
#biom convert -i merged_taxon_biome.tsv -o merged_taxon_biome.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

biomot = import_biom('0_merged_OTU/merged_taxon_biome.biom')
metadata = import_qiime_sample_data('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00007/mapping_final_combined.tsv')
phylo_object = merge_phyloseq(biomot, metadata)
phylo_object@sam_data$eRAP_ID<-as.factor(as.character(phylo_object@sam_data$eRAP_ID))

#Plot Different species richness metrics (various alpha diversity metrics)
p<-plot_richness(phylo_object, x = "combo",color = "combo", measures="Shannon")
p+geom_boxplot(data = p$data, aes(x = combo, y = value,color = combo), alpha = 0.1)

tax_df<-as.data.frame(phylo_object@tax_table)
tax_glom(phylo_object,taxrank = "Rank7")


phylo_object_filt<-filter_taxa(phylo_object, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


total = median(sample_sums(phylo_object_filt))
standf = function(x, t=total) round(t * (x / sum(x)))
phylo_object_normalized = transform_sample_counts(phylo_object_filt, standf)

otu_table<-phylo_object_normalized@otu_table
sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])


phylo_object_CV_filtered = filter_taxa(phylo_object_normalized, function(x) sd(x)/mean(x) > 3.0, TRUE)

otu_table<-phylo_object_CV_filtered@otu_table
sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])

phylo_object.ord <- ordinate(phylo_object, "NMDS", "bray", weighted=TRUE)
p2 = plot_ordination(phylo_object, phylo_object.ord, type="samples", color="eRAP_ID", shape="combo")+
  stat_ellipse()
p2

phylo_object_normalized.ord <- ordinate(phylo_object_normalized, "NMDS", "bray", weighted=TRUE)
p2 = plot_ordination(phylo_object_normalized, phylo_object_normalized.ord, type="samples", color="eRAP_ID", shape="combo")+
  stat_ellipse()
p2

p2 = plot_ordination(phylo_object, phylo_object.ord, type="samples", color="combo", shape="ABX_admin_24hrprior_sample_collection_bool")+
  stat_ellipse()+ scale_color_manual(values=c("blue", "red", "orange",'dark green'))
#p3<-p2 + geom_polygon(aes(fill=combo)) + geom_point(size=5) + ggtitle("samples")
png('nmds_clusters.png',height = 2000,width = 3000,res=300)
print(p2)
dev.off()

p4 = plot_ordination(phylo_object_normalized, phylo_object_normalized.ord, type="split", color="combo", shape="ABX_admin_24hrprior_sample_collection_bool", label="combo", title="split") 
p4

ggplotly(p2)

plotly.offline.plot(p3, filename='name.html')

saveWidget(
  widget=as_widget(p3),   
  'plot.tmp',
  selfcontained=F);

p<-plot_richness(phylo_object, x = "combo",color = "combo", measures="Shannon")
p+geom_boxplot(data = p$data, aes(x = combo, y = value,color = combo), alpha = 0.1)


table(phylo_object@sam_data$combo)

#phy=qza_to_phyloseq("0_merged_OTU/table.qza", "1_aligned_OTU/rooted-tree.qza", "4_taxons/taxonomy.qza","mapping_final_combined.tsv")