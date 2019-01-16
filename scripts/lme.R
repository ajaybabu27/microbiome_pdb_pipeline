#Linear Mixed Effect Modeling
library(lme4)
library(MuMIn)
library(nlme)
library(ggplot2)
library(htmlwidgets)
library(qiime2R)
library(tidyverse)
library(reshape2)
library('GoodmanKruskal')
library(plotly)
library("phyloseq")

setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00008')

set.seed(100)

mapping_file<-read.table(file='mapping_final_combined.tsv',sep='\t',header = T,comment.char = '')

#Shannon metrics (raw)
file_ls <- as.character(unzip("2_core-metrics-results/shannon_vector.qza", list = TRUE)$Name)
shannon_stats=unz(filename = file_ls[grep("alpha-diversity",file_ls)], description = "2_core-metrics-results/shannon_vector.qza")
shannon_stats_df_raw=read.table(shannon_stats,sep='\t',header = T)

mapping_file_shannon<-merge(mapping_file,shannon_stats_df_raw,by=1)
mapping_file_shannon[is.na(mapping_file_shannon$eRAP_ID),'eRAP_ID']<-0
#mapping_file_shannon<-mapping_file_shannon[!is.na(mapping_file_shannon$caseERAPID),]
mapping_file_shannon$eRAP_ID<-as.factor(mapping_file_shannon$eRAP_ID)

#Load PCR QC
PCR_qc<-read.table('pcr_qc.csv',sep=',',header = T)
mapping_file_shannon<-merge(mapping_file_shannon,PCR_qc,by=1)


#Assess Multicollinearity

library(GGally)

mapping_file_shannon$hc_bool<-F
mapping_file_shannon[grepl('HC',mapping_file_shannon$X.SampleID),'hc_bool']<-T

mapping_file_shannon_sub<-mapping_file_shannon[mapping_file_shannon$bacteria_reads<100000,]

ggpairs(mapping_file_shannon_sub[,c('bacteria_reads','input','PCR2','shannon')],lower = list(mapping = aes(color = mapping_file_shannon_sub$hc_bool)))

ggpairs(mapping_file_shannon_sub[,c('bacteria_reads','input','PCR2','shannon')],lower = list(mapping = aes(color = mapping_file_shannon_sub$hc_bool)))

independent=c('bacteria_reads','input','CDItestResult','CaseControlAnnot','ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool', 
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'met_int_admin_1wprior_sample_collection_bool',
              'Systematic_Plate_ID')


library(usdm)
predictor_df<-mapping_file_shannon[,independent]
vif_raw=as.data.frame(vif(predictor_df))
vif_raw$VIF<-round(vif_raw$VIF,2)
write.table(file='vif_raw.tsv',vif_raw,sep='\t',row.names = F)


independent=c('CDItestResult','CaseControlAnnot','ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool', 
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'met_int_admin_1wprior_sample_collection_bool',
              'Systematic_Plate_ID','illumina_runid')

predictor_df<-mapping_file_shannon[,independent]
GKmatrix1 <- GKtauDataframe(predictor_df)
plot(GKmatrix1)

independent=c('ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool', 
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'met_int_admin_1wprior_sample_collection_bool',
              'Systematic_Plate_ID','illumina_runid')

predictor_df<-mapping_file_shannon[,independent]
GKmatrix1 <- GKtauDataframe(predictor_df)
png('categorical_removed',height = 5000,width = 7000,res=600)
plot(GKmatrix1)
dev.off()


independent=c('bacteria_reads','input','ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool', 
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'met_int_admin_1wprior_sample_collection_bool',
              'Systematic_Plate_ID')

predictor_df<-mapping_file_shannon[,independent]
vif_filt=as.data.frame(vif(predictor_df))
vif_filt$VIF<-round(vif_filt$VIF,2)
write.table(file='vif_filt.tsv',vif_filt,sep='\t',row.names = F)
#mapping_file_shannon$bacteria_reads_cor<-mapping_file_shannon$bacteria_reads*1.5

cdi_biom.model<-lmer(shannon ~ bacteria_reads+CDItestResult + CaseControlAnnot + 
                       ABX_admin_1wprior_sample_collection_bool + input+PCR2+
                       vanco_oral_admin_1wprior_sample_collection_bool+
                       vanco_int_admin_1wprior_sample_collection_bool+
                       met_oral_admin_1wprior_sample_collection_bool+
                       met_int_admin_1wprior_sample_collection_bool+
                       (1|Systematic_Plate_ID)+
                       (1|eRAP_ID), data=mapping_file_shannon, REML=F)


df_R2_full_raw<-NULL
df_R2_reduced_raw<-NULL
cdi_biom.model<-lmer(shannon ~ bacteria_reads+
                       ABX_admin_24hrprior_sample_collection_bool +
                       ABX_admin_1wprior_sample_collection_bool +
                       vanco_oral_admin_1wprior_sample_collection_bool+
                       met_oral_admin_1wprior_sample_collection_bool+
                       vanco_int_admin_1wprior_sample_collection_bool+
                       met_int_admin_1wprior_sample_collection_bool+
                       (1|Systematic_Plate_ID)+
                       (1|eRAP_ID), data=mapping_file_shannon, REML=F)




full_R2f<-r.squaredGLMM(cdi_biom.model)[1]
full_R2m<-r.squaredGLMM(cdi_biom.model)[2]-r.squaredGLMM(cdi_biom.model)[1]
full_R2<-r.squaredGLMM(cdi_biom.model)[2]

df_R2_full_raw<-rbind(df_R2_full_raw,c(full_R2,full_R2f,full_R2m))
colnames(df_R2_full_raw)<-c('R2','R2-fixed','R2-mixed')

dependent="shannon ~ "
independent=c('bacteria_reads',
              'ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool',
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
                'met_int_admin_1wprior_sample_collection_bool',
              'eRAP_ID','Systematic_Plate_ID')

for (y in independent){
  print(y)
  print(table(mapping_file_shannon[,y]))
  formula_right=paste(independent[! independent %in% c(y,'Systematic_Plate_ID',
                                                       'eRAP_ID')],collapse = '+')  
  if(y=='Systematic_Plate_ID'){
    formula_right=paste(formula_right,'(1|eRAP_ID)',sep='+')
  }
  else if(y=='eRAP_ID'){
    formula_right=paste(formula_right,'(1|Systematic_Plate_ID)',sep='+')
  }
  else{
    formula_right=paste(formula_right,'(1|Systematic_Plate_ID)','(1|eRAP_ID)',sep='+')
  }
  formula=as.formula(paste(dependent,formula_right,sep=''))
  #print(formula)
  cdi_biom.model.reduced<-lmer(formula, data=mapping_file_shannon, REML=F)
  
  reduced_R2<-r.squaredGLMM(cdi_biom.model.reduced)[2]
  cont=full_R2-reduced_R2
  pval<-anova(cdi_biom.model.reduced,cdi_biom.model)$`Pr(>Chisq)`[2]
  
  df_R2_reduced_raw<-rbind(df_R2_reduced_raw,c(y,reduced_R2,cont,pval))
  
}

colnames(df_R2_reduced_raw)<-c('cov_removed','reduced_R2','contribution','pval')

df_R2_reduced_raw<-as.data.frame(df_R2_reduced_raw)
df_R2_reduced_raw$reduced_R2<-as.numeric(as.character(df_R2_reduced_raw$reduced_R2))
df_R2_reduced_raw$contribution<-as.numeric(as.character(df_R2_reduced_raw$contribution))
df_R2_reduced_raw$pval<-as.numeric(as.character(df_R2_reduced_raw$pval))
df_R2_reduced_raw$pval<-round(df_R2_reduced_raw$pval,2)
df_R2_reduced_raw[,c(2:3)]<-round(100*df_R2_reduced_raw[,c(2:3)],2)

write.table(file='loo.tsv',df_R2_reduced_raw,row.names = F,sep='\t')

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

vif.mer(cdi_biom.model)


summary(cdi_biom.model)
library(car)
vif(cdi_biom.model)

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


#iterate through rarefied shannon values
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
###              OTU based mixed effect lm                           ###
########################################################################

tax_df_filt<-as.data.frame(phylo_object_phylum@otu_table)

tax_table

tax_df_filt<-merge(tax_df_filt,phylo_object_phylum@tax_table[,c('Rank2')],by=0)
rownames(tax_df_filt)<-tax_df_filt$Rank2
tax_df_filt[,c('Row.names','Rank2')]<-NULL

taxa_joined<-as.data.frame(apply(phylo_object@tax_table, 1, paste, collapse="_"))
colnames(taxa_joined)[1]<-'full_taxa'
tax_df_filt<-merge(tax_df_filt,taxa_joined,by=0)
tax_df_filt$Row.names<-NULL
tax_df_filt<-aggregate(. ~ full_taxa, tax_df_filt, sum)
rownames(tax_df_filt)<-tax_df_filt$full_taxa
tax_df_filt$full_taxa<-NULL

rownames(tax_df_filt)<-gsub("\\[", "", rownames(tax_df_filt))
rownames(tax_df_filt)<-gsub("\\]", "", rownames(tax_df_filt))
rownames(tax_df_filt)<-gsub("-", "_", rownames(tax_df_filt))

tax_df_filt_t<-t(tax_df_filt)

mapping_file_shannon_otus<-merge(mapping_file_shannon,tax_df_filt_t,by.x=1,by.y=0)

df_R2_full_otu_normalized<-NULL
df_R2_reduced_otu_normalized<-NULL

for (dependent_var in rownames(tax_df_filt)){
  
  dependent=paste(dependent_var," ~ ",sep='')
  independent=c('bacteria_reads',
                'ABX_admin_24hrprior_sample_collection_bool',
                'ABX_admin_1wprior_sample_collection_bool',
                'vanco_oral_admin_1wprior_sample_collection_bool',
                'met_oral_admin_1wprior_sample_collection_bool',
                'vanco_int_admin_1wprior_sample_collection_bool',
                'met_int_admin_1wprior_sample_collection_bool',
                'eRAP_ID','Systematic_Plate_ID')
  formula_right=paste(independent[! independent %in% c('Systematic_Plate_ID',
                                                       'eRAP_ID')],collapse = '+')
  
  formula_right=paste(formula_right,'(1|Systematic_Plate_ID)','(1|eRAP_ID)',sep='+')
  formula=as.formula(paste(dependent,formula_right,sep=''))
  
  #print(formula)
  cdi_biom.model<-lmer(formula, data=mapping_file_shannon_otus, REML=F)
  
  
  
  
  full_R2f<-r.squaredGLMM(cdi_biom.model)[1]
  full_R2m<-r.squaredGLMM(cdi_biom.model)[2]-r.squaredGLMM(cdi_biom.model)[1]
  full_R2<-r.squaredGLMM(cdi_biom.model)[2]
  
  df_R2_full_otu_normalized<-rbind(df_R2_full_otu_normalized,c(dependent_var,full_R2,full_R2f,full_R2m))
 
  for (y in independent){
    print(y)
    print(table(mapping_file_shannon_otus[,y]))
    formula_right=paste(independent[! independent %in% c(y,'Systematic_Plate_ID',
                                                         'eRAP_ID')],collapse = '+')  
    if(y=='Systematic_Plate_ID'){
      formula_right=paste(formula_right,'(1|eRAP_ID)',sep='+')
    }
    else if(y=='eRAP_ID'){
      formula_right=paste(formula_right,'(1|Systematic_Plate_ID)',sep='+')
    }
    else{
      formula_right=paste(formula_right,'(1|Systematic_Plate_ID)','(1|eRAP_ID)',sep='+')
    }
    formula=as.formula(paste(dependent,formula_right,sep=''))
    #print(formula)
    cdi_biom.model.reduced<-lmer(formula, data=mapping_file_shannon_otus, REML=F)
    
    reduced_R2<-r.squaredGLMM(cdi_biom.model.reduced)[2]
    cont=full_R2-reduced_R2
    pval<-anova(cdi_biom.model.reduced,cdi_biom.model)$`Pr(>Chisq)`[2]
    
    df_R2_reduced_otu_normalized<-rbind(df_R2_reduced_otu_normalized,c(dependent_var,y,reduced_R2,cont,pval))
    
  }
  
  
  
  
}
colnames(df_R2_full_otu_normalized)<-c('Dependent_Var','R2','R2-fixed','R2-mixed')

df_R2_full_otu_normalized<-as.data.frame(df_R2_full_otu_normalized)
df_R2_full_otu_normalized[,c(2:4)]<-apply(df_R2_full_otu_normalized[,c(2:4)],2,as.character)
df_R2_full_otu_normalized[,c(2:4)]<-apply(df_R2_full_otu_normalized[,c(2:4)],2,as.numeric)
df_R2_full_otu_normalized[,c(2:4)]<-round(100*df_R2_full_otu_normalized[,c(2:4)],2)

df_R2_reduced_otu_normalized<-as.data.frame(df_R2_reduced_otu_normalized)
df_R2_reduced_otu_normalized[,c(3:4)]<-apply(df_R2_reduced_otu_normalized[,c(3:4)],2,as.character)
df_R2_reduced_otu_normalized[,c(3:4)]<-apply(df_R2_reduced_otu_normalized[,c(3:4)],2,as.numeric)
df_R2_reduced_otu_normalized[,c(3:4)]<-round(100*df_R2_reduced_otu_normalized[,c(3:4)],2)
  
colnames(df_R2_reduced_otu_normalized)<-c('Dependent_Var','covariate','R2-reduced%','Cont_per','pval')

png('R2_OTU_dist_norm_phylum_raw.png',width=4000,height=2000,res=300)
a<-ggplot()+
  geom_violin(data=df_R2_reduced_otu_normalized,aes(covariate,Cont_per))+
geom_jitter(data=df_R2_reduced_otu_normalized,aes(covariate,Cont_per,color=df_R2_reduced_otu_normalized$Dependent_Var), position = position_jitter(width = 0.05))+
  geom_jitter(data=df_R2_reduced_raw,aes(cov_removed,contribution), position = position_jitter(width = 0.05),shape=10,size=5)+
        theme(axis.text.x = element_text(angle = 65, hjust = 1))+ theme(legend.position="none")#+
  #coord_cartesian(ylim = c(0, 25))
print(a)
dev.off()
df_R2_reduced_otu_normalized$firmicute_bool<-df_R2_reduced_otu_normalized$Dependent_Var%in%c("Firmicutes")
a<-ggplot()+geom_jitter(data=df_R2_reduced_otu_normalized,aes(covariate,Cont_per,color=df_R2_reduced_otu_normalized$Dependent_Var,shape=df_R2_reduced_otu_normalized$firmicute_bool), position = position_jitter(width = 0.05))+
  geom_jitter(data=df_R2_reduced_raw,aes(cov_removed,contribution), position = position_jitter(width = 0.05),shape=10,size=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="none")
#coord_cartesian(ylim = c(0, 25)) 

b<-ggplotly(a)
htmlwidgets::saveWidget(b, "all_otus_variance_partition_normalized.html")

########################################################################
###The following code Will be moved to seperate script once finalized###
########################################################################



setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00009')

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

#file_ls <- as.character(unzip("/home/ajay/Desktop/minerva/sc/orga/scratch/kumara22/test_microbiome_analysis/0_merged_OTU/table.qza", list = TRUE)$Name)
#biom_file=unz(filename = paste(file_ls[1],'/data/feature-table.biom',sep=''),description = "table.qza")

#phyloloseq Analysis
setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00009')

#merging mapping files
mapping_file1<-read.table(file='/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434/all_samples_QC/mapping_final.tsv',sep='\t',header = T,comment.char = '')
mapping_file1<-mapping_file1[mapping_file1$CaseControlAnnot!="healthy control" & mapping_file1$SampleType=="Stool",]
mapping_file1$illumina_runid='H434'

mapping_file2<-read.table(file='/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00395/all_samples_QC/mapping_final.tsv',sep='\t',header = T,comment.char = '')
mapping_file2<-mapping_file2[mapping_file2$CaseControlAnnot!="healthy control" & mapping_file2$SampleType=="Stool",]

mapping_file3<-read.table(file='/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00396/all_samples_QC/mapping_final.tsv',sep='\t',header = T,comment.char = '')
mapping_file3<-mapping_file3[mapping_file3$CaseControlAnnot!="healthy control" & mapping_file3$SampleType=="Stool",]

mapping_file_merged<-rbind(mapping_file1,mapping_file2,mapping_file3)
mapping_file_merged<-mapping_file_merged[mapping_file_merged$bacteria_reads>500,]

write.table(mapping_file_merged,file='/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00008/0_merged_OTU/mapping_final_combined.tsv',
            row.names=F,sep='\t',quote = F)

mapping_file_merged<-read.table('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00009/mapping_final_combined.tsv',
                                sep='\t',header=T)

SVs<-read_qza("0_merged_OTU/table.qza")
biom_tsv<-as.data.frame(SVs$data)

taxon<-read_qza("4_taxons/taxonomy.qza")
taxonomy_tsv<-taxon$data
taxonomy_tsv$Confidence<-NULL
colnames(taxonomy_tsv)<-c('Feature.ID','taxonomy')

biom_tsv=merge(biom_tsv,taxonomy_tsv,by.x = 0, by.y = 1)

write.table(biom_tsv,file = '0_merged_OTU/merged_taxon_biome.tsv',sep='\t',row.names = F,quote = F)

#bash command
#biom convert -i merged_taxon_biome.tsv -o merged_taxon_biome.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

biomot = import_biom('0_merged_OTU/merged_taxon_biome.biom')
metadata = import_qiime_sample_data('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00009/mapping_final_combined.tsv')
phylo_object = merge_phyloseq(biomot, metadata)
phylo_object@sam_data$eRAP_ID<-as.factor(as.character(phylo_object@sam_data$eRAP_ID))
phylo_object@sam_data[is.na(phylo_object@sam_data$cdi_order),'cdi_order']<-0
phylo_object@sam_data$cdi_order<-as.factor(as.character(phylo_object@sam_data$cdi_order))

phylo_object@sam_data$combo_cdi_order<-interaction(phylo_object@sam_data$combo,phylo_object@sam_data$cdi_order)

phylo_object@sam_data$combo_cdi_order<-factor(as.character(phylo_object@sam_data$combo_cdi_order),levels=c('CD Pos Case samples.1',
                                                                                                           'CD Pos Case samples.2',
                                                                                                           'CD Pos Case samples.3',
                                                                                                           'CD Pos Case samples.4',
                                                                                                           'CD Pos Case samples.5',
                                                                                                           'CD Pos Case samples.0',
                                                                                                           'CD Neg Case samples.0',
                                                                                                           'Matched Control Samples.0',
                                                                                                           'Healthy control Samples.0'))


#Plot Different species richness metrics (various alpha diversity metrics)

give.n <- function(x){
  return(c(y = -0.2, label = length(x)))
}

png('alpha_diversity_cdi_order.png',width = 2500,height = 3000,res=300)
p<-plot_richness(phylo_object, x = "combo_cdi_order",color = "combo_cdi_order", measures="Shannon")
p+geom_boxplot(data = p$data, aes(x = combo_cdi_order, y = value,color = combo_cdi_order), alpha = 0.1)+ 
  stat_summary(fun.data = give.n, geom = "text")
dev.off()

rank_names(phylo_object)

phylo_object_phylum<-tax_glom(phylo_object,taxrank = "Rank2")
tax_df<-as.data.frame(phylo_object_phylum@otu_table)

otu_table<-phylo_object@otu_table
sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])

phylo_object_filt<-filter_taxa(phylo_object, function(x) sum(x > 3) > (0.01*length(x)), TRUE)
tax_df_filt<-as.data.frame(phylo_object_filt@otu_table)
tax_df_filt_missing<-tax_df[!rownames(tax_df)%in%rownames(tax_df_filt),]
max(apply(tax_df_filt_missing,1,max))

total = median(sample_sums(phylo_object_filt))
standf = function(x, t=total) round(t * (x / sum(x)))
phylo_object_normalized = transform_sample_counts(phylo_object_filt, standf)

otu_table<-phylo_object_normalized@otu_table
sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])
phylo_object_normalized@sam_data$num_otus<-apply(otu_table,2,function(x) sum(x>10000))


phylo_object_CV_filtered = filter_taxa(phylo_object_normalized, function(x) sd(x)/mean(x) > 3.0, TRUE)

otu_table<-phylo_object_CV_filtered@otu_table
sparsity=sum(otu_table == 0)/(dim(otu_table)[1]*dim(otu_table)[2])

phylo_object.ord <- ordinate(phylo_object, "NMDS", "bray", weighted=TRUE)
p2 = plot_ordination(phylo_object, phylo_object.ord, type="samples", color="eRAP_ID", shape="combo")+
  stat_ellipse()
p2

mapping_file_shannon$bacteria_reads<-quantile(mapping_file_shannon$bacteria_reads)
phylo_object_normalized.ord <- ordinate(phylo_object_normalized, "NMDS", "bray", weighted=TRUE)
p2 = plot_ordination(phylo_object_normalized, phylo_object_normalized.ord, type="samples", color="CDItestResult", shape="CaseControlAnnot")+
  aes(text  = paste("  <i>Patient ID:</i> ",  eRAP_ID, "<br>",
                                "  <i>bacteria_reads:</i> ", bacteria_reads, "<br>",
                                "  <i>num_otus_above_1000:</i> ", num_otus, "<br>",
                                sep=""))
ggplotly(p2)
interactive_df<-as.data.frame(phylo_object_normalized.ord$points)
interactive_df<-merge(interactive_df,mapping_file_shannon[,c(independent)])

independent=c('bacteria_reads','CDItestResult','CaseControlAnnot','ABX_admin_24hrprior_sample_collection_bool',
              'ABX_admin_1wprior_sample_collection_bool', 
              'vanco_oral_admin_1wprior_sample_collection_bool',
              'vanco_int_admin_1wprior_sample_collection_bool',
              'met_oral_admin_1wprior_sample_collection_bool',
              'met_int_admin_1wprior_sample_collection_bool',
              'eRAP_ID','Systematic_Plate_ID','illumina_runid')

for(ind in independent){
print(ind)  
  if(ind=='eRAP_ID'){
    plot_ordination(phylo_object_normalized , phylo_object_normalized.ord, type="samples", color=ind,title="Raw")+
      stat_ellipse()
    ggsave(paste(ind,'_Raw_phylum.pdf',sep=''),width=10,height=8)
  }
  else if(is.factor(mapping_file_shannon[,ind]) | is.logical(mapping_file_shannon[,ind])){
  #pdf(paste(y,'.pdf',sep=''),width=10,height=8)
plot_ordination(phylo_object , phylo_object.ord, type="samples", color=ind,title="Raw")+
  stat_ellipse()+ scale_color_manual(values=c("blue", "red", "orange",'dark green','purple','brown','magenta','black'))
    ggsave(paste(ind,'_Raw_phylum.pdf',sep=''),width=10,height=8)
  }  else{
    plot_ordination(phylo_object , phylo_object.ord, type="samples", color=ind,title="Raw")
    ggsave(paste(ind,'_Raw_phylum.pdf',sep=''),width=10,height=8)
      
  }

#dev.off()
}

plot_ordination(phylo_object, phylo_object.ord, type="samples", color="cdi_order")+
  stat_ellipse()#+ scale_color_manual(values=c("blue", "red", "orange",'dark green'))


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