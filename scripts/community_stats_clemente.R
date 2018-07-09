setwd('/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434/all_samples_QC/mc_QC/clemente/')

summary_df=read.table(file='summary.tsv',sep='\t',header=T)

row.names(summary_df)<-summary_df$Sample
summary_df$Sample<-NULL
theoretical_df=read.table(file='/home/ajay/Desktop/minerva/sc/orga/projects/InfectiousDisease/reference-db/microbial_community_standards/jose_mc_16s_abundance.csv',sep='\t',header=T,row.names = 1)
summary_df<-summary_df[,c(rownames(theoretical_df),'total')]
summary_df<-summary_df[grep('M',rownames(summary_df)),] #for run3
summary_df$total_mapped_reads=rowSums(summary_df[,c(1:15)]) 
summary_df$primary_mapped_per<-(summary_df$total_mapped_reads/summary_df$total)*100
summary_df_per<-(summary_df[,c(1:15)]/summary_df$total_mapped_reads)*100

theoretical_df[,c(1:5)]=100*theoretical_df[,c(1:5)]
theoretical_df_t=t(theoretical_df[,c('M1','M6','M13','M14','M15')])

jose_summary_df<-rbind(summary_df_per,theoretical_df_t)
summary_df_per_t<-t(jose_summary_df)
class(summary_df_per_t)<-"numeric"
#write.table(summary_df_per_t,file='clemente_corr_plot_df.tsv',sep='\t',col.names = NA)
jose_summary_df$expt_cond<-row.names(jose_summary_df)
jose_summary_df$expt_cond<-as.factor(jose_summary_df$expt_cond)

cor_tab2<-cor(summary_df_per_t)
cor_tab2<-as.data.frame(cor_tab2)

cor_tab<-as.matrix(dist(jose_summary_df))
cor_tab<-as.data.frame(cor_tab)

cor_tab_plot_vals<-data.frame(expt_cond=c("M1","M1.1","M1.2","M1.L", "M6","M6.L", "M13","M13.L"),
                                  cor_value=c('0_1',paste(round(cor_tab['M1','M1.1'],2),round(cor_tab2['M1','M1.1'],2),sep='_'),
                                              paste(round(cor_tab['M1','M1.2'],2),round(cor_tab2['M1','M1.2'],2),sep='_'),
                                              paste(round(cor_tab['M1','M1.L'],2),round(cor_tab2['M1','M1.L'],2),sep='_'),
                                              '0_1',paste(round(cor_tab['M6','M6.L'],2),round(cor_tab2['M6','M6.L'],2),sep='_'),
                                              '0_1',paste(round(cor_tab['M13','M13.L'],2),round(cor_tab2['M13','M13.L'],2),sep='_')))





library(reshape2)
library(ggplot2)
summary_df_per_melt<-melt(jose_summary_df)
summary_df_per_melt<-merge(summary_df_per_melt,cor_tab_plot_vals,by = 'expt_cond',all.x = T)
summary_df_per_melt<-merge(summary_df_per_melt,summary_df[,c('total','primary_mapped_per')],by.x = 1,by.y = 0,all.x=T)
summary_df_per_melt$mapped_per<-round(summary_df_per_melt$primary_mapped_per, 2)
summary_df_per_melt$total<-paste(round(summary_df_per_melt$total/1000,2),'K',sep='')
summary_df_per_melt<-merge(summary_df_per_melt,theoretical_df[,c('gram','GC.genome')],by.x='variable',by.y=0,all.x=T)
summary_df_per_melt$var_fill<-paste(summary_df_per_melt$variable,'(',summary_df_per_melt$gram,')',summary_df_per_melt$GC.genome)
summary_df_per_melt$GC.genome<-as.character(summary_df_per_melt$GC.genome)
summary_df_per_melt$GC.genome<-substr(summary_df_per_melt$GC.genome, 1, nchar(summary_df_per_melt$GC.genome)-1)
summary_df_per_melt$GC.genome<-as.numeric(summary_df_per_melt$GC.genome)
summary_df_per_melt$gram<-factor(as.character(summary_df_per_melt$gram),levels=c('+','-'))
summary_df_per_melt<-summary_df_per_melt[order(summary_df_per_melt$gram,summary_df_per_melt$GC.genome),]
summary_df_per_melt$var_fill<-factor(summary_df_per_melt$var_fill,levels=as.character(unique(summary_df_per_melt$var_fill)))

summary_df_per_melt<-summary_df_per_melt[!grepl('M14|M15',summary_df_per_melt$expt_cond),]

ggplot(data=summary_df_per_melt, aes(x = summary_df_per_melt$expt_cond, y = summary_df_per_melt$value,
                                     fill=summary_df_per_melt$var_fill,label=signif(summary_df_per_melt$value, digits = 2))) +
  geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=-0.1,size=40))+
  geom_text(size = 16, position = position_stack(vjust = 0.5))+
  geom_text(data=summary_df_per_melt,aes(x=summary_df_per_melt$expt_cond,y=100,label=summary_df_per_melt$cor_value),hjust=-0.10,vjust=-0.1,angle=65,size=20)+
  labs(x='Samples',y='% Abundance',fill='Organism (Gram +/-) GC%',size=40)+
  theme(legend.text=element_text(size=40),plot.title = element_text(size=40),legend.key.size = unit(5,"line"),
        legend.title=element_text(size=50),
        axis.text=element_text(size=40),
        axis.title=element_text(size=35,face="bold"))+
  expand_limits(y = max(100 * 1.05))+
  geom_point(aes(x = summary_df_per_melt$expt_cond,y=summary_df_per_melt$primary_mapped_per),color='white',size=5,show.legend=F)
ggsave("run7_summary_mc_clemente_stacked.pdf", width = 126, height = 126,unit='cm',dpi=200)
ggsave("run7_summary_mc_clemente_stacked.png", width = 126, height = 126,unit='cm',dpi=200)

#write.table(summary_df_per_melt,file='mc_melt.tsv',sep='\t',row.names = F)

library(GGally)
library(scales)
out_data<-as.data.frame(summary_df_per_t)
colnames(out_data)[1]<-"M13_illumina_1A"
colnames(out_data)[3]<-"M1_illumina_1A"
colnames(out_data)[5]<-"M6_illumina_1A"

out_data_cle=read.table('../microany00005/clemente_corr_plot_df.tsv',sep='\t',header = T,row.names = 1)
out_data_run2=read.table('../microany00003/clemente_corr_plot_df.tsv',sep='\t',header = T,row.names = 1)
colnames(out_data_run2)<-c(paste('run2',colnames(out_data_run2),sep="_"))
out_data<-merge(out_data,out_data_cle[,c(1:5)],by=0)
out_data<-merge(out_data,out_data_run2[,c(1:5)],by.x=1,by.y=0)

out_data$color<-out_data$Row.names
out_data$Row.names<-NULL

out_data<-merge(out_data,theoretical_df[,c(6,7)],by.x='color',by.y=0)
out_data$GC.genome<-as.character(out_data$GC.genome)
out_data$GC.genome<-substr(out_data$GC.genome, 1, nchar(out_data$GC.genome)-1)
out_data$GC.genome<-as.numeric(out_data$GC.genome)
out_data$gram<-factor(as.character(out_data$gram),levels=c('+','-'))
out_data$color<-factor(out_data$color,levels = out_data$color)
out_data<-out_data[order(out_data$gram,out_data$GC.genome),]
out_data$color<-as.character(out_data$color)
out_data$color<-factor(out_data$color,levels=unique(out_data$color))

#write.table(out_data,file='corr_plot.tsv',sep='\t',row.names = F)

combo_plot<-function(p,p1,p2){
  g2 <- ggplotGrob(p2)
  colors <- g2$grobs[[6]]$children[[3]]$gp$fill
  # Change background color to tiles in the upper triangular matrix of plots 
  idx <- 1
  for (k1 in 1:(p-1)) {
    for (k2 in (k1+1):p) {
      plt <- getPlot(p1,k1,k2) +
        theme(panel.background = element_rect(fill = colors[idx], color="white"),
              panel.grid.major = element_line(color=colors[idx])
        )
      p1 <- putPlot(p1,plt,k1,k2)
      idx <- idx+1
    }
  }
  return(p1)
}

#Run 2
out_data[2:16]<-log2(out_data[2:16])
col_sel<-c(5,15,7)
col_sel<-colnames(out_data)[grepl('M15',colnames(out_data))]
col_sel[2:3]<-col_sel[3:2]


#RUn 4_5
out_data[2:22]<-log2(out_data[2:22])

col_sel<-colnames(out_data)[grepl('M6',colnames(out_data))]
col_sel=c("M6_illumina_1A","M6_IDT_1A","run2_M6_1A_1_1ng","cle_M6_1A","M6")
col_sel<-c(4,5,21,16,8)

out_data[2:21]=log2(out_data[2:21])
col_sel<-colnames(out_data)[grepl('M6',colnames(out_data)) & !grepl('cle',colnames(out_data))]
col_sel<-c(2,3,4,20,7)
col_sel[2:3]<-col_sel[3:2]
p1<-ggpairs(out_data,axisLabels='internal',lower=list(mapping = aes(colour = color,size=GC.genome,shape=gram)),columns = col_sel)
p2<-ggcorr(out_data[,c(col_sel)],label_round = 2,label = T,label_color = "black")
p<-5
print(combo_plot(p,p1,p2))
ggsave("run7_summary_mc_M1_corr_plot_dna_conc.pdf", width = 50, height = 30,unit='cm',dpi=200)
ggsave("run7_summary_mc_M1_corr_plot_dna_conc.png", width = 50, height = 30,unit='cm',dpi=200)

