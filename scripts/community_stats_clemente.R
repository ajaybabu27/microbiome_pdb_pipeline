#Author: Ajay
#Date: 10/10/2017
#Description: Generate abundance charts for Clemente Lab Microbial Community libraries

args = commandArgs(trailingOnly=TRUE)

#Sets the directory containing abundance tables of Clemente MC libraries as working directory. 
setwd(args[1])

#Read abundance file
summary_df=read.table(file='summary.tsv',sep='\t',header=T)

row.names(summary_df)<-summary_df$Sample
summary_df$Sample<-NULL
theoretical_df=read.table(file='/sc/orga/projects/InfectiousDisease/reference-db/microbial_community_standards/jose_mc_16s_abundance.csv',sep='\t',header=T,row.names = 1)
summary_df<-summary_df[,c(rownames(theoretical_df),'total')]
summary_df<-summary_df[grep('M',rownames(summary_df)),] #for run3
summary_df$total_mapped_reads=rowSums(summary_df[,c(1:15)]) 
summary_df$primary_mapped_per<-(summary_df$total_mapped_reads/summary_df$total)*100
summary_df_per<-(summary_df[,c(1:15)]/summary_df$total_mapped_reads)*100 # Convert to %

theoretical_df[,c(1:5)]=100*theoretical_df[,c(1:5)]
theoretical_df_t=t(theoretical_df[,c('M1','M6','M13','M14','M15')])

jose_summary_df<-rbind(summary_df_per,theoretical_df_t)
summary_df_per_t<-t(jose_summary_df)
class(summary_df_per_t)<-"numeric"
jose_summary_df$expt_cond<-row.names(jose_summary_df)
jose_summary_df$expt_cond<-as.factor(jose_summary_df$expt_cond)

#Calculate correlation
cor_tab<-cor(summary_df_per_t)
cor_tab<-as.data.frame(cor_tab)

#Calculate euclidean distance
dist_tab<-as.matrix(dist(jose_summary_df))
dist_tab<-as.data.frame(dist_tab)

cor_tab_plot_vals<-data.frame(expt_cond=c("M1","M1.1","M1.2","M1.L", "M6","M6.L", "M13","M13.L"),
                                  cor_value=c('0_1',paste(round(dist_tab['M1','M1.1'],2),round(cor_tab['M1','M1.1'],2),sep='_'),
                                              paste(round(dist_tab['M1','M1.2'],2),round(cor_tab['M1','M1.2'],2),sep='_'),
                                              paste(round(dist_tab['M1','M1.L'],2),round(cor_tab['M1','M1.L'],2),sep='_'),
                                              '0_1',paste(round(dist_tab['M6','M6.L'],2),round(cor_tab['M6','M6.L'],2),sep='_'),
                                              '0_1',paste(round(dist_tab['M13','M13.L'],2),round(cor_tab['M13','M13.L'],2),sep='_')))

								
#Draw stacked abundance plot comparing samples with theoretical distribution. 								
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

ggsave("summary_mc_clemente_stacked.pdf", width = 126, height = 126,unit='cm',dpi=200)
ggsave("summary_mc_clemente_stacked.png", width = 126, height = 126,unit='cm',dpi=200)

#Draw pairwise scatter plot comparing samples with theoretical distribution. 		
library(GGally)
library(scales)
out_data<-as.data.frame(summary_df_per_t)
out_data$color<-rownames(out_data)

out_data<-merge(out_data,theoretical_df[,c(6,7)],by=0)
out_data$GC.genome<-as.character(out_data$GC.genome)
out_data$GC.genome<-substr(out_data$GC.genome, 1, nchar(out_data$GC.genome)-1)
out_data$GC.genome<-as.numeric(out_data$GC.genome)
out_data$gram<-factor(as.character(out_data$gram),levels=c('+','-'))
out_data$color<-factor(out_data$color,levels = out_data$color)
out_data<-out_data[order(out_data$gram,out_data$GC.genome),]
out_data$color<-as.character(out_data$color)
out_data$color<-factor(out_data$color,levels=unique(out_data$color))


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

out_data[2:11]=log2(out_data[2:11])

for (sample_group in c('M1','M6','M13')){

	if (sample_group=='M1'){
	  col_sel<-colnames(out_data)[grepl('M1\\.',colnames(out_data))]
	  col_sel<-c(col_sel,'M1')
	  
	}
	
	else{
		col_sel<-colnames(out_data)[grepl(sample_group,colnames(out_data))]
	}
	
	p1<-ggpairs(out_data,axisLabels='internal',lower=list(mapping = aes(colour = color,size=GC.genome,shape=gram)),columns = col_sel)
	p2<-ggcorr(out_data[,c(col_sel)],label_round = 2,label = T,label_color = "black")
	p<-length(col_sel)
	
	pdf(paste("summary_mc_",sample_group,"_corr_plot_dna_conc.pdf",sep=''), width=15, height=15,onefile=F)
	print(combo_plot(p,p1,p2))
	dev.off()


}

#Draw edit distance charts
library(reshape2)
library(ggplot2)
library(scales)

test <- read.table("../PhiX/PhiXQC_edit_dist.txt",sep='\t',header = T)

test_m_temp<-melt(test,na.rm = T)

xpos <- c(Inf,Inf)
ypos <- c(Inf,Inf)
mean_char=paste('Mean=',as.character(round(mean(test_m_temp$value),digits=2)))
sd_char=paste('SD=',as.character(round(sd(test_m_temp$value),digits = 2)))
annotateText <- c(mean_char,sd_char)
hjustvar<-c(1,1) 
vjustvar<-c(1,2.5)

pdf("clemente_edit_dist_perhist.pdf", onefile = TRUE,width=15, height=5)
p<-ggplot(test_m_temp, aes(x = factor(test_m_temp$value))) +  
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  scale_y_continuous(labels = percent_format()) + geom_vline(xintercept = quantile(test_m_temp$value,.90)+1,colour='red')+
  labs(x='edit_distance',y='percentage_distribution',title="PhiX Edit Distance (Pre-QC)")+
  geom_text(data=as.data.frame(annotateText),aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))+
  scale_x_discrete(limits=as.character(seq(0,10)))

print(p)
#Plot percentage histogram plots for edit distances
library(Rmisc)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)

test <- read.csv("clemente_postqc_edit_dist.txt",sep='\t',header = T)

test_m<-melt(test,na.rm = T)

list_plot=list()

for (x in levels(test_m$variable)){
  
  test_m_temp<-test_m[test_m$variable==x,]
  
  xpos <- c(Inf,Inf)
  ypos <- c(Inf,Inf)
  mean_char=paste('Mean=',as.character(round(mean(test_m_temp$value),digits=2)))
  sd_char=paste('SD=',as.character(round(sd(test_m_temp$value),digits = 2)))
  annotateText <- c(mean_char,sd_char)
  hjustvar<-c(1,1) 
  vjustvar<-c(1,2.5)
  
  
  p<-ggplot(test_m_temp, aes(x = factor(test_m_temp$value))) +  
    geom_bar(aes(y = (..count..)/sum(..count..))) + 
    scale_y_continuous(labels = percent_format()) + geom_vline(xintercept = quantile(test_m_temp$value,.90)+1,colour='red')+
    labs(x='edit_distance',y='percentage_distribution',title=x)+
    geom_text(data=as.data.frame(annotateText),aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))+
    scale_x_discrete(limits=as.character(seq(0,10)))
  
  print(p)
  
  list_plot[[x]]<-p
  
  
  
  #do.call("grid.arrange", p)  
  
}
dev.off()