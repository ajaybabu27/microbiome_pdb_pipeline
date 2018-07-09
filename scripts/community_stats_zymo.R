args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

summary_df=read.table(file='summary.tsv',sep='\t',header=T)

row.names(summary_df)<-summary_df$Sample
summary_df$Sample<-NULL
theoretical_df=read.table(file='/sc/orga/projects/InfectiousDisease/reference-db/microbial_community_standards/zymo_mc_16s_abundance.csv',sep=',',header=T,row.names = 1)
summary_df<-summary_df[,c(rownames(theoretical_df),'total')]
summary_df$total_mapped_reads=rowSums(summary_df[,c(1:8)]) 
summary_df$mapped_per<-(summary_df$total_mapped_reads/summary_df$total)*100

summary_df_per<-(summary_df[,c(1:8)]/summary_df$total_mapped_reads)*100

summary_df_per<-rbind(summary_df_per,theoretical_df$zymo_theoretical)
row.names(summary_df_per)[nrow(summary_df)+1]<-'Theoretical'

summary_df_per_t<-t(summary_df_per)
summary_df_per$expt_cond<-row.names(summary_df_per)

cor_tab<-as.matrix(dist(summary_df_per))
cor_tab<-as.data.frame(cor_tab)
cor_tab<-cor_tab[nrow(summary_df)+1,]
cor_tab<-t(cor_tab)
cor_tab<-as.data.frame(cor_tab)
cor_tab$Theoretical<-round(cor_tab$Theoretical, 2)

cor_tab2<-cor(summary_df_per_t)
cor_tab2<-as.data.frame(cor_tab2)
cor_tab2<-cor_tab2[nrow(summary_df)+1,]
cor_tab2<-t(cor_tab2)
cor_tab2<-as.data.frame(cor_tab2)
cor_tab2$Theoretical<-round(cor_tab2$Theoretical, 2)

cor_tab$Theoretical<-paste(cor_tab$Theoretical,cor_tab2$Theoretical,sep="_")


library(reshape2)
library(ggplot2)
summary_df_per_melt<-melt(summary_df_per)
summary_df_per_melt<-merge(summary_df_per_melt,cor_tab,by.x = 1,by.y = 0)

summary_df_per_melt<-merge(summary_df_per_melt,summary_df[,c('total','mapped_per')],by.x = 1,by.y = 0,all.x=T)
summary_df_per_melt$mapped_per<-round(summary_df_per_melt$mapped_per, 2)
summary_df_per_melt$total<-paste(round(summary_df_per_melt$total/1000,2),'K',sep='')
summary_df_per_melt<-merge(summary_df_per_melt,theoretical_df[,c('gram','GC.genome')],by.x='variable',by.y=0,all.x=T)
summary_df_per_melt$var_fill<-paste(summary_df_per_melt$variable,'(',summary_df_per_melt$gram,')',summary_df_per_melt$GC.genome)
summary_df_per_melt$GC.genome<-as.character(summary_df_per_melt$GC.genome)
summary_df_per_melt$GC.genome<-substr(summary_df_per_melt$GC.genome, 1, nchar(summary_df_per_melt$GC.genome)-1)
summary_df_per_melt$GC.genome<-as.numeric(summary_df_per_melt$GC.genome)
summary_df_per_melt$gram<-factor(as.character(summary_df_per_melt$gram),levels=c('+','-'))
summary_df_per_melt<-summary_df_per_melt[order(summary_df_per_melt$gram,summary_df_per_melt$GC.genome),]
summary_df_per_melt$var_fill<-factor(summary_df_per_melt$var_fill,levels=as.character(unique(summary_df_per_melt$var_fill)))

ggplot(data=summary_df_per_melt, aes(x = summary_df_per_melt$expt_cond, y = summary_df_per_melt$value,
                                     fill=summary_df_per_melt$var_fill,label=signif(summary_df_per_melt$value, digits = 2))) +
  geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=-0.1,size=40))+
  geom_text(size = 16, position = position_stack(vjust = 0.5))+
  geom_text(data=summary_df_per_melt,aes(x=summary_df_per_melt$expt_cond,y=100,label=summary_df_per_melt$Theoretical),vjust=-0.05,hjust=-0.1,angle=65,size=20)+
  labs(x='Samples',y='% Abundance',fill='Organism (Gram +/-) GC%',size=40)+
  theme(legend.text=element_text(size=40),plot.title = element_text(size=40),legend.key.size = unit(5,"line"),
        legend.title=element_text(size=50),
        axis.text=element_text(size=40),
        axis.title=element_text(size=35,face="bold"))+
  expand_limits(y = max(100 * 1.06))+
   geom_point(aes(x = summary_df_per_melt$expt_cond,y=summary_df_per_melt$mapped_per),color='white',size=10,show.legend=F)
ggsave("summary_mc_zymo_stacked.pdf", width = 126, height = 126,unit='cm',dpi=200)
ggsave("summary_mc_zymo_stacked.png", width = 126, height = 126,unit='cm',dpi=200)

library(GGally)
library(scales)
out_data<-as.data.frame(summary_df_per_t)
out_data$color<-rownames(out_data)
out_data<-merge(out_data,theoretical_df[,c(2,3)],by.x='color',by.y=0)
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
      )+scale_x_continuous(breaks = seq(0, 30, by = 6)) +
      scale_y_continuous(breaks = seq(0, 30, by = 6))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}
return(p1)
}

p1<-ggpairs(out_data,axisLabels='internal',lower=list(mapping = aes(colour = color,size=GC.genome,shape=gram)),columns = c(2:4))
p2<-ggcorr(out_data[,c(2:4)],label_round = 2,label = T,label_color = "black")
p<-3

#print(combo_plot(p,p1,p2))

ggsave("summary_mc_zymo_corr_plot_zymobar.pdf", width = 50, height = 30,unit='cm',dpi=200)
ggsave("summary_mc_zymo_corr_plot_zymobar.png", width = 50, height = 30,unit='cm',dpi=200)

library(reshape2)
library(ggplot2)
library(scales)

test <- read.table("../phiX/PhiXQC_edit_dist.txt",sep='\t',header = T)

test_m_temp<-melt(test,na.rm = T)

xpos <- c(Inf,Inf)
ypos <- c(Inf,Inf)
mean_char=paste('Mean=',as.character(round(mean(test_m_temp$value),digits=2)))
sd_char=paste('SD=',as.character(round(sd(test_m_temp$value),digits = 2)))
annotateText <- c(mean_char,sd_char)
hjustvar<-c(1,1) 
vjustvar<-c(1,2.5)

pdf("zymo_edit_dist_perhist.pdf", onefile = TRUE,width=15, height=5)
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

test <- read.csv("zymo_postqc_edit_dist.txt",sep='\t',header = T)

test_m<-melt(test,na.rm = T)

list_plot=list()


#for (x in zymo_order){
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

