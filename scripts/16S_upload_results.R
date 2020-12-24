library(getopt)
args = matrix(c('work_dir'  , 'w', 2, "character", "Working directory",
                'illum_rid' , 'r', 2, "character", "Illumina RunID",
                'dbhost'    , 'x', 2, "character", "Database host",
                'dbname'    , 'n', 2, "character", "Database name",
                'dbuser'    , 'u', 2, "character", "Database username",
                'dbpass'    , 'p', 2, "character", "Database password",
                'output'    , 'o', 2, "character", "Output file",
                'help'      , 'h', 0, "logical",   "Brief help message"
), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$dbhost)    ) { opt$dbhost    = "data1"    }
if ( is.null(opt$dbname)    ) { opt$dbname    = "vanbah01_pathogens" }
if ( is.null(opt$dbuser)    ) { opt$dbuser    = "pathogendb_rw"      }

#############
# LIBRARIES #
#############

suppressMessages(library(RMySQL));
suppressMessages(library(DBI));


mydb = dbConnect(MySQL(), user=opt$dbuser, password=opt$dbpass, dbname=opt$dbname, host=opt$dbhost);

e_rs = dbSendQuery(mydb, "SELECT extract_ID, concentration FROM tExtracts where extract_ID LIKE '%.1%';")
extracts_table             = fetch(e_rs, n=-1);
extracts_table$extract_mod<-gsub("\\..*","",extracts_table$extract_ID)


RunID=opt$illum_rid

setwd(opt$work_dir)

file_ls <- as.character(unzip("all_samples_QC/dada2/denoising_stats.qza", list = TRUE)$Name)
dada_stats=unz(filename = file_ls[grep("stats.tsv",file_ls)], description = "all_samples_QC/dada2/denoising_stats.qza")
dada_stats_df=read.csv(dada_stats,sep='\t',comment.char = "#")
dada_stats_df$bad_reads_per=(((dada_stats_df$input-dada_stats_df$denoised))/dada_stats_df$input)*100
dada_stats_df$unmerged_per=(((dada_stats_df$input-dada_stats_df$merged)/dada_stats_df$input)*100)-dada_stats_df$bad_reads_per
dada_stats_df$chimeric_per=(((dada_stats_df$input-dada_stats_df$non.chimeric)/dada_stats_df$input)*100)-(dada_stats_df$bad_reads_per+dada_stats_df$unmerged_per)

mapfile = "all_samples_QC/mapping_final.tsv"
map_df<-read.csv(mapfile,header = T,sep='\t')

all_stats_df=merge(dada_stats_df,map_df,by=1)
all_stats_df[,c(13:17)][is.na(all_stats_df[,c(13:17)])]<-0

all_stats_df$bacteria_reads_per=(all_stats_df$bacteria_reads/all_stats_df$input)*100
all_stats_df$viral_reads_per=(all_stats_df$viral_reads/all_stats_df$input)*100
#all_stats_df$cdiff_reads_per=(all_stats_df$cdiff_reads/all_stats_df$input)*100
#all_stats_df$bacteria_reads_per=all_stats_df$bacteria_reads_per-all_stats_df$cdiff_reads_per
all_stats_df$human_reads_per=(all_stats_df$human_reads/all_stats_df$input)*100
all_stats_df$unclassified_reads_per=(all_stats_df$unclassified_reads/all_stats_df$input)*100


run_qual_stats_per<-all_stats_df[,c(1,2,7:9,43:46)]
run_qual_stats_per$unclassified_reads_per=100-rowSums(all_stats_df[,c(7:9,43:45)])
#run_qual_stats_per[,c(6:10)]<-100*(run_qual_stats_per[,c(6:10)]/run_qual_stats_per$input)
run_qual_stats_per$total_per<-rowSums(run_qual_stats_per[,c(3:9)],na.rm=TRUE)
#run_qual_stats_per$unclassified_reads<-0
run_qual_stats_per[,-c(1,2)]<-round(run_qual_stats_per[,-c(1,2)],2)
#run_qual_stats_per$unclassified_reads<-100-rowSums(run_qual_stats_per[,-c(1,2,6)],na.rm = T)
#run_qual_stats_per[run_qual_stats_per$unclassified_reads<0,'unclassified_reads']<-0
#run_qual_stats_per[126,'unclassified_reads']<-0
run_qual_stats_per$total_per<-NULL
run_qual_stats_per$input_label<-paste(as.character(run_qual_stats_per$input/1000),'K',sep='')

run_qual_stats_per$sample.id<-as.character(run_qual_stats_per$sample.id)
sample_order<-run_qual_stats_per[order(run_qual_stats_per$input),'sample.id']
run_qual_stats_per$sample.id<-factor(run_qual_stats_per$sample.id,levels=sample_order)

library(reshape2)
run_qual_stats_per$input<-as.factor(as.character(run_qual_stats_per$input))

run_qual_stats_per_m<-melt(run_qual_stats_per)
run_qual_stats_per_m$value[is.na(run_qual_stats_per_m$value)]<-0

library(ggplot2)

png('all_samples_QC/run_quality_QC_plot.png',width = 8000,height = 2000,res=200)
ggplot(data=run_qual_stats_per_m, aes(x = run_qual_stats_per_m$sample.id, y = run_qual_stats_per_m$value,
                                      fill=run_qual_stats_per_m$variable,
                                     label=signif(run_qual_stats_per_m$value, digits = 2))) +
  geom_bar(stat='identity')+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=-0.1,size=10))+
  geom_text(data=run_qual_stats_per_m,aes(x=run_qual_stats_per_m$sample.id,y=100,label=run_qual_stats_per_m$input_label),vjust=-0.05,hjust=-0.1,angle=65,size=5)+
  labs(x='Samples',y='% Abundance')+
  expand_limits(y = max(100 * 1.06))
dev.off()

#DNA concentration charts. Need to develop more for future. 
if(F){
map_df$sample_id_mod<-gsub("\\..*","",map_df$X.SampleID)
map_df<-merge(map_df,extracts_table,by.x=35,by.y=3)

ggplot(data=map_df,aes(x=concentration,y=bacteria_reads))+geom_point()
ggplot(data=map_df,aes(x=concentration,y=bacteria_reads))+geom_point()+coord_cartesian(ylim = c(0, 10000)) 
}

mydb = dbConnect(MySQL(), user=opt$dbuser, password=opt$dbpass, dbname=opt$dbname, host=opt$dbhost);

for (r in 1:nrow(all_stats_df)){
  sample_name=as.character(all_stats_df[r,'sample.id'])
  if(grepl("^[[:alpha:]]{2}\\d{5}",sample_name)){
    sample_name<-gsub("\\..*","",sample_name)
    extract_id=paste(sample_name,'.1A',sep='')
    microbiome_id=paste(sample_name,'.',RunID,sep='')
    #print(extract_id)
  }
  
  else{
    extract_id=sample_name
    microbiome_id=sample_name
  }
 
  sequencing_type='paired_end_16S_V4'
  read_length='250'
  platform='illumina'
  unmerged_per=round(all_stats_df[r,'unmerged_per'],2)
  chimeric_per=round(all_stats_df[r,'chimeric_per'],2)
  bad_reads_per=round(all_stats_df[r,'bad_reads_per'],2)
  read_count=all_stats_df[r,'input']
  unclassified_reads_per=round(all_stats_df[r,'unclassified_reads_per'],2)
  bacteria_reads_per=round(all_stats_df[r,'bacteria_reads_per'],2)
  viral_reads_per=round(all_stats_df[r,'viral_reads_per'],2)
  #cdiff_reads_per=round(all_stats_df[r,'cdiff_reads_per'],2)
  human_reads_per=round(all_stats_df[r,'human_reads_per'],2)
  
  
  dbSendQuery(mydb, paste("INSERT INTO tMicrobiomes ",
                          "(microbiome_ID,extract_ID,sequencing_type,run_ID,platform,read_length,read_count,percent_unmerged_reads_16S,percent_chimeric_reads_16S,percent_bad_reads,percent_unclassified_reads,percent_bacteria,percent_virus,percent_human) ",
                          "VALUES('",microbiome_id,"','",extract_id,"','",sequencing_type,"','",RunID,"','",platform,"','",read_length,"','",read_count,"','",unmerged_per,"','",chimeric_per,
                          "','",bad_reads_per,"','",unclassified_reads_per,"','",bacteria_reads_per,"','",viral_reads_per,"','",human_reads_per,"') ON DUPLICATE KEY UPDATE microbiome_ID = VALUES (microbiome_ID)",sep=""))
                          
 
  
  #dbSendQuery(mydb, paste("INSERT INTO tMicrobiome_linker ",
  #                        "( microbiome_ID,analysis_ID) ",
  #                        "VALUES('",microbiome_id,"','",analysis_id,"')",sep=''))
}

#Compute summary stats for run
total_nonbarcoded_reads<-sum(all_stats_df$input)
phix_data<-read.table(file='all_samples_QC/mc_QC/PhiX/PhiXQC.sorted.seq',sep='')
phix_reads<-phix_data$V1[1]
total_barcoded_reads<-phix_data$V1[2]
total_reads<-total_barcoded_reads+total_nonbarcoded_reads
phix_percentage<-round((phix_reads/(total_reads))*100,2)

dbSendQuery(mydb, paste("INSERT INTO tMicrobiome_runs ",
                        "(run_ID,total_reads,barcoded_reads,non_barcoded_reads,phix_reads,phix_per) ",
                        "VALUES('",RunID,"','",total_reads,"','",total_barcoded_reads,"','",total_nonbarcoded_reads,"','",phix_reads,"','",phix_percentage,"') ON DUPLICATE KEY UPDATE run_ID = VALUES (run_ID)",sep = ''))

