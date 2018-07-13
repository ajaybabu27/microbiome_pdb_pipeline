library(getopt)
library(RMySQL)
library(DBI)

args = matrix(c('work_dir'  , 'w', 2, "character", "Working directory",
                'dbhost'    , 'x', 2, "character", "Database host",
                'dbname'    , 'n', 2, "character", "Database name",
                'dbuser'    , 'u', 2, "character", "Database username",
                'dbpass'    , 'p', 2, "character", "Database password",
                'output'    , 'o', 2, "character", "Output file",
                'help'      , 'h', 0, "logical",   "Brief help message"
), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$dbhost)    ) { opt$dbhost    = "data1.hpc.mssm.edu"    }
if ( is.null(opt$dbname)    ) { opt$dbname    = "vanbah01_pathogens" }
if ( is.null(opt$dbuser)    ) { opt$dbuser    = "pathogendb_rw"      }
if ( is.null(opt$output)    ) { opt$output    = "tmp"      }

# Set up database connection
mydb = dbConnect(MySQL(), user=opt$dbuser, password=opt$dbpass, dbname=opt$dbname, host=opt$dbhost)

#Load mapping file
setwd(opt$work_dir)
raw_mapping_file<-read.table('mapping.tsv',header = T,sep='\t',comment.char = "",check.names = F)
raw_mapping_file$SampleID<-gsub("\\..*","",raw_mapping_file$`#SampleID`)

a_rs = dbSendQuery(mydb, "select * from tCdiffProjectSamples");
cdi_cohort_table = fetch(a_rs, n=-1)
cdi_cohort_table$case_eRAP_ID[is.na(cdi_cohort_table$case_eRAP_ID)]<-cdi_cohort_table$eRAP_ID[is.na(cdi_cohort_table$case_eRAP_ID)]
cdi_cohort_table$cdi_test_quikchek[grepl('Pos',cdi_cohort_table$cdi_test_quikchek)]<-'CD Pos'
cdi_cohort_table$cdi_test_quikchek[grepl('Neg',cdi_cohort_table$cdi_test_quikchek)]<-'CD Neg'


cdi_cohort_table$cdi_union[cdi_cohort_table$cdi_test_PCR=='Positive' | cdi_cohort_table$cdi_test_quikchek=='CD Pos']<-'CD Pos' 
cdi_cohort_table$cdi_union[cdi_cohort_table$subject_type=='control']<-'CD Neg'
cdi_cohort_table$cdi_union[cdi_cohort_table$cdi_test_PCR=='Negative' | cdi_cohort_table$cdi_test_quikchek=='CD Neg']<-'CD Neg' 

cdi_cohort_table2<-cdi_cohort_table[,c('specimen_ID','eRAP_ID','case_eRAP_ID','cdi_union','subject_type','specimen_sampling_date')]
colnames(cdi_cohort_table2)<-c('specimen_ID','eRAP_ID','caseERAPID','CDItestResult','CaseControlAnnot','SamplingDate')

mapping_file<-merge(raw_mapping_file,cdi_cohort_table2,by.x='SampleID',by.y='specimen_ID',all.x=T)
mapping_file$SampleID<-NULL
mapping_file$SampleType<-'Stool'
mapping_file$`#SampleID`<-as.character(mapping_file$`#SampleID`)
mapping_file$SampleType[grepl('Empty|M1|M6|M13|Zymo|HC|Neg',mapping_file$`#SampleID`)]<-mapping_file$`#SampleID`[grepl('Empty|M1|M6|M13|Zymo|HC|Neg',mapping_file$`#SampleID`)]


mapping_file$CaseControlAnnot[grepl('HC',mapping_file$`#SampleID`)]<-'healthy control'
mapping_file$CDItestResult[grepl('HC',mapping_file$`#SampleID`)]<-'CD Neg'

mapping_file$combo<-interaction(mapping_file$CDItestResult,mapping_file$CaseControlAnnot)

mapping_file$combo=factor(as.character(mapping_file$combo),levels = c("CD Pos.case","CD Neg.case","CD Neg.control","CD Neg.healthy control"),
                                                 labels=c("CD Pos Case samples","CD Neg Case samples","Matched Control Samples","Healthy control Samples"))


mapping_file$SampleNameDate<-interaction(mapping_file$`#SampleID`, mapping_file$SamplingDate)

#write.table(file='mapping_final.tsv',mapping_file,row.names = F,sep='\t')


#ABX Info addition

b_rs = dbSendQuery(mydb, "SELECT tpa.eRAP_ID,tpa.order_start_date,tpa.order_end_date, tabxl.synonym,tpa.route FROM `tPatientAntibiotics` tpa, tABX_lookup tabxl where tpa.med_ID=tabxl.medication_code");
abx_table = fetch(b_rs, n=-1)

mapping_file$ABX_admin_on_sample_collection<-NA
mapping_file$ABX_admin_24hrprior_sample_collection<-NA

for (x in 1:nrow(mapping_file)){
erap_id=mapping_file[x,'eRAP_ID']
if (is.na(erap_id)){
  next
}
sample_date=as.Date(mapping_file[x,'SamplingDate'])
abx_subset=abx_table[abx_table$eRAP_ID==erap_id,]
abx_subset$order_start_date<-as.Date(abx_subset$order_start_date)
abx_subset$order_end_date<-as.Date(abx_subset$order_end_date)
abx_list=c()
abx_list_24=c()
for (y in 1:nrow(abx_subset)){
  if ( (sample_date >= abx_subset[y,'order_start_date'] & sample_date <= abx_subset[y,'order_end_date'] )  ){
    
    abx_list=c(abx_list,abx_subset[y,'synonym'])        
    if((sample_date - abx_subset[y,'order_start_date'])>=1){
      abx_list_24=c(abx_list_24,abx_subset[y,'synonym'])
    }
    
  }
  
  if ((sample_date - abx_subset[y,'order_end_date'])==1){
    abx_list_24=c(abx_list_24,abx_subset[y,'synonym'])
  }
  
  }

  if(length(abx_list)>0){
    mapping_file[x,'ABX_admin_on_sample_collection']<-paste(unique(abx_list), collapse = ';')
  }

  if(length(abx_list_24)>0){
    mapping_file[x,'ABX_admin_24hrprior_sample_collection']<-paste(unique(abx_list_24), collapse = ';')
  }
}

write.table(file='mapping_final.tsv',mapping_file,row.names = F,sep='\t')

