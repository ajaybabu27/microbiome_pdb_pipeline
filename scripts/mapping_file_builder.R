#Author: Ajay
#Date: 10/10/2017
#Description: Generate mapping file/meta data file for the 16S run.

library(getopt)
library(RMySQL)
library(DBI)

args = matrix(c('work_dir'  , 'w', 2, "character", "Working directory",
                'run_id'    , 'r', 2, "character", "Illumina Run ID" ,
                'dbhost'    , 'x', 2, "character", "Database host",
                'dbname'    , 'n', 2, "character", "Database name",
                'dbuser'    , 'u', 2, "character", "Database username",
                'dbpass'    , 'p', 2, "character", "Database password",
                'help'      , 'h', 0, "logical",   "Brief help message"
), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$dbhost)    ) { opt$dbhost    = "data1"    }
if ( is.null(opt$dbname)    ) { opt$dbname    = "vanbah01_pathogens" }
if ( is.null(opt$dbuser)    ) { opt$dbuser    = "pathogendb_rw"      }


# Set up database connection
mydb = dbConnect(MySQL(), user=opt$dbuser, password=opt$dbpass, dbname=opt$dbname, host=opt$dbhost)

#Load skeleton mapping file
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
mapping_file$SampleType[grepl('Empty|M1|M6|M13|Zymo|Neg',mapping_file$`#SampleID`)]<-mapping_file$`#SampleID`[grepl('Empty|M1|M6|M13|Zymo|Neg',mapping_file$`#SampleID`)]


mapping_file$CaseControlAnnot[grepl('HC',mapping_file$`#SampleID`)]<-'healthy control'
mapping_file$CDItestResult[grepl('HC',mapping_file$`#SampleID`)]<-'CD Neg'



mapping_file$combo<-interaction(mapping_file$CDItestResult,mapping_file$CaseControlAnnot)

mapping_file$combo=factor(as.character(mapping_file$combo),levels = c("CD Pos.case","CD Neg.case","CD Neg.control","CD Neg.healthy control"),
                                                 labels=c("CD Pos Case samples","CD Neg Case samples","Matched Control Samples","Healthy control Samples"))

mapping_file$SampleNameDate<-interaction(mapping_file$`#SampleID`, mapping_file$SamplingDate)

#Add CDI infection order
mapping_file$cdi_order<-NA
for(pat in unique(mapping_file$caseERAPID)){
  temp_cdi_order_df<-mapping_file[mapping_file$eRAP_ID==pat & grepl("CD",mapping_file$`#SampleID`) ,]
  if(nrow(temp_cdi_order_df)==0 | is.na(pat)){
    print('skipped')
    next
  }
  else{
    temp_cdi_order_df<-temp_cdi_order_df[order(temp_cdi_order_df$SamplingDate),]
    temp_cdi_order_df$cdi_order<-seq(1,nrow(temp_cdi_order_df))
    mapping_file$cdi_order <- replace(mapping_file$cdi_order, mapping_file$`#SampleID` %in% temp_cdi_order_df$`#SampleID`, temp_cdi_order_df$cdi_order)
  }
    
}

#ABX Info addition

b_rs = dbSendQuery(mydb, "SELECT tpa.eRAP_ID,tpa.order_start_date,tpa.order_end_date, tabxl.synonym,tpa.route FROM `tPatientAntibiotics` tpa, tABX_lookup tabxl where tpa.med_ID=tabxl.medication_code");
abx_table = fetch(b_rs, n=-1)

mapping_file$ABX_admin_on_sample_collection<-NA
mapping_file$ABX_admin_24hrprior_sample_collection<-NA
mapping_file$ABX_admin_1wprior_sample_collection<-NA
mapping_file$vanco_oral_admin_1wprior_sample_collection<-NA
mapping_file$vanco_int_admin_1wprior_sample_collection<-NA
mapping_file$met_oral_admin_1wprior_sample_collection<-NA
mapping_file$met_int_admin_1wprior_sample_collection<-NA

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
	abx_list_1w=c()
	abx_list_1w_vanco_oral=c()
	abx_list_1w_met_oral=c()
	abx_list_1w_vanco_int=c()
	abx_list_1w_met_int=c()
	if (nrow(abx_subset)>0){
  	for (y in 1:nrow(abx_subset)){
  		if ( (sample_date >= abx_subset[y,'order_start_date'] & sample_date <= abx_subset[y,'order_end_date'] )  ){
  		
  			abx_list=c(abx_list,abx_subset[y,'synonym'])        
  			if((sample_date - abx_subset[y,'order_start_date'])>=1){
  				abx_list_24=c(abx_list_24,abx_subset[y,'synonym'])
  				abx_list_1w=c(abx_list_1w,abx_subset[y,'synonym'])
  				if(abx_subset[y,'synonym']=='VANCOMYCIN'){
  				  if(abx_subset[y,'route']=='oral'){
  				    abx_list_1w_vanco_oral=c(abx_list_1w_vanco_oral,abx_subset[y,'synonym'])			  
  				  }
  				  else if(abx_subset[y,'route']=='intravenous'){
  				    abx_list_1w_vanco_int=c(abx_list_1w_vanco_int,abx_subset[y,'synonym'])			  
  				  }
  				}
  				if(abx_subset[y,'synonym']=='METRONIDAZOLE'){
  				  if(abx_subset[y,'route']=='oral'){
  				    abx_list_1w_met_oral=c(abx_list_1w_met_oral,abx_subset[y,'synonym'])			  
  				  }
  				  else if(abx_subset[y,'route']=='intravenous'){
  				    abx_list_1w_met_int=c(abx_list_1w_met_int,abx_subset[y,'synonym'])			  
  				  }
  				}
  			
  			}
  		
  		}
  	  
  		if ((sample_date - abx_subset[y,'order_end_date'])==1){
  			abx_list_24=c(abx_list_24,abx_subset[y,'synonym'])
  		}
  	  
  	  if ((sample_date - abx_subset[y,'order_end_date'])<=7){
  	    abx_list_1w=c(abx_list_1w,abx_subset[y,'synonym'])
  	    
  	    if(abx_subset[y,'synonym']=='VANCOMYCIN'){
  	      if(abx_subset[y,'route']=='oral'){
  	        abx_list_1w_vanco_oral=c(abx_list_1w_vanco_oral,abx_subset[y,'synonym'])			  
  	      }
  	      else if(abx_subset[y,'route']=='intravenous'){
  	        abx_list_1w_vanco_int=c(abx_list_1w_vanco_int,abx_subset[y,'synonym'])			  
  	      }
  	    }
  	    if(abx_subset[y,'synonym']=='METRONIDAZOLE'){
  	      if(abx_subset[y,'route']=='oral'){
  	        abx_list_1w_met_oral=c(abx_list_1w_met_oral,abx_subset[y,'synonym'])			  
  	      }
  	      else if(abx_subset[y,'route']=='intravenous'){
  	        abx_list_1w_met_int=c(abx_list_1w_met_int,abx_subset[y,'synonym'])			  
  	      }
  	    }
  	    
  	    
  	  }
  	  
  	 }
  }
	if(length(abx_list)>0){
		mapping_file[x,'ABX_admin_on_sample_collection']<-paste(unique(abx_list), collapse = ';')
	}

	if(length(abx_list_24)>0){
		mapping_file[x,'ABX_admin_24hrprior_sample_collection']<-paste(unique(abx_list_24), collapse = ';')
	}
	
	if(length(abx_list_1w)>0){
	  mapping_file[x,'ABX_admin_1wprior_sample_collection']<-paste(unique(abx_list_1w), collapse = ';')
	}
	
	if(length(abx_list_1w_vanco_int)>0){
	  mapping_file[x,'vanco_int_admin_1wprior_sample_collection']<-paste(unique(abx_list_1w_vanco_int), collapse = ';')
	}
	if(length(abx_list_1w_vanco_oral)>0){
	  mapping_file[x,'vanco_oral_admin_1wprior_sample_collection']<-paste(unique(abx_list_1w_vanco_oral), collapse = ';')
	}
	if(length(abx_list_1w_met_int)>0){
	  mapping_file[x,'met_int_admin_1wprior_sample_collection']<-paste(unique(abx_list_1w_met_int), collapse = ';')
	}
	if(length(abx_list_1w_met_oral)>0){
	  mapping_file[x,'met_oral_admin_1wprior_sample_collection']<-paste(unique(abx_list_1w_met_oral), collapse = ';')
	}
}

mapping_file$ABX_admin_on_sample_collection_bool<-!is.na(mapping_file$ABX_admin_on_sample_collection)
mapping_file$ABX_admin_24hrprior_sample_collection_bool<-!is.na(mapping_file$ABX_admin_24hrprior_sample_collection)
mapping_file$ABX_admin_1wprior_sample_collection_bool<-!is.na(mapping_file$ABX_admin_1wprior_sample_collection)
mapping_file$vanco_int_admin_1wprior_sample_collection_bool<-!is.na(mapping_file$vanco_int_admin_1wprior_sample_collection)
mapping_file$met_int_admin_1wprior_sample_collection_bool<-!is.na(mapping_file$met_int_admin_1wprior_sample_collection)
mapping_file$vanco_oral_admin_1wprior_sample_collection_bool<-!is.na(mapping_file$vanco_oral_admin_1wprior_sample_collection)
mapping_file$met_oral_admin_1wprior_sample_collection_bool<-!is.na(mapping_file$met_oral_admin_1wprior_sample_collection)

#Add plate ID

r_rs = dbSendQuery(mydb, paste("SELECT * FROM `tIlluminaCoreSubmissions` where Sequence_Run_ID='",opt$run_id,"'",sep=''));
runs_table = fetch(r_rs, n=-1)
plates=sort(runs_table$Sequence_Plate_ID)

p_rs = dbSendQuery(mydb, paste("SELECT * FROM `tNGS_Plates` where Systematic_Plate_ID in ('",plates[1],"','",plates[2],"')",sep=''))
plates_table = fetch(p_rs, n=-1)
plates_table$Extract_Name<-gsub("-","\\.",plates_table$Extract_Name)
plates_table$Extract_Name<-gsub(" ","\\.",plates_table$Extract_Name)
plates_table[grepl("HC",plates_table$Extract_Name),'Extract_Name']<-'HC'
plates_table[plates_table$Extract_Name %in% c("HC",'M1','Zymo','Neg') & plates_table$Systematic_Plate_ID==plates[1],'Extract_Name']<-paste(plates_table[plates_table$Extract_Name %in% c("HC",'M1','Zymo','Neg') & plates_table$Systematic_Plate_ID==plates[1],'Extract_Name'],'.1',sep='')
plates_table[plates_table$Extract_Name %in% c("HC",'M1','Zymo','Neg') & plates_table$Systematic_Plate_ID==plates[2],'Extract_Name']<-paste(plates_table[plates_table$Extract_Name %in% c("HC",'M1','Zymo','Neg') & plates_table$Systematic_Plate_ID==plates[2],'Extract_Name'],'.2',sep='')
mapping_file$SampleID<-mapping_file$`#SampleID`
mapping_file[grepl("RC",mapping_file$SampleID),'SampleID']<-gsub("\\..*","",mapping_file[grepl("RC",mapping_file$`#SampleID`),'SampleID'])

plates_table$Extract_Name=gsub(".1A","",plates_table$Extract_Name)

mapping_file=merge(mapping_file,plates_table[,c('Extract_Name','Systematic_Plate_ID')],by.x='SampleID',by.y='Extract_Name',all.x=T)
mapping_file$SampleID<-NULL

mapping_file$illumina_runid<-opt$run_id

write.table(file='mapping_final.tsv',mapping_file,row.names = F,sep='\t')


#Add patient Age, BMI, WBC and Creatnine levels

#e_rs=dbSendQuery(mydb, "SELECT * FROM `tPatientEncounter`");



