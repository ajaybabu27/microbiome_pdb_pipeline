# 16S Pipeline 

Detailed description of tasks explaining workflows from QC to various analysis steps. 

## Flow chart of 16S Pipeline 

![16s_rake_pipeline.png](https://github.com/ajaybabu27/microbiome_pdb_pipeline/blob/master/docs/16s_rake_pipeline.jpg)


## Requirements

### Folder structure for directory containing FASTQ raw reads 

The first task create_manifest_file requires the folder structure described below:
```
FASTQ_FOLDER
├── Sample_ID
│   ├── 0_Raw 
│       ├── Sample_ID_R1.fastq.gz 
│       └── Sample_ID_R2.fastq.gz 
├── Sample_ID 2
├── Sample_ID 3
├── Sample_ID N
├── Undetermined_S0_L001_R1_001.fastq.gz (Non-barcoded forward reads) 
├── Undetermined_S0_L001_R2_001.fastq.gz (Non-barcoded reverse reads) 
```
The task "run_MC_QC" currently runs analysis for Zymogen and Clemente Lab Microbial communities. The folders associated with these libraries are
expected to be named in the following manner:
```


```
Note that 

   
