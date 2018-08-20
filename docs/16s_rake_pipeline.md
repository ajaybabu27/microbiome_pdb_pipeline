# 16S Pipeline 

Detailed description of tasks explaining workflows from QC to various analysis steps. 

## Flow chart of 16S Pipeline 

Flowchart describing various rake tasks involved in the key steps of the pipeline.

![16s_rake_pipeline.png](https://github.com/ajaybabu27/microbiome_pdb_pipeline/blob/master/docs/16s_rake_pipeline.jpg)
For more information on the tools used in the tasks mentioned in the flowchart, look through ![Qiime documentation] (https://docs.qiime2.org/2018.6/)
An excellent tutorial on a typical ![Qiime 16S workflow] (https://docs.qiime2.org/2018.6/tutorials/moving-pictures/)


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
The task "run_MC_QC" currently runs analysis for Zymogen and Clemente Lab Microbial communities (MC). The folders associated with these libraries are
expected to be named in the following manner:
```
Zymo* (e.g. Zymo-1,Zymo-2)
M* (e.g. M1,M6,M13)
```
Note that these libraries start with letters Z (Zymo MC) and M (Celemente MC) and it is important that the other sample IDs don't start with these letters. 

The 16S FASTQ files are stored in the following location after parsing to organize the above mentioned format:
```
/sc/orga/projects/InfectiousDisease/microbiome-output/samples/$Run_ID
```
The above location also contains the initial QC file. 
!!!Warning!!! Do not use this folder as QC output until unless you want to represent the results in PathogenDB front end website. All the files in this folders are directly linked to the website. 


### Resource Requirements

Most of the tasks require a minimum of 12 cores and memory of 40 GB at the minimum.


### User Input Requirements

Most of the tasks are automated except for one task `create_postQC_biome_file` which requires the user to input the minimum number of reads needed for filtering low quality libraries. This threshold should be set based on the QC analysis results. 

## Tasks

The tasks for 16S pipeline are run in the following order:

QC Tasks: 
1. `create_manifest_file` 
2. `run_dada2` 
3. `run_kraken` 
4. `run_MC_QC` 
5. `create_postQC_biome_file` 

Analysis Tasks: 
1. `prepare_Qiime_analysis` 
2. `run_Qiime_analysis` 





   
