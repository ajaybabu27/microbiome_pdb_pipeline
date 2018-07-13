#Author: Ajay
#Date: 07/02/2018
#Description: Write out Qiime post qc fasta sequences per sample after performing QC. Takes the QC directory as input argument. 

import os
import sys
import pandas as pd
from Bio import SeqIO
import copy

working_directory=sys.argv[1]

#Store fasta sequences in Biopython object. 
record_dict = SeqIO.to_dict(SeqIO.parse(working_directory+"/all_samples_QC/dada2/dna-sequences.fasta", "fasta"))

#Store biome file in pandas object.
tsv_file=open(working_directory+'/all_samples_QC/dada2/feature-table.tsv','r')
df=pd.read_csv(tsv_file, header=1,sep='\t')

#Store feature IDs (fasta file header names)
rep_seq_list=df['#OTU ID']

#Iterate through each column (sample names) in pandas dataframe to print out per sample fasta file.
for column in df:
	
	if column=='#OTU ID':
		continue

	if len([i for i in ['M1','M6','M13','HC','Neg','Zymo'] if column.startswith(i)])>0:
		sample_folder=column.replace('.','-')
	else:
		sample_folder=column.replace('.','_')
		
	with open(working_directory+'/'+sample_folder+'/1_'+column+'.postqc.fasta','w') as fasta_out:
		for ind in list(df[column][df[column]>0].index.values):
			
			id_orig=copy.deepcopy(record_dict[df['#OTU ID'][ind]].id.split(" ")[0])					
			sequence=copy.deepcopy(record_dict[df['#OTU ID'][ind]])
			for count, copy_num in enumerate(range(0,int(df[column][ind]))):

				sequence.id=id_orig+"_"+str(count+1)
				
				SeqIO.write(sequence, fasta_out, "fasta")
			
