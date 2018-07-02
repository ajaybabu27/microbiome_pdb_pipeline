#Author: Ajay
#Date: 07/02/2018
#Description: Write out fasta sequences per sample after performing QC

import os
import sys
import pandas as pd
from Bio import SeqIO

working_directory=sys.argv[1]


record_dict = SeqIO.to_dict(SeqIO.parse(working_directory+"/all_samples_QC/dada2/dna-sequences.fasta", "fasta"))

file=open(working_directory+'/all_samples_QC/dada2/feature-table.tsv','r')
df=pd.read_csv(file, header=1,sep='\t')

rep_seq_list=df['#OTU ID']

for column in df:
	
	if column=='#OTU ID':
		continue

	if len([i for i in ['M1','M6','M13','HC','Neg','Zymo'] if column.startswith(i)])>0:
		sample_folder=column.replace('.','-')
	else:
		sample_folder=column.replace('.','_')
		
	with open(working_directory+'/'+sample_folder+'/1_'+column+'.postqc.fasta','w') as fasta_out:
		for ind in list(df[column][df[column]>0].index.values):
			id_orig=record_dict[df['#OTU ID'][ind]].id.split(" ")[0]
			sequence=record_dict[df['#OTU ID'][ind]]
			for count, copy_num in enumerate(range(0,int(df[column][ind]))):

				sequence.id=id_orig+"_"+str(count+1)
				
				SeqIO.write(sequence, fasta_out, "fasta")
			
