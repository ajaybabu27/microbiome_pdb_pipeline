#Author: Ajay
#Date: 10/10/2017
#Description: Concatenates files with reads mapped per contig in a microbial community library sample. Takes a folder containing read stats files. 

import os  
import sys
dict_samples={}

working_directory=sys.argv[1]

#Iterate through files in directory and find files ending with ".sorted.seq"
for fn in os.listdir(working_directory):
     
	if fn.endswith(".seq") :
		sample=fn.split('.sorted.seq')[0]

		dict_samples[sample]={}	#Declare dictionary to store values in all files. 
		with open(working_directory+'/'+fn,'r') as bwa_out:
			for line in bwa_out:
				line2=line.rstrip()
				line3=line2.rsplit('\t')
				org=line3[1].rsplit('_')
				org="_".join(org[:2])

			if org in dict_samples[sample].keys():
				dict_samples[sample][org]=dict_samples[sample][org]+float(line3[0])
			else:	
				dict_samples[sample][org]=float(line3[0])

sample_list=sorted(dict_samples.iterkeys())
organism_list=sorted(dict_samples[sample_list[0]].iterkeys())

#Iterate through the dictionary to print the values in a unified summary file. 
with open(working_directory+'/summary.tsv','w') as bwa_out_summ:

	bwa_out_summ.write("Sample"+'\t'+ '\t'.join(organism_list)+"\n")
	for sample in sample_list:
		print_list=[]
		for org in organism_list:
			try:
				print_list.append(str(dict_samples[sample][org]))
			except KeyError:
				print_list.append('0')
				bwa_out_summ.write(sample+'\t'+ '\t'.join(print_list)+"\n")
		
