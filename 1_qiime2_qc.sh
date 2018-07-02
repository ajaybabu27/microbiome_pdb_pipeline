#!/bin/bash

#Author: Ajay
#Description: Qiime2 QC workflow
#Date: 04/27/2018

#module purge
#module load qiime2

working_directory=/sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434
analysis_directory=/sc/orga/projects/InfectiousDisease/microbiome-output/analyses/microany00007_q2

export FONTCONFIG_PATH=/etc/fonts
export LC_ALL=aa_DJ.utf8
export LANG=C.UTF-8


mkdir -p $analysis_directory


#echo -e "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tunclassified_reads\tbacteria_reads\tviral_reads\tcdiff_reads\thuman_reads" > $analysis_directory/mapping.tsv


#echo -e sample-id,absolute-filepath,direction > $analysis_directory/manifest.csv

for sample in $working_directory/*/ ; do
	
	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)
	
	if [[ $sample_id_mod =~ ^[A-Z]{2}[0-9]{5}.* ]]; then
		sample_id_mod2=(`echo $sample_id_mod | sed -e 's/\.1//g' | awk '{print $1"_1A"}' `) 		
	else
		sample_id_mod2=$sample_id_mod
	fi
	
	


	#kraken --threads 36 -db /tmp/db_custom_06192018 --output $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.txt" --fasta-input $analysis_directory/2_dada2/filtered_seq/$sample_id_mod.fasta

	#kraken-report --db /tmp/db_custom_06192018 $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.txt" > $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" && \
#ktImportTaxonomy $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" -t 5 -m 3 -s 0 -o $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.kronaQC.html"

	#unclassified_reads=(`grep "unclassified" $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`) 
	#bacteria_reads=(`grep "Bacteria" $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`) 
	#viral_reads=(`grep "Viruses" $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	#cdiff_reads=(`grep "Clostridioides difficile" $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	#human_reads=(`grep "Homo sapiens" $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)

	#forward_lib=$sample"0_Raw"/$sample_id"_"R1.fastq.gz
	#reverse_lib=$sample"0_Raw"/$sample_id"_"R2.fastq.gz

	#echo -e $sample_id_mod","$forward_lib",forward" >> $analysis_directory/manifest.csv
	#echo -e $sample_id_mod","$reverse_lib",reverse" >> $analysis_directory/manifest.csv

	#echo -e $sample_id_mod"\t\t\t\t"$unclassified_reads"\t"$bacteria_reads"\t"$viral_reads"\t"$cdiff_reads"\t"$human_reads >> $analysis_directory/mapping.tsv

	#ln -s $analysis_directory/2_dada2/kraken_out/$sample_id_mod"_kraken_db_custom_out.kronaQC.html" /hpc/users/kumara22/www/kraken_qc/$sample_id_mod2"_H434".qc.krona.html

	

done

: <<'END'
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $analysis_directory/manifest.csv \
  --output-path $analysis_directory/paired-end-demux.qza \
  --source-format PairedEndFastqManifestPhred33

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $analysis_directory/paired-end-demux.qza \
  --p-cores 8 \
  --p-front-f GTGCCAGCMGCCGCGGTAA \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --o-trimmed-sequences $analysis_directory/paired-end-demux-trimmed.qza

qiime demux summarize \
  --i-data $analysis_directory/paired-end-demux-trimmed.qza \
  --o-visualization $analysis_directory/paired-end-demux-trimmed.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $analysis_directory/paired-end-demux-trimmed.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads 36  \
  --output-dir $analysis_directory/2_dada2 

qiime tools export \
  table.qza \
  --output-dir $analysis_directory/2_dada2

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

qiime feature-table summarize \
  --i-table $analysis_directory/2_dada2/table.qza \
  --o-visualization $analysis_directory/2_dada2/table.qzv \
  --m-sample-metadata-file $analysis_directory/mapping_final.tsv

qiime feature-table tabulate-seqs \
  --i-data $analysis_directory/2_dada2/representative_sequences.qza \
  --o-visualization $analysis_directory/2_dada2/representative_sequences.qzv

qiime alignment mafft \
  --i-sequences $analysis_directory/2_dada2/representative_sequences.qza \
  --o-alignment $analysis_directory/2_dada2/aligned-representative_sequences.qza \
  --p-n-threads 36

qiime alignment mask \
  --i-alignment $analysis_directory/2_dada2/aligned-representative_sequences.qza \
  --o-masked-alignment $analysis_directory/2_dada2/masked-aligned-representative_sequences.qza

qiime phylogeny fasttree \
  --i-alignment $analysis_directory/2_dada2/masked-aligned-representative_sequences.qza \
  --o-tree $analysis_directory/2_dada2/unrooted-tree.qza \
  --p-n-threads 36

qiime phylogeny midpoint-root \
  --i-tree $analysis_directory/2_dada2/unrooted-tree.qza \
  --o-rooted-tree $analysis_directory/2_dada2/rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $analysis_directory/2_dada2/rooted-tree.qza \
  --i-table $analysis_directory/2_dada2/table.qza \
  --p-sampling-depth 4000 \
  --m-metadata-file $analysis_directory/mapping_final.tsv \
  --output-dir $analysis_directory/4_core-metrics-results \
  --p-n-jobs 36

qiime feature-classifier classify-sklearn \
  --i-classifier $analysis_directory/ref_classifier/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads $analysis_directory/2_dada2/representative_sequences.qza \
  --o-classification $analysis_directory/5_taxons/taxonomy.qza \
  --p-n-jobs 36

qiime metadata tabulate \
  --m-input-file $analysis_directory/5_taxons/taxonomy.qza \
  --o-visualization $analysis_directory/5_taxons/taxonomy.qzv
  
qiime taxa barplot \
  --i-table $analysis_directory/2_dada2/table.qza \
  --i-taxonomy $analysis_directory/5_taxons/taxonomy.qza \
  --m-metadata-file $analysis_directory/mapping_final.tsv \
  --o-visualization $analysis_directory/5_taxons/taxa-bar-plots.qzv

#Filter based on sample
qiime feature-table filter-samples \
  --i-table $analysis_directory/2_dada2/table.qza \
  --m-metadata-file $analysis_directory/mapping_final.tsv \
  --p-where "caseERAPID='191235'" \
  --o-filtered-table $analysis_directory/191235/191235.qza

qiime diversity core-metrics \
  --i-table $analysis_directory/191235/191235.qza \
  --p-sampling-depth 30000 \
  --m-metadata-file $analysis_directory/mapping_final_191235.tsv \
  --output-dir $analysis_directory/191235/191235_core-metrics-results \
  --p-n-jobs 36

qiime emperor plot \
  --i-pcoa $analysis_directory/191235/191235_core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file $analysis_directory/mapping_final_191235.tsv \
  --p-custom-axes numDays \
  --o-visualization $analysis_directory/191235/191235_core-metrics-results/bray-curtis-emperor-DaysSinceFirstCaseCD.qzv
END
