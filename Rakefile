require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
require 'bundler/setup'
require 'rspec/core/rake_task'
include Colors

task :default => :check

LSF = LSFClient.new
LSF.disable! if ENV['LSF_DISABLED']  # Run everything locally if set (useful for debugging)

REPO_DIR = File.dirname(__FILE__)


#######
# Other environment variables that may be set by the user for specific tasks (see README.md)
#######

RUN_ID = ENV['RUN_ID']
FASTQ_DIR = ENV['FASTQ_DIR']
QC_DIR = ENV['QC_DIR']
ANALYSIS_DIR = ENV['ANALYSIS_DIR']
READS_THRESHOLD= ENV['READS_THRESHOLD']

task :env do
     
  sc_orga_scratch = "/sc/orga/scratch/#{ENV['USER']}"
  
  ENV['TMP'] ||= Dir.exists?(sc_orga_scratch) ? sc_orga_scratch : "/tmp"
  
  ENV['PERL5LIB'] ||= "/usr/bin/perl5.10.1"

end

file "#{REPO_DIR}/scripts/env.sh" => "#{REPO_DIR}/scripts/example.env.sh" do
  cp "#{REPO_DIR}/scripts/example.env.sh", "#{REPO_DIR}/scripts/env.sh"
end

ENV_ERROR = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."

desc "Checks environment variables and requirements before running tasks"
task :check => [:env, "#{REPO_DIR}/scripts/env.sh"] do

  mkdir_p ENV['TMP'] or abort "FATAL: set TMP to a directory that can store scratch files"
 
end

# ========================
# = create_manifest_file =
# ========================

desc "creates manifest file for Qiime2"
task :create_manifest_file => [:check, "#{QC_DIR}/#{RUN_ID}/all_samples_QC/manifest.csv"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/manifest.csv" do |t|

 system <<-SH or abort

    sample_directory=#{FASTQ_DIR}/#{RUN_ID}
    output_directory=#{QC_DIR}/#{RUN_ID}
    
    mkdir -p $output_directory/all_samples_QC
	
    echo -e sample-id,absolute-filepath,direction > $output_directory/all_samples_QC/manifest.csv
    
    for sample in $sample_directory/*/ ; do	

	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)

	if [[ $sample_id = 'all_samples_QC' ]]; then

		continue

	fi

	
	forward_lib=$sample"0_Raw"/$sample_id"_"R1.fastq.gz
    reverse_lib=$sample"0_Raw"/$sample_id"_"R2.fastq.gz

    echo -e $sample_id_mod","$forward_lib",forward" >> $output_directory/all_samples_QC/manifest.csv
    echo -e $sample_id_mod","$reverse_lib",reverse" >> $output_directory/all_samples_QC/manifest.csv
	
    done     
  
    
 SH

end

# =============
# = run_dada2 =
# =============

desc "Perform Qiime2 QC steps"
task :run_dada2 => [:check, "#{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza" => "#{QC_DIR}/#{RUN_ID}/all_samples_QC/manifest.csv" do |t|

 system <<-SH or abort
    
    module purge
    module load qiime2/2018.4

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path #{QC_DIR}/#{RUN_ID}/all_samples_QC/manifest.csv \
    --output-path #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux.qza \
    --source-format PairedEndFastqManifestPhred33

    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux.qza \
    --p-cores 12 \
    --p-front-f GTGCCAGCMGCCGCGGTAA \
    --p-front-r GGACTACHVGGGTWTCTAAT \
    --o-trimmed-sequences #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux-trimmed.qza

    qiime demux summarize \
    --i-data #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux-trimmed.qza \
    --o-visualization #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux-trimmed.qzv

    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs #{QC_DIR}/#{RUN_ID}/all_samples_QC/paired-end-demux-trimmed.qza \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-n-threads 0  \
    --output-dir #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2   
    
 SH

end

# ==============
# = run_kraken =
# ==============


desc "Export sequences that pass Qiime QC and run Kraken for contaminant analysis"
task :run_kraken => [:check,"#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping_final.tsv"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping_final.tsv" => "#{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza" do |t|

 system <<-SH or abort

  module purge
  module load qiime2/2018.4
  #module load Krona
  module load kraken/2.0.7
  

  qiime tools export \
  #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza \
  --output-dir #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2

  qiime tools export \
  #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/representative_sequences.qza \
  --output-dir #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2

  biom convert -i #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/feature-table.biom -o #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/feature-table.tsv --to-tsv

  module purge
  module load python py_packages
  python #{REPO_DIR}/scripts/post_qc_reads_export.py #{QC_DIR}/#{RUN_ID}

  num_cores=(`grep -c ^processor /proc/cpuinfo`)

  echo -e "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tunclassified_reads\tbacteria_reads\tviral_reads\tcdiff_reads\thuman_reads" > #{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv

  sample_directory=#{QC_DIR}/#{RUN_ID}
  for sample in $sample_directory/*/ ; do

	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)

	if [[ $sample_id = 'all_samples_QC' ]]; then

		continue

	fi
	
	mkdir -p $sample"2_kraken"
	
	kraken2 --threads 12 -db /sc/orga/projects/InfectiousDisease/reference-db/kraken/kraken2_db_standard --output $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.txt" --report $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" $sample"1_"$sample_id_mod.postqc.fasta

	ktImportTaxonomy $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" -t 5 -m 3 -s 0 -o $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.kronaQC.html"

	unclassified_reads=(`grep "unclassified" $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	bacteria_reads=(`grep "Bacteria" $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	viral_reads=(`grep "Viruses" $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	cdiff_reads=(`grep "Clostridioides difficile" $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	human_reads=(`grep "Homo sapiens" $sample"2_"kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)

	echo -e $sample_id_mod"\t\t\t\t"$unclassified_reads"\t"$bacteria_reads"\t"$viral_reads"\t"$cdiff_reads"\t"$human_reads >> #{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv

   done

  #process PhiX reads

  mkdir -p #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX

  kraken2 --threads 12 -db /sc/orga/projects/InfectiousDisease/reference-db/kraken/kraken2_db_standard --output #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX/PhiX_kraken_db_custom_out.txt --report #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX/PhiXQC_kraken_db_custom_out.QC.kreport --paired #{QC_DIR}/#{RUN_ID}/Undetermined_S0_L001_R1_001.fastq.gz #{QC_DIR}/#{RUN_ID}/Undetermined_S0_L001_R2_001.fastq.gz
    
  ktImportTaxonomy #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX/PhiXQC_kraken_db_custom_out.QC.kreport -t 5 -m 3 -s 0 -o #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX/PhiXQC_kraken_db_custom_out.kronaQC.html

  module purge
  module load R
     
  #Output metadata file
  Rscript #{REPO_DIR}/scripts/mapping_file_builder.R -w #{QC_DIR}/#{RUN_ID}/all_samples_QC -p #{ENV['PDBPASS']} -r #{RUN_ID}



 SH

end

# =============
# = run_MC_QC =
# =============

desc "Calculate abundance and edit distances by comparing expected vs oberved MC measures"
task :run_MC_QC => [:check,"#{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_qc/zymo/zymo_edit_dist_perhist.pdf"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_qc/zymo/zymo_edit_dist_perhist.pdf" => "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping_final.tsv" do |t|

 system <<-SH or abort
  module purge
  module load bwa
  module load samtools/1.7
  module load R
  module load bamtools

  sample_directory=#{QC_DIR}/#{RUN_ID}
  mkdir -p sample_directory/all_samples_QC/mc_QC/zymo
  mkdir -p sample_directory/all_samples_QC/mc_QC/clemente

  # Process PhiX reads
  
  output_directory=#{ENV['TMP']}/mc_out/#{RUN_ID}/PhiX  
  #{REPO_DIR}/scripts/phiX.sh #{QC_DIR}/#{RUN_ID} $output_directory #{REPO_DIR}
  
  cp $output_directory/PhiXQC_edit_dist.txt $output_directory/PhiXQC.sorted.seq #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/PhiX

  # Process Zymo and Celemente MC samples
  mkdir -p #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/zymo
  mkdir -p #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/clemente

  for sample in $sample_directory/*/ ; do
  	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)
	query_fasta=$sample"1_"$sample_id_mod".postqc.fasta"
		
	if [[ $sample_id =~ ^M.* ]]; then
	
		ref_fasta_file=/sc/orga/projects/InfectiousDisease/reference-db/microbial_community_standards/jose_mc.fasta
		output_directory=#{ENV['TMP']}/mc_out/#{RUN_ID}/clemente
		mkdir -p $output_directory
		{REPO_DIR}/scripts/bwa_mapper.sh $output_directory $sample_id_mod $query_fasta $ref_fasta_file		
		
	elif [[ $sample_id =~ ^Z.* ]]; then
		ref_fasta_file=/sc/orga/projects/InfectiousDisease/reference-db/microbial_community_standards/zymo_mc.fasta
		output_directory=#{ENV['TMP']}/mc_out/#{RUN_ID}/zymo
		mkdir -p $output_directory
		{REPO_DIR}/scripts/bwa_mapper.sh $output_directory $sample_id_mod $query_fasta $ref_fasta_file		

	fi
	
  done

  output_directory=#{ENV['TMP']}/mc_out/#{RUN_ID}/clemente
  python #{REPO_DIR}/scripts/bwa_mapper_parse.py $output_directory
  paste $output_directory/*edit* > #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/clemente/clemente_postqc_edit_dist.txt
  cp $output_directory/summary.tsv #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/clemente/summary.tsv
  Rscript #{REPO_DIR}/scripts/community_stats_clemente.R #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/clemente

  output_directory=#{ENV['TMP']}/mc_out/#{RUN_ID}/zymo
  python #{REPO_DIR}/scripts/bwa_mapper_parse.py $output_directory
  paste $output_directory/*edit* > #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/zymo/zymo_postqc_edit_dist.txt
  cp $output_directory/summary.tsv #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/zymo/summary.tsv
  Rscript #{REPO_DIR}/scripts/community_stats_zymo.R #{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_QC/zymo

 SH

end

# ==================
# = create_QC_page =
# ==================

#desc "Create a landing page for all QC results"

#Add Qiime QC and MC QC charts
#Filter biom table to only include QC passed samples. Remove MC and Negative and Empty samples. i.e. only include stool samples
#Link out to Biom File, Metadata file

# =================
# = upload_QC_pdb =
# =================

#desc "push QC info to PathogenDB tables"

#task :run_MC_QC => [:check,"#{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_qc/zymo/zymo_edit_dist_perhist.pdf"]
#file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mc_qc/zymo/zymo_edit_dist_perhist.pdf" => "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv" do |t| 

# system <<-SH or abort

#  module purge
#  module load R

#  Rscript #{REPO_DIR}/scripts/16S_upload_results.R -w #{QC_DIR}/#{RUN_ID} -r #{RUN_ID} 

# SH

#end

# ============================
# = create_postQC_biome_file =
# ============================

#desc "Create a post QC biome file based on QC analysis"
task :create_postQC_biome_file => [:check, "#{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/filtered_representative_sequences.qza"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/filtered_representative_sequences.qza" do |t|
  
  system <<-SH or abort

    module purge
    module load qiime2/2018.4
    
    mkdir -p #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out
    
    qiime feature-table filter-samples \
    --i-table #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza \
    --p-min-frequency #{READS_THRESHOLD} \
    --m-metadata-file #{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping_final.tsv \
    --p-where "SampleType='Stool'" \
    --o-filtered-table #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/table_rem_bad_lib.qza

	qiime feature-table filter-features \
	--i-table #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/table_rem_bad_lib.qza \
	--p-min-frequency 1 \
	--o-filtered-table #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza

	qiime feature-table filter-seqs \
	--i-data #{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/representative_sequences.qza \
	--i-table #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
	--o-filtered-data #{QC_DIR}/#{RUN_ID}/all_samples_QC/final_biom_out/filtered_representative_sequences.qza

  SH

end

# ==========================
# = prepare_Qiime_analysis =
# ==========================

desc "Prepare analysis run"
task :prepare_Qiime_analysis => [:check, "#{ANALYSIS_DIR}/5_taxons/taxonomy.qzv"]
file "#{ANALYSIS_DIR}/5_taxons/taxonomy.qzv" do |t|

#Merege BIOM tables and metatables.

end

# ======================
# = run_Qiime_analysis =
# ======================

desc "Calculate Alpha Diversity, Beta Diversity and Taxon classification for singe/multiple Runs"
task :run_Qiime_analysis => ["#{ANALYSIS_DIR}/4_taxons/taxa-bar-plots.qzv"]
file "#{ANALYSIS_DIR}/4_taxons/taxa-bar-plots.qzv" do |t|

system <<-SH or abort

  module purge
  module load qiime2/2018.4

  #Construct phylogeny for diversity analyses

  mkdir -p #{ANALYSIS_DIR}/1_aligned_OTU

  qiime alignment mafft \
    --i-sequences #{ANALYSIS_DIR}/0_merged_OTU/representative_sequences.qza \
    --o-alignment #{ANALYSIS_DIR}/1_aligned_OTU/aligned-representative_sequences.qza \
    --p-n-threads 12

  qiime alignment mask \
    --i-alignment #{ANALYSIS_DIR}/1_aligned_OTU/aligned-representative_sequences.qza \
    --o-masked-alignment #{ANALYSIS_DIR}/1_aligned_OTU/masked-aligned-representative_sequences.qza

  qiime phylogeny fasttree \
    --i-alignment #{ANALYSIS_DIR}/1_aligned_OTU/masked-aligned-representative_sequences.qza \
    --o-tree #{ANALYSIS_DIR}/1_aligned_OTU/unrooted-tree.qza \
    --p-n-threads -1

  qiime phylogeny midpoint-root \
    --i-tree #{ANALYSIS_DIR}/1_aligned_OTU/unrooted-tree.qza \
    --o-rooted-tree #{ANALYSIS_DIR}/1_aligned_OTU/rooted-tree.qza

  #Perform diversity analyses
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny #{ANALYSIS_DIR}/1_aligned_OTU/rooted-tree.qza \
    --i-table #{ANALYSIS_DIR}/0_merged_OTU/table.qza \
    --p-sampling-depth #{READS_THRESHOLD} \
    --m-metadata-file #{ANALYSIS_DIR}/mapping_final_combined.tsv \
    --output-dir #{ANALYSIS_DIR}/2_core-metrics-results \
    --p-n-jobs -1
  
  mkdir -p #{ANALYSIS_DIR}/3_alpha_diversity
  qiime diversity alpha-rarefaction \
    --i-table #{ANALYSIS_DIR}/0_merged_OTU/table.qza \
    --i-phylogeny #{ANALYSIS_DIR}/1_aligned_OTU/rooted-tree.qza \
    --p-max-depth #{READS_THRESHOLD} \
    --m-metadata-file #{ANALYSIS_DIR}/mapping_final_combined.tsv \
    --o-visualization #{ANALYSIS_DIR}/3_alpha_diversity/alpha-rarefaction.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity #{ANALYSIS_DIR}/2_core-metrics-results/shannon_vector.qza \
    --m-metadata-file #{ANALYSIS_DIR}/mapping_final_combined.tsv \
    --o-visualization #{ANALYSIS_DIR}/3_alpha_diversity/shannon-group-significance.qzv

  #Perform taxon classification
  mkdir -p #{ANALYSIS_DIR}/4_taxons
  qiime feature-classifier classify-sklearn \
    --i-classifier /sc/orga/projects/InfectiousDisease/reference-db/gg_13_8_otus/gg-13-8-99-515-806-nb-classifier.qza \
    --i-reads #{ANALYSIS_DIR}/0_merged_OTU/representative_sequences.qza \
    --o-classification #{ANALYSIS_DIR}/4_taxons/taxonomy.qza \
    --p-n-jobs -1

  qiime metadata tabulate \
    --m-input-file #{ANALYSIS_DIR}/4_taxons/taxonomy.qza \
    --o-visualization #{ANALYSIS_DIR}/4_taxons/taxonomy.qzv
  
  qiime taxa barplot \
    --i-table #{ANALYSIS_DIR}/0_merged_OTU/table.qza \
    --i-taxonomy #{ANALYSIS_DIR}/4_taxons/taxonomy.qza \
    --m-metadata-file #{ANALYSIS_DIR}/mapping_final_combined.tsv \
    --o-visualization #{ANALYSIS_DIR}/4_taxons/taxa-bar-plots.qzv


SH

end
# ==============================
# = create_Qiime_analysis_page =
# ==============================

desc "Create landing page for all Qiime analyses" 

# =======================
# = upload_analysis_pdb =
# =======================

desc "push analysis results to PDB tables" 


