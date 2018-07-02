require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
require 'bundler/setup'
require 'rspec/core/rake_task'
include Colors

task :default => :create_manifest_file

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

file "#{REPO_DIR}/scripts/env.sh" => "#{REPO_DIR}/scripts/example.env.sh" do
  cp "#{REPO_DIR}/scripts/example.env.sh", "#{REPO_DIR}/scripts/env.sh"
end

ENV_ERROR = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."

# ========================
# = create_manifest_file =
# ========================

desc "creates manifest file for Qiime2"
task :create_manifest_file 
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/manifest.csv" do |t|

 system <<-SH or abort

    sample_directory=#{FASTQ_DIR}/#{RUN_ID}
    
    mkdir -p $sample_directory/all_samples_QC
	
    echo -e sample-id,absolute-filepath,direction > $sample_directory/all_samples_QC/manifest.csv
    
    for sample in $sample_directory/*/ ; do	

	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)

	if [[ $sample_id = 'all_samples_QC' ]]; then

		continue

	fi

	
	forward_lib=$sample"0_Raw"/$sample_id"_"R1.fastq.gz
        reverse_lib=$sample"0_Raw"/$sample_id"_"R2.fastq.gz

        echo -e $sample_id_mod","$forward_lib",forward" >> $sample_directory/all_samples_QC/manifest.csv
        echo -e $sample_id_mod","$reverse_lib",reverse" >> $sample_directory/all_samples_QC/manifest.csv
	
    done     
  
    
 SH

end

# =============
# = run_dada2 =
# =============

desc "Perform Qiime2 QC steps and export sequences that pass QC"
task :run_dada2 => ["#{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza"]
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
    --p-cores 8 \
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


desc "Perform Qiime2 QC steps and export sequences that pass QC"
task :run_kraken => ["#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv"]
file "#{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv" => "#{QC_DIR}/#{RUN_ID}/all_samples_QC/dada2/table.qza" do |t|

 system <<-SH or abort

  module purge
  module load qiime2/2018.4
    
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

  sample_directory=#{QC_DIR}/#{RUN_ID}

  mkdir -p /tmp/db_custom_06192018
  cp /sc/orga/projects/InfectiousDisease/reference-db/kraken/db_custom_06192018/* /tmp/db_custom_06192018

  echo -e "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tunclassified_reads\tbacteria_reads\tviral_reads\tcdiff_reads\thuman_reads" > #{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv

  for sample in $sample_directory/*/ ; do

	sample_id=(`basename $sample`)
	sample_id_mod=(`echo $sample_id | tr '_' '.' | tr '-' '.' | awk '{gsub("_1", "");print}' `)

	if [[ $sample_id = 'all_samples_QC' ]]; then

		continue

	fi
	
	mkdir -p $sample/2_kraken	
	
	kraken --threads 12 -db /tmp/db_custom_06192018 --output $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.txt" --fasta-input $sample/1_$sample_id_mod.fasta

	kraken-report --db /tmp/db_custom_06192018 $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.txt" > $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" && \
ktImportTaxonomy $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" -t 5 -m 3 -s 0 -o $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.kronaQC.html"

	unclassified_reads=(`grep "unclassified" $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`) 
	bacteria_reads=(`grep "Bacteria" $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`) 
	viral_reads=(`grep "Viruses" $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	cdiff_reads=(`grep "Clostridioides difficile" $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)
	human_reads=(`grep "Homo sapiens" $sample/2_kraken/$sample_id_mod"_kraken_db_custom_out.QC.kreport" | cut -d$'\t' -f 2`)

	echo -e $sample_id_mod"\t\t\t\t"$unclassified_reads"\t"$bacteria_reads"\t"$viral_reads"\t"$cdiff_reads"\t"$human_reads >> #{QC_DIR}/#{RUN_ID}/all_samples_QC/mapping.tsv

  
 SH

end

# =================
# = upload_QC_pdb =
# =================



