export FONTCONFIG_PATH=/etc/fonts
export LC_ALL=aa_DJ.utf8
export LANG=C.UTF-8

module purge
module load anaconda3
source activate qiime2-2019.1

trpth='/sc/arion/projects/InfectiousDisease/microbiome-output/analyses/microany00020'
mkdir -p $trpth/0_merged_OTU

qiime feature-table merge \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/H434/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00396/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00504/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00543/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00544/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00705/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01222/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01246/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01247/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01282/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01334/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01419/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01420/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --i-tables /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01540/all_samples_QC/final_biom_out/table_rem_bad_lib_zero_features.qza \
  --o-merged-table $trpth/0_merged_OTU/table.qza

qiime feature-table merge-seqs \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/H434/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00396/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00504/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00543/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00544/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00705/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01222/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01246/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01247/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01282/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01334/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01419/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01420/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --i-data /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01540/all_samples_QC/final_biom_out/filtered_representative_sequences.qza \
  --o-merged-data $trpth/0_merged_OTU/representative_sequences.qza

tab[1]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/H434/all_samples_QC/mapping_final.tsv'
tab[2]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00396/all_samples_QC/mapping_final.tsv'
tab[3]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00504/all_samples_QC/mapping_final.tsv'
tab[4]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00543/all_samples_QC/mapping_final.tsv'
tab[5]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00544/all_samples_QC/mapping_final.tsv'
tab[6]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD00705/all_samples_QC/mapping_final.tsv'
tab[7]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01222/all_samples_QC/mapping_final.tsv'
tab[8]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01246/all_samples_QC/mapping_final.tsv'
tab[9]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01247/all_samples_QC/mapping_final.tsv'
tab[10]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01282/all_samples_QC/mapping_final.tsv'
tab[11]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01334/all_samples_QC/mapping_final.tsv'
tab[12]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01419/all_samples_QC/mapping_final.tsv'
tab[13]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01420/all_samples_QC/mapping_final.tsv'
tab[14]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539/all_samples_QC/mapping_final.tsv'
tab[15]='/sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01540/all_samples_QC/mapping_final.tsv'

otpth_tab='mapping_final_combined.tsv'

# run script
# ----------
touch "$trpth"/"$otpth_tab"
head -n 1  "${tab[1]}" > "$trpth"/"$otpth_tab"
for ((i=1;i<=15;i++)); do
   tail -n +2 "${tab[$i]}" >> "$trpth"/"$otpth_tab"
done
