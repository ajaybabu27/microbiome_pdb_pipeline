#Author: Ajay
#Date: 10/10/2017
#Description: Generate read and edit distance stats for various microbial commmunity libraries. 

output_directory=$1
sample_id=$2
query_fasta=$3
fasta_file=$4

mkdir -p $output_directory

#Map post qc microbial community sample fasta file to appropriate reference genome.
bwa mem -t 12 $fasta_file $query_fasta > $output_directory/$sample_id.sam 2> $output_directory/$sample_id.log

# convert sam to bam
samtools view -@ 12 -bS $output_directory/$sample_id.sam -o $output_directory/$sample_id.bam
samtools sort $output_directory/$sample_id.bam > $output_directory/$sample_id.sorted.bam

#Filter for primary mapped reads with more than 200bp aligned. 
total_lib_count=(`grep ">" $query_fasta | wc -l`) #change to postQC
samtools view -b -F s $output_directory/$sample_id.sorted.bam > $output_directory/$sample_id.sorted.primary.bam
samtools view -h $output_directory/$sample_id.sorted.primary.bam  | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > 200 or /^@/' | samtools view -bS - > $output_directory/$sample_id.sorted.primary.200.bam
samtools index $output_directory/$sample_id.sorted.primary.200.bam
samtools idxstats $output_directory/$sample_id.sorted.primary.200.bam | awk '{ print $3 "\t" $1}' | head -n -1 > $output_directory/$sample_id".sorted.seq" 

#compute edit distances and downstream stats form the filtered BAM file.
echo $sample_id
echo postqc"_"$sample_id > $output_directory/$sample_id"_"edit_dist.txt
samtools view -h $output_directory/$sample_id.sorted.primary.200.bam | awk '{ if (NR>11) {print $12} }' | cut -d':' -f3 >> $output_directory/$sample_id"_"edit_dist.txt

num_multimapped_alignments=(`samtools flagstat $output_directory/$sample_id.sorted.primary.200.bam | grep supplementary | cut -d ' ' -f1`)
mean_ed=(`samtools view -h $output_directory/$sample_id.sorted.primary.200.bam | awk '{ if (NR>11) {print $12} }' | cut -d':' -f3 | R --slave -e 'x <- scan(file="stdin",quiet=TRUE); mean(x)' | cut -d ' ' -f2`)
sd_ed=(`samtools view -h $output_directory/$sample_id.sorted.primary.200.bam | awk '{ if (NR>11) {print $12} }' | cut -d':' -f3 | R --slave -e 'x <- scan(file="stdin",quiet=TRUE); sd(x)' | cut -d ' ' -f2`)

echo -e $total_lib_count"\t"total >> $output_directory/$sample_id".sorted.seq"
echo -e $num_multimapped_alignments"\t"multimapped_alignments >> $output_directory/$sample_id".sorted.seq"
echo -e $mean_ed"\t"mean_edit_distance >> $output_directory/$sample_id".sorted.seq"
echo -e $sd_ed"\t"sd_edit_distance >> $output_directory/$sample_id".sorted.seq"


