#Author: Ajay
#Date : 10/10/2017
#Description: Get read mapping stats for PhiX spike in control

module load R
module load fastqc
module load seqtk

phiX_reference=/sc/orga/projects/InfectiousDisease/reference-db/phiX_control
working_directory=$1
output_directory=$2
repo_directory=$3

#PhiX forward and reverse libraries
forward_lib=$working_directory/Undetermined_S0_L001_R1_001.fastq.gz
reverse_lib=$working_directory/Undetermined_S0_L001_R2_001.fastq.gz

mkdir -p $output_directory

#PhiX reference genome
fasta_file=$phiX_reference/PhiX-reference-genome.fasta
sample_id=PhiXQC


num_cores=(`grep -c ^processor /proc/cpuinfo`)

#Map fastq files to reference. 
bwa mem -t $num_cores $fasta_file $forward_lib $reverse_lib > $output_directory/$sample_id.sam 2> $output_directory/$sample_id.log

echo preqc"_"$sample_id > $output_directory/$sample_id"_"edit_dist.txt
cat $output_directory/$sample_id.sam | awk '{ if (NR>11) {print $12} }' | cut -d':' -f3 >> $output_directory/$sample_id"_"edit_dist.txt

samtools view -@ $num_cores -bS $output_directory/$sample_id.sam -o $output_directory/$sample_id.bam 

samtools sort -@ $num_cores $output_directory/$sample_id.bam > $output_directory/$sample_id.sorted.bam

total_lib_count=(`gzip -cd $forward_lib | wc -l | awk '{print $1/4}'`)

#Filter reads to output primary mapped paired reads
$repo_directory/scripts/reads_per_seq.sh $output_directory/$sample_id.sorted.bam

#Get unmapped reads
samtools view -@ $num_cores -b -f 4 $output_directory/$sample_id.sorted.bam > $output_directory/$sample_id.unmapped.bam

samtools bam2fq $output_directory/$sample_id.unmapped.bam | seqtk seq -A > $output_directory/umapped.fa

#Extract and print out edit distance stats.
num_multimapped_alignments=(`samtools flagstat $output_directory/$sample_id.sam | grep supplementary | cut -d ' ' -f1`)
mean_ed=(`cat $output_directory/$sample_id.sam | awk '{ print $12 }' | cut -d':' -f3 | R --slave -e 'x <- scan(file="stdin",quiet=TRUE); mean(x)' | cut -d ' ' -f2`)
sd_ed=(`cat $output_directory/$sample_id.sam | awk '{ print $12 }' | cut -d':' -f3 | R --slave -e 'x <- scan(file="stdin",quiet=TRUE); sd(x)' | cut -d ' ' -f2`)
	
echo -e $total_lib_count"\t"total >> $output_directory/$sample_id".sorted.seq"
echo -e $num_multimapped_alignments"\t"multimapped_alignments >> $output_directory/$sample_id".sorted.seq"
echo -e $mean_ed"\t"mean_edit_distance >> $output_directory/$sample_id".sorted.seq"
echo -e $sd_ed"\t"sd_edit_distance >> $output_directory/$sample_id".sorted.seq"

