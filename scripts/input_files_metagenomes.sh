cd /sc/arion/projects/InfectiousDisease/illumina-fastq/TD01539
mkdir -p /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539

for samples in */* ; do
	#samples=${samples%?}

	IFS='/' read -r folder sample_file <<< "$samples"
	IFS='_' read -r id1 id2 id3 id4 <<< "$sample_file"
	IFS='-' read -r id expt expt2 <<< "$id1"
	
	if [[ $expt == '' ]] ;
	then		
		sample_name=$id

	elif [[ $id == 'HC' ]] ;
	then
		sample_name=$id"_"$expt"_"$expt2

	elif [[ $id == 'RC'* ]]	;
	then 
		sample_name=$id	
		
	else
		sample_name=$id"_"$expt
		
	fi

		
	
	mkdir -p /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539/$sample_name/0_Raw
		
	cp $samples /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539/$sample_name/0_Raw/$sample_name"_"$id3.fastq.gz
	
	
done

#cp /sc/arion/projects/InfectiousDisease/illumina-fastq/TD01539/* /sc/arion/projects/InfectiousDisease/microbiome-output/samples/TD01539

: <<'END'

cd /sc/orga/projects/InfectiousDisease/illumina-fastq/TD01118
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01118

for samples in */* ; do
	#samples=${samples%?}

	IFS='/' read -r folder sample_file <<< "$samples"
	IFS='_' read -r id1 id2 id3 id4 id5 <<< "$sample_file"
	#IFS='-' read -r id expt <<< "$id1"
	
#	if [[ $expt == '' ]] ;
#	then		
		sample_name=$id1
	
#	else

#		sample_name=$id"_"$expt
		
#	fi

		
	
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01118/$sample_name/0_Raw
		
	cp $samples /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01118/$sample_name/0_Raw/$sample_name"_"$id3.fastq.gz
	
	
done



cd /sc/orga/projects/InfectiousDisease/illumina-fastq/TD01004
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01004

for samples in */ ; do
	samples=${samples%?}
	
	IFS='_' read -r id1 id2 <<< "$samples"
	IFS='-' read -r id expt <<< "$id1"
	
	if [[ $expt == '' ]] ;
	then		
		sample_name=$id
	
	else

		sample_name=$id"_"$expt
		
	fi

		
	echo $sample_name
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01004/$sample_name/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01004/$sample_name/0_Raw/$sample_name"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01004/$sample_name/0_Raw/$sample_name"_"R2.fastq.gz
	
done

cp /sc/orga/projects/InfectiousDisease/illumina-fastq/TD01004/* /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD01004



cd /sc/orga/projects/InfectiousDisease/illumina-fastq/TD00908
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00908

for samples in */* ; do
	#samples=${samples%?}

	IFS='/' read -r folder sample_file <<< "$samples"
	IFS='_' read -r id1 id2 id3 id4 id5 <<< "$sample_file"
	IFS='-' read -r id expt <<< "$id1"
	
	if [[ $expt == '' ]] ;
	then		
		sample_name=$id
	
	else

		sample_name=$id"_"$expt
		
	fi

		
	
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00908/$sample_name/0_Raw
		
	cp $samples /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00908/$sample_name/0_Raw/$sample_name"_"$id4.fastq.gz
	
	
done

cp /sc/orga/projects/InfectiousDisease/illumina-fastq/TD00908/* /sc/orga/projects/InfectiousDisease/microbiome-output/samples/TD00908

cd /sc/orga/projects/InfectiousDisease/illumina-fastq/H434
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434

for samples in * ; do
	echo $samples	
	IFS='_' read -r rem id expt <<< "$samples"
	
	if [[ $expt == '' ]] ;
	then
		
		sample_name=$id
	else
		sample_name=$id"_"$expt
	fi
	echo $sample_name/0_Raw
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434/$sample_name/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434/$sample_name/0_Raw/$sample_name"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H434/$sample_name/0_Raw/$sample_name"_"R2.fastq.gz
	
done


cd /sc/orga/projects/InfectiousDisease/illumina-fastq/H399
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H399

for samples in * ; do
		
	# Rename sample names
	
	if [[ $samples == *Zymo-Q* ]] ;
	then
		sample_name='ZymoExtract_Q_1A'

	elif [[ $samples == *Zymo* ]] ;
	then
		sample_name='ZymoDNA_1A'
	else
		sample_name=$(echo `basename $samples | cut -d"_" -f 2` )"_1A"
		
	fi
		
	echo $samples
	echo $sample_name/0_Raw
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H399/$sample_name/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H399/$sample_name/0_Raw/$sample_name"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/H399/$sample_name/0_Raw/$sample_name"_"R2.fastq.gz
	
done



cd /sc/orga/projects/InfectiousDisease/illumina-fastq/G665
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G665

for samples in * ; do
		
	# Rename sample names
	
	if [[ $samples == *Zymo_DNA* ]] ;
	then
		sample_name='ZymoDNA_IDT_1A'

	elif [[ $samples == *Zymo_extrct* ]] ;
	then
		sample_name='ZymoExtract_IDT_1A'
	else
		sample_name=$(echo `basename $samples | cut -d"_" -f 2 | cut -d"-" -f 1` )"_IDT_1A"
		if [[ $samples == *_3* ]] ;
		then
			sample_name=$sample_name"_IDT_clonal"
		fi		
	fi
		
	echo $samples
	echo $sample_name/0_Raw
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G665/$sample_name/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G665/$sample_name/0_Raw/$sample_name"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G665/$sample_name/0_Raw/$sample_name"_"R2.fastq.gz
	
done



cd /sc/orga/projects/InfectiousDisease/illumina-fastq/G664/QC_J259.G664_HVB-_16S_samples.PE.Library_Type_Not_Found.Pipe_Line_Settings_Not_Found.Species_Not_Found/Raw/Library_Not_Found.IlluminaHiSeq2500.Library_Type_Not_Found
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G664

for samples in * ; do
		
	# Rename sample names
	
	if [[ $samples == *Zymo_DNA* ]] ;
	then
		sample_name='ZymoDNA_1A'

	elif [[ $samples == *Zymo_extrct* ]] ;
	then
		sample_name='ZymoExtract_1A'
	else
		sample_name=$(echo `basename $samples | cut -d"_" -f 2 | cut -d"-" -f 1` )"_1A"
		if [[ $samples == *_3* ]] ;
		then
			sample_name=$sample_name"_clonal"
		fi		
	fi
		
	#echo $samples
	#echo $sample_name/0_Raw
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G664/$sample_name/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G664/$sample_name/0_Raw/$sample_name"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G664/$sample_name/0_Raw/$sample_name"_"R2.fastq.gz

done



cd /sc/orga/projects/InfectiousDisease/illumina-fastq/G539
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G539

for samples in * ; do
	# Rename IDT and illumina primer samples
	if [[ $samples == *IDT* ]] ;
	then
		primer='IDT'
	else
		primer='illumina'		
	fi
	
	# Rename sample names
	
	if [[ $samples == *Zymo_DNA* ]] ;
	then
		sample_name='ZymoDNA_1A'
	
	else
		sample_name=$(echo `basename $samples | cut -d"_" -f 2 | cut -d"-" -f 1` )"_1A"		
	fi
		
	#echo $samples
	#echo $sample_name"_"$primer/0_Raw
	mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G539/$sample_name"_"$primer/0_Raw
	
	cp $samples/*R1*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G539/$sample_name"_"$primer/0_Raw/$sample_name"_"$primer"_"R1.fastq.gz
	cp $samples/*R2*.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G539/$sample_name"_"$primer/0_Raw/$sample_name"_"$primer"_"R2.fastq.gz

	
done


cd /sc/orga/projects/InfectiousDisease/illumina-fastq/G358
mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G358

for samples in * ; do
	# Rename duplicates
	if [[ $samples == *dupl* ]] ;
	then
		dupl=2
	else
		dupl=1		
	fi
	
	# Rename sample names
	if [[ $samples == *CD* ]] ;
	then
		sample_name='CD01204_1A'
	elif [[ $samples == *Zymo_DNA* ]] ;
	then
		sample_name='ZymoDNA_1A'
	elif [[ $samples == *Zymo_MoBio* ]] ;
	then
		sample_name='ZymoMoBio_1A'
	else
		sample_name=$(echo `basename $samples | cut -d"_" -f 1` )"_1A"		
	fi
	
	# Rename experimental variables
	if [[ $samples == *0_5_ng* ]] ;
	then
		expt_cond='0_5ng'
	elif [[ $samples == *_1_ng* || $samples == *_1ng*  ]] ;
	then
		expt_cond='1ng'
	elif [[ $samples == *_7_5_ng* ]] ;
	then
		expt_cond='7_5ng'
	elif [[ $samples == *_2_5_ng* ]] ;
	then
		expt_cond='2_5ng'
	elif [[ $samples == *_5_ng* ]] ;
	then
		expt_cond='5ng'	
	
	elif [[ $samples == *_10_ng* ]] ;
	then
		expt_cond='10ng'
	fi
	echo $samples
	echo $sample_name"_"$dupl"_"$expt_cond/0_Raw
	#mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G358/$sample_name"_"$dupl"_"$expt_cond/0_Raw
	
	#cp $samples/*_1.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G358/$sample_name"_"$dupl"_"$expt_cond/0_Raw/$sample_name"_"$dupl"_"$expt_cond"_"R1.fastq.gz
	#cp $samples/*_2.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/G358/$sample_name"_"$dupl"_"$expt_cond/0_Raw/$sample_name"_"$dupl"_"$expt_cond"_"R2.fastq.gz

	
done


cd /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/QC_M476.F932_HVB_16S_C_diff.PE.Custom_Capture.NULL.Species_Not_Found/Raw/DNA.IlluminaHiSeq2500.Custom_Capture
cp /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/QC_M476.F932_HVB_16S_C_diff.PE.Custom_Capture.NULL.Species_Not_Found/* /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932

for samples in * ; do
		
	#for x in `find $samples -name "*.fastq.gz" | sort`; do

	if [[ $samples == *dup* ]] ;
	then
		
		sample_name=$(echo `basename $samples | cut -d"_" -f 3` )"_1A_2"
		expt_condition=$(echo `basename $samples | cut -d"_" -f 4` )

		if [[ $expt_condition == 62 ]] ;
		then	
			expt_condition='62_5C'

		elif [[ $expt_condition == 0 ]] ;
		then 
			expt_condition='0_1ng'

		elif [[ $expt_condition == 1 ]] ;
		then 
			expt_condition='1ng'

		elif [[ $expt_condition == 10 ]] ;
		then 
			expt_condition='10ng'		

		
		fi

		mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw
		
		cp $samples/*R1_001.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw/$sample_name"_"$expt_condition"_"R1.fastq.gz
		cp $samples/*R2_001.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw/$sample_name"_"$expt_condition"_"R2.fastq.gz


	

	else
		
		sample_name=$(echo `basename $samples | cut -d"_" -f 2` )"_1A_1"
		expt_condition=$(echo `basename $samples | cut -d"_" -f 3` )

		if [[ $expt_condition == 62 ]] ;
		then	
			expt_condition='62_5C'

		elif [[ $expt_condition == 0 ]] ;
		then 
			expt_condition='0_1ng'

		elif [[ $expt_condition == 1 ]] ;
		then 
			expt_condition='1ng'

		elif [[ $expt_condition == 10 ]] ;
		then 
			expt_condition='10ng'		

		elif [ $sample_name == "CD01190" ] || [ $sample_name == "RC01151" ] ;
		then
			expt_condition='55C'
		fi

		mkdir -p /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw
		
		cp $samples/*R1_001.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw/$sample_name"_"$expt_condition"_"R1.fastq.gz
		cp $samples/*R2_001.fastq.gz /sc/orga/projects/InfectiousDisease/microbiome-output/samples/F932/$sample_name"_"$expt_condition/0_Raw/$sample_name"_"$expt_condition"_"R2.fastq.gz
		
		
	fi
	
	done


module load pigz




cd /sc/orga/projects/InfectiousDisease/cdiff_metagenomics


for run_id in 2016-05-06_D495 2016-06-27_D966 2016-06-27_D967 2016-06-27_D968 2016-08-23_E305 2016-09-12_E306 2016-09-23_E562 2016-09-26_E547  ; do 
	
	run_id_exact=$(echo $run_id | cut -d"_" -f2)
	mkdir -p /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact
	cp $run_id/* /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact
	
	for x in `find $run_id -name "*.fastq" -o -name "*.fastq.gz" | sort`; do
		sample_name=$(echo `basename $x | cut -d"_" -f 1` )"_1A"
		lib_name=$(echo `basename $x | cut -d"_" -f 1` )"_1A_"$(echo `basename $x | cut -d"_" -f 4` )
		mkdir -p /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name/0_Raw

		if [ ${x: -3} == ".gz" ]; then
			
			cp $x /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name/0_Raw/$lib_name.fastq.gz &		

		else

			submitjob 1 -c 12 -m 50 -q private -g bashia02 pigz -c -p 12 -k $x \> /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name/0_Raw/$lib_name.fastq.gz		
		fi	
		
	done		
		
done



for run_id in 2016-09-15_E575; do #2016-05-06_D495 2016-09-15_E575 2016-06-27_D966 2016-06-27_D967 2016-06-27_D968 2016-08-23_E305 2016-09-12_E306 2016-09-23_E562 2016-09-26_E547
	
	run_id_exact=$(echo $run_id | cut -d"_" -f2)
	mkdir -p /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact
	cp $run_id/* /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact
	
	for x in `find $run_id/mapping -name "*.fastq" -o -name "*.fastq.gz" | sort`; do

		sample_name=$(echo `basename $x | cut -d"_" -f 1` )"_1A"
		lib_name=$(echo `basename $x | cut -d"_" -f 1` )"_1A_XT_R"$(echo `basename $x | cut -d"_" -f 2 | cut -d"." -f 1`)
		mkdir -p /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name"_XT"/0_Raw

		if [ ${x: -3} == ".gz" ]; then
			
			cp $x /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name"_XT"/0_Raw/$lib_name.fastq.gz &		

		else

			submitjob 1 -c 12 -m 50 -q private -g bashia02 pigz -c -p 12 -k $x \> /sc/orga/scratch/kumara22/microbiome-output/samples/$run_id_exact/$sample_name/0_Raw/$lib_name.fastq.gz		
		fi	
		
	done		
		
done

END
