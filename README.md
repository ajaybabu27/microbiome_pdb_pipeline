# Microbiome PDB Pipeline
Computational workflows for analyzing 16S and metagenomics data 
## Requirements

Pipeline currently works on Minerva HPC environment. With little tweaking it should be portable in any environment of choice. Users need access to Minerva and InfectiousDisease project. Users also need access to PathogenDB database to query and upload meta and final data points for further query. NOTE: This pipeline only contains components for the 16S workflow and the metagenomics workflow will be added in the future. 
## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure your environment first using the included script:

    $ cp scripts/example.env.sh scripts/env.sh
    $ $EDITOR scripts/env.sh    

All the default variables should work for any Minerva user with appropriate permissions. Then, you can source the script into your shell and install required gems locally into the `vendor/bundle` directory as follows:

    $ source scripts/env.sh
    $ bundle install --deployment

When this is complete, you should be able to run `rake` to kick off the pipeline as follows. However, first read **[Environment variables](#required-environment-variables)** below, as certain tasks require more variables to be set before being invoked.

    $ rake -T                    # list the available tasks
    $ rake $TASK_NAME            # run the task named $TASK_NAME
    $ FOO="bar" rake $TASK_NAME  # run $TASK_NAME with FOO set to "bar"

When firing up the pipeline in a new shell, **remember to always `source scripts/env.sh` before running `rake`.**

## Workflows

There are two main workflows,
1) [16S Workflow](https://github.com/ajaybabu27/microbiome_pdb_pipeline/blob/master/docs/16s_rake_pipeline.md)
2) Metagenomics Workflow

Click on the links above to read about more details about the tasks that are part of these workflows. 

### Required environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ RUN_ID=H434 rake create_manifest

If a required environment variable isn't present when a task is run and there is no default value, rake will abort with an error message.

Variable             | Required by                                             | Default | Purpose
---------------------|---------------------------------------------------------|---------|-----------------------------------
`RUN_ID`             | all QC tasks                                            | (none)  | Illumina RUN ID 
`FASTQ_DIR`          | `create_manifest_file`,`run_kraken`,`run_MC_QC`         | (none)  | Directory containing the raw fastq files 
`QC_DIR`             | all QC tasks                                            | (none)  | Directory were all the QC output files will be saved
`READS_THRESHOLD`    | `create_postQC_biome_file` & all analysis tasks         | (none)  | Set number of reads for sample filtering threshold and other rarefaction analysis
`READ_IDS`           | `prepare_Qiime_analysis`								   | (none)  | Illumina RUN IDs for merging runs. Enter as Comma Seperated Values. 
`ANALYSIS_DIR`       | all analysis tasks                                      | (none)  | Directory were all the files generated from various analyses will be stored


### Running as a `bsub` task

You may also want to run the pipeline as a non-interactive job on the cluster.  The benefit of this approach is that you can reserve specific resources in advance to decrease the likelihood of the job running out of memory or exceeding other system limits.  For this, the `scripts/example.microbiome_pdb_wrapper` should be copied, modified as appropriate, and then can be submitted with `bsub` as in the following example:

    $ bsub -R 'rusage[mem=4000] span[hosts=1]' -P acc_InfectiousDisease -W "24:00" \
            -L /bin/bash -q premium -n 12 -J H434 \
            -o "%J.stdout" -eo "%J.stderr" \
			microbiome_pdb_wrapper \
            RUN_ID=H434 \
            FASTQ_DIR=/sc/arion/project/InfectiousDisease/microbiome_output/samples \
	    QC_DIR=/sc/arion/project/InfectiousDisease/microbiome_output/samples \
            run_mc_qc
