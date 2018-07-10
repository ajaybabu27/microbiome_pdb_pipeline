# Microbiome PDB Pipeline
Computational workflows for analyzing 16S and metagenomics data 
## Requirements

Pipeline currently works on Minerva HPC environment. With little tweaking it should be portable in any environment of choice. Users need access to Minerva and InfectiousDisease project. Users also need access to PathogenDB database to upload meta and final data points for further query. 
## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure your environment first using the included script:

    $ cp scripts/example.env.sh scripts/env.sh
    $ $EDITOR scripts/env.sh    

All the default variables should work for any Minerva user with appropriate permissions. Then, you can source the script into your shell and install required gems locally into the `vendor/bundle` directory as follows:

    $ source scripts/env.sh
    $ bundle install --deployment

When this is complete, you should be able to run `rake` to kick off the pipeline as follows. However, first read **[Environment variables](#required-environment-variables)** below, as certain tasks require more variables to be set before being invoked.  A description of the typical sequence for assembling a genome is described below in **[Tasks](#tasks)**.

    $ rake -T                    # list the available tasks
    $ rake $TASK_NAME            # run the task named $TASK_NAME
    $ FOO="bar" rake $TASK_NAME  # run $TASK_NAME with FOO set to "bar"

When firing up the pipeline in a new shell, **remember to always `source scripts/env.sh` before running `rake`.**

##Workflows

There are two main workflows,
1) 16S Workflow
2) Metagenomics Workflow

Click on the links above to read about more details about the tasks that are part of these workflows. 

### Required environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ RUN_ID=H434 rake create_manifest

If a required environment variable isn't present when a task is run and there is no default value, rake will abort with an error message.

Variable             | Required by                                             | Default | Purpose
---------------------|---------------------------------------------------------|---------|-----------------------------------
`RUN_ID`             | all tasks                                               | (none)  | Illumina RUN ID
`QC_DIR`        |   | (none)  | The ID of the job on the SMRTPortal with your reads.
