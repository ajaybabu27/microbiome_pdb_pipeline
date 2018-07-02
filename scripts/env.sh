#!/bin/bash

module unload ruby
module load ruby

# Defaults will probably work for these

export TMP="/sc/orga/scratch/$USER/tmp"
export SHARED_DIR="/sc/orga/scratch/$USER/shared_dir"
export CLUSTER="LSF_PSP"

# If running from interactive1/interactive2, need to run requests through internal HTTP proxy
export HTTP_PROXY="http://proxy.mgmt.hpc.mssm.edu:8123"
