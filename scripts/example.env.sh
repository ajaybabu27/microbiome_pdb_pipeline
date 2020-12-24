#!/bin/bash

# Enter the directory where the local REPO is stored
export REPO_DIR=""

# Load required modules
module load ruby/2.2.0

#Fill in the read/write privilege PDB password here
export PDBPASS=""

# Setting up environments for Qiime2
export FONTCONFIG_PATH=/etc/fonts
export LC_ALL=aa_DJ.utf8
export LANG=C.UTF-8
