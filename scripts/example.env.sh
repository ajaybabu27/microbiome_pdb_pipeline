# Must be configured by the end user
export REPO_DIR=""

# Load required modules
module load ruby/2.2.0

# Defaults will probably work for these
export TMP="/sc/orga/scratch/$USER/tmp"

#Fill in the PDB password here
export PDBPASS=""

# If running from interactive1/interactive2, need to run requests through internal HTTP proxy
export HTTP_PROXY="http://proxy.mgmt.hpc.mssm.edu:8123"
