#!/bin/bash
#SBATCH -J GeoMx_Pipeline				# A single job name for the array
#SBATCH -n 5                                            # Number of cores
#SBATCH -N 1                                            # All cores on one machine
#SBATCH --mem 64000                                     # in MB
#SBATCH -t 0-20:00                                      # Maximum execution time (D-HH:MM)
#SBATCH --qos bbdefault					# Request castles
#SBATCH --mail-type=FAIL,END                            # Tell me when
#SBATCH --account=croftap-stia-atac       # Bear account


# This is a slurm job submission script for the production of file types that can be loaded onto the GeoMx

# Important head of each script
set -e # Fail the script on the first error
module purge; module load bluebear # For reproducibility. Get the default apps

## Load your required Apps
module load bear-apps/2021b
module load GeoMxNGSPipeline/2.3.3.10 # Required package to run script

## Run the actual script ========================================================
### CHANGE THE NUMBER OF THREADS IN THE CONFIG.INI FILE TO 8
geomxngspipeline --in=/rds/projects/c/croftap-visium-manuscript-01/Geomix_synovium/fastqs/Elizabeth_Clay_UOB451525_PlateALL \
--out=/rds/projects/c/croftap-visium-manuscript-01/Geomix_synovium/DCC_files_plateALL \
--ini=/rds/projects/c/croftap-visium-manuscript-01/Geomix_synovium/Annie_2_20231213T1200_REVERSE/Annie_2_20231213T1200_GNP_config.ini \
--save-interim-files=true \
--threads=5 \

## Change "NAME_OF_INDIVIDUAL" to the name of the folder in which you are running
## Change "NAME_OF_PROJECT" to the name of the project, usually prefixed by "GB"
## Change "NAME_OF_CONFIG_FILE" to the name of the configuration file (".ini" file type) which was downloaded from the GeoMx
## None of the files here require double inverted commas surrounding them

