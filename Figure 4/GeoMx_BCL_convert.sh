#!/bin/bash
#SBATCH -J BCL_convert				# A single job name for the array
#SBATCH -n 8                                            # Number of cores
#SBATCH -N 1                                            # All cores on one machine
#SBATCH --mem 64G                                     # in MB
#SBATCH --qos castles				# Request Castles Node
#SBATCH -t 0-05:00                                      # Maximum execution time (D-HH:MM)
#SBATCH --mail-type=FAIL,END                            # Tell me when
#SBATCH --account=mcmurraj-bta-geomx-data-storage       # Bear account


# This is a slurm job submission script for conversion of BCL files to FASTQ

# Important head of each script
set -e # Fail the script on the first error
module purge; module load bluebear # For reproducibility. Get the default apps

## Load your required Apps
module load bear-apps/2022a
module load bcl-convert/4.0.3-2el7.x86_64

## Run the actual script ========================================================
#### ENSURE THAT NO-LANE-SPLITTING IS SET TO FALSE IN THE SAMPLE SHEET ####
bcl-convert --bcl-input-directory /rds/projects/m/mcmurraj-bta-geomx-data-storage/Count_Data/Illumina/GB500923-AH \
--output-directory /rds/projects/m/mcmurraj-bta-geomx-data-storage/Count_Data/Illumina/GB500923-AH/Data/Intensities/BaseCalls/GB500923-AH_FASTQ \
--sample-sheet /rds/projects/m/mcmurraj-bta-geomx-data-storage/Count_Data/Illumina/GB500923-AH/SampleSheet.csv
