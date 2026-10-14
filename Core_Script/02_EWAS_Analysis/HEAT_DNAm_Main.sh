#!/bin/bash

#SBATCH --job-name=HEAT_DNAm_MAIN # Set a Job name
#SBATCH --partition=long # Set a partition(time) or queue to submit to
#SBATCH --account=generic # Set a user account to submit as
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yu.zhao@isglobal.org # Where to send mail events
#SBATCH --ntasks=1 # Run a single tasks
#SBATCH --cpus-per-task=8 
#SBATCH --mem=100gb # Job memory request
#SBATCH --output=~/010_EWAS_temp_YZ/results/log/test_%j.log # Standard output and error log ~ change to your own path

# Load the module environment suitable for the job.
source ~/.bashrc
source activate heat_Pla_DNAm #name of your configured environment with installed required packages

# And finally run the job calling R as you wold do using system R installation
R CMD BATCH ~/010_EWAS_temp_YZ/script/EWAS_Analysis/HEAT_DNAm_MAIN.R # ~ change to your own path