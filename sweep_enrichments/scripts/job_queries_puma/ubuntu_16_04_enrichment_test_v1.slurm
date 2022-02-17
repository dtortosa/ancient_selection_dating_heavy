#!/bin/bash
 
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=sweep_enrichment_met_genes
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=denard
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=dftortosa@email.arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### OPTIONAL. Select the cores bought by David.
#SBATCH --qos=user_qos_denard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=25
### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses, like the iHS calculation with hapbin per population. 
#SBATCH --cpus-per-task=1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=30gb
### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH.
### SBATCH --mem-per-cpu=12gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=20:00:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed. Recommendation from Puma is use "module --ignore-cache load "singularity"" instead of "module load singularity". It is also possible your cache file is out-of-date.
module load singularity
### change to your script’s directory
cd /xdisk/denard/dftortosa/climate_adaptation_met_genes/
### Run your work
singularity run ubuntu_16_04_enrichment_test_v1.sif