#!/bin/bash
 
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=hap_map_calc
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=denard
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
#SBATCH --mail-type=end
### Optional. Specify email address to use for notification
#SBATCH --mail-user=dftortosa@email.arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### OPTIONAL. Select the cores bought by David.
#SBATCH --qos=user_qos_denard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=115
### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287)
#SBATCH --cpus-per-task=1
### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
###SBATCH --mem=400gb
### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory).
#SBATCH --mem-per-cpu=5gb
### set the constraint for high memory nodes in case you use a lot of memory per node
###SBATCH --constraint=hi_mem
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=240:00:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed. Recommendation from Puma is use "module --ignore-cache load "singularity"" instead of "module load singularity". It is also possible your cache file is out-of-date.
module load singularity
### change to your script’s directory
cd /xdisk/denard/dftortosa/climate_adaptation_met_genes/hg38_mig/
### Run your work
singularity run 01_ubuntu_20_04_hg38_mig_hap_map.sif
