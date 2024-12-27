#!/bin/bash
chmod +x ./01_models_benchmark.py
sbatch ./slurm_files/01_models_benchmark_elastic_net.slurm;
sbatch ./slurm_files/01_models_benchmark_random_forest.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets.slurm;

n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}')
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel
#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1
