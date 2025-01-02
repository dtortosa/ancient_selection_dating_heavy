#!/bin/bash
chmod +x ./01_models_benchmark.py
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_0.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_1.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_2.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_3.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_4.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_5.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_6.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_7.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_8.slurm;
#sbatch ./slurm_files/01_models_benchmark_elastic_net_cv_partition_9.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_0.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_1.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_2.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_3.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_4.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_5.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_6.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_7.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_8.slurm;
#sbatch ./slurm_files/01_models_benchmark_random_forest_cv_partition_9.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_0.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_1.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_2.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_3.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_4.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_5.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_6.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_7.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_8.slurm;
sbatch ./slurm_files/01_models_benchmark_xgboost_cv_partition_9.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_0.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_1.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_2.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_3.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_4.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_5.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_6.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_7.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_8.slurm;
sbatch ./slurm_files/01_models_benchmark_neural_nets_cv_partition_9.slurm;

n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}')
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel
#chmod +x ./00b_master_bash_script.sh; ./00b_master_bash_script.sh > 00b_master_bash_script.out 2>&1
