#!/bin/bash
chmod +x ../02_ihs_modeling_across_pops.py
sbatch ACBD.slurm;
sbatch ASWD.slurm;
sbatch BEBD.slurm;
sbatch CDXD.slurm;
sbatch CEUD.slurm;
sbatch CHBD.slurm;
sbatch CHSD.slurm;
sbatch CLMD.slurm;
sbatch ESND.slurm;
sbatch FIND.slurm;
sbatch GBRD.slurm;
sbatch GIHD.slurm;
sbatch GWDD.slurm;
sbatch IBSD.slurm;
sbatch ITUD.slurm;
sbatch JPTD.slurm;
sbatch KHVD.slurm;
sbatch LWKD.slurm;
sbatch MSLD.slurm;
sbatch MXLD.slurm;
sbatch PELD.slurm;
sbatch PJLD.slurm;
sbatch PURD.slurm;
sbatch STUD.slurm;
sbatch TSID.slurm;
sbatch YRID.slurm;
n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}');
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel
#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1