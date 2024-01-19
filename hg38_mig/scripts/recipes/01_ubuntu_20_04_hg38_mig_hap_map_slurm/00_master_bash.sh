#!/bin/bash
chmod +x ../01d_hap_map_calcs.py
qsub GBR.slurm; qsub FIN.slurm; qsub CHS.slurm; qsub PUR.slurm; qsub CDX.slurm; qsub CLM.slurm; qsub IBS.slurm; qsub PEL.slurm; qsub PJL.slurm; qsub KHV.slurm; qsub ACB.slurm; qsub GWD.slurm; qsub ESN.slurm; qsub BEB.slurm; qsub MSL.slurm; qsub STU.slurm; qsub ITU.slurm; qsub CEU.slurm; qsub YRI.slurm; qsub CHB.slurm; qsub JPT.slurm; qsub LWK.slurm; qsub ASW.slurm; qsub MXL.slurm; qsub TSI.slurm; qsub GIH.slurm; 
n_jobs=$(squeue -u dftortosa | awk '{if(NR!=1){count++}}END{print count}')
echo 'WE HAVE ' $n_jobs ' jobs'
#to stop all jobs #squeue -u dftortosa | awk '{if(NR!=1){print $1}}' | xargs qdel
#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1