#!/bin/bash 
	#to run this script: chmod +x script.sh; ./script.sh
	#!/bin/sh does not work with my terminal en msi of David.
	#if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
	#you can save the output and the errors
		#./bash_optuna_processes.sh > bash_optuna_processes.out #only output
		#./bash_optuna_processes.sh 2> error.out #only error
		#./bash_optuna_processes.sh > bash_optuna_processes.out 2> error.out #both in different files
		#./bash_optuna_processes.sh > bash_optuna_processes.out 2>&1 #both in the same file
		#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/

#run script for calculating the distance of genes to the interest genes of each selective pressure
python3.9 ./scripts/00a_pressure_distance_calc_v1.py --selective_pressure="bat" --n_cores=5 > ./scripts/00a_pressure_distance_calc_v1_bat.out 2>&1
python3.9 ./scripts/00a_pressure_distance_calc_v1.py --selective_pressure="smt" --n_cores=5 > ./scripts/00a_pressure_distance_calc_v1_smt.out 2>&1
