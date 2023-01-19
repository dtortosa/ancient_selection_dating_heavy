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

#run several scripts each one runnning an optuna process with different seeds but all using the same relational database
#seeds changed in the second run
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=8765 --n_trials=20000 &
sleep 30s
	#we need delay the next process in order to leave time the first one to create
	#the database and avoid errors because two processes are creating the same db
	#https://stackoverflow.com/questions/49944364/how-to-run-two-commands-but-with-a-delay-on-the-second-command-without-stopping
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=4572 --n_trials=20000 &
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=5234 --n_trials=20000 &
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=0976 --n_trials=20000 &
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=3413 --n_trials=20000