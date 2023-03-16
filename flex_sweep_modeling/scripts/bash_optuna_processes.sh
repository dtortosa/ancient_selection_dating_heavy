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
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=14744 --n_trials=20000 &
#--optuna_seed=2508
#--optuna_seed=
#--optuna_seed=8765
#--optuna_seed=0432
sleep 30s
	#we need delay the next process in order to leave time the first one to create
	#the database and avoid errors because two processes are creating the same db
	#https://stackoverflow.com/questions/49944364/how-to-run-two-commands-but-with-a-delay-on-the-second-command-without-stopping
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=33956 --n_trials=20000 &
#--optuna_seed=8364
#--optuna_seed=0327
#--optuna_seed=4572
#--optuna_seed=9456
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=61584 --n_trials=20000 &
#--optuna_seed=2094
#--optuna_seed=8869
#--optuna_seed=5234
#--optuna_seed=0467
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=61242 --n_trials=20000 &
#--optuna_seed=5813
#--optuna_seed=5391
#--optuna_seed=0976
#--optuna_seed=1624
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=43316 --n_trials=20000
#--optuna_seed=0636
#--optuna_seed=0139
#--optuna_seed=3413
#--optuna_seed=8563