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
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=75558 --n_trials=20000 &
#--optuna_seed=87918
#--optuna_seed=61581
#--optuna_seed=30914
#--optuna_seed=14744
#--optuna_seed=2508
#--optuna_seed=
#--optuna_seed=8765
#--optuna_seed=0432
sleep 30s
	#we need delay the next process in order to leave time the first one to create
	#the database and avoid errors because two processes are creating the same db
	#https://stackoverflow.com/questions/49944364/how-to-run-two-commands-but-with-a-delay-on-the-second-command-without-stopping
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=80669 --n_trials=20000 &
#--optuna_seed=64345
#--optuna_seed=65582
#--optuna_seed=93557
#--optuna_seed=33956
#--optuna_seed=8364
#--optuna_seed=0327
#--optuna_seed=4572
#--optuna_seed=9456
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=78154 --n_trials=20000 &
#--optuna_seed=52844
#--optuna_seed=71673
#--optuna_seed=96611
#--optuna_seed=61584
#--optuna_seed=2094
#--optuna_seed=8869
#--optuna_seed=5234
#--optuna_seed=0467
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=66568 --n_trials=20000 &
#--optuna_seed=35814
#--optuna_seed=24335
#--optuna_seed=38039
#--optuna_seed=61242
#--optuna_seed=5813
#--optuna_seed=5391
#--optuna_seed=0976
#--optuna_seed=1624
python3.9 /opt/scripts/flex_sweep_continuous_modeling_yoruba.py --optuna_seed=95550 --n_trials=20000
#--optuna_seed=84470
#--optuna_seed=85753
#--optuna_seed=29454
#--optuna_seed=43316
#--optuna_seed=0636
#--optuna_seed=0139
#--optuna_seed=3413
#--optuna_seed=8563