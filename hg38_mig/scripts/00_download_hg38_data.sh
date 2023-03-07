#!/bin/bash 
	#to run this script: chmod +x script.sh; ./script.sh
	#!/bin/sh does not work with my terminal en msi of David.
	#if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
	#you can save the output and the errors
		#./master_inputs_slim.sh > master_inputs_slim.out #only output
		#./master_inputs_slim.sh 2> error.out #only error
		#./master_inputs_slim.sh > master_inputs_slim.out 2> error.out #both in different files
		#./master_inputs_slim.sh > master_inputs_slim.out 2>&1 #both in the same file
		#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/




######################################################################################
############################### MASTER SCRIPT FOR FIGURES ############################
######################################################################################

#Master script for copying and making figures and tables of main text and supplementary. We run all the scripts previously created to obtain tables, figures and finally compile the tex files

#save the working directory
path_working_dir="/home/dftortosa/singularity/ihs_deep_learning/slim_simulations/scripts"

#set the working directory
cd $path_working_dir