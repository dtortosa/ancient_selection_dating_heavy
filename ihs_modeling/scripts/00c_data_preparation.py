#!/usr/bin/env python3
# coding: utf-8
    #to run this script: chmod +x script.py; ./script.py
    #!/bin/sh does not work with my terminal en msi of David.
    #if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
    #you can save the output and the errors
        #./script.py > script.out #only output
        #./script.py 2> error.out #only error
        #./script.py > script.out 2> error.out #both in different files
        #./script.py > script.out 2>&1 #both in the same file
        #https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/



########################################################
######## COMPARE DIFFERENT MODELS IN PREDICTION ########
########################################################

#This script will prepare the input data for 01_models_benchmark.py




##############################
#region INITIAL STEPS ########
##############################


########################################
# define function to print text nicely #
########################################

#text="checking function to print nicely"
#header=1
def print_text(text, header=2):
    if header==1:
        print("\n#######################################\n#######################################")
        print(text)
        print("#######################################\n#######################################")
    elif header==2:
        print("\n###### " + text + " ######")
    elif header==3:
        print("\n## " + text + " ##")
    elif header==4:
        print("\n# " + text + " #")
print_text("checking function to print nicely: header 1", header=1)
print_text("checking function to print nicely: header 2", header=2)
print_text("checking function to print nicely: header 3", header=3)
print_text("checking function to print nicely: header 4", header=4)



########################################
# define function to run bash commands #
########################################

#create a wrapper for subprocess.run in order to define a set of arguments and avoid typing them each time. We will ensure that we are using bash always and not sh.
from subprocess import run, PIPE
#command="ls"
def run_bash(command, return_value=False):

    #run the command
    complete_process = run(
        command, 
        shell=True,
        executable="/bin/bash", 
        stdout=PIPE,
        stderr=PIPE, 
        text=True)
    #we have to use popen in order to ensure we use bash, os.system does not allow that
        #shell=True to execute the command through the shell. This means that the command is passed as a string, and shell-specific features, such as wildcard expansion and variable substitution, can be used.
            #THIS IS DANGEROUS IF UNTRUSTED DATA
        #executable="/bin/bash" to ensure the use of bash instead of sh
        #stdout=PIPE to capture the output into an python object. You can also capture the error doing stderr=PIPE. stdout and stderr are the standard output and error
            #you could also use capture_output=True to capture both stdout and stderr
        #text=True will return the stdout and stderr as string, otherwise as bytes
            #https://www.datacamp.com/tutorial/python-subprocess
            #https://docs.python.org/3/library/subprocess.html#subprocess.run

    #this generated a CompletedProcess instance where you can get
        #args: The command and arguments that were run.
        #returncode: The return code of the subprocess.
        #stdout: The standard output of the subprocess, as a bytes object.
        #stderr: The standard error of the subprocess, as a bytes object.

    #if the return code is 0, i.e., success 
        #https://askubuntu.com/a/892605
    if complete_process.returncode==0:

        #if stderr is empty
        if complete_process.stderr=="":

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
        else:

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #print the standard error without stopping
            print(complete_process.stderr)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
    else:
        #print the standard error and stop
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print_text("check behaviour run_bash", header=1)
print_text("see working directory", header=2)
run_bash("pwd")
print_text("list files/folders there", header=2)
run_bash("ls")



##########
# set WD #
##########
print_text("Set the working directory", header=1)
import os
os.chdir("/home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating_heavy_analyses/dating_climate_adaptation/ihs_modeling/")

# endregion






##################################
# region data preparation ########
##################################
print_text("data preparation", header=1)
print_text("Set the window size selected for those variables that are calculated within windows centered around genes", header=2)
gene_window_size = "1000kb"
    #see MDR paper about window selection and why we focused in 1000kb gene windows
print(f"The selected window size is {gene_window_size}")



print_text("load and clean input data", header=2)
print_text("load mean iHS and number of data points", header=3)
#ihs is the response
response="mean_ihs_" + gene_window_size
print(f"The response is {response}")
import pandas as pd
ihs_yoruba_mean = pd.read_csv( \
    "../../../method_deep_heavy_analyses/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows/YRID_mean_ihs_gene_windows_final_v1.txt.gz", \
    sep="\t", \
    header=0, \
    low_memory=False, \
    compression="gzip").loc[:,["gene_id", response]]
print(ihs_yoruba_mean)
ihs_yoruba_n = pd.read_csv( \
    "../../../method_deep_heavy_analyses/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows/YRID_n_ihs_gene_windows_final_v1.txt.gz", \
    sep="\t", \
    header=0, \
    low_memory=False, \
    compression="gzip").loc[:,["gene_id", "n_ihs_"+gene_window_size]]
print(ihs_yoruba_n)


print_text("load flex-sweep probability and most predictors into pandas", header=3)
#we will not use flexsweep just the predictors
import pandas as pd
exclude_columns = ["predicted_class", "prob(sweep)", "number_thermogenic_1000kb", "number_vips_1000kb"]
main_predictors = pd.read_csv( \
    "../flex_sweep_modeling/data/flex_sweep_closest_window_center.txt.gz", \
    sep=",", \
    header=0, \
    low_memory=False, \
    compression="gzip", \
    usecols=lambda column: column not in exclude_columns
)
    #usecols is a parameter of the pd.read_csv function that specifies which columns to read from the CSV file.
    #lambda column: column not in exclude_columns is a lambda function that takes a column name as input and returns True if the column name is not in the exclude_columns list, and False otherwise.
    #When pd.read_csv is called with this lambda function as the usecols parameter, it will include only the columns for which the lambda function returns True.
#check we have removed the flexsweep columns
if(len([col for col in exclude_columns if col in main_predictors.columns]) !=0):
    raise ValueError("ERROR: FALSE! WE HAVE NOT EXCLUDED THE COLUMNS WE WANTED TO EXCLUDE")
else:
    print(main_predictors)


print_text("load BAT and SMT data, only 1% percentile for now", header=3)
p_value_percentile = 1
#the 1% means we have selected those genes in top 1% of p-value. In the case of the BAT, this is the BAT connectome used in the paper. So we have the distance of each coding gene to the closest gene inside the top 1% of the BAT and SMT connectomes.
bat_distance = pd.read_csv( \
    "../flex_sweep_modeling/data/bat_distance/bat_data_percentile_" + str(p_value_percentile) + ".tsv",
    sep='\t', 
    header=0, 
    low_memory=False)
smt_distance = pd.read_csv( \
    "../flex_sweep_modeling/data/smt_distance/smt_data_percentile_" + str(p_value_percentile) + ".tsv",
    sep='\t', 
    header=0, 
    low_memory=False)
print(bat_distance)
print(smt_distance)


print_text("merge all data.frames", header=3)
print_text("make list with all DFs", header=4)
list_dataframes = [ihs_yoruba_mean, ihs_yoruba_n, main_predictors, bat_distance, smt_distance]
print(list_dataframes)


print_text("merge all of them with reduce", header=4)
from functools import reduce
final_data_yoruba = reduce( \
    lambda x, y: \
        pd.merge( \
            left=x, \
            right=y, \
            how="left", \
            on="gene_id"), \
    list_dataframes \
)
        #reduce applies a function of two arguments cumulatively to the items of a sequence, from left to right, so as to reduce the sequence to a single value. For example, reduce(lambda x, y: x+y, [1, 2, 3, 4, 5]) calculates ((((1+2)+3)+4)+5).
        #therefore, in the first iteration, we merge the first DF of the list (left) with the second (right). The resulting DF (left) is then merged with the third DF (right), and so on....
        #the first DF is going to be the original dataset with flex-sweep probability and most of predictors, therefore using left we ensure we use only those rows with gene_id in the original DF.
            #in any case, we should not have NA for distance to BAT and SMT as these have been calculated using the original gene coordinate file, which was in turn used to calculate many of our predictors.
print(final_data_yoruba)

print_text("count the number of NANs in distance BAT and SMT genes", header=4)
n_nans_bat = sum(final_data_yoruba["bat_distance_percentile_1"].isna())
n_nans_smt = sum(final_data_yoruba["smt_distance_percentile_1"].isna())
if(n_nans_bat <= (final_data_yoruba.shape[0]*0.03)) & (n_nans_smt <= (final_data_yoruba.shape[0]*0.03)):
    print(f"We have {n_nans_bat} NANs for the distance to the closest BAT gene")
    print(f"We have {n_nans_smt} NANs for the distance to the closest SMT gene")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 3% OF NANs FOR DISANCE TO INTEREST GENES")

print_text("remove NANs", header=4)
final_data_yoruba = final_data_yoruba.dropna()
print(final_data_yoruba)

print_text("save the data", header=4)
final_data_yoruba.to_csv( \
    "./data/YRID_modeling_dataset_v1.tsv.gz", \
    sep="\t", \
    compression="gzip", \
    index=False \
)

# endregion





##########
# FINISH #
##########
print_text("finish", header=1)
#chmod +x ./scripts/00c_data_preparation.py; ./scripts/00c_data_preparation.py > ./scripts/00c_data_preparation.out 2>&1




