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



####################################
######## CREATE SLURM FILES ########
####################################

# This script generates the slurm files to run the models benchmarking in the cluster.



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



###################################################
# create slurm files iterating across model names #
###################################################
print_text("create slurm files iterating across model names", header=1)
print_text("list of models", header=2)
list_models = ["elastic_net", "random_forest", "xgboost", "neural_nets"]
print(list_models)

print_text("loop", header=2)
#model=list_models[0]
for model in list_models:

    #print the model
    print_text(model, header=3)

    #create a string
    print_text("create the string", header=3)
    string_to_write = " \
        #!/bin/bash \n \
        # -------------------------------------------------------------- \n \
        ### PART 1: Requests resources to run your job. \n \
        # -------------------------------------------------------------- \n \
        ###info about slurm commands: https://slurm.schedmd.com/pdfs/summary.pdf \n \
        ### Optional. Set the job name \n \
        #SBATCH --job-name=01_models_benchmark_" + model + " \n \
        ### Optional. Set the output filename. \n \
        ### SLURM reads %x as the job name and %j as the job ID \n \
        #SBATCH --output=%x-%j.out \n \
        #SBATCH --error=%x-%j.err \n \
        ### Optional. Request email when job begins and ends \n \
        ### SBATCH --mail-type=ALL \n \
        ### Optional. Specify email address to use for notification \n \
        ### SBATCH --mail-user=dftortosa@gmail.com \n \
        ### REQUIRED. Set the partition for your job. \n \
        #SBATCH --partition=albaicin \n \
        ### REQUIRED. Set the number of cores and nodes that will be used for this job. It seems that each node has 56 cores, so you have to adjust accordingly \n \
        #SBATCH --nodes=1 \n \
        #SBATCH --ntasks=50 \n \
        #SBATCH --ntasks-per-node=50 \n \
        ### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287) \n \
        #SBATCH --cpus-per-task=1 \n \
        ### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits \n \
        ###SBATCH --mem=400gb \n \
        ### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory). \n \
        #SBATCH --mem-per-cpu=3gb \n \
        ### set the constraint for high memory nodes in case you use a lot of memory per node. Normal nodes have a 512Gb limit. \n \
        ###SBATCH --constraint=hi_mem \n \
        ### REQUIRED. Specify the time required for this job, hhh:mm:ss \n \
        #SBATCH --time=72:00:00 \n \
        # -------------------------------------------------------------- \n \
        ### PART 2: Executes bash commands to run your job \n \
        # -------------------------------------------------------------- \n \
        module load singularity \n \
        ### change to your script directory \n \
        cd /home/UGR002/dsalazar/climahealth/ihs_modeling \n \
            #/home is the stable directory, while scratch is where results of analyses can be stored temporary, as stuff gets removed after 20 days \n \
        ### Run your work \n \
        singularity exec ./containers/01_models_benchmark.sif ./scripts/01_models_benchmark.py --model_name='" + model + "' > ./scripts/01_models_benchmark_" + model + ".out 2>&1"

    print_text("remove spaces at the beginning", header=3)
    string_to_write = string_to_write.replace("         ", "")

    print_text("Open the file in write mode", header=3)
    with open("./scripts/slurm_files/01_models_benchmark_" + model + ".slurm", 'w') as file:
        # Write the string to the file
        file.write(string_to_write)

    print_text("end of the model: " + model, header=3)



###########################
# prepare the bash script #
###########################
print_text("create bash file", header=1)

print_text("open bash string", header=2)
bash_string="#!/bin/bash\nchmod +x ./01_models_benchmark.py\n"

print_text("add the models", header=2)
for model in list_models:
    bash_string=bash_string+"sbatch ./slurm_files/01_models_benchmark_"+model+".slurm;\n" 

print_text("add check", header=2)
bash_string = bash_string+"\nn_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}')\n"
bash_string = bash_string+"echo 'WE HAVE ' $n_jobs ' jobs'\n"
bash_string = bash_string+"#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel\n"
bash_string = bash_string+"#chmod +x ./00b_master_bash_script.sh; ./00b_master_bash_script.sh > 00b_master_bash_script.out 2>&1\n" 

print_text("see bash script", header=2)
print(bash_string)

print_text("Open the file in write mode", header=2)
with open("./scripts/00b_master_bash_script.sh", 'w') as file:
    file.write(bash_string)



###########################
# check we have all files #
###########################
print_text("check we have all files", header=1)
print_text("we should have 1 file per model", header=2)
n_files=run_bash(" \
    cd ./scripts/slurm_files; \
    n_files=$(ls -1 | wc -l); \
    echo $n_files", return_value=True \
).strip()
if(n_files!=str(len(list_models))):
    raise ValueError("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL models")


#check we have all batches in the master bash script
sbatch_appearance = run_bash(" \
    grep \
        'sbatch' \
        --only-matching \
        ./scripts/00b_master_bash_script.sh", return_value=True \
).strip().split("\n")
    #--only-matching: show only nonempty parts of lines that match
        #we get only "sbatch" as many times as it appears even if it is repeated in the same row
        #https://stackoverflow.com/a/3249761
if(len(sbatch_appearance)!=len(list_models)):
    raise ValueError("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL models")



##########
# FINISH #
##########
print_text("finish", header=1)
#chmod +x ./scripts/00a_models_benchmark_slurm_files.py; ./scripts/00a_models_benchmark_slurm_files.py > ./scripts/00a_models_benchmark_slurm_files.out 2>&1
