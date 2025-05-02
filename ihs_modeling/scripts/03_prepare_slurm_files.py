#!/usr/bin/env python3.9
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

#this script creates slurm files.




#######################################
# region INITIAL ANOTATIONS AND STEPS #
#######################################


###########
# imports #
###########

from itertools import product
import textwrap
import os
import numpy as np

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

# endregion






#############################
# region create slurm files #
#############################
print_text("create slurm files", header=1)
print_text("preparations", header=2)
print_text("open folder", header=3)
run_bash("\
    mkdir \
        -p \
        ./scripts/slurm_files_across_pops/; \
")

print_text("get a list of all populations for which we have average iHS across gene windows", header=3)
print_text("list of populations", header=4)
list_ihs_files = os.listdir("../../../method_deep_heavy_analyses/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows/")
    #this is the folder where we have the mean iHS files for all populations. These were generated for the MDR paper

print_text("extract the POP name for each file and then get the unique names", header=4)
#i=list_ihs_files[0]
list_pops = np.unique([i.split("_")[0] for i in list_ihs_files])


print_text("create a slurm file for each combination", header=2)
#pop=list_pops[0]
for pop in list_pops:

    print_text(f"starting slurm file for {pop}", header=3)
    slurm_file_content = f"""
    #!/bin/bash
    # --------------------------------------------------------------
    ### PART 1: Requests resources to run your job.
    # --------------------------------------------------------------
    ###info about slurm commands: https://slurm.schedmd.com/pdfs/summary.pdf
    ### Optional. Set the job name
    #SBATCH --job-name=dnn_ihs_modeling_{pop}
    ### Optional. Set the output filename.
    ### SLURM reads %x as the job name and %j as the job ID
    #SBATCH --output=%x-%j.out
    #SBATCH --error=%x-%j.err
    ### Optional. Request email when job begins and ends
    ### SBATCH --mail-type=ALL
    ### Optional. Specify email address to use for notification
    ### SBATCH --mail-user=dftortosa@gmail.com
    ### REQUIRED. Set the partition for your job.
    #SBATCH --partition=albaicin
    ### REQUIRED. Set the number of cores and nodes that will be used for this job. It seems that each node has 56 cores, so you have to adjust accordingly
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --ntasks-per-node=1
    ### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287)
    #SBATCH --cpus-per-task=15
    ### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
    ###SBATCH --mem=400gb
    ### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory).
    #SBATCH --mem-per-cpu=15gb
    ### set the constraint for high memory nodes in case you use a lot of memory per node. Normal nodes have a 512Gb limit.
    ###SBATCH --constraint=hi_mem
    ### REQUIRED. Specify the time required for this job, hhh:mm:ss
    #SBATCH --time=24:00:00

    
    # --------------------------------------------------------------
    ### PART 2: Executes bash commands to run your job
    # --------------------------------------------------------------
    module load singularity
    ### change to your directory
    cd /home/UGR002/dsalazar/climahealth/ihs_modeling
        #/home is the stable directory, while scratch is where results of analyses can be stored temporary, as stuff gets removed after 20 days
    ### Run your work
    singularity exec ./containers/03_explore_selected_model_class.sif ./scripts/02_ihs_modeling_across_pops.py --pop_name={pop} --n_iterations=10 > ./02_ihs_modeling_across_pops_{pop}.out 2>&1
    """

    #remove the first empty line and the spaces at the beginning of the lines
    slurm_file_content_update = textwrap.dedent("\n".join(slurm_file_content.splitlines()[1:]))
        #slurm_file_content.splitlines() Splits the multi-line string slurm_file_content into a list of individual lines. Indeed, "slurm_file_content.splitlines()[0]" is the first empty line, while "slurm_file_content.splitlines()[1]" is the line with the shebang.
        #so we select all lines except the first one (i.e., the empty one)
        #Joins the remaining lines back into a single string, with each line separated by a newline character (\n).
        #using "textwrap.dedent", take the result and remove any common leading whitespace (including tabs) from all lines in the string.

    # Save the SLURM file
    with open(f"./scripts/slurm_files_across_pops/{pop}.slurm", "w") as slurm_file:
        slurm_file.write(slurm_file_content_update)

# endregion







#############################
# region MASTER BASH SCRIPT #
#############################
print_text("create master bash script", header=1)
print_text("set output file", header=2)
output_file = "./scripts/slurm_files_across_pops/master_bash_script.sh"

print_text("Iterate over the combinations using also the index of each combination", header=2)
#idx=0; pop=list_pops[0]
for idx, pop in enumerate(list_pops):
    
    #create the line to add
    line = f"sbatch {pop}.slurm;\n"

    #write or append to the file
    if idx == 0:
        
        #for the first item, create the file and write the line
        with open(output_file, "w") as file:
            file.write("#!/bin/bash\n") #Add the shebang at the top and
            file.write("chmod +x ../02_ihs_modeling_across_pops.py\n") #make the main script executable
            file.write(line) #add the line of the dataset - phenotype combination
    else:

        #append the line of the selected dataset - phenotype combination, as it is not the first one and the file has been already created
        with open(output_file, "a") as file:
            file.write(line)

#finish the file adding the last line
with open(output_file, "a") as file:
    file.write("n_jobs=$(squeue -u dsalazar | awk '{if(NR!=1){count++}}END{print count}');\necho 'WE HAVE ' $n_jobs ' jobs'\n#to stop all jobs #squeue -u dsalazar | awk '{if(NR!=1){print $1}}' | xargs scancel\n#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1")

# endregion






print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating_heavy_analyses/dating_climate_adaptation/ihs_modeling/
#chmod +x ./scripts/03_prepare_slurm_files.py
#singularity exec ./containers/03_explore_selected_model_class.sif ./scripts/03_prepare_slurm_files.py > ./03_prepare_slurm_files.out 2>&1
#grep -Ei 'error|false|fail' ./03_prepare_slurm_files.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
