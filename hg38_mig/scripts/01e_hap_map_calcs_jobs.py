#!/usr/bin/env python3.10
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



#################################################
######## CREATE SLURM FILES FOR EACH POP ########
#################################################


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
        elif ("[W::bcf_calc_ac] Incorrect number of AC fields at" in complete_process.stderr):

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #print the standard error without stopping
            print("THIS WARNING SHOULD BE ONLY IN DUMMY DATA NOT IN 1KGP DATA as AC field should have only 1 value per line in the data. " + complete_process.stderr)

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
print("\n#######################################\n#######################################")
print("see working directory")
print("#######################################\n#######################################")
run_bash("pwd")
print("\n#######################################\n#######################################")
print("list files/folders there")
print("#######################################\n#######################################")
run_bash("ls")




#############
# pops prep #
#############
print_text("Preparate pedigree data", header=1)
#load original 2504 unrelated samples from phase 3. This includes sample IDs and pop/superpop codes. This is the definitive ped we will use both for pop codes and sample IDs
#Data:
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#Readme
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220
import pandas as pd
original_unrel_ped = pd.read_csv(
    "data/pedigrees/integrated_call_samples_v3.20130502.ALL.panel.txt", 
    sep="\t", 
    header=0, 
    low_memory=False)

#load pedigree of the latest version of the phased data that has sample IDs, sex and parents but no pop names. We will use this to father/mother IDs
#Downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
samples_pedigree = pd.read_csv(
    "data/pedigrees/1kGP.3202_samples.pedigree_info.txt", 
    sep=" ", 
    header=0, 
    low_memory=False)

#load also a pedigree present in the main directory of the high coverage data. This has sample and pop IDs, but parents and sex are different with respect to the pedigree of the new sample. We will use this to compare pop codes with the original ped.
#downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
samples_pedigree_pop = pd.read_csv(
    "data/pedigrees/20130606_g1k_3202_samples_ped_population.txt", 
    sep=" ", 
    header=0, 
    low_memory=False)

#load the final pedigree we will use to select samples per pop
unrelated_samples = pd.read_csv(
    "data/pedigrees/unrelated_samples.tsv", 
    sep="\t", 
    header=0, 
    low_memory=False)

#FOR DETAILS ABOUT THE SELECTION OF PEDIGREES AND THE REMOVAL OF RELATED SAMPLES, LOOK AT "01a_selecting_pedegree.py"




###################
# calculate files #
###################
print_text("Calculate the files", header=1)

#set the number of cores
n_cores=7

#get pop names
pop_names=unrelated_samples["pop"].unique()
if(len(pop_names)!=26):
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE ALL POPS BEFORE CREATING SLURM FILES")

#open folder to save slurm files
run_bash("mkdir -p ./scripts/recipes/01_ubuntu_20_04_hg38_mig_hap_map_slurm/")

#run loop to create the slurm files
#pop="GBR"
for pop in pop_names:

    with open("./scripts/recipes/01_ubuntu_20_04_hg38_mig_hap_map_slurm/"+pop+".slurm", "w") as f:
        f.write("""\
#!/bin/bash 
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=hap_map_calc_""" + pop + """
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=denard
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
#SBATCH --mail-type=end
### Optional. Specify email address to use for notification
#SBATCH --mail-user=dftortosa@email.arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### OPTIONAL. Select the cores bought by David.
#SBATCH --qos=user_qos_denard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=""" + str(n_cores) + """
### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses. For example, you run 22 independent processes for each chromosome and then you also parallelize the calculations inside each chromosome (https://stackoverflow.com/a/51141287)
#SBATCH --cpus-per-task=1
### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
###SBATCH --mem=400gb
### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH. If you get a .core file, this usually means that the program fails because it asked for too much memory, so it creates a record of the working memory at the time that can be used for debugging. MPI jobs will usually create a core file for each task. You should increase memory limits (https://researchcomputing.princeton.edu/support/knowledge-base/memory).
#SBATCH --mem-per-cpu=10gb
### set the constraint for high memory nodes in case you use a lot of memory per node
###SBATCH --constraint=hi_mem
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=90:00:00

# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### change to your script’s directory
cd /xdisk/denard/dftortosa/climate_adaptation_met_genes/hg38_mig/
### create a folder to save the global outputs
mkdir -p ./scripts/01_hap_map_calcs_outputs/global_outputs
### Run your work
singularity exec 01_ubuntu_20_04_hg38_mig_hap_map.sif ./scripts/01d_hap_map_calcs.py --pop='""" + pop + """' --n_cores=""" + str(n_cores) +  """ > ./scripts/01_hap_map_calcs_outputs/global_outputs/""" + pop + """.out 2>&1
""")
    #singularity exec lets you to use a container to run a program outside of the container.
    #we can run a python script if we previously do "chmod +x script.py"
    #then add arguments for the script, in this case, the population name and the number of cores we want
    #save the output in a file
    #https://stackoverflow.com/a/49516807

#now create the a short bash script in order to do qsub for each slurm file
with open("./scripts/recipes/01_ubuntu_20_04_hg38_mig_hap_map_slurm/00_master_bash.sh", "w") as f:
    #open bash script
    f.write("#!/bin/bash\n")
    #give rights to run the python script
    f.write("chmod +x ../01d_hap_map_calcs.py\n")
    #send slurm jobs
    #pop="GBR"
    for pop in pop_names:
        f.write("qsub " + pop + ".slurm; ")
    #count the number of jobs on the queue
    f.write("\nn_jobs=$(squeue -u dftortosa | awk '{if(NR!=1){count++}}END{print count}')\n")
    f.write("echo 'WE HAVE ' $n_jobs ' jobs'\n")
    f.write("#to stop all jobs #squeue -u dftortosa | awk '{if(NR!=1){print $1}}' | xargs qdel\n")
    f.write("#chmod +x 00_master_bash.sh; ./00_master_bash.sh > 00_master_bash.out 2>&1")


###########################
# check we have all files #
###########################
print_text("check we have all files", header=1)

#we should have 1 file per pop (26) plus the mast bash script
n_files=int(run_bash("\
    cd ./scripts/recipes/01_ubuntu_20_04_hg38_mig_hap_map_slurm/; \
    n_files=$(ls | awk 'END{print NR}'); \
    echo $n_files; \
    ", return_value=True).strip())
if(n_files!=len(pop_names)+1):
    raise ValueError("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL POPULATIONS")

#check we have all pops in the master bash script
qsub_appearance = run_bash(" \
    grep \
        'qsub' \
        --only-matching \
        ./scripts/recipes/01_ubuntu_20_04_hg38_mig_hap_map_slurm/00_master_bash.sh", return_value=True).strip().split("\n")
    #--only-matching: show only nonempty parts of lines that match
        #we get only "qsub" as many times as it appears even if it is repeated in the same row
        #https://stackoverflow.com/a/3249761
    #remove the final new line (\n) with strip and the rest with split
if(len(qsub_appearance)!=len(pop_names)):
    raise ValueError("ERROR! FALSE! WE HAVE NOT OBTAINED JOB FILES FOR ALL POPULATIONS")

##########
# FINISH #
##########
print_text("FINISH", header=1)
#chmod +x ./scripts/01e_hap_map_calcs_jobs.py; ./scripts/01e_hap_map_calcs_jobs.py > ./scripts/01e_hap_map_calcs_jobs.out 2>&1
