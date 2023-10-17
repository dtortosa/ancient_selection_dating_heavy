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




########################################################################
######## CHECK THE OLD AND NEW VCF FILES AFTER VEP ARE THE SAME ########
########################################################################

#We have run again VEP installing it this time with the option "f" in INSTALL.pl, meaning that fastas have been downloaded during installation. These are not the ancestral fastas, see next lines:

#By pointing VEP to a FASTA file (or directory containing several files), it is possible to retrieve reference sequence locally when using --cache or --offline. This enables VEP to retrieve HGVS notations (--hgvs), check the reference sequence given in input data (--check_ref), and construct transcript models from a GFF or GTF file without accessing a database

#IMPORTANT:
    #FASTA files can be set up using the installer; files set up using the installer are automatically detected by VEP when using --cache or --offline; you should not need to use --fasta to manually specify them.
    #Therefore, even if we do not use --fasta, the fasta files could be still used by VEP because we use both --cache and --offline, so better to download them just in case so VEP can use them is required. The whole container weights 1.7GB, so it acceptable.

#see next link for more info
    #https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta




##################
#### Starting ####
##################


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




######################
# do the comparisons #
######################
print_text("do the comparisons", header=1)
print_text("run function across chromosomes", header=2)
print_text("define function", header=3)
#chrom=1
def check_per_chrom(chrom):

    print_text("Starting chromosome " + str(chrom), header=4)
    #convert the chrom number to string
    chrom_str = str(chrom)

    print_text("use cmp to compare byte by byte both files and then obtain the exit status", header=4)
    exit_status_check = run_bash(" \
        cmp \
            --silent \
            <( \
                bcftools view \
                    --no-header \
                    ./results/00_vep_vcf_files/chr" + chrom_str + "/1kGP_high_coverage_Illumina.chr" + chrom_str + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz \
            ) \
            <( \
                bcftools view \
                    --no-header \
                    ./results/00_vep_vcf_files_old/chr" + chrom_str + "/1kGP_high_coverage_Illumina.chr" + chrom_str + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz \
            ); \
        echo $?", return_value=True).strip()
    #we are using bash substitution
        #Piping the stdout of a command into the stdin of another is a powerful technique. But, what if you need to pipe the stdout of multiple commands? This is where process substitution comes in.
            #You have to use "<(command_list)"
            #Process substitution uses /dev/fd/<n> files to send the results of the process(es) within parentheses to another process. [1]
            #Caution 
                #There is no space between the the "<" or ">" and the parentheses. Space there would give an error message.
            #https://tldp.org/LDP/abs/html/process-sub.html
        #It seems that "/dev/fd" includes stdin and stdout.
            #For example, the reason /dev/stdin and friends exist is that sometimes a program requires a file name, but you want to tell it to use a file that's already open (a pipe, for instance). So you can pass /dev/stdin to tell the program to read its standard input.
            #https://askubuntu.com/a/183243
        #we can use two process substitutions for the two inputs of cmp. We could even use subprocess for one input and use an saved file as another input
            #cmp --silent <(echo "hola mundo") <(echo "hola mundo"); echo $?
                #exist status -> 0
            #cmp --silent <(echo "hola mundo") <(echo "hola mundos"); echo $?
                #exist status -> 1
            #echo "hola mundo">file_1.txt; cmp --silent file_1.txt <(echo "hola mundo"); echo $?; rm file_1.txt
            #https://stackoverflow.com/a/48672774/12772630

    print_text("convert to boolean", header=4)
    test_bool = [True if exit_status_check=="0" else False][0]
    print(test_bool)

    print("save the result")
    return test_bool


print_text("create list with all chromosomes", header=3)
print("get chromosome names")
chromosomes = [str(i) for i in range(1, 23, 1)]

print("we are going to analyze 22 chromosomes?")
print((len(chromosomes) == 22))

print("See them")
print(chromosomes)


print_text("run it in parallel", header=3)
print_text("open the pool", header=4)
import multiprocessing as mp
pool = mp.Pool(8)

print_text("run function across pandas rows", header=4)
bool_results = pool.map(check_per_chrom, chromosomes)

print_text("close the pool", header=4)
pool.close()


print_text("now check that we have true for all chromosomes", header=3)
if (len(bool_results)==22) & (sum(bool_results) == len(bool_results)):
    print("We have all TRUE, so we can just remove the previous VCF files because they are identical to the new ones")
    print("YES! GOOD TO GO!")
else:
    raise ValueError("ERROR! FALSE! We have at least one chromosome for which the old and the new VCF files (after adding the download of fasta ('f' option in INSTALL.pl for VEP) are not the same")


###remove files???? too much space...