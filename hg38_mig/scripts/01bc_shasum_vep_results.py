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



print_text("run loop across chromosomes", header=2)
print_text("open empty list to save booleans from comparisons", header=3)
bool_results=list()


print_text("run the loop", header=3)
#chrom=1
for chrom in range(1,23,1):

    print_text("Starting chromosome " + str(chrom), header=4)
    #convert the chrom number to string
    chrom_str = str(chrom)

    print("calculate sh256 digest of the new and old VCF files")
    digest_new = run_bash(" \
        bcftools view \
            --no-header \
            ./results/00_vep_vcf_files/chr" + chrom_str + "/1kGP_high_coverage_Illumina.chr" + chrom_str + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -n 10 | \
        shasum --algorithm 256", return_value=True).strip()
    digest_old = run_bash(" \
        bcftools view \
            --no-header \
            ./results/00_vep_vcf_files_old/chr" + chrom_str + "/1kGP_high_coverage_Illumina.chr" + chrom_str + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -n 10 | \
        shasum --algorithm 256", return_value=True).strip()
        #you have several hash functions available, but you should use sha256
            #Hash functions are one-way functions that takes a message as its argument. This message can be of any size, but the function always returns a fixed size hash. A hash is considered impossible (within the bounds of practicality) to reverse and to find two different messages with the same hash (called a collision).
        #Differences between some hash functions
            #MD5 was invented in the early 1990s and is considered flawed and obsolete by now.
            #SHA1 was also developed in the early 1990s. It is considered stronger than MD5, but not strong enough. Its use is currently being withdrawn from the digital signature on X.509 digital certificates.
            #SHA256 is the currently recommended hash function.
        #https://unix.stackexchange.com/a/260519
        #sha256 is very resistant to collision, i.e., two different files generating the same digest, at least at the moment of writing this code. Therefore, if two files have the same digest, we can consider them as the same, unaltered.

    print("check both digests are the same")
    test_bool = digest_new == digest_old
    print(test_bool)

    print("save the result")
    bool_results.append(test_bool)


print_text("now check that we have true for all chromosomes", header=3)
if(sum(bool_results) == len(bool_results)):
    print("We have all TRUE, so we can just remove the previous VCF files because they are identical to the new ones")
    run_bash(" \
        rm \
            --recursive \
            --force \
            ./results/00_vep_vcf_files_old")
else:
    raise ValueError("ERROR! FALSE! We have at least one chromosome for which the old and the new VCF files (after adding the download of fasta ('f' option in INSTALL.pl for VEP) are not the same")
