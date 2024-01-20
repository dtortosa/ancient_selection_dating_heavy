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



###################################
######## CHECK THE OUTPUTS ########
###################################


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



######
# wd #
######

#import os
#print(os.getcwd())

#you can also use the "!" character, but only whithin ipython, it does not work if you run the script as a program
#message="hola mundo"
#!echo {message}
#!pwd
    #https://jakevdp.github.io/PythonDataScienceHandbook/01.05-ipython-and-shell-commands.html



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



########################################################
#### Do some checks after analyzing all chromosomes ####
########################################################
print_text("Do some checks after analyzing all chromosomes", header=1)
print_text("create array with all combinations of pops and chromosomes", header=2)
print_text("get pop and chromosome names", header=3)
pop_names=unrelated_samples["pop"].unique()
#pop_names=["GBR", "PUR", "PEL"]
chromosomes = [i for i in range(1, 23, 1)]
print("we are going to analyze 26 pops and 22 chromosomes?")
print((len(pop_names) == 26) & (len(chromosomes) == 22))
print("See them")
print(pop_names)
print(chromosomes)


print_text("get all the combinations but first make a dummy example", header=3)
import itertools

print_text("dummy example to get all possible combinations of two lists", header=4)
dummy_x = ["marbella", "cuzco", "granada"]
dummy_y = [1, 2, 3]
#product get all possible combinations between the two lists
dumm_combinations = [x+"_"+str(y) for x in dummy_x for y in dummy_y]
print(dumm_combinations)
    #first for each each value of X, and then for each value of Y, combine X and Y, so combine X1 with Y1, X1 with Y2, .... X2 with Y1, X2 with Y2 and so on...
    #y has to be converted to string with it is integer
print("Do we have all dummy combinations?")
print(len(dumm_combinations) == len(dummy_x)*len(dummy_y))

print_text("get all combinations from the actual pops and chromosomes", header=4)
full_combinations_pop_chroms = [pop+"_"+str(chrom) for pop in pop_names for chrom in chromosomes]

print("Do we have all combinations of chromosomes and populations?")
print(len(full_combinations_pop_chroms) == len(pop_names) * len(chromosomes))

print("is this equivalent to itertools.product?")
print(\
    full_combinations_pop_chroms == \
    [i[0]+"_"+str(i[1]) for i in list(itertools.product(pop_names, chromosomes))])
    #itertools.product gives a tuple for each combination, so you can extract both elements and bind them with join.
        #https://stackoverflow.com/a/34032549/12772630
        #https://docs.python.org/3/library/itertools.html#itertools.product


print_text("check we do NOT have any errors in the output files of all chromosomes*pops combinations, also calculate the percentage of SNPs lost", header=2)
print_text("run loop across chromosomes*pops combinations", header=3)
snps_lost_percentage = []
snps_remaining_list = []
#combination=full_combinations_pop_chroms[0]
for combination in full_combinations_pop_chroms:
    print_text("Doing combination " + combination, header=4)

    print_text("split the combination name", header=4)
    comb_pop = combination.split("_")[0]
    comb_chrom = combination.split("_")[1]

    print_text("count number of cases with 'error' or 'false' in the output file", header=4)
    count_error_false = run_bash(" \
        grep \
            'error|false' \
            --extended-regexp \
            --ignore-case \
            --count \
            ./scripts/01_hap_map_calcs_outputs/" + comb_pop + "/chr" + comb_chrom + "_" + comb_pop + ".out || \
        [[ $? == 1 ]]", return_value=True).strip()
        #using grep, look for 
            #"error" OR "false" using "|". In order to avoid escaping the symbol, i.e., "\|", we need to use the flag "--extended-regexp"
                #"error" include any string combing after like "errorS", "errores", etc... If the string "error" is present alone or in combination with other strings, you will get a hit
            #ignore the case, so "Error" and "ERROR" are also included
            #get the count, not the rows matching
        #if the exit code of grep indicates error run the code after "||". This is the function of "||", run the code at the right only if the code at the left failed
            #check if the exit code ("$?") is 1, if so, this will give 0 as exist status (i.e., no error and True) and give the count, which is zero, as stdout. As explained below, if the count is zero (stdout=0), the exit code is 1 in grep:
                #0: no error and one or more lines were selected
                #1: no error but no lines were selected
                #>1: an error occurred
            #This is different from other programs where exit code equals to 1 is error, and we coded that accordingly in run_bash.
            #Because of this, in this particular case, we add an additional line in case grep gives non-zero exist status, and avoid error if the exit status is 1.
                #If grep gives "1" as exit status because the string is not present in the file, we run [[ $? == 1 ]]. This will give "0" as exist status if the previous exit status was "1", while maintaining the previous stdout, i.e., the "count=0" because grep did not find the string in the file.
                #If the exit status is >1 and thus, there is an error, this will give "1" as exist status and run_bash will fail, so we are not hiding errors. 
                #If the exist status is 0, "||" avoids running the conditional (I have checked looking for "ancestral"), so we are good.
            #https://unix.stackexchange.com/a/427598
            #https://pubs.opengroup.org/onlinepubs/9699919799/utilities/grep.html#tag_20_55_14
        #https://linuxize.com/post/grep-multiple-patterns/

    print_text("check the count of problematic cases is zero", header=4)
    if count_error_false == "0":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE ERROR/FALSE IN THE OUTPUT FILE OF COMBINATION " + combination)

    print_text("check we have the row of FINISH", header=4)
    check_finish = run_bash(" \
        grep \
            '## FINISH ##' \
            --count \
            ./scripts/01_hap_map_calcs_outputs/" + comb_pop + "/chr" + comb_chrom + "_" + comb_pop + ".out || \
        [[ $? == 1 ]]", return_value=True).strip()
        #check we have the output indicating finish with a specific way, in uppercase and separated by "#", so we avoid the flag "--ignore-case". Ask for the count.
        #as in the previous case, we add a line to avoid errors if the count is zero. We will deal with the lack of FINISH in the next line indicating also the chromosome name for which we found an error
    if check_finish == "1":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! THE SCRIPT HAS NOT FINISHED FOR COMBINATION " + combination)

    print_text("check we have do NOT have the warning that should be only present in the dummy example", header=4)
    check_warning_to_avoid = run_bash(" \
        grep \
            'THIS WARNING SHOULD BE ONLY IN DUMMY DATA NOT IN 1KGP DATA as AC field should have only 1 value per line in the data' \
            --count \
            ./scripts/01_hap_map_calcs_outputs/" + comb_pop + "/chr" + comb_chrom + "_" + comb_pop + ".out || \
        [[ $? == 1 ]]", return_value=True).strip()
        #check we have the output indicating finish with a specific way, in uppercase and separated by "#", so we avoid the flag "--ignore-case". Ask for the count.
        #as in the previous case, we add a line to avoid errors if the count is zero. We will deal with the lack of FINISH in the next line indicating also the chromosome name for which we found an error
    if check_warning_to_avoid == "0":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! THE SCRIPT FOR COMBINATION " + combination + "HAS A WARNING THAT SHOULD NOT APPEAR HERE")

    print_text("extract the number of SNPs lost due to the filtering", header=4)
    row_results = run_bash(" \
        grep \
            'During filtering, we have lost' \
            ./scripts/01_hap_map_calcs_outputs/" + comb_pop + "/chr" + comb_chrom + "_" + comb_pop + ".out", return_value=True).strip()
        #look for the row including these results

    print_text("split the row, extract the numbers and calculate the percentage of lost snps and the number of remaining snps", header=4)
    row_results_split = row_results.split(" ")
    percent_lost=int(row_results_split[5])/int(row_results_split[8])*100
    if(percent_lost>=90):
        raise ValueError("ERROR! FALSE! THE PERCENTAGE OF SNPS LOST IS 90 OR HIGHER FOR COMBINATION " + combination)
    remaining_snps=int(row_results_split[8])-int(row_results_split[5])

    print_text("check we have the correct number of SNPs and samples in the hap and the map files", header=4)
    print("calculate the number of rows (SNPs) and columns (samples*2) in the hap file")
    hap_shape=run_bash("\
        awk \
            'BEGIN{ \
                FS=OFS=\" \" \
            }END{ \
                print NR,NF\
            }' \
            <( \
                gunzip \
                    --stdout \
                    ./results/03_hap_map_files/" + comb_pop + "/chr" + comb_chrom + "/chr" + comb_chrom + "_" + comb_pop + "_IMPUTE2.hap.gz \
            )", return_value=True).strip()
    hap_n_rows=int(hap_shape.split(" ")[0])
    hap_n_columns=int(hap_shape.split(" ")[1])
    print("calculate the number of rows (SNPs) and columns (chrom, id, gen pos, physical pos) in the map file")
    map_shape=run_bash("\
        awk \
            'BEGIN{ \
                FS=OFS=\" \" \
            }END{ \
                print NR,NF\
            }' \
            <( \
                gunzip \
                    --stdout \
                    ./results/03_hap_map_files/" + comb_pop + "/chr" + comb_chrom + "/chr" + comb_chrom + "_" + comb_pop + "_selscan.map.gz \
            )", return_value=True).strip()
    map_n_rows=int(map_shape.split(" ")[0])
    map_n_columns=int(map_shape.split(" ")[1])
    print("the number of columns of hap should be the number of samples of the pop times 2. It should be 4 in am")
    if((unrelated_samples.loc[unrelated_samples["pop"]==comb_pop,:].shape[0]!=hap_n_columns/2) | (map_n_columns!=4)):
        raise ValueError("ERROR! FALSE! THE NUMBER OF COLUMNS OF THE HAP FILE IS NOT THE NUMBER OF SAMPLES TIMES 2 OR 4 IN MAP FOR COMBINATION " + combination)
    if((hap_n_rows!=remaining_snps) | (map_n_rows != remaining_snps)):
        raise ValueError("ERROR! FALSE! THE NUMBER OF ROWS OF THE HAP AND MAP FILES IS NOT THE SAME THAN THE NUMBRE OF REMAINING SNPS CALCULATED AT THE END FOR COMBINATION " + combination)

    print_text("save", header=4)
    snps_lost_percentage.append(percent_lost)
    snps_remaining_list.append(remaining_snps)


print_text("see the percentiles of percentage of SNPs lost across all chromosomes and populations", header=3)
print("Do we have calculated all combinations")
print(len(snps_lost_percentage)==len(full_combinations_pop_chroms))
print(len(snps_remaining_list)==len(full_combinations_pop_chroms))

print("calculate the percentiles of the number of snps lost across combinations")
import numpy as np
for i in [0.025, 0.1,0.25,0.4,0.5,0.6,0.75,0.9, 0.975]:
    print("Percentile " + str(i) + "%: " + str(np.quantile(snps_lost_percentage, i)))
print("max percentage of lost "+str(np.max(snps_lost_percentage)))
print("This is combination " + full_combinations_pop_chroms[np.where(snps_lost_percentage==np.max(snps_lost_percentage))[0][0]])

print("calculate the percentiles of the number of snps left across combinations")
for i in [0.025, 0.1,0.25,0.4,0.5,0.6,0.75,0.9, 0.975]:
    print("Percentile " + str(i) + "%: " + str(np.quantile(snps_remaining_list, i)))
print("min number of snps left "+str(np.min(snps_remaining_list)))
print("This is combination " + full_combinations_pop_chroms[np.where(snps_remaining_list==np.min(snps_remaining_list))[0][0]])


print_text("check we do NOT have any errors in the global output of each pop, also check that we reached the end of the script", header=2)
print_text("run loop across pops combinations", header=3)
#selected_pop=pop_names[0]
for selected_pop in pop_names:
    print_text("Doing pop " + selected_pop, header=4)

    print_text("count number of cases with 'error' or 'false' in the output file", header=4)
    count_error_false = run_bash(" \
        grep \
            'error|false' \
            --extended-regexp \
            --ignore-case \
            --count \
            ./scripts/01_hap_map_calcs_outputs/global_outputs/" + selected_pop + ".out || \
        [[ $? == 1 ]]", return_value=True).strip()
        #using grep, look for 
            #"error" OR "false" using "|". In order to avoid escaping the symbol, i.e., "\|", we need to use the flag "--extended-regexp"
                #"error" include any string combing after like "errorS", "errores", etc... If the string "error" is present alone or in combination with other strings, you will get a hit
            #ignore the case, so "Error" and "ERROR" are also included
            #get the count, not the rows matching
        #if the exit code of grep indicates error run the code after "||". This is the function of "||", run the code at the right only if the code at the left failed
            #check if the exit code ("$?") is 1, if so, this will give 0 as exist status (i.e., no error and True) and give the count, which is zero, as stdout. As explained below, if the count is zero (stdout=0), the exit code is 1 in grep:
                #0: no error and one or more lines were selected
                #1: no error but no lines were selected
                #>1: an error occurred
            #This is different from other programs where exit code equals to 1 is error, and we coded that accordingly in run_bash.
            #Because of this, in this particular case, we add an additional line in case grep gives non-zero exist status, and avoid error if the exit status is 1.
                #If grep gives "1" as exit status because the string is not present in the file, we run [[ $? == 1 ]]. This will give "0" as exist status if the previous exit status was "1", while maintaining the previous stdout, i.e., the "count=0" because grep did not find the string in the file.
                #If the exit status is >1 and thus, there is an error, this will give "1" as exist status and run_bash will fail, so we are not hiding errors. 
                #If the exist status is 0, "||" avoids running the conditional (I have checked looking for "ancestral"), so we are good.
            #https://unix.stackexchange.com/a/427598
            #https://pubs.opengroup.org/onlinepubs/9699919799/utilities/grep.html#tag_20_55_14
        #https://linuxize.com/post/grep-multiple-patterns/

    print_text("check the count of problematic cases is zero", header=4)
    if count_error_false == "0":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE ERROR/FALSE IN THE GLOBAL OUTPUT FILE OF POP " + selected_pop)

    print_text("check we have the row of Next steps", header=4)
    check_finish = run_bash(" \
        grep \
            'Next steps' \
            --count \
            ./scripts/01_hap_map_calcs_outputs/global_outputs/" + selected_pop + ".out || \
        [[ $? == 1 ]]", return_value=True).strip()
        #check we have the output indicating finish with a specific way, in uppercase and separated by "#", so we avoid the flag "--ignore-case". Ask for the count.
        #as in the previous case, we add a line to avoid errors if the count is zero. We will deal with the lack of FINISH in the next line indicating also the chromosome name for which we found an error
    if check_finish == "1":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! THE SCRIPT HAS NOT FINISHED FOR POP " + selected_pop)



print_text("FINISH", header=1)
#chmod +x ./scripts/01f_hap_map_calcs_check_outputs.py; ./scripts/01f_hap_map_calcs_check_outputs.py > ./scripts/01f_hap_map_calcs_check_outputs.out 2>&1
