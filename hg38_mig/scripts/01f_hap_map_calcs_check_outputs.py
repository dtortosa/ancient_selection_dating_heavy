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
print_text("define function", header=3)
#combination=full_combinations_pop_chroms[0]
def heavy_checks(combination):
    print_text("Doing combination " + combination, header=4)

    print_text("split the combination name", header=4)
    comb_pop = combination.split("_")[0]
    comb_chrom = combination.split("_")[1]

    print_text("count number of problematic cases in the output file", header=4)
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

    print_text("check we have the correct number of SNPs and samples in the hap, map and vcf files. Also check we have selected the correct ID samples in the VCF file", header=4)
    print("extract the IDs of the samples included in the VCF file and the rows without counting the header")
    id_samples_n_row_raw=run_bash("\
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\" \
            }{ \
                if($0 ~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/){ \
                    print $0; \
                    last_row_header=NR \
                } \
            }END{ \
                print NR-last_row_header \
            }' \
            <( \
                gunzip \
                    --stdout \
                    ./results/01_cleaned_vep_vcf_files/" + comb_pop + "/chr" + comb_chrom + "/1kGP_high_coverage_Illumina.chr" + comb_chrom + ".filtered.SNV_phased_panel.vep.anc_up." + comb_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz \
            )", return_value=True).strip()
        #we are just getting the last row of the header with the sample IDs and getting the number of rows of the VCF file after subtracting the number of rows of the header, i.e., getting the number of SNPs.
    vcf_n_rows=int(id_samples_n_row_raw.split("\n")[1])
    id_samples_raw=id_samples_n_row_raw.split("\n")[0]
        #the row with IDs and the number of rows are separated by "\n".
    id_samples=id_samples_raw.split("\t")[9:]
        #within the last row of the header, i.e., the one with sample IDs, the columns are separated with tabs
        #We are assuming here that sample ID start at position 9 (counting from zero). If, for any reason, you analyze a different VCF file format, then we will not be selecting all sample IDs and the check will fail, so no problem. We would come here to check what is going on.
    print("check the IDs obtained from the last version of the VCF file are the same than those of the pop file")
    if(id_samples!=unrelated_samples.loc[unrelated_samples["pop"]==comb_pop, "sample"].tolist()):
        raise ValueError("ERROR! FALSE! WE HAVE NOT SELECTED THE CORRECT SAMPLE IDS FOR COMBINATION " + combination)

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
    print("the number of rows of the hap, map and VCF files should be the same than the one indicated as the remaining number of snps")
    if((hap_n_rows!=remaining_snps) | (map_n_rows != remaining_snps) | (vcf_n_rows != remaining_snps)):
        raise ValueError("ERROR! FALSE! THE NUMBER OF ROWS OF THE HAP AND MAP FILES IS NOT THE SAME THAN THE NUMBRE OF REMAINING SNPS CALCULATED AT THE END FOR COMBINATION " + combination)

    print_text("return the combination, the number of samples (i.e., half hap columns), remaining snps and the percentage of snps lost", header=4)
    return([combination, hap_n_columns/2, remaining_snps, percent_lost])


print_text("run function", header=3)
#heavy_checks(full_combinations_pop_chroms[0])
import multiprocessing as mp
pool = mp.Pool(20)
results_map = pool.map( \
    func=heavy_checks, \
    iterable=full_combinations_pop_chroms)
        #Apply `func` to each element in `iterable`, collecting the results in a list that is returned
        #if you need more arguments, you can use "partial" to fix the value of other argument and only iterate across the combination
pool.close()
print_text("see first 10 gene ids", header=4)
print(results_map[0:10])

print_text("convert to DF, set column names and check we do have all the combinations", header=4)
print("create final DF")
results_df=pd.DataFrame( \
    results_map, \
    columns=["combination", "n_samples", "remaining_snps", "percent_lost"])
print(results_df)
print("create columns for pop and chrom separately")
#x=results_df.iloc[0,:]
pop_column=results_df.apply(lambda x: x["combination"].split("_")[0], axis=1)
chrom_column=results_df.apply(lambda x: x["combination"].split("_")[1], axis=1)
    #in each row, select the combination and split it, selecting the pop and the chrom, respectively
results_df.insert(1, "pop", pop_column)
results_df.insert(2, "chrom", chrom_column)
    #insert both series in a specific location
        #https://stackoverflow.com/a/18674915
if(not results_df["combination"].equals(results_df["pop"]+"_"+results_df["chrom"])):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM GENERATING THE POP-CHROM COLUMNS")
if(results_df.columns.tolist()!=["combination", "pop", "chrom", "n_samples", "remaining_snps", "percent_lost"]):
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM GENERATING THE POP-CHROM COLUMNS")
print("check no NA and all combinations have been analyzed")
if(not results_df.dropna().equals(results_df)):
    raise ValueError("ERROR! FALSE! WE HAVE NOT ANALYZED ALL COMBINATIONS")
if(results_df["combination"].unique().tolist() != full_combinations_pop_chroms):
    raise ValueError("ERROR! FALSE! WE HAVE NOT ANALYZED ALL COMBINATIONS")
if(results_df.shape[0]!=26*22):
    raise ValueError("ERROR! FALSE! WE HAVE NOT ANALYZED ALL COMBINATIONS")
print("check that we have the same number of sample for the combinations of the same pop")
#x=results_df.loc[results_df["pop"]=="GBR",:]
if(sum(results_df.groupby(["pop"]).apply(lambda x: len(x["n_samples"].unique())==1))!= len(results_df["pop"].unique())):
    raise ValueError("ERROR! FALSE! THE CHROMS OF THE SAME POP DO NOT HAVE THE SAME NUMBER OF SAMPLES")
    #for each group of rows of the same pop, extract the number of samples, select unique cases and count them, check it is 1. the total number of Trues should be the same than unique number of pops
        #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.groupby.html
print("check we have 22 chromosomes in each population")
#x=results_df.loc[results_df["pop"]=="GBR",:]
if(sum(results_df.groupby(["pop"]).apply(lambda x: x["chrom"].values.tolist()==[str(item) for item in range(1, 23, 1)]))!= len(results_df["pop"].unique())):
    raise ValueError("ERROR! FALSE! WE DO NOT HAVE FROM 1 TO 22 CHROMOSOMES IN ALL POPS")

print_text("calculate percentiles per pop", header=4)
#x=results_df.loc[results_df["pop"]=="GBR",:]
results_df_quantiles=results_df.groupby(["pop"])[["percent_lost", "remaining_snps"]].quantile([0.025, 0.5, 0.975])
    #you get two columns with the percentiles of each pop one after other
results_df_quantiles.index.names=["pop", "quantiles"]
    #add the name of the multiindex
        #one index is population
        #the second column of the index is the quantile
    #these two indexes will be columns when writing
results_df_quantiles=results_df_quantiles.reset_index()
    #convert the index as columns to do operatiosn later on
        #https://stackoverflow.com/a/20461206
if(results_df_quantiles.shape[0]!=26*3):
    raise ValueError("ERROR! FALSE! PROBLEM IN PERCENTILE FILE. we should have 3 QUANTILES FOR EACH OF THE 26 POPS")

print_text("plot the number of snps lost and retained across pops and chromosomes", header=4)
print("create an array with the colors for each chromosome")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.lines import Line2D
#import seaborn as sns
#palette = sns.color_palette(None, len(results_df["chrom"].unique()))
    #https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
palette = cm.rainbow(np.linspace(0, 1, len(results_df["chrom"].unique())))
    #Return evenly spaced numbers over a specified interval
        #we want 22 numbers (one per chromosome) between 0 and 1
    #then create a colormap object based on lookup tables using these linear segments.
        #The lookup table is generated using linear interpolation for each primary color, with the 0-1 domain divided into any number of segments
    #the result is 22 arrays with 4 numbers in each one, I guess indicating the color 
print("get the pop names sorted by the median number of remaining snps")
pops_sort=results_df_quantiles[results_df_quantiles["quantiles"]==0.5].sort_values("remaining_snps", ascending=False)["pop"]
    #first the pops with more snps
print("plot chroms of each pop as a loop")
fig, axs = plt.subplots(2)
    #https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
#plot the data of each population (i.e., pandas series) separately, first pop 1, then pop 2, so we can set a different shape of each population and also apply the same color pattern for chromosomes across all of them
    #https://stackoverflow.com/a/12236808
#pop_index=1; pop=pops_sort.iloc[1,] 
for pop_index, pop in enumerate(pops_sort):
    
    #select the data (raw, not median) of the selected pop
    subset_pop=results_df.loc[results_df["pop"]==pop]
    
    #caculate the position in X axis for the datapoints (i.e., value per chromosome) for the selected pop
    xtics_pos=np.arange(1,23,1)+(22*pop_index)
        #for the first pop, pop_index=0, thus we sum 0 leading to 1...22
        #for the second, pop_index=1, thus we sum 22, leading to 23, 24...
        #and so on...

    #get a list the possible markers
    list_markers=list(Line2D.markers.keys())[:-4]
        #we use a dict with the symbols as keys
        #avoid the last 4 cases which are None
            #https://stackoverflow.com/a/59648414
            #https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers

    #plot the remaining and lost snps in each position and disable xticks
    axs[0].scatter(x=xtics_pos, y=subset_pop["remaining_snps"], s=20, c=palette, marker=list_markers[pop_index])
        #https://stackoverflow.com/a/12236808
    axs[1].scatter(x=xtics_pos, y=subset_pop["percent_lost"], s=20, c=palette, marker=list_markers[pop_index])
        #s: marker size
        #c: marker colors
            #it can be an array of values to be colormapped
        #marker: shape to plot

    #disable xticks
    axs[0].set_xticks([])
    axs[1].set_xticks([])
        #https://www.geeksforgeeks.org/how-to-hide-axis-text-ticks-or-tick-labels-in-matplotlib/

    #set titles of X and Y axes for each plot
    axs[0].set_ylabel("Number of SNPs retained", fontsize=12.5)
    axs[1].set_ylabel("Percentage of SNPs lost", fontsize=12.5)
    axs[1].set_xlabel("Population", fontsize=12.5)

    #add the title of the whole plot
    fig.suptitle(t="Results of filtering", x=0.5, y=0.92, fontsize=17)
        #https://stackoverflow.com/a/7066293

#plot the population names sorted by the number fo remaining snps putting them in the middle of the points in each case
axs[1].set_xticks(ticks=[i for i in range(12, results_df.shape[0]+1, 22)], labels=pops_sort, fontsize=12.5)
    #we create a range staring at 19 every 22 chromosomes, to cover the chromosomes of each pop. Se we are jumping from position the middle of chromosomes of one population to the next one

#get a list of the chromosomes
list_chrom=[chrom for chrom in range(1, len(results_df["chrom"].unique())+1, 1)]

#create a list with the objects indicating the color of each chromosome in the legend
list_color_points=list()
#index_chrom=0; chrom=list_chrom[0]
for index_chrom, chrom in enumerate(list_chrom):

    #selected point color
    selected_point_color=Line2D([0], [0], marker="o", color="white", markerfacecolor=palette[index_chrom], label=chrom, markersize=5)
        #xdata, ydata
            #coordinates of the line
        #marker: shape, you can select a line or a dot
        #color: color of the line (a line is always plotted indepently of the shape, if you do not want it, use "white")
        #markerfacecolor: color of the shape
        #label to be shown
        #size of the marker shape
            #https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
            #https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D

    #append
    list_color_points.append(selected_point_color)

#plot the legend with each chromosome
axs[0].legend(labels=list_chrom, handles=list_color_points, bbox_to_anchor=(1.12, 1.05), fontsize=15, markerscale=1.5, title="Chromosome", title_fontsize=12.5)
    #we use handle to explicitly listing the artists in the legend
        #https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
    #bbox_to_anchor for having the legend outside the plot
        #https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
    #fontsize is the size of the labels
    #markerscale is the size of the marker with respect of the size of the marker in the plot

#set size of the plot
plt.gcf().set_size_inches(w=18, h=9)
    #width and height in inches
        #https://stackoverflow.com/a/72199494

#save
plt.savefig(fname="./results/03_hap_map_files/snps_lost_retain.png", dpi=200)
plt.close()

print_text("save the tables", header=4)
results_df.to_csv( \
    "./results/03_hap_map_files/snps_lost_retain.tsv", \
    sep="\t", \
    header=True, \
    index=False)
results_df_quantiles.to_csv( \
    "./results/03_hap_map_files/snps_lost_retain_quantiles.tsv", \
    sep="\t", \
    header=True, \
    index=False)


print_text("see the percentiles of percentage of SNPs lost across all chromosomes and populations", header=3)
print_text("calculate the percentiles of the number of snps lost across combinations", header=4)
import numpy as np
for i in [0.025, 0.1,0.25,0.4,0.5,0.6,0.75,0.9, 0.975]:
    print("Percentile " + str(i) + "%: " + str(results_df["percent_lost"].quantile(i)))
print("max percentage of lost "+str(results_df["percent_lost"].max()))
print("This is combination " + results_df.loc[results_df["percent_lost"]==results_df["percent_lost"].max(), "combination"].values)

print_text("calculate the percentiles of the number of snps left across combinations", header=4)
for i in [0.025, 0.1,0.25,0.4,0.5,0.6,0.75,0.9, 0.975]:
    print("Percentile " + str(i) + "%: " + str(results_df["remaining_snps"].quantile(i)))
print("min number of snps left "+str(results_df["remaining_snps"].min()))
print("This is combination " + results_df.loc[results_df["remaining_snps"]==results_df["remaining_snps"].min(), "combination"].values)



print_text("check we do NOT have any errors in the global output of each pop, also check that we reached the end of the script", header=2)
print_text("run loop across pops combinations", header=3)
#selected_pop=pop_names[0]
for selected_pop in pop_names:
    print_text("Doing pop " + selected_pop, header=4)

    print_text("count number of problematic cases in the output file", header=4)
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
