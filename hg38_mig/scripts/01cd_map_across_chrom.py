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



######################################################
######## CALCULATE HAP AND MAP FILES FOR HG38 ########
######################################################

#We are going to migrate to hg38 as the new release of 1KGP matches this genome reference version (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38). This script is going to take the VCF files for SNPs in the hg38 high coverage data from 1000GP and obtain hap and map files that will be used as input by flex sweep.

#We will use as input the VCF files that include the ancestral allele as a new field, which was obtained with VEP.

#this script will calculate only the map files across chromosomes. In a next step, we will calculate map and hap files per pop. In this way we avoid to calculate genetic positions of SNPs several times, just one per SNP.

#For more information about the download of the pedigrees and the VCF files with the original 1KGP data, see 01a_selecting_pedegree.py and 01b_vep_ancestral.py, respectively.



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



###############
# folder prep #
###############
print_text("folder prep", header=1)

#save name of path including input vcf files
input_vcfs_path = "results/00_vep_vcf_files"

#create folders to save the results
run_bash(" \
    mkdir \
        -p ./scripts/00b_map_calcs_outputs; \
    mkdir \
        -p ./results/00b_map_files")
    #-p: no error if the folder already exists, make parent directories as needed



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



#################
# bcftools prep #
#################
print_text("Prepare bcftool", header=1)

#see version of bcftools
print("\n#######################################\n#######################################")
print("see bcftools version")
print("#######################################\n#######################################")
run_bash("bcftools -v")
    #bcftools cheatsheet
        #https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b#file-bcftools-cheat-sheet



##############################################################
#### function to calculate the genetic position of a SNP ####
##############################################################
print_text("function to calculate the genetic position of a SNP", header=1)
print_text("chr define the function", header=2)
#TO DEBUG, look for snp_map_raw inside master_processor
#selected_snp_id=snp_map_raw.iloc[31967]["id"] #snp with cM value just in its position for chromosome 1
#selected_snp_id=snp_map_raw.iloc[33877]["id"] #snp with cM values at both sides for chromosome 1
#selected_snp_id=snp_map_raw.iloc[0]["id"] #snp without cM values around for chromosome 1
import numpy as np
def gen_pos(selected_snp_id):

    #extract the number of the selected chromosome
    selected_chromosome=selected_snp_id.split(":")[0].split("chr")[1]
    
    #select the global variables of the SNP and deCODE maps for the selected chromosome
    snp_map_raw = eval("snp_map_raw_chr" + selected_chromosome)
    decode2019_map_subset = eval("decode2019_map_subset_chr" + selected_chromosome)
        #these maps are created in the parent function (see below) for the selected chromosome and then converted to global variables that can be accessed by any other function
        #we save them with the original name of these maps within this function, but given we are inside the child function and we are not converting them to global variables, then this would not affect to the original variables with the same name. See dummy example:
            #a=1
            #print("See the global 'a' variable just created: " + str(a))
            #def eso():
            #    a=eval("a")
            #    a=a+1
            #    return(a)
            #print("See 'a' modified within the function, i.e., local: a=" + str(eso()))
            #print("See again 'a' outside of the function, i.e., global: a=" + str(a))
            #the global variable "a" is NOT modified by the function eso()

    #extract the row in the raw map for the selected SNP
    selected_snp_row = snp_map_raw.loc[snp_map_raw["id"] == selected_snp_id,:]

    #check we have the correct chromosome
    check_0 = selected_snp_row["chr"].unique()[0] == "chr"+str(selected_chromosome)

    #extract position of the selected snp
    selected_snp_physical_pos = selected_snp_row["pos"].iloc[0]

    #extract old ID (this follows VFP file format so we can use them to filter it)
    selected_snp_old_id = selected_snp_row["id_old"].iloc[0]

    #select those deCODE intervals that are at least 1 MB close to the selected SNP: I have search for genetic position data around each SNP of each population, 1MB at each side. I think remember you told me that if we do not find data points at both side of the SNP, we can safely remove it. In that way, we include areas with low recombination (possible haplotypes) but not areas with a lot of missing data.
    decode2019_map_subset_around_snp = decode2019_map_subset.loc[\
        (decode2019_map_subset["end"] >= (selected_snp_physical_pos - 1000000)) & \
        (decode2019_map_subset["end"] <= (selected_snp_physical_pos + 1000000)), :]
            #we are only interested in the END coordinate because the cM data of each interval came from the end of the interval. Indeed, the start coordinate of the next interval is the same of the end of the previous one. Therefore, we focus on end coordinate.
            #Note that we are using 1000000 directly. If a SNP is at 1000001, then 1000001-1000000=1, length(1:1000001) is equal to 1000001, which is not exactly 1MB, but this is only 1 base of difference. This is not important.
            #Also note that for SNPs between base 0 and 1000kb, the difference between SNP position and 1000kb will be negative, but this is OK:
                #If a SNP is before base 1000kb, then there is less than 1000kb bases to look for decode intervals before the SNP, reaching base 0.
                #therefore, having a negative value would mean the same than just look up to zero (there are not decode intervals below zero).
                    #SNP at position 500kb
                        #500kb+1000kb=1500kb
                        #500kb-1000kb=-500kb
                        #there is no enough space at the left of the SNP to look for decode intervals up to 1000kb, so we have to reach 0, which is the same than looking for values equal or higher than a negative given that no decode interval has a negative position.
                #this will be the case until a SNP in position 1001kb, as 1001kb-1000kb would be 1, so we do not look for decode intervals below 1.
                #As we move foward from 1000kb, the lower limit of the window starts moving away from 1.
                #indeed, using the absolute value would not work
                    #SNP at position 500kb
                        #1000kb+500kb=1500kb
                        #1000kb-500kb=500kb
                        #the lower limit cannot be 500kb, when the SNP is at 500kb. We would automatically lose this SNP.
                        #you would need -500 to 1500kb.
                #one concern about this is that for some SNPs we are looking for decode intervals in a smaller region, but SNPs below 1MB are not frequent. 
                    #For example, in chromosome 1, only 1339 out of 800K are at a coordinate below 1000kb. Therefore, this does not seem to be a problem. I have not tested it, but I guess the same would go for SNPs close to the end of the chromosome, this would be a small proportion of the total number of SNPs.
                    #More important, there are NO decode interval below base 500kb, so we will discard any SNP before that base because no cM value will be available to the left in order to interpolate. Therefore, the importance of this issue is very limited.       

    #select those intervals before and after the selected SNP
    intervals_lower_end = decode2019_map_subset_around_snp.loc[(decode2019_map_subset_around_snp["end"] < selected_snp_physical_pos), :]
    intervals_upper_end = decode2019_map_subset_around_snp.loc[(decode2019_map_subset_around_snp["end"] > selected_snp_physical_pos), :]

    #select those decode intervals with the same position than the selected SNP
    interval_same_pos = decode2019_map_subset_around_snp.loc[decode2019_map_subset_around_snp["end"] == selected_snp_physical_pos, :]

    #if we have deCODE intervals 1MB around the selected SNP, i.e., we have intervals at both sides, intervals ending before and after the selected SNP OR we have a deCODE interval ending exactly at the SNP. In the second case if you have cM value in the exact position of the selected SNP, then you do not need intervals at both sides.
    if (intervals_lower_end.shape[0]>0) & (intervals_upper_end.shape[0]>0) | (interval_same_pos.shape[0]>0): 
        #the first condition do not need equal because in the next condition (after "|") we consider the option of equal coordinate between window extreme and deCODE end interval.

        #checks
        check_1= \
            ( \
                (decode2019_map_subset_around_snp["end"] >= (selected_snp_physical_pos - 10**6)) & \
                (decode2019_map_subset_around_snp["end"] <= (selected_snp_physical_pos + 10**6)) \
            ).sum() == decode2019_map_subset_around_snp.shape[0]
        
        #if we dot NOT have an interval with an end coordinate exactly similar to the selected SNP
        if (interval_same_pos.shape[0] == 0):

            #check
            check_2 = (intervals_lower_end["end"] < selected_snp_physical_pos).sum() == intervals_lower_end.shape[0]
            check_3 = (intervals_upper_end["end"] > selected_snp_physical_pos).sum() == intervals_upper_end.shape[0]

            #from the intervals below the extreme window, select the biggest and hence closest to the extreme window   
            lowest_interval = intervals_lower_end.loc[intervals_lower_end["end"] == max(intervals_lower_end["end"]),:]
                #we cannot have two cases with the same value because the intervals are not overlapped and they are also in increasing order, the coordinate of an interval is bigger than the previous one.
                    #This was checked in decode_conversion_hg19_v3.R

            #from the intervals above the extreme window, select the smallest and hence closest to the extreme window
            highest_interval = intervals_upper_end.loc[intervals_upper_end["end"] == min(intervals_upper_end["end"]),:] 
                #we cannot have two cases with the same value because the intervals are not overlapped and they are also in increasing order, the coordinate of an interval is bigger than the previous one.
                    #This was checked in decode_conversion_hg19_v3.R

            #check that the end coordinate with lowest difference respect the SNP is the selected in the previous step both for the lower and higher intervals
            check_4a = intervals_lower_end.loc[ \
                np.abs(intervals_lower_end["end"]-selected_snp_physical_pos) == \
                np.min(np.abs(intervals_lower_end["end"]-selected_snp_physical_pos)), \
                "end"].to_list() == lowest_interval["end"].to_list()
            check_4b = intervals_upper_end.loc[\
                np.abs(intervals_upper_end["end"]-selected_snp_physical_pos) == \
                np.min(np.abs(intervals_upper_end["end"]-selected_snp_physical_pos)), \
                "end"].to_list() == highest_interval["end"].to_list()

            #stop if we have more than 1 closest interval in each side 
            if (lowest_interval.shape[0] > 1) | (highest_interval.shape[0] > 1):
                raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: we have more than 1 closest interval in each side")


            ##calculate cM value of the snp
            #extract the centimorgan of each the closest deCODE intervals to the SNP
            left_cM = lowest_interval["cM"].to_numpy()[0]
            right_cM = highest_interval["cM"].to_numpy()[0]

            #calculate the distance from each interval to the SNP
            distance_left_end = (selected_snp_physical_pos - lowest_interval["end"]).to_numpy()[0]
            distance_right_end = (highest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]
                #We do not need to include both extremes, we want the distance from one point to another. Imagine the window begins at 1 and ends at 3. Including both extremes, the size of the window is 3, you have 3 bases. However, the distance from the point 1 to 3 is 2 (3-1=2). We want the distance between two points with centiMorgan values.

            #check that calculating the distance with abs and changing order gives the same result
            check_5a = distance_left_end == abs(lowest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]
            check_5b = distance_right_end == abs(highest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]

            #check that the sum of the physical distance of each deCODE end point to the SNP is the same than the total distance between the deCODE end points
            check_6 = distance_left_end + distance_right_end == np.abs(lowest_interval["end"].to_numpy()[0] - highest_interval["end"].to_numpy()[0])

            #calculate the genetic distance using the formula of David
            genetic_distance = left_cM + (right_cM - left_cM) * distance_left_end / (distance_left_end + distance_right_end)
                #Explanation of David: In that case, you have to find the genetic map position of a single SNP (or in general a position of the genome). To assign a position to each SNP, you can use the two genetic map positions (end points as you described) left and right form the SNP. You can then consider that the genetic position increases linearly between the two left and right positions. For example, if a SNP is between a genetic position on the left at 100 cM, and a genetic position on the right at 102 cM, and the SNP is located 20 kb from the left genetic position but 80kb from the right position, then the SNP will be located at genetic position: 100 cM + (102 cM-100 cM) * 20kb / (20kb+80kb) = 100.4 cM. Of course, if the SNP is right on the coordinate of an end point, then just use the genetic map position directly for that SNP.
                #My explanation: What David is doing is 100 + ((102-100)*20)/(20+80). This gives exactly 100.4. David is using the rule of three (https://en.wikipedia.org/wiki/Cross-multiplication#Rule_of_Three). You have three points, A, B and C. If the physical distance distance A-C is 100 kb (20+80) and the genetic distance between these points is 2 cM (102-100) , what would be the genetic distance between A-B if these points are separated by 20 kb? ((102-100 cM) * 20 kb) / (20+80 kb); ((2 cM) * 20 kb) / (100 kb); ((2 cM) * 20 kb) / (100 kb); (40 cM * kb) / 100 kb; 0.4 cM. 0.4 is the genetic distance between A and B. Now we can sum 0.4 to the genetic position of A, to get the genetic position of B in the genome. 100 cM + 0.4 cM = 100.4 cM.
                #If the point for which you calculate the genetic distance is exactly in the middle of the two points with 100 and 102 cM of genetic distance, the resulting genetic distance would be exactly in the middle, i.e., 101: 100 + ((102-100)*50)/(50+50). 50 is the physical distance between cM point and the point of interest. 
                #This method assumes that relationship between genetic distance and physical distance between two points is lineal and stable, so you can estimate the genetic distance based on the physical distance in the genomic region encompassed by these points. Note that you are using points that are at least 1MB close to the point under study, therefore, we are estimating the genetic distance using the relationship between physical and genetic distance in a specific genomic region, not the whole genome.
                #see figure 31 for further details.
        else:

            #if not and hence we have an deCODE interval exactly in the SNP position

            #stop if we have more than 1 interval with the exact position of the selected SNP
            if (interval_same_pos.shape[0] > 1):
                raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM HERE: We have more than two deCODE intervals with the exact physical position of the selected SNP.")

            #save the genetic position
            genetic_distance = interval_same_pos["cM"].to_numpy()[0]

            #set NA for the rest of results. They are not needed.
            check_2 = np.nan
            check_3 = np.nan                
            check_4a = np.nan
            check_4b = np.nan
            check_5a = np.nan
            check_5b = np.nan
            check_6 = np.nan
            left_cM = np.nan
            right_cM = np.nan
            distance_left_end = np.nan
            distance_right_end = np.nan
    else:

        #if not, and hence we cannot calculate the cM of SNP
        genetic_distance = np.nan

        #set NA for the rest of results. They are not needed.
        check_1 = np.nan
        check_2 = np.nan
        check_3 = np.nan
        check_4a = np.nan
        check_4b = np.nan
        check_5a = np.nan
        check_5b = np.nan
        check_6 = np.nan
        left_cM = np.nan
        right_cM = np.nan
        distance_left_end = np.nan
        distance_right_end = np.nan

    #save results
    return(tuple([selected_chromosome, selected_snp_id, selected_snp_old_id, selected_snp_physical_pos, check_0, check_1, check_2, check_3, check_4a, check_4b, check_5a, check_5b, check_6, genetic_distance, left_cM, right_cM, distance_left_end, distance_right_end]))





######################################################
#### function to calculate map files for ALL SNPs ####
######################################################
print_text("function to calculate map files for ALL SNPs", header=1)
#selected_chromosome="1"; debugging=True; debug_file_size=50000
def master_processor(selected_chromosome, debugging=False, debug_file_size=None):

    #redirect standard output ONLY when running for production
    if debugging == False:
        import sys
        original_stdout = sys.stdout
            #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
        file_output = open("./scripts/00b_map_calcs_outputs/chr" + selected_chromosome + ".out", "w")
            #open a file where stdout will be saved
        sys.stdout = file_output
            #redirect stdout to that file
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/
            #https://stackoverflow.com/a/23838153



    print_text("Initial operations", header=2)
    print_text("Create a new folder for the selected chromosome", header=3)
    run_bash(" \
        mkdir \
            -p \
            ./results/00b_map_files/chr" + selected_chromosome + "/")

    print_text("chr" + selected_chromosome + ": see VCF file version", header=3)
    vcf_version = run_bash(" \
        bcftools head \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        grep -i '^##fileformat'", return_value=True).strip()
            #use bcftools to see the header and then select the row starting with ##fileformat to see the version.
            #https://www.htslib.org/howtos/headers.html
    print("This script assumes VCFv4.2 (see script for explanations about this format). Do we have that VCF version?")
    if vcf_version == "##fileformat=VCFv4.2":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WITH THE VCF VERSION")

    #The Variant Call Format Specification v4.2 is the one used in the 1000GP high coverage
        #https://samtools.github.io/hts-specs/VCFv4.2.pdf

    #Fixed fields of this format v4.2
        #CHROM - chromosome: An identifier from the reference genome or an angle-bracketed ID String (“<ID>”) pointing to a contig in the assembly file (cf. the ##assembly line in the header). All entries for a specific CHROM should form a contiguous block within the VCF file. (String, no whitespace permitted, Required).
        #POS - position: The reference position, WITH THE 1ST BASE HAVING POSITION 1. Positions are sorted numerically, in increasing order, within each reference sequence CHROM. It is permitted to have multiple records with the same POS. Telomeres are indicated by using positions 0 or N+1, where N is the length of the corresponding chromosome or contig. (Integer, Required).
            #THIS IS 1 BASED.
        #ID - identifier: Semicolon-separated list of unique identifiers where available. If this is a dbSNP variant it is encouraged to use the rs number(s). NO IDENTIFIER SHOULD BE PRESENT IN MORE THAN ONE DATA RECORD. If there is no identifier available, then the missing value should be used. (String, no whitespace or semicolons permitted).
        #REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String. See Variant Call Format Specification v4.2 for further details.
        #ALT - alternate base(s): Comma separated list of alternate non-reference alleles. These alleles do not have to be called in any of the samples. See Variant Call Format Specification v4.2 for further details.
        # QUAL - quality: Phred-scaled quality score for the assertion made in ALT. ee Variant Call Format Specification v4.2 for further details.
        #FILTER - filter status: PASS if this position has passed all filters, i.e., a call is made at this position. See Variant Call Format Specification v4.2 for further details.
            #Our datasets was previously filter selecting only those variants with PASS, see readme.
                #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf
        #INFO - additional information: (String, no whitespace, semicolons, or equals-signs permitted; commas are permitted only as delimiters for lists of values) INFO fields are encoded as a semicolon-separated series of short keys with optional values in the format: <key>=<data>[,data]. If no keys are present, the missing value must be used. Arbitrary keys are permitted, although the following sub-fields are reserved (albeit optional):
            #NOTE ABOUT GENOTYPES USED TO CALCULATE INFO DATA:
                #First, note, that INFO is fixed, you cannot modify, so whatever you do ni the VCF like subseting individuals or filtering variants, INFO wal remain without out changes. 
                #For example, if you look the number of alleles after subseting samples, it is going to give you the number of alleles considering the whole sample, not just for the subset of samples. I have checked this.
                #Second, I understand that the called genotypes are those included in my VCF file before subseting or doing anything and primary data should be original data before the filtering that generated the VCF I have received from 1000 genomes project.
            #AA: ancestral allele
            #AC: allele count in genotypes, for each ALT allele, in the same order as listed.
            #AF: allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not from called genotypes.
            #AN: total number of alleles in called genotypes
            #BQ: RMS base quality at this position
            #CIGAR: cigar string describing how to align an alternate allele to the reference allele
            #DB: dbSNP membership
            #DP: combined depth across samples, e.g. DP=154
            #EN: end position of the variant described in this record (for use with symbolic alleles)
            #H2: membership in hapmap2
            #H3: membership in hapmap3
            #MQ: RMS mapping quality, e.g. MQ=52
            #MQ0: Number of MAPQ == 0 reads covering this record
            #NS: Number of samples with data
            #SB: strand bias at this position
            #SOMATIC: indicates that the record is a somatic mutation, for cancer genomics
            #VALIDATED: validated by follow-up experiment
            #1000G:  membership in 1000 Genomes
            #You can add more fields to INFO, that can be specific for your study if you want to add information per variant.
        #Genotype fields
            #If genotype information is present, then the same types of data must be present for all samples. First a FORMAT field is given specifying the data types and order (colon-separated alphanumeric String). This is followed by one field per sample, with the colon-separated data in this field corresponding to the types specified in the format. The first sub-field must always be the genotype (GT) if it is present. There are no required sub-fields. 1KGP only has GT, if other fields present, you would have 0|0.... al genotypes, the : number : number :...
                #GT : genotype, encoded as allele values separated by either of / or | (UNPHASED AND PHASED RESPECTIVELY). The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc. For haploid calls, e.g. on Y, male nonpseudoautosomal X, or mitochondria, only one allele value should be given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, ‘.’ should be specified for each missing allele in the GT field (for example ‘./.’ for a diploid genotype and ‘.’ for haploid genotype). The meanings of the separators are as follows (see the PS field below for more details on incorporating phasing information into the genotypes):
                    #/ : genotype unphased
                    #| : genotype phased


    print_text("chr " + selected_chromosome + ": see first 10 samples", header=3)
    run_bash(" \
        bcftools query \
            --list-samples \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -10")
        #https://samtools.github.io/bcftools/howtos/query.html


    print_text("chr " + selected_chromosome + ": the number of samples is equal to the number of samples considering unrelated individuals and trios-duos?", header=3)
    run_bash(" \
        n_samples=$( \
            bcftools query \
                --list-samples \
                " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
                wc -l); \
        if [[ $n_samples -eq " + str(samples_pedigree.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #list the samples and count them
        #if the number of samples is equal to the number of samples in the pedigree, perfect

    print_text("chr " + selected_chromosome + ": show the variant type, ID, chromosome, position, alleles and frequency for the first snps", header=3)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF\n' \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
            head -5")
        #select the format of the query indicating the columns you want to show per SNP.
            #you can include data from INFO
            #end with \n to have different lines per SNPs


    print_text("chr " + selected_chromosome + ": see genotypes of first samples", header=3)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -2")


    print_text("Select the input VCF for genetic map calculation. Use only a subset of the data if we are on debugging mode", header=3)
    if debugging==True:
        run_bash(" \
            bcftools view \
                " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
            head -n " + str(debug_file_size) + " > ./results/00b_map_files/chr" + selected_chromosome + "/00b_map_calc_debug_subset_chr" + selected_chromosome + ".vcf; \
            ls -l")
        input_vcf_file_map_calc="./results/00b_map_files/chr" + selected_chromosome + "/00b_map_calc_debug_subset_chr" + selected_chromosome + ".vcf"
    else:
        run_bash(" \
            cp \
                " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz \
                ./results/00b_map_files/chr" + selected_chromosome)
        input_vcf_file_map_calc="./results/00b_map_files/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz"
    print(input_vcf_file_map_calc)
        #--types
            #comma-separated list of variant types to select. Site is selected if any of the ALT alleles is of the type requested. Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags like INFO/INDEL or INFO/VT. Use --include to select based on INFO tags.
            #therefore, I understand that they look if REF/ALT has 1 base or several. If you have "AACCCC", you have an indel. You have "A", then you have a single nucleotide polymorphism.
        #https://samtools.github.io/bcftools/bcftools.html

    print_text("check we have only SNPs", header=4)
    print("count the number of rows of the VCF when we exclude SNPs")
    n_rows_no_snps = run_bash(" \
        bcftools view \
            --no-header \
            --drop-genotypes \
            --exclude-types snps \
            " + input_vcf_file_map_calc + " | \
        awk \
            'END{print NR}'", return_value=True).strip()
        #--exclude-types
            #comma-separated list of variant types to exclude. Site is excluded if any of the ALT alleles is of the type requested. Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags like INFO/INDEL or INFO/VT. Use --exclude to exclude based on INFO tags.
                #therefore, I understand that they look if REF/ALT has 1 base or several. If you have "AACCCC", you have an indel. You have "A", then you have a single nucleotide polymorphism.
        #https://samtools.github.io/bcftools/bcftools.html
    print("check that the number of rows with non-SNP variants is actually zero")
    if n_rows_no_snps=="0":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE STILL HAVE NON-SNPS AFTER FILTERING OUT INDELS, ETC...")


    print_text("extract the position (index) of several columns in the VCF file, so we can be sure we are selecting these columns in later steps", header=3)
    print_text("obtain the position of these columns using awk", header=4)
    indexes_chrom_pos = run_bash(" \
        bcftools view \
            --header \
            " + input_vcf_file_map_calc + " | \
        awk \
            'BEGIN{FS=\"\t\"}; \
            END{ \
                for(i=1;i<=NF;i++){ \
                    if($i == \"#CHROM\"){ \
                        chrom_index=i \
                    }; \
                    if($i == \"POS\"){ \
                        pos_index=i \
                    }; \
                    if($i == \"ID\"){ \
                        id_index=i \
                    }; \
                    if($i == \"REF\"){ \
                        ref_index=i \
                    }; \
                    if($i == \"ALT\"){ \
                        alt_index=i \
                    }; \
                    if($i == \"FILTER\"){ \
                        filter_index=i \
                    }; \
                    if($i == \"INFO\"){ \
                        info_index=i \
                    }; \
                    if($i == \"FORMAT\"){ \
                        format_index=i \
                    }; \
                    if($i ~/^HG/ || $i ~/^NA/){ \
                        n_samples++ \
                    } \
                }; \
                n_fields=NF; \
                printf \"n_fields=%s,n_samples=%s,chrom=%s,pos=%s,id=%s,ref=%s,alt=%s,filter=%s,info=%s,format=%s\", n_fields, n_samples, chrom_index, pos_index, id_index, ref_index, alt_index, filter_index, info_index, format_index \
            }'", return_value=True).strip()
        #get the header of the VCF file after extracting AA from CSQ, removing CSQ and select SNPs, just like we are going to do when we replace lower for upper case in the next line
        #open the header with AWK
            #when you reach the last line, which includes the headers
                #run loop across fields, i.e., headers
                    #if the header is CHROM or POS then save the index of the field as a new variable, chrom_index and pos_index, respectively.
                        #"i" is the index, like 1, 2, 3... because of this, when you do $i is like you are doing $1.
                        #if we save "i", we are saving the index, the number
                        #https://unix.stackexchange.com/a/616495
                    #also look for REF, ALT, FILTER because we will use these columns later
                    #if the header starts with HG or NA, add 1 to the count of n_samples, because this is a GT column for a given sample
                        #"~" let you use regular expression
                        #"/.../" is a regular expression to match text that meet condition
                        #"^" text that starts with...
                        #https://unix.stackexchange.com/a/72763
                #save the number of fields
            #then print the number of fields, the number of sample and the index of both CHROM and POS
                #https://www.gnu.org/software/gawk/manual/html_node/Printf-Examples.html
    print("extract the numbers")
    n_fields = indexes_chrom_pos.split(",")[0].replace("n_fields=", "")
    n_samples = indexes_chrom_pos.split(",")[1].replace("n_samples=", "")
    index_chrom = indexes_chrom_pos.split(",")[2].replace("chrom=", "")
    index_pos = indexes_chrom_pos.split(",")[3].replace("pos=", "")
    index_id = indexes_chrom_pos.split(",")[4].replace("id=", "")
    index_ref = indexes_chrom_pos.split(",")[5].replace("ref=", "")
    index_alt = indexes_chrom_pos.split(",")[6].replace("alt=", "")
    index_filter = indexes_chrom_pos.split(",")[7].replace("filter=", "")
    index_info = indexes_chrom_pos.split(",")[8].replace("info=", "")
    index_format = indexes_chrom_pos.split(",")[9].replace("format=", "")
    print("total number of fields: " + n_fields)
    print("number of samples: " + n_samples)
    print("index of column CHROM: " + index_chrom)
    print("index of column POS: " + index_pos)
    print("index of column ID: " + index_id)
    print("index of column REF: " + index_ref)
    print("index of column ALT: " + index_alt)
    print("index of column FILTER: " + index_filter)
    print("index of column INFO: " + index_info)
    print("index of column FORMAT: " + index_format)
    print("the total number of fields minus the number of samples should be 9. The number of fixed fields in VCF v4.2 is 8 and then FORMAT, which in our case only has GT, thus we should have 9 fields. Also, the index of CHROM, POS, REF, ALT and FILTER should be 1, 2, 4, 5 and 7, respectively")
    if (int(n_fields)-int(n_samples) == 9) & (index_chrom=="1") & (index_pos=="2") & (index_id=="3") & (index_ref=="4") & (index_alt=="5") & (index_filter=="7") & (index_info=="8") & (index_format=="9"):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WITH THE FIELDS OF THE VCF FILE BEFORE CONVERTING TO UPPER CASE ANCESTRAL ALLELES")




    ##########################################
    # calculate map file within selected pop #
    ##########################################
    print_text("calculate map file within selected pop", header=2)
        #we are going to calculate the map of each population directly using the SNPs in its vcf file FILTERED WITHIN population. Then, when we know for which of the SNPs we have genetic position, we can further filter the VCF file and convert to hap.")
        #an alternative would be just take all the SNPs in the raw VCF file and calculate their genetic position.
            #it should be ok regarding the ID of the SNPs because REF/ALT are split when using multiallelics, and these fields are NOT switched based on the frequency of the SNPs in specific subsets.
        #I am going for the first option just to be completely sure I am using the SNPs (and positions) of the selected population and because I have the script almost ready.
            #If it is too slow this option, think about the other one.

    #Instructions david
        #I understand that SNPs without genetic position are NO useful for any summary statistic, right? So I can safely remove these SNPs from the VCF and hap files right?
            #ok
        #format of ID
            #in the map file can I use the format "CHROM:POS_REF_ALT" for the ID? 
            #I think remember that the map files I originally got from you in the previous project (before decode2019 conversion) used as ID just the physical position. Not sure if there is any specific reason for doing that.
            #not asked, by irrelevant question
        #data format decode map
            #in the decode map, they say clearly that the data is aligned to hg38.
            #they do not specify if the coordinates are 1 or 0-based, but I understand these are 1-based from what they say: "Begin (start point position of interval in GRCh38 coordinates) and End (end point position of interval in GRCh38 coordinates)"
                #I asked to Bjarni Vilhjálmur Halldórsson
                    #I would like to confirm that the three genetic maps published as supplementary files (Data S1-S3, i.e., maternal, paternal and average maps) have 1-based coordinates. From the description, I understand this is the case, but I want to double check that with you to ensure these are not 0-based.
                #He answered
                    #This shouldn't be the case, but let us know if you think there is a problem.
                #I asked David to interpret what Halldórsson said, and he told me:
                    #About the deCode map, he is saying that the coordinates are 1-based, we can proceed.
            #1KGP data is also aligned to hg38 and is 1-based.
            #Therefore I can just use the position of the SNPs in 1KGP to calculate their genetic position in the decode map, right?
                #David said that there is not problem if the map and hap files are in the same format. 
                #that is the case because the map file is calculated with decode map, and hap file with 1KGP VCF file, and as I said, I have no reason to think that the decode map is not 1-based.



    print_text("create the raw genetic map", header=3)
    print_text("extract snps from the cleaned VCF file", header=4)
    run_bash(" \
        bcftools view \
            --no-header \
            ./" + input_vcf_file_map_calc + " | \
        awk \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\" \"; \
                index_chrom=" + index_chrom + "; \
                index_pos=" + index_pos + "; \
                index_id=" + index_id + "; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + " \
            }{ \
                print $index_chrom, $index_id, $index_pos, $index_ref, $index_alt \
            }' | \
        gzip --force > ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_raw.map.gz; \
        gunzip \
            --stdout \
            ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_raw.map.gz | head -20")
        #load the VCF file without header
        #select the columns of interest with AWK
            #we want chromosome, position, REF/ALT to compare with the positions in the raw hap file (see below)
            #output "space" separated to meet salescan requirements
        #compress
        #take a look at it

    #Note about the format of the positions
    #pos in VCF files v4.2 is 1-based according to the specification file (this is the format of 1KGP data). Therefore, we have here 1-based coordinates.
        #POS - position: The reference position, with the 1st base having position 1. Positions are sorted numerically, in increasing order, within each reference sequence CHROM. It is permitted to have multiple records with the same POS. Telomeres are indicated by using positions 0 or N+1, where N is the length of the corresponding chromosome or contig. (Integer, Required)
            #https://samtools.github.io/hts-specs/VCFv4.2.pdf

    #required format according to hapbin
        #The map files (--map) should be in the same format as used by Selscan with one row per variant and four space-separated columns specifiying 
            #chromosome, 
            #locus ID, 
            #genetic position
            #physical position.
        #therefore, we still need to 
            #include the genetic position
            #remove the allele names


    print_text("load the raw_map file", header=4)
    snp_map_raw = pd.read_csv(\
        "./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_raw.map.gz", \
        sep=" ", \
        header=None, \
        low_memory=False)
    print(snp_map_raw)


    print_text("chr " + selected_chromosome + " : check we have the correct number of SNPs in the map file loaded in python", header=4)
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./" + input_vcf_file_map_calc + " | \
            wc -l); \
        if [[ $n_snps -eq " + str(snp_map_raw.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 


    print_text("rename the columns", header=4)
    snp_map_raw = snp_map_raw.rename( \
        {0: "chr", 1: "id_old", 2: "pos", 3: "ref", 4: "alt"}, \
        axis=1)
        #we name ID as old because this is the ID coming from the VCF file, which we need for selecting those variants in the VCF file with genetic position. The final ID will be in plink format, see below.
        #use a dict with old and new column names. indicated we are renaming columns (axis=1)
    print("see map file with renamed columns")
    print(snp_map_raw)


    print_text("chr " + selected_chromosome + ": create a new ID variable following plink 2.0 format and check it was correctly created", header=4)
    #the original ID column is in the format "CHR:POS:REF:ALT" but the problem is that some SNPs have their position indicated in the ID is shifted. 1KGP authors separated multiallelic SNPs in different lines with bcftools norm and then shifted their position so they could be phased, combing back to the original position afterwards (see next line). I guess during that process, they updated the IDs using chrom and the shifted position was used in the ID. Therefore, even if POS comes back to the original position, the ID remains with the shifted position. I have checked several multiallelic SNPs, and they have all the same issue with the position in the ID.
        #From README ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf"): "SHAPEIT2 does not handle multiallelic variant phasing. To phase both biallelic and multiallelic variants we first split the multiallelics into separate rows while left-aligning and normalizing INDELs using bcftools norm tool (Li, 2011). Next, we shifted the position of multiallelic variants (2nd, 3rd, etc ALT alleles) by 1 or more bp (depending on how many ALT alleles there are at a given position) to ensure a unique start position for all variants, which is required for SHAPEIT2. We shifted the positions back to the original ones after phasing".
    #This is not very convenient because when later we create the hap files, SNPs will have the plink format for ID to avoid strand flips (CHROM:POS_REF_ALT), so we are going to have different IDs between hap and map files, making more difficult to do checks.
    #Therefore, we are going to update the ID of each SNP using plink format, and ensuring in this way SNPs will be names the same in hap and map files.
    #bcftools convert to hap format will use "CHROM:POS_REF_ALT" for the IDs, so we need to follow that format.
    snp_map_raw["id"] = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
    print("check")
    check_id = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
        #make a series combining chromosome, pos, ref and alt, and using the corresponding separators
    print(check_id.equals(snp_map_raw["id"]))
        #check it is identical to id
    print(snp_map_raw)
    print("chr " + selected_chromosome + ": remove the ref/alt columns as we have this information already included in the ID")
    snp_map_raw = snp_map_raw.drop(["ref", "alt"], axis=1)
    print(snp_map_raw)


    print_text("load and explore the decode2019 map", header=3)
    #I know that the original 2019 decode map is alligned to hg38. Also, I assume that the decode2019 map is 1-based because they do not specify is 0-based. I assume that if you say anything, base 1 in your coordinates is base 1 in the genome. I assume this is the default.
        #Data S3.genetic.map.final.sexavg.gor.gz:
            #average genetic map computed from the paternal and maternal genetic maps, which were in turn computed from the paternal and maternal crossover, respectively. The data columns are as follows: Chr (chromosome), Begin (start point position of interval in GRCh38 coordinates), End (end point position of interval in GRCh38 coordinates), cMperMb (recombination rate in interval), cM (centiMorgan location of END POINT of interval)
            #Page 85 of "aau1043-halldorsson-sm-revision1.pdf"
        #Importantly, we have confirmation about this from the first author of the paper:
            #I sent an email to the first author of the deCODE 2019 map paper asking whether the coordinates are 1- or 0-based, but I got an answer a bit difficult to understand for me.
            #I asked this: "I would like to confirm that the three genetic maps published as supplementary files (Data S1-S3, i.e., maternal, paternal and average maps) have 1-based coordinates. From the description, I understand this is the case, but I want to double check that with you to ensure these are not 0-based."
            #He replied: "This shouldn't be the case, but let us know if you think there is a problem."
            #David answered: "About the deCode map, he is saying that the coordinates are 1-based, we can proceed."

    #1KGP is aligned to hg38 (see paper) and coordinates are 1-based as VCF format 4.2 has 1-based coordinates.

    #therefore, we have the same position format in both datasets, so we can just use the decode 2019 map to calculate the genetic position of each SNP.

 
    print_text("chr " + selected_chromosome + ": see first lines of the Data S3 of decode paper, which is the sex average map (see above). The file has header", header=4)
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        head -20")

    print_text("chr " + selected_chromosome + ": remove the header and save", header=4)
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\"; \
                selected_chromosome=\"chr" + selected_chromosome + "\"; \
                header=\"yes\"; \
            }{ \
                if($0 ~ /Chr\tBegin\tEnd\tcMperMb\tcM/){print $0; header=\"no\"}; \
                if(header == \"no\" && $1==selected_chromosome){print $0} \
            }' | \
        gzip \
            --force > ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_decode_map.tsv.gz; \
        gunzip \
            --stdout \
            ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_decode_map.tsv.gz | \
        head -5")
        #decompress the map
        #load into awk
            #begin
                #using tabs as delimiter for inputs and outputs
                #also set the variable header as yes
                #add the chromosome name as a variable
            #if the row is not the column names, header remains as "yes"
                #use regex so we can use "\t" as a pattern not look for
            #but if we reach the column names, print them and then set header as "no", so the next lines with positions will be printed
            #only print if header="no" and the interval belongs to the selected chromosome.
    
    print_text("chr " + selected_chromosome + ": load decode 2019 map into python", header=4)
    decode2019_map = pd.read_csv(\
        "./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_decode_map.tsv.gz", \
        sep="\t", \
        header=0, \
        low_memory=False)

    print_text("chr " + selected_chromosome + ": load the whole decode map to check we have selected the correct SNPs", header=4)
    decode2019_map_check = pd.read_csv(\
        "./data/decode_2019/aau1043_datas3.gz", \
        sep="\t", \
        header=0, \
        low_memory=False, \
        skiprows=[i for i in range(0,7,1)])
        #skip the first rows with the header
    print("select the rows of the selected chromosome, reset the index and then check the resulting DF is identical to the decode map generated with awk")
    print(decode2019_map_check \
        .loc[decode2019_map_check["Chr"] == "chr"+selected_chromosome, :] \
        .reset_index(drop=True) \
        .equals(decode2019_map))
    print("remove the whole decode map")
    del(decode2019_map_check)
    import gc
    gc.collect()

    print_text("rename columns in lower case", header=4)
    decode2019_map=decode2019_map.rename({"Chr":"chr", "Begin":"begin", "End":"end", "cMperMb":"cM_Mb", "cM":"cM"}, axis=1)
    print(decode2019_map)
        #Average genetic map computed from the paternal and maternal genetic maps.
        #The data columns are as follows:
        #Chr (chromosome)
        #Begin (start point position of interval in GRCh38 coordinates)
        #End (end point position of interval in GRCh38 coordinates)
        #cMperMb (recombination rate in interval)
        #cM (centiMorgan location of END POINT of interval)

        #I did a lot of checks on this map regarding overlapping of the intervals etc, check recomb_v3.R in method_deep paper for further details.

    print_text("chr " + selected_chromosome + ": subset decode map for the selected chromosome. We already did that with awk, but are doing it again just in case", header=4)
    decode2019_map_subset = decode2019_map.loc[decode2019_map["chr"] == "chr"+str(selected_chromosome),:]
    print(decode2019_map_subset)
    print("Do we selected the correct chromosome?")
    unique_chrom_decode_subset = decode2019_map_subset["chr"].unique()
    if (len(unique_chrom_decode_subset)==1) & (unique_chrom_decode_subset=="chr"+str(selected_chromosome)) & (decode2019_map_subset.equals(decode2019_map)):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: The subset of deCODE map does not have the correct chromosome")

    print_text("remove the full decode map from python", header=4)
    del(decode2019_map)
    import gc
    gc.collect()


    print_text("perform the calculation of the genetic position of each SNP", header=3)
    print_text("chr " + selected_chromosome + ": Convert the raw SNP map and the deCODE map for the select chromosome to global variable", header=4)
    globals()["snp_map_raw_chr"+selected_chromosome] = snp_map_raw
    globals()["decode2019_map_subset_chr"+selected_chromosome] = decode2019_map_subset

    print_text("chr " + selected_chromosome + ": Run the function on just one snp", header=4)
    print(gen_pos(snp_map_raw.iloc[10]["id"]))

    print_text("chr " + selected_chromosome + ": run function across SNPs using a pool for parallelization", header=4)
    #set the pool
    import concurrent.futures as mp
    nested_executor=mp.ProcessPoolExecutor(max_workers=5)
        #The ProcessPoolExecutor class is an Executor subclass that uses a pool of processes to execute calls asynchronously. ProcessPoolExecutor uses the multiprocessing module, which allows it to side-step the Global Interpreter Lock but also means that only picklable objects can be executed and returned.
            #https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ProcessPoolExecutor
            #max_workers
                #the number of processes
                #An Executor subclass that executes calls asynchronously using a pool of at most max_workers processes. If max_workers is None or not given, it will default to the number of processors on the machine. If max_workers is less than or equal to 0, then a ValueError will be raised.
            #mp_context
                #mp_context can be a multiprocessing context or None. It will be used to launch the workers. If mp_context is None or not given, the default multiprocessing context is used.
                #we are using the default context on POSIX except macOS, which is "fork".
                    #"fork": The parent process uses os.fork() to fork the Python interpreter. The child process, when it begins, is effectively identical to the parent process. All resources of the parent are inherited by the child process. Note that safely forking a multithreaded process is problematic.
                    #set_start_method() of multiprocessing can set the type of context
                        #from multiprocessing import set_start_method
                        #mp_context=set_start_method("spawn")
                    #Users of macOS or users of libc or malloc implementations other than those typically found in glibc to date are among those already more likely to experience deadlocks running such code.
                    #https://docs.python.org/3/library/os.html#os.fork
                    #https://docs.python.org/3/library/multiprocessing.html#multiprocessing-start-methods
                #when forking, the child inherits everything (global variables, modules...) except threads that are already running in the parent. In contrast, when using spawn, you get a fresh python interpreter. The problem with the former option is that the child could try to write something in the same place than the parent is also writing. I do not fully understand why, but it seems the file can be lock the first time and then the child cannot open it. This is exactly correct, but you can get an idea.
                    #https://pythonspeed.com/articles/python-multiprocessing/
                    #https://britishgeologicalsurvey.github.io/science/python-forking-vs-spawn/
                    #https://stackoverflow.com/questions/64095876/multiprocessing-fork-vs-spawn
                #This is a deadlock, resources shared are locked, and none of the parallel processes can move forward. For example: Suppose process P1 is waiting for a resource R1 currently being used by process P2. Meanwhile, P2 is waiting for resource R2 that's being used by P1. Neither process is able to proceed.
                    #https://devopedia.org/deadlock
                #I do not think we have a problem in my case, because the initial fork of master processor generates independent processes that are working on different chromosomes. There is one file, the deCODE map, that is used for each of the different processes, but just read, then a subset of that file for the corresponding selected chromosome is created, meaning each process writes in different files. Therefore, I do not see a problem here.
                #In any case, I am doing checks to check whether all SNPs of the corresponding chromosome have been analyzed and also compare for the first few first thousand lines with map, which does not use parallelization.
    #run the genetic position function in parallel within the pool
    #IF DEBUGGING, CHECK gen_pos AT THIS POINT
    final_genetic_pos_list = nested_executor.map(gen_pos, snp_map_raw["id"])
    final_genetic_pos = list(final_genetic_pos_list)
        #the first line is lazy, so we need to actually ask for a list with the results, i,e., a list of tuples as the function return a tuple for each SNP id
    #close the pool
    nested_executor.shutdown(wait=True)
        #wait: If True then shutdown will not return until all running futures have finished executing and the resources used by the

    print_text("convert the tuple to DF and add the column names", header=4)
    final_genetic_pos_df = pd.DataFrame(final_genetic_pos, columns=["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "selected_snp_physical_pos", "check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6", "genetic_distance", "left_cM", "right_cM", "distance_left_end", "distance_right_end"])
    print("see results:")
    print(final_genetic_pos_df)

    print_text("calculate a few thousand rows with the original 'map' function and compare with the previous results to check we did not mess with the parallelization", header=4)
    #randomly select SNP IDs
    if debugging == True:
        n_snps_to_check = 10000
    elif debugging == False:
        n_snps_to_check = 50000
    ids_to_check = snp_map_raw.sample(n=n_snps_to_check, replace=False)["id"]
        #Return a random sample of items from an axis of object
        #we set replacement as False, so the same row cannot be select two times
        #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sample.html
    #calculate the genetic position of these SNPs using map
    final_genetic_pos_2 = list(map(gen_pos, ids_to_check))
    #convert the result to pandas DF
    final_genetic_pos_2_df = pd.DataFrame(final_genetic_pos_2, columns=["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "selected_snp_physical_pos", "check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6", "genetic_distance", "left_cM", "right_cM", "distance_left_end", "distance_right_end"])
    #select the randomly selected SNPs from the previous DF with results, which has been obtained using the new parallelization method
    results_new_method = final_genetic_pos_df \
        .iloc[np.where(final_genetic_pos_df["selected_snp_id"].isin(ids_to_check))[0],:] \
        .set_index("selected_snp_id", drop=False) \
        .sort_index() \
        .reset_index(drop=True) \
        .fillna(0)
        #select rows where the ID is included in the Series with the IDs to be checked
        #convert the SNP IDs as index, but without removing it as a column
            #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.set_index.html
        #sort the rows using the new index, i.e., SNP ID
        #remove the index, not including it as a new column
        #fill NANs with zero
    #format the results obtained with map so they can be compared with the previous results obtained with the new parallelization method
    results_old_method = final_genetic_pos_2_df \
        .set_index("selected_snp_id", drop=False) \
        .sort_index() \
        .reset_index(drop=True) \
        .fillna(0)
        #like before, set SNP ID as index, reorder using it, and then remove the index, finish converting NANs to zero
    #check whether each value is NOT the same between both approaches. Sum 
    check_old_new_parallel_method = (results_new_method != results_old_method).apply(sum, axis=0) == 0
        #axis=0 so if performs the operation considering all rows of a given column, i.e., number of True cases per column
        #the total number of Trues per column should be zero.
    print(check_old_new_parallel_method)
    #the check should be True for all columns, if not raise error
    if (check_old_new_parallel_method.shape[0] != check_old_new_parallel_method.sum()):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM: The genetic position calculation for 50K random rows is not the same between 'map' and the new parallelization method.")

    print_text("delete the global variables created for the selected chromosome", header=4)
    del globals()["snp_map_raw_chr"+selected_chromosome]
    del globals()["decode2019_map_subset_chr"+selected_chromosome]
        #https://stackoverflow.com/a/47717307

    print_text("chr " + selected_chromosome + ": all checks of genetic position calculation are True?", header=4)
    checks_across_snps_pos = final_genetic_pos_df[["check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6"]].all(axis=0, skipna=True)
        #Exclude NA/null values. If the entire row/column is NA and skipna is True, then the result will be True, as for an empty row/column. If skipna is False, then NA are treated as True, because these are not equal to zero.
    print(checks_across_snps_pos)
        #important:
            #all() does not consider nan, so if you have nan and the rest True, the output is True.
            #this is ok for us, because we use nan in some conditions.
    if (checks_across_snps_pos.shape[0] != checks_across_snps_pos.sum()):
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM: There is an error in one of the checks for the calculation of the genetic position of SNPs")

    print_text("chr " + selected_chromosome + ": check we have the correct number of SNPs in the calculation of genetic position", header=4)
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./" + input_vcf_file_map_calc + " | \
            awk 'END{print NR}'); \
        if [[ $n_snps -eq " + str(final_genetic_pos_df.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    print_text("chr " + selected_chromosome + ": check that we have the exact same snps than in the raw map", header=4)
    list_check_all_snps = list()
    list_check_all_snps.append(np.array_equal(
        snp_map_raw["chr"].to_numpy(),
        ("chr" + final_genetic_pos_df["selected_chromosome"]).to_numpy()))
    list_check_all_snps.append(np.array_equal(
        snp_map_raw["id"].to_numpy(),
        final_genetic_pos_df["selected_snp_id"].to_numpy()))
    list_check_all_snps.append(np.array_equal(
        snp_map_raw["id_old"].to_numpy(),
        final_genetic_pos_df["selected_snp_old_id"].to_numpy()))
    list_check_all_snps.append(np.array_equal(
        snp_map_raw["pos"].to_numpy(),
        final_genetic_pos_df["selected_snp_physical_pos"].to_numpy()))
    print(list_check_all_snps)
    if False in list_check_all_snps:
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM: There is an error when comparing the raw SNP map and the new one created with deCODE data")
    
    print_text("recalculate genetic distance of each SNP", header=4)
    #but we exclude SNPs that have cM exactly in their position in the decode map or SNPs without decode data around
    final_genetic_pos_df_check_dist_calc = final_genetic_pos_df.loc[\
        (~final_genetic_pos_df["genetic_distance"].isna()) & \
        (~final_genetic_pos_df["left_cM"].isna()), :]
    #calculate the physical distance between the end of the two deCODE ranges, that is, the distance from the closest interval to the SNP of the left PLUS the distance from the SNP to the closest deCODE interval from the right
    phys_distance_decode_intervals = final_genetic_pos_df_check_dist_calc["distance_left_end"] + final_genetic_pos_df_check_dist_calc["distance_right_end"]

    #calculate the genetic distance between the end of the two deCODE interval closest to the SNP (on both sides, right and left)
    gen_distance_decode_intervals = final_genetic_pos_df_check_dist_calc["right_cM"] - final_genetic_pos_df_check_dist_calc["left_cM"]

    #extract the physical distance from the SNP to the closest deCODE interval to the left
    phy_distance_left_decode = final_genetic_pos_df_check_dist_calc["distance_left_end"]

    #calculate the genetic distance: If the physical distance between the end of deCODE intervals (phys_distance_decode_intervals) corresponds with a known genetic distance (gen_distance_decode_intervals), the physical distance from the closest deCODE intervals from the left to the selected SNP (phy_distance_left_decode) would correspond with X; thus X = (gen_distance_decode_intervals*phy_distance_left_decode)/phys_distance_decode_intervals X is the genetic distance from the closest deCODE interval from the left to the SNP If you sum this to the genetic position of that closest deCODE interval from the left (final_genetic_pos$left_cM), you would have the genetic distance of the selected SNP. The genetic position of that interval gives the cM value until that point, and you just calculated the rest of cM increase until the SNP 
    new_genetic_distance = ((gen_distance_decode_intervals * phy_distance_left_decode) / phys_distance_decode_intervals) + final_genetic_pos_df_check_dist_calc["left_cM"]
        #see figure 31 and the calculation of genetic distance in "recomb_calc" function for the full explanation using the words of David and also my interpretation.

    print("chr " + selected_chromosome + ": compare the new genetic distance and the distance previously calculated. It is ok not having True here if the next check is ok")
    raw_check_gen_dis = new_genetic_distance == final_genetic_pos_df_check_dist_calc["genetic_distance"]
    if (raw_check_gen_dis.shape[0] != raw_check_gen_dis.sum()):
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM: There is an error in the calculation of the genetic distance of each SNP")

    print("chr " + selected_chromosome + ": check that the cases with NA for last checks but with genetic distance are the cases of SNPs in a position exactly with deCODe data")
    cases_gen_pos_no_last_checks = final_genetic_pos_df.loc[\
        (final_genetic_pos_df["check_4a"].isna()) & \
        (~final_genetic_pos_df["genetic_distance"].isna()), "selected_snp_id"]
    #extract the ID of snps with a position that have deCODE genetic position
    snps_with_decode_data = final_genetic_pos_df.loc[\
        final_genetic_pos_df["selected_snp_physical_pos"].isin(decode2019_map_subset["end"]), "selected_snp_id"]
    #print the check
    print(cases_gen_pos_no_last_checks.equals(snps_with_decode_data))

    print_text("check whether the SNPs with duplicated physical position are different alleles of the same multiallelic SNP", header=4)
    print("first check that no SNP has a duplicated ID")
    print(final_genetic_pos_df["selected_snp_id"].duplicated(keep=False).sum()==0)
    print(final_genetic_pos_df["selected_snp_old_id"].duplicated(keep=False).sum()==0)
        #keep=False: Mark all duplicates as True

    print("check that we only have the selected chromosome")
    print(final_genetic_pos_df["selected_chromosome"].unique() == selected_chromosome)

    print("select those SNPs with a duplicated physical position")
    snps_duplicated_position = final_genetic_pos_df.loc[final_genetic_pos_df["selected_snp_physical_pos"].duplicated(keep=False), :]
    print(snps_duplicated_position)

    print("define a function to check that each of the duplicated SNPs are indeed multiallelic SNPs")
    snps_duplicated_position_unique_positions = snps_duplicated_position["selected_snp_physical_pos"].unique()
        #get the list of physical positions through which we have to iterate
    #selected_position=snps_duplicated_position_unique_positions[0]
    def check_dup_pos_fun(selected_position):
        
        #select the SNPs with the selected physical position
        selected_snps_pos = snps_duplicated_position.loc[snps_duplicated_position["selected_snp_physical_pos"] == selected_position, :]
        
        #for each of the rows, i.e., SNPs, with the selected physical position, 
        old_ids_split = selected_snps_pos.apply(lambda x: x["selected_snp_old_id"].split(":"), axis=1)
            #select the old ID and then split it by ":"

        #check that the chromosome is the same in both SNPs 
        check_0 = old_ids_split.iloc[0][0] == old_ids_split.iloc[1][0]

        #check that the physical position just differ in a few bases.
        diff_position = np.abs(int(old_ids_split.iloc[0][1]) - int(old_ids_split.iloc[1][1]))
        check_1 = (diff_position >= 1 and diff_position < 15)
            #The position saved in the old ID seems to be the one used after separating multiallelic SNP shifting the position of their alleles in order to do phasing (SHAPEIT2 does not handle multiallelic SNPs). They put back the alleles to the original positions after phasing, but maybe they did not change the IDs because these SNPs have different position in the ID but their actual position is the same.
            #Note that we could have more than 3 multiallelic SNPs, so the shift of positions could be larger than 1, but I do not expect the difference to be very high, i.e., more than 15 bases.
            #"we shifted the position of multiallelic variants (2nd, 3rd, etc ALT alleles) by 1 or more bp (depending on how many ALT alleles there are at a given position) to ensure a unique start position for all variants, which is required for SHAPEIT2"
                #https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf

        #check that the REFs are the same but the ALTs are different, i.e., the expected behaviour of the alleles of the same multiallelic variant
        check_2 = old_ids_split.iloc[0][2] == old_ids_split.iloc[1][2]
        check_3 = old_ids_split.iloc[0][3] != old_ids_split.iloc[1][3]

        #return the results as a tuple
        return(tuple([selected_position, check_0, check_1, check_2, check_3]))

    print("run the function with just one position")
    print(check_dup_pos_fun(snps_duplicated_position_unique_positions[0]))

    print("run the function across all duplicated positions")
    check_dup_pos = list(map(check_dup_pos_fun, snps_duplicated_position_unique_positions))

    print("convert the resulting list of tuples into a DF")
    check_dup_pos_df = pd.DataFrame(check_dup_pos, columns=["selected_position", "check_0","check_1","check_2","check_3"])

    print("check we have analyzed all the duplicated positions")
    print(np.array_equal(check_dup_pos_df["selected_position"].to_numpy(), snps_duplicated_position_unique_positions))

    print("check we do NOT have NA")
    #x=check_dup_pos_df["check_0"]
    #x=check_dup_pos_df["selected_position"]
    print(check_dup_pos_df.apply(lambda x: (x==True).all(skipna=True) if x.name!="selected_position" else "No need to check", axis=0))
        #in each column, check whether each row is equals to True, then check that the result is True for ALL rows of the column.
            #skipna=True
                #If the entire row/column is NA and skipna is True, then the result will be True, as for an empty row/column.
                #no problem because we have previously checked that no NA is present. In addition, np.nan==True gives false, so we should have False for NA cases before doing all().
        #do this only for the columns that are NOT the "selected_position". If the column is the "selected_position", then just print that no check is required.


    print_text("chr " + selected_chromosome + ": prepare final map file", header=3)
    print_text("subset only the columns for map files", header=4)
    final_genetic_pos_map_file = final_genetic_pos_df[["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "genetic_distance", "selected_snp_physical_pos"]]        
    #add chrom
    chrom_column = final_genetic_pos_map_file.pop("selected_chromosome")
    final_genetic_pos_map_file.insert(0, "selected_chromosome", "chr"+chrom_column)
    del(chrom_column)
        #pop the chromosome column
        #then insert it as the first column (indicating the corresponding name) but adding "chr"
    print(final_genetic_pos_map_file)
        #save the chromosome, ID, genetic position and physical position. This is the format expected by hapbin
            #https://github.com/evotools/hapbin

    print_text("select those rows of the map for which the SNP have genetic position", header=4)
    final_genetic_pos_map_file_pruned = final_genetic_pos_map_file.loc[~final_genetic_pos_map_file["genetic_distance"].isna(), :]

    print("see the total number of SNPs and the number of SNPs removed due to the lack of genetic position")
    print(final_genetic_pos_map_file.shape[0])
    print(final_genetic_pos_map_file.shape[0] - final_genetic_pos_map_file_pruned.shape[0])

    print_text("select old_id from the map that have genetic position. this ID is the original retained from the VCF file, so we can use it to subset the VCF file", header=4)
    snps_id_with_gen_pos = final_genetic_pos_map_file_pruned.loc[:, "selected_snp_old_id"]

    print("check")
    print(final_genetic_pos_map_file_pruned["selected_snp_old_id"].equals(snps_id_with_gen_pos))


    print_text("save the names in a txt file", header=4)
    with open(r"./results/00b_map_files/chr" + selected_chromosome + "/list_snps_with_gen_pos_" + selected_chromosome + ".tsv", "w") as fp:
        fp.write("\n".join(snps_id_with_gen_pos))
            #each name in a different line so we have to add "\n" to the name
            #https://pynative.com/python-write-list-to-file/
        fp.write("\n")
            #add empty line at the end

    print("compress")
    run_bash(" \
        gzip \
            --force \
            ./results/00b_map_files/chr" + selected_chromosome + "/list_snps_with_gen_pos_" + selected_chromosome + ".tsv")

    print("take a look to the file")
    run_bash(" \
        gunzip \
            --stdout \
            ./results/00b_map_files/chr" + selected_chromosome + "/list_snps_with_gen_pos_" + selected_chromosome + ".tsv.gz | \
        awk \
            '{if(NR<20){print $0}}'")

    print_text("plot the genetic vs physical distance vs recombination rate to check whether there is a correlation", header=4)
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1)
    #make more vertical space between subplots
    fig.subplots_adjust(hspace=0.6)
            #https://stackoverflow.com/a/5159405/12772630
    fig.suptitle("Chromosome " + selected_chromosome)
    ax1.set_title("Genetic vs physical distance of SNPs")
    ax1.scatter(x=final_genetic_pos_map_file_pruned["selected_snp_physical_pos"], y=final_genetic_pos_map_file_pruned["genetic_distance"], s=0.1)
    ax1.set_xlabel("Pair-base position")
    ax1.set_ylabel("Genetic distance")
    ax2.set_title("Recombination rate vs physical position of deCODE intervals")
    decode2019_map_subset_plot = decode2019_map_subset.loc[ \
        (decode2019_map_subset["end"]>=np.min(final_genetic_pos_map_file_pruned["selected_snp_physical_pos"])) & \
        (decode2019_map_subset["end"]<=np.max(final_genetic_pos_map_file_pruned["selected_snp_physical_pos"])),:]
        #get only intervals overlapped with the SNPs
    ax2.scatter(x=decode2019_map_subset_plot["end"], y=decode2019_map_subset_plot["cM_Mb"], s=0.5)
    ax2.set_xlabel("Pair-base end position")
    ax2.set_ylabel("Recombination rate (cM/Mb)")
    plt.savefig( \
        fname="./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_plot_genetic_vs_phyiscal_distance.png")
    plt.close()

    print_text("see final map and save", header=4)
    print_text("remove old ID as we have already filtered the VCF file and check", header=4)
    final_genetic_pos_map_file_pruned = final_genetic_pos_map_file_pruned.drop(["selected_snp_old_id"], axis=1)
    print(final_genetic_pos_map_file_pruned.columns == ["selected_chromosome", "selected_snp_id", "genetic_distance", "selected_snp_physical_pos"])
    print(final_genetic_pos_map_file_pruned)
        #required format according to hapbin
            #The map files (--map) should be in the same format as used by Selscan with one row per variant and four space-separated columns specifiying 
                #chromosome, 
                #locus ID, 
                #genetic position
                #physical position.
    final_genetic_pos_map_file_pruned.to_csv(\
        "./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_selscan.map.gz", \
        sep=" ", \
        header=False, \
        index=False)

    print_text("check we have the correct number of rows and columns in the map file", header=4)
    run_bash("\
        n_rows=$( \
            gunzip \
                --stdout \
                ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NR}'); \
        n_cols=$( \
            gunzip \
                --stdout \
                ./results/00b_map_files/chr" + selected_chromosome + "/chr" + selected_chromosome + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NF}'); \
        if [[ $n_cols -eq 4 && $n_rows -eq " + str(final_genetic_pos_map_file_pruned.shape[0]) + " ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress map file to stdout and then calculate the number of rows (NR) and fields (NF). Do that only after the whole file has been read (END)
                #https://www.gnu.org/software/gawk/manual/html_node/Using-BEGIN_002fEND.html
            #the number of columns (fields) should be 4 following salescan format, while the number of rows should be equal to the number of snps we have in the map file loaded in python, which was indeed used to write this .map file.


    print_text("remove files not required anymore", header=3)
    run_bash(" \
        rm " + input_vcf_file_map_calc + "; \
        ls -l ./results/00b_map_files/chr" + selected_chromosome)


    print_text("FINISH", header=3)


    print_text("restore sys.stdout using the previously saved reference to it. This is useful if you intend to use stdout for other things only required if we are in production, as we changed stdout only in that case", header=3)
    if debugging==False:
        sys.stdout = original_stdout
        file_output.close()
            #redirect stdout to the original stsdout and close the file where the output has been saved
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/
            #https://stackoverflow.com/a/23838153



#####################
#### paralellize ####
#####################
print_text("paralellize", header=1)
print_text("create list with all chromosomes", header=2)
print_text("get chromosome names", header=3)
chromosomes = [str(i) for i in range(1, 23, 1)]

print_text("we are going to analyze 22 chromosomes?", header=3)
print((len(chromosomes) == 22))

print_text("See them", header=3)
print(chromosomes)


print_text("run parallel analyses", header=2)
print_text("open the pool", header=3)
import concurrent.futures as mp
master_executor = mp.ProcessPoolExecutor(max_workers=22)
    #we are using 22 processes, 1 per chromosome. Then, a nested executor will be run inside master_processor within each chromosome. That nested executor will use 5 processes, i.e., we can have up to 5 processes working on the calculation of the genetic position of each chromosome. Therefore 22*5 makes a total of 110 processes.
    #see the nested parallelization in master_processor for details about concurrent.futures

print_text("run function across chromosomes", header=3)
master_executor.map(master_processor, chromosomes)
    #in the nested executor, we obtained a list of tuples as result. In order to get that, I had to force nested_executor.map to do the list using list(nested_executor.map())
    #I guess in this case is not lazy because we are not getting an object as output but just run the code of master_processor and save the corresponding outputs as files.

print_text("close the pool", header=3)
master_executor.shutdown(wait=True)
    #wait: If True then shutdown will not return until all running futures have finished executing and the resources used by the




########################################################
#### Do some checks after analyzing all chromosomes ####
########################################################
print_text("Do some checks after analyzing all chromosomes", header=1)
print_text("Number of SNPs lost due to the lack of genetic position", header=2)
print_text("empty lists to save the counts", header=3)
count_snps=[]
count_snps_lost = []


print_text("run loop across chromosomes", header=3)
#chrom=1
for chrom in chromosomes:
    print("Doing chromosome " + str(chrom))

    print("from the output file of the selected chromosome, extract the row with the counts of high and low confidence ancestral alleles")
    row_results = run_bash(" \
        grep \
            --after-context=2 \
            'see the total number of SNPs and the number of SNPs removed due to the lack of genetic position' \
            ./scripts/00b_map_calcs_outputs/chr" + str(chrom) + ".out", return_value=True).strip()
        #--after-context: number of lines to print after the match
            #we want also the next two lines as these include the total number of SNPs and the SNPs lost due to the lack of genetic position

    print("see the row")
    print(row_results)

    print("split the row")
    row_results_split = row_results.split("\n")
    print(row_results_split)

    print("check that we only have one row, and it can be split in three parts with comma")
    if ("\n" not in row_results_split) & (len(row_results_split)==3):
        print("YES! GOOD TO GO!")
    else: 
        raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs lost due to lack of genetic position")

    print("append to each list the corresponding count")
    count_snps.append(row_results_split[1])
    count_snps_lost.append(row_results_split[2])


print_text("process the results", header=3)
print_text("convert to int each count if the count is NOt zero", header=4)
count_snps = [int(x) if x!="" else 0 for x in count_snps]
count_snps_lost = [int(x) if x!="" else 0 for x in count_snps_lost]

print_text("check we have all chromosomes", header=4)
#selected_list=[count_snps, count_snps_lost][0]
check_count = [len(selected_list) == 22 for selected_list in [count_snps, count_snps_lost]]
if sum(check_count) == len(check_count):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs lost due to the lack of genetic position")

print_text("make the sum across chromosomes", header=4)
count_snps_total = sum(count_snps)
count_snps_lost_total = sum(count_snps_lost)

print_text("see the sums", header=4)
print("Total number of SNPs is: " + str(count_snps_total))
print("Total number of SNPs lost due to lack of genetic position: " + str(count_snps_lost_total))

print_text("calculate the percentage of SNPs lost due to the lack of genetic position across all chromosomes", header=4)
print(count_snps_lost_total/count_snps_total*100)



print_text("now check we do NOT have any errors in the output files of all chromosomes", header=2)
print_text("run loop across chromosomes", header=3)
#chrom=1
for chrom in chromosomes:
    print_text("Doing chromosome " + str(chrom), header=4)

    print_text("count number of cases with 'error' or 'false' in the output file", header=4)
    count_error_false = run_bash(" \
        grep \
            'error|false' \
            --extended-regexp \
            --ignore-case \
            --count \
            ./scripts/00b_map_calcs_outputs/chr" + str(chrom) + ".out || \
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
        raise ValueError("ERROR! FALSE! WE HAVE ERROR/FALSE IN THE OUTPUT FILE OF CHROMOSOME NUMBER " + str(chrom))

    print_text("check we have the row of FINISH", header=4)
    check_finish = run_bash(" \
        grep \
            '## FINISH ##' \
            --count \
            ./scripts/00b_map_calcs_outputs/chr" + str(chrom) + ".out  || \
        [[ $? == 1 ]]", return_value=True).strip()
        #check we have the output indicating finish with a specific way, in uppercase and separated by "#", so we avoid the flag "--ignore-case". Ask for the count.
        #as in the previous case, we add a line to avoid errors if the count is zero. We will deal with the lack of FINISH in the next line indicating also the chromosome name for which we found an error
    if check_finish == "1":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! THE SCRIPT HAS NOT FINISHED FOR CHROMOSOME NUMBER " + str(chrom))




####################
#### Next steps ####
####################
print_text("Next steps", header=1)
#add nested parallelism
    #check that you are using executors correctly
    #check that you get the same than with mp
    #run in the HPC
#check the script from "por aqui"
#build again container with the final script and run it in the HPC
#after 33 hours running, only chromosomes 19, 21 and 22 where done, which are the shortest ones. We at least will need 120 hours (5 days)
#check how many SNPs we lose due to the lack of genetic position
#check the plots of physical vs genetic distance. We should see higher increases in genetic distance in the SNPs when recombination increases in the deCODE intervals
    #show david abrupt changes in chromosome 1... they correlate with increases in recombination rate, so we should be fine
    #also ask him about the accesibility mask of the next step. We already losing a lot of SNPs with the filters...
