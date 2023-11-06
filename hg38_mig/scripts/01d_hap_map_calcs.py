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

#To calculate hap files per pop, we will use as input the VCF files that include the ancestral allele as a new field, which was obtained with VEP. We will use the map files calculates across chromosomes to calculate map files per pop.

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

#save name of path including input vcf files
input_vcfs_path = "results/00_vep_vcf_files"

#create folders to save the results
run_bash(" \
    mkdir \
        -p ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep; \
    mkdir \
        -p ./results/01_cleaned_vep_vcf_files; \
    mkdir \
        -p ./results/02_hap_map_files_raw; \
    mkdir \
        -p ./results/03_hap_map_files; \
    mkdir \
        -p ./scripts/01_hap_map_calcs_outputs")
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



####################################################
# check bcftools's behaviour with a dummy vcf file including ancestral alleles calculated with VEP #
####################################################
print_text("Check bcftools's behaviour with a dummy vcf file including ancestral alleles calculated with VEP", header=1)

print_text("explore the dummy VCF file and do some operations", header=2)
#
print_text("see the dummy variants", header=3)
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #We have variants with different characteristics. In some cases, the AN, AC fields are not correct, but we use this to check whether the different commands we use can correctly update these fields.
            #Biallelic SNPs:
                #SNP rs6054247 chr20 14280 T A 5 0 0.2 AA: . AA_upcase: . GTs: 0|0 0|0 1|1
                #SNP rs6054248 chr20 14290 C A 5 0 0.2 AA: . AA_upcase: . GTs: 0|0 0|0 1|.
                #SNP rs6054249 chr20 14300 C A 5 1 0.2 AA: . AA_upcase: . GTs: 0|0 0|0 1|.
                #SNP rs6054250 chr20 14310 C A 3 0 0 AA: N AA_upcase: N GTs: 1|0 .|. .|.
                #SNP rs6054251 chr20 14320 C A 6 2 0.333 AA: - AA_upcase: - GTs: 0|0 0|0 1|1
                #SNP rs6054252 chr20 14350 C A 4 2 0.5 AA: c AA_upcase: C GTs: 1|0 1|0 .|.
                #SNP rs6054257 chr20 14370 G A 6 3 0.5 AA: c AA_upcase: C GTs: 1|0 1|1 0|0
                #SNP rs6054255 chr20 14371 G C 6 3 0.667 AA: . AA_upcase: . GTs: 0|1 1|1 1|0
                #SNP rs6040351 chr20 17330 T A 6 2 0.333 AA: . AA_upcase: . GTs: 0|0 1|0 1|0
            #Multiallelic SNPs:
                #SNP rs6040355 chr20 1110696 A G,T 6 2,2 0.333,0.333 AA: G,G,G,G,G,G,G,G,G,G,G,G,G,G AA_upcase: G,G,G,G,G,G,G,G,G,G,G,G,G,G GTs: 1|1 2|2 0|0
                #SNP rs6040356 chr20 1110697 A C,T 6 2,1 0.333,0 AA: C,C,C,C,C,C,C,C,C,C,C,C,C,C AA_upcase: C,C,C,C,C,C,C,C,C,C,C,C,C,C GTs: 1|1 0|0 0|0
                #SNP rs6040357 chr20 1110698 A G,T 6 3,3 0.5,0.5 AA: A,A,A,A,A,A,A,A,A,A,A,A,A,A AA_upcase: A,A,A,A,A,A,A,A,A,A,A,A,A,A GTs: 1|1 2|2 1|2
                #SNP rs6040358 chr20 1110699 G A,T 6 2,0 0.333,0 AA: G,G,G,G,G,G,G,G,G,G,G,G,G,G AA_upcase: G,G,G,G,G,G,G,G,G,G,G,G,G,G GTs: 1|1 0|0 .|.
                #SNP rs6040359 chr20 1110700 A G,T 4 0,3 0,0.6 AA: T,T,T,T,T,T,T,T,T,T,T,T,T,T AA_upcase: T,T,T,T,T,T,T,T,T,T,T,T,T,T GTs: 2|2 2|2 2|.
            #Exact duplicate SNPs, i.e., pos, chr, REF, alt
                #SNP rs6040360 chr20 1110701 A G 6 3 0.5 AA: G,G,G,G,G,G,G AA_upcase: G,G,G,G,G,G,G GTs: 1|1 0|0 0|1
                #SNP rs6040360_copy chr20 1110701 A G 6 3 0.5 AA: G,G,G,G,G,G,G AA_upcase: G,G,G,G,G,G,G GTs: 1|1 1|0 0|1
            #SNP with unphased data
                #SNP rs6040361 chr20 1110702 A G 6 3 0.5 AA: C,C,C,C,C,C,C AA_upcase: C,C,C,C,C,C,C GTs: 1|1 0/0 0|1
            #microsatellite (indel)
                #INDEL microsat1 chr20 1110703 GTC G,GTCT . . . AA: TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT AA_upcase: TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT,TT GTs: 0/1 0/2 1/1


#
print_text("split multiallelic SNPs ", header=3)
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
    #split multiallelic SNPs in different lines
        #--multiallelic -snps
            #split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+).
            #https://samtools.github.io/bcftools/bcftools.html#norm
    #output
        #"Lines   total/split/realigned/skipped: 18/5/0/0"
            #as we are splitting multiallelics and we have 5 multi, 5 lines are split.
    #when a multiallelic snp (e.g., REF=A, ALT=G,T) is split in several lines, these lines have the same position and REF allele, but different ALT. This was done in 1000 genomes data. They also add 1 to the position of each new line to do the phasing, but then they put all lines back in the same position.
    #Both lines of the same SNP have the same ancestral allele, which makes sense.
    #In these new lines, the ALT is now only one, but the AC still has two fields, i.e., count for both ALTs in both lines. the 1000 genomes project SNPs that are multiallelic have their AC field just with one value so I guess they used +fill-tags to update these fields (see below).
            #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf 
        #First scenario: rs6040355
            #the first line shows the genotypes for A and G, so a sample that is G|G would be 1|1, but a T|T sample would be 0|0, so A and T gets 0 in this line. This is ok, because the next line shows genotypes for A and T, and that T|T sample will be 1|1.
                #This case is multiallelic, so we should remove both lines. We can just combine these two lines with norm --multiallelic, as it looks for snps with same pos and similar REF to merge (see when --multiallelic +snps is first used in this script).
        #Second scenario: rs6040356
            #we could also have a case where the second ALT does not exist in our subpop, then the first line would be 0s (for A) and 1s (for G), while the second one would be all zeros, because T (1) is not present.
                #We need to retain the first row but not the second. This can be solved removing monomorphic so the second row with all 0 is removed, then pass duplicates filter, which should not affect the first line, as its "sister" row has been removed. When then we combine lines with multiallelic +snps, this remaining line should not have any other line with the same position and REF/ALT.
        #Third scenario: rs6040357
            #we could also have a case where the first line is all 1 (G) except for the cases that are T (0), zero in this case, not A (REF). In the second row, all 1 (T), except those that are G (0) in this case, not A (REF). Therefore, we do not have the REF in the subset, and only the two ALT alleles.
                #We have several ALT alleles, and no REF, so we should remove these SNPs, all lines. This case would be solved just by merging these lines with --multiallelic and then filter out those with number of ALT>2.
        #Fourth scenario: rs6040358
            #we could also have a scenario similar where the second ALT is not present in the subset, and we also have missing data.
            #the second row is for the second ALT and is "0|0 0|0 .|.", thus should be removed with our approach that uses AC=AN and AC=0 to remove monomorphic.
        #Fifth scenario: rs6040359
            #the first ALT is not present, indeed all samples are homozigous for the second ALT (T) except the last that has one T and then missing.
            #the first row (0|0 0|0 .|.) and the second (1|1 1|1 1|.) should be removed. Our approach will consder both as monomorphic as in the first case AC=0 and in the second AC=AN.

#
print_text("select only two samples to check whether AC, AN and AF changes", header=3)
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #the order you use in --samples is the order that will follow the genotypes. So if you --samples NA00002,NA00001, you get first the genotype of NA00002 and then those of NA00001. In other words, the genotype columns are ordered based on the sample ID you used as input.
        #you can see how AN and AC are updated but not AF. 
            #For example, rs6040357 has three copies of the two ALT alleles (AC=3,3) and a total number of 6 alleles (AN=6), but after removing the third sample (1|2), AC becomes 2,2 and AN becomes 4.
                #Also note that both lines are updated with the same allele count. They both have 2,2. Therefore, the two lines are still connected and bcftools consider them as part of the same SNP.
            #Similarly, rs6040351 has AN=6, 2 copies of the allele and frequency of 2/6=0.33, but after removing the third sample (1|0), AN becomes 4, AC is 1, BUT AF remains 0.333, when it should be 0.25.
        #Also note that REF/ALT are NOT updated
            #rs6054251 has "0|0 0|0 1|1", being REF=C and ALT=A. Therefore, when the third sample is filtered out, the SNP has "0|0 0|0", there is only C. 
            #However, the REF and ALT columns remain the same.
        #Indeed, there is a command (--no-update) to avoid (re)calculating INFO fields for the subset, currently including only INFO/AC and INFO/AN. Therefore, they are clearly saying that AN and AC are the only INFO fields updated for the subset.
            #You can update other fields using bcftools +fill
        #Also note that --multiallelic -snps does NOT update AC and AN after splitting multiallelic SNPs, so they still have two allele counts even they are separated in different lines, being these counts connected (see above)
        #SUMMARY: AC and AN are updated after subsetting, but not AF.

#
print_text("select SNPs", header=3)
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools view \
        --types snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #--types/--exclude-types LIST
            #Select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other

#
print_text("see monomorphic snps: rs6054249 (0|0 0|0 1|.) is not considered as monomorphic which is right, given that the last genotype not only has missing but also an ALT allele. rs6054248 (0|0 0|0 1|.) is not included even having AC erronously set as zero, thanks we update AC field with +fill-tags. The cases considered as monomorphic are those with all REF (e.g., rs6040358) or all ALT (e.g., rs6040359) irrespectively if they have missing. For example, rs6040359 has '1|1 1|1 1|.' as GT and it is considered monomorphic despite having a missing allele in the last genotype", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #first update the AN and AC fields so we have the current number of non-ref alleles for each SNP and we avoid two allele counts for multiallelic because they are now separated in different lines (--multiallelic -snps) as found in the 1KGP data. 
            #remember that AC is "allele count in genotypes, for each ALT allele" while AN is "total number of alleles in called genotypes".
            #therefore
                #it is IMPORTANT that you do NOT have multiallelic SNPs, so we do not have INFO/AC with several numbers
                    #--samples does NOT remove the allele count of the other ALT! see below
                #also key to have this updated after selecting samples, so we have the number of alleles in the subset
            #see below for further details about "+fill-tags"
        #then include those variants whose allele count
            #is equal to the total number of alleles, i.e., all alleles are ALT
            #is equal to zero, i.e., there is no ALT so all alleles are REF 
        #https://www.biostars.org/p/360620/

#
print_text("see monomorphic snps after subseting two first samples: rs6054249 is now considered monomorphic because it has 0|0 0|0 1|., being 1|. removed after filtering out the third sample. We get a warning after subsetting because --sampls updated AC, but for splitted multiallelic snps, we have several allele counts in each line, so the first time it sees this (rs6040355), print warning and no more. But AC is correctly updated anyway, leaving just one count, thus breaking the connection between the allele counts of the lines of a given multiallelic snp", header=3)
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")

#
print_text("exlcude monomorphic snps", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")

#
print_text("check that subseting and the updating AN, AC, AF fields does not change REF/AF: I select only the third sample. For rs6054247, the genotype is 0|0 0|0 1|1, being the third sample the only one with ALT. Thus, after subseting, the genotype is 1|1. AN and AC are both 2 now, while AF is 1. In other words, the ALT is the most frequent, but this alleles remains as the ALT. REF and ALT are not switched", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --samples NA00003 | \
    bcftools +fill-tags \
        -- --tags AN,AC,AF | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")

#
print_text("see what missing definition uses bcftools in stats: Stats considers as missing both '.|.' and '0|.', resulting in 5 missing for the third sample (3 cases with .|. and 2 cases with 1|.). This makes sense because we are looking for the proportion of genotypes with missing data, '0|.' has missing data besides the allele '0'. Therefore, it should be counted as missing. We can consider this using GT='mis' (see below). Note that we are not going to filter by missingness anyways because the 1KGP data is already filtered for 5% missingness and David said that we should stick to that", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools stats \
        --samples -")
        #produce stats per sample (--sample) considering all samples ("-")

#
print_text("filter by frequency of missing. This approach counts the number of missing genotypes with GT='mis' ('.|.', '.', '0|.') and divided by the number of samples, so you get the proportion of samples with missing genotypes. The strength of this approach over calculating the tag F_MISSING (+fill-tags) is that you are counting anything with '.', so even '1|.' would be considered missing. To me this is more correct because we are filtering by genotype missingness, i.e., the proportion of genotypes with missing. Therefore, '1|.' has missing even though it has one allele. It is different in previous filters when we want to check if a SNP is monomorphic. In that situation, having 0|0 1|. is not monomorphic if we count the last allele in the genotype with missing. At the same time the genotype has missing (it should be count for missingness) and has an ALT allele (it is not monomorphic), we should consider both aspects. Note that we are not going to filter by missingness anyways because the 1KGP data is already filtered for 5% missingness and David said that we should stick to that", header=3)
print("#### SNPs with 1/3 or more of missing: we have snps with 1 out of 3 missing genotypes (e.g., rs6040358) and variants with 2 out of 3 missing (rs6054250) ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES >= 1/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
print("#### SNPs with 2/3 of missing: We get a variant with 2 out of 3 missing (rs6054250) ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES = 2/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
print("#### select variants with less than 1/3 of missing: We have lost all variants with missing because we only have variants 2/3 and 1/3 of missing ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 1/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #Lines of multiallelic rs6040358 have been removed because the last sample has a missing genotype. Note that rs6040358 has not been excluded as monomorphic as it has 1 and 0, but because it has missing.
print("#### select variants with less than 0.05 of missing: We have lost all variants with missing because we only have variants 2/3 and 1/3 of missing ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #1/2 would be 0.5, i.e., 50%. Therefore, 5% would be 0.05 (i.e., 0.5/10).

#
print_text("remove now the exact duplicates: it does not matter if we remove duplicates before or after subsetting samples. Subseting does not update chrom, pos and REF/ALT, which are the fields used by '--rm-dup exact' to remove exact duplicates. In an hypothetically scenario where we have two duplicated SNPs and one of them is monomorphic, we would remove first the monomorphic due to the previous filters and then retain the second. If you use --rm-dup before, you would select the SNP that appears first, that could be the monomorphic or the other one. Removing mono first force to select always the second, while removing duplicates before makes possible to select one or the other. I do not think this is relevant", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. As you can see for rs6040360, it does not matter the ID or the genotypes, this command only targets chr, pos, ref, and alt. 
        #bcftools norm --rm-dup exact selects only the first appearance and remove the next. Consequently, rs6040360_copy is filtered out, retaining rs6040360.
            #https://github.com/samtools/bcftools/issues/1089
        #https://samtools.github.io/bcftools/bcftools.html#norm

#
print_text("combine lines of multialleic snps in just one line", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #by combining the different lines of each multiallelic SNP, we update the ALT column, which now includes several ALT alleles. Therefore, we can filter out these SNPs with --max-alleles.
        #--multiallelic +snps
            #split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps
                #https://samtools.github.io/bcftools/bcftools.html#norm
        #this combines snps with the same position and at least equal REF or ALT
            #the normal scenario is to have same REF and different ALTs, so it just add the different ALTs in the ALT column.
                #this will be removed when filtering by max-alleles.
            #I have checked that it also merge when REF is different but ALT is the same, in that case it takes the REF of the first SNP.
                #We should NOT have two SNPs with the same position, same ALT but different REF.
                #According to the specification file, they say for REF "Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted". I understand that multiple bases are possible for indels but together, not separated by comma! Indeed, they specify that in ALT: "Comma separated list of alternate non-reference alleles".
                #In addition, if you have a line with two REF alleles separated by comma, bcftools considers this variant like OHTER instead of SNP! so by definition it cannot be a SNP.
                #According to the 2022 1KGP project, "BCFtools v1.9 (Li, 2011) was used to split multiallelic variants into multiple rows". Therefore, it is NOT possible that they could load SNPs with 2 REFs into bcftools, so there is no two SNPs with the same pos, ALT and different REF. It is not possible.
                #In summary, we should not worry about this scenario because it is unlikely to be present in our current data.
            #If merges when REF and ALT are the same.
                #If we have SNPs with the same REF and ALT and position, these will be removed with dup (see above), so no problem.
            #when REF and ALT are different, no merging is done and get an error.
                #if REF and ALT are different, we get error and we will see.
        #some genotypes appear as not phased within the new fused multiallelic SNPs, but that is not a problem because these will be removed anyway.

#
print_text("check what happens if we remove the third sample and then rs6040355 have no REF anymore, only two ALTs: Despite losing the REF, both lines of rs6040355 are merged into a multiallelic SNP, thus it can be removed with --max-alleles 2. This makes sense because --multiallelics +snps looks for snps with the same position and same REF or ALT columns", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")

#
print_text("now we can remove those SNPs with more than 1 allele in ALT. rs6040358 is now included because the second ALT is not present, only REF and ALT, having no missing (it has missing for the last sample that is filtered out)", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #--min-alleles/--max-alleles INT  
            #Minimum/maximum number of alleles listed in REF and ALT
        #This looks at REF/ALT, so we need these columns updated in order to use this. If you just do this when the multiallelic SNPs are splitted across different lines, then you would not remove them, because each one is presented like a biallelic SNP.
        #As explained above, it seems we cannot have more than 1 allele in REF, but we can have several allele in the ALT column, thus we are indeed filtering here those SNPs with a given number of ALT alleles.

#   
print_text("select only phased data", header=3)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #rs6040361 has "1|1 0/0 0|1" as GT, thus the second sample does not have phased data.
        #it is removed after applying the phased filter.
        #We separate --phased from --max/min-allele filter to be completely sure that first we remove multiallelic snps, and then we remove those snps with unphased data for any sample.
        #I guess it is not mandatory to do it in this order because if a SNP that is multiallelic has also unphased data, will be removed anyways, and biallelic and phased SNPs will not be affected. 
        #The order is more important in the case of --samples, as you could remove a SNP that has unphased data for some samples, but none of these samples are included in the subset, so you are removing the SNP while it has complete phased data for the samples you selected.

#   
print_text("select only SNPs in regions with high accessibility", header=3)
#general info
    #The 1000 Genomes Project created what they defined as accessibility masks for the pilot phase, phase one and phase three of the Project. Some other studies have similar files.
    #In phase three of the 1000 Genomes Project, using the pilot criteria 95.9% of the genome was found to be accessible. For the stricter mask created during phase three, 76.9% was found to be accessible. A detailed description of the accessibility masks created during phase three, the final phase of the Project, can be found in section 9.2 of the supplementary material for the main publication. The percentages quoted are for non-N bases.
    #While the above was generated on GRCh37, similar files were created on GRCh38 for the reanalysis of the 1000 Genomes Project data on GRCh38. HGSVC2 also have files listing regions of the genome that were not analysed.
        #https://www.internationalgenome.org/faq/are-there-any-genomic-regions-that-have-not-been-studied/
#directory
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/README.accessible_genome_mask.20160622
    #This directory contains GRCh38 genomic masks that help identify regions of the genome that are more or less accessible to next generation sequencing methods using short reads.
    #Each mask is summarized in a fasta file (one fasta file per chromosome), coded as follows:
        #N - the base is an N in the reference genome GRCh37
            #I guess these are regions not found in hg19
        #L - depth of coverage is much lower than average
        #H - depth of coverage is much higher than average
        #Z - too many reads with zero mapping quality overlap this position
        #Q - the average mapping quality at the position is too low
        #P - the base passed all filters
            #these are the regions accepted
        #0 - an overlapping base was never observed in aligned reads
    #Regions marked as L, H, Z or Q are less accessible to short reads. Although they can still be analyzed they are more prone to false positives. 
    #The calculation of accessibility is based on the alignment of whole genome low coverage data of 2691 phase3 samples to GRCh38.
    #The masks are useful for (a) comparing accessibility using current technologies to accessibility in the pilot project, and (b) population genetic analysis (such as estimates of mutation rate) that must focus on genomic regions with very low false positive and false negative rates.
    #there are two masks
        #the pilot masks that uses the definition of the pilot project and it is less stringent, accepting 95% of the genome. In other words, this mask was created during the pilot, the definition, and then applied to phase 3.
        #the strict mask that uses a more stringent definition and I think this is based on phase 3, accepting only 76% of the genome. This mask was defined during phase 3.
    #IMPORTANT:
        #We are going to use the bed file that contains the regions marked as "passed" (P) across all chromosomes. Therefore, we avoid LHZQ regions, which are more prone to false positives (see above).
#mask data here
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/
    #we will use the bed files with the regions has to be retained
        #"In addition to masked fasta files, a bed file of all passed sites can be found in this directory."

print_text("create a dummy mask", header=4)
run_bash(" \
    cd ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/; \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}' \
        <<< \"chr20\t14369\t14370\tpilot\nchr20\t1110695\t1110700\tpilot\nchr20\t1110690\t1110696\tpilot\nchr1\t14369\t17330\tpilot\" \
    > dummy_pilot_mask.bed; \
    cat dummy_pilot_mask.bed")
    #we are not using echo to create the file from a string, because "echo" changes a lot between versions. In some of them you need to use "-e" to be able to use slashs (e.g., \t), while in others not.
        #https://stackoverflow.com/a/8161194/12772630
    #we are directly using a string as input to awk using "<<<". We use a here-string:
        #<<< is known as here-string. Instead of typing in text, you give a pre-made string of text to a program. For example, with such program as bc we can do bc <<< 5*4 to just get output for that specific case, no need to run bc interactively. Think of it as the equivalent of echo '5*4' | bc. Here-strings in bash are implemented via temporary files, usually in the format /tmp/sh-thd.<random string>.
            #https://askubuntu.com/a/678919
            #https://stackoverflow.com/a/10959179/12772630
    #load a string with columns separated as tabs ("\t") and several lines indicated with "\n". Use as field delimiter both for input and output "\t". Print all fields in each row.

print_text("compress the dummy (mask) bed file to match what we will do with the real data. bcftools --targets-file can take a .bed.gz file as input (see below)", header=4)
run_bash(" \
    gzip \
        --force \
        --keep \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed")
    #--force: force overwrite of output file and compress links
    #--keep: keep (don't delete) input files

print_text("apply the dummy mask", header=4)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #select variants inside interest regions. We can do this after applying the rest of filters because if some of the previous filters consider position, it is the same position. For example, you remove exact duplicates, i.e., SNPs with the same chromosome, position and REF/ALT. Therefore, if we have two rows with the same position, they will be both remove if they are outside the interest regions, it does not matter if remove them together or first the duplicate and then remaining non-duplicate row. The same goes for multiallelic....
        #Also, no info field will be modified, we are just removing complete rows, i.e., SNPs, not samples or genotypes.
            #--targets-file
                #Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index.
                    #Therefore, we do not need to compress the VCF file and create and index that bcftools can use to randomly access different positions.
                #Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). 
                    #The columns of the tab-delimited file can contain either positions (two-column format: CHROM, POS) or intervals (three-column format: CHROM, BEG, END), but not both. Positions are 1-based and inclusive. 
                    #The columns of the tab-delimited BED file are also CHROM, POS and END (trailing columns [like pilot/strict] are ignored), BUT COORDINATES ARE 0-BASED, HALF-OPEN. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive).
                        #THIS IS OUR CASE.
                        #0-based means that position 0 in the bed file is base 1 in a chromosome, while position 1 in the bed file is base 2 in the chromosome.
                        #Half-open interval means that one of the extremes, the upper bound in this case, is not included.
                        #https://en.wikipedia.org/wiki/BED_(file_format)
                    #Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed.
                    #Note that sequence names must match exactly, "chr20" is not the same as "20".
                        #we have chrXX in both the bed and the VCF files
                    #Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE.
                        #I have checked that adding intervals of chromosomes not present in the VCF file does NOT change anything. Therefore, bcftools is not including variants of one chromosome because that coordinate in other chromosome is accepted.
                #Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets.
                #Note that sequence names must match exactly, "chr20" is not the same as "20".
                #Note that -t cannot be used in combination with -T.
                #Unlike -r, targets can be prefixed with "^" to request logical complement. For example, "^X,Y,MT" indicates that sequences X, Y and MT should be skipped.
                #Another difference is that using --regions, the existence of overlapping regions within FILE can result in duplicated out of order positions in the output.
                    #In other words, the resulting order of the vcf is compromised, but this does not seem to be the case for --targets.
                    #I have also check this in the dummy file with two overlapped intervals, seeing that there is no duplication of SNPs because of this.
                    #https://github.com/samtools/bcftools/issues/57
                #Yet another difference between the -t/-T and -r/-R is that -r/-R checks for proper overlaps and considers both POS and the end position of an indel, while -t/-T considers the POS coordinate only (by default; see also --regions-overlap and --targets-overlap).
                    #This is the default behaviour that can be changed using --targets-overlap
            #--targets-overlap: not used this time.
                #This option controls how overlapping records are determined: 
                    #set to pos or 0 if the VCF record has to have POS inside a region (this corresponds to the default behavior of -t/-T); 
                        #This option does NOT consider INDELS with POS at the end of a region because their end coordinate is outside. It starts just at the end of the regions, thus at least 1 base is outside.
                    #set to record or 1 if also overlapping records with POS outside a region should be included (this is the default behavior of -r/-R, and includes indels with POS at the end of a region, which are technically outside the region); 
                    #or set to variant or 2 to include only true overlapping variation (compare the full VCF representation "TA>T-" vs the true sequence variation "A>-").
            #Summary:
                #If you have INDELS, you should use --targets-overlap records so the INDEL is included even if it ends outside of the interval.
                #--targets is slower as it does not use index, but I avoid the indexing. I can just use it on my vcf file.
        #the filtering works
            #chr20:14369-14370
                #This interval starts in the chromosome 20 at position 14370 (not 14369! BED files are 0-based so we start at 0!) and ends and position 14371, being the latter not included.
                #rs6054257 (chr20:14370) falls fully within the interval (14370-14371), but rs6054255 (chr20:14371) is in the upper bound and it is not included.
            #chr20:1110695-1110700
                #rs6040356 (chr20:1110697) and rs6040358 (chr20:1110699) fully within.
                #rs6040360 (chr20:1110701) is in the last base of an interval (1110700 in BED is 1110701 in genome), out.
            #chr20:1110690-1110696
                #this interval is overlapped with the previous one, but there is no impact on the results. I do not see duplicated positions, so it seems this method is not sensitive to that.
            #chr1:14369-17330
                #rs6040351 (chr20:17330) does not fall within this or the previous intervals. Note that this interval is in chr1 but this variant is in chr20. This means that --targets correctly selects only intervals of the corresponding chromosome.

#
print_text("check you can filter the BED file by chromosome using awk", header=4)
run_bash(" \
    gunzip -c ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    awk \
        -F '\t' \
        '{if ($1 == \"chr20\") print $0}' \
        > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask_chr20.bed; \
    gzip --force ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask_chr20.bed; \
    gunzip -c ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask_chr20.bed.gz")
        #decompress the dummy bed file and sent it to stdout
        #process it with awk
            #-F is the delimiter, "\t" in this case
            #for each row, if the first field (column) has a value of "chr20", print all fields and save into a file.
            #compress the file even if a compressed file with the same name exists (--force).
                #https://unix.stackexchange.com/questions/399560/using-awk-to-select-rows-with-specific-value-in-specific-column
                #https://stackoverflow.com/questions/2961635/using-awk-to-print-all-columns-from-the-nth-to-the-last

#
print_text("check the cleaning of the BED file with awk", header=4)
run_bash(" \
    uniq_chrom=$(\
        gunzip -c ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask_chr20.bed.gz | \
        awk \
            -F '\t' \
            '{print $1}' | \
        uniq); \
    if [[ $uniq_chrom == 'chr20' ]];then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")
        #decompress the BED file and send it to stdout
        #process with awk
            #delimiter (-F) is tab
            #print all rows for the first column (chromosome name)
        #select the uniq cases
        #save into a variable
        #if the variable is "chr20", perfect

#
print_text("use +fill-tags to update INFO fields", header=3)
print_text("before using +fill-tags to update INFO fields, see the sample 1 and 2 to check AF is not updated", header=4)
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #As expected, after the subset of samples AN and AC are updated considering only the selected samples, while AF is not. For example, rs6040351 has only 1 ALT allele after selecting the two first samples, so AN=4, AC=1 but AF is still 0.333 instead of 0.25 (1/4).

print_text("update AF with fill-tags", header=4)
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools +fill-tags \
        -- --tags AN,AC,AF | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #we use +fill-tags to update fields that are not updated after subsetting like frequency of the alternative allele and create some additional fields.
            #+fill-tags has updated AF, so for example, rs6040355 has 0.5 for the frequency of both ALTs, while before it was 0.333. Similarly, the frequency of the ALT is 0.25 in rs6040351, which is correct.
        #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed due to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields.

#
print_text("update and create more fields with fill-tags", header=4)
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #we have create several fields that are included the 1KGDP vcf files or are related:
            #INFO/AN: Total number of alleles in called genotypes
            #INFO/AC: Allele count in genotypes
                #According to vcf specification file v2: 
                    #allele count in genotypes, for each ALT allele, in the same order as listed
            #INFO/AC_Hom: Allele counts in homozygous genotypes
            #INFO/AC_Het: Allele counts in heterozygous genotypes
            #INFO/AF: Allele frequency from FORMAT/GT or AC,AN if FORMAT/GT is not present
            #INFO/MAF: Frequency of the second most common allele
            #INFO/ExcHet: Test excess heterozygosity; 1=good, 0=bad
            #INFO/HWE: HWE test (PMID:15789306); 1=good, 0=bad
            #INFO/NS: Number of samples with data
        #This works after selecting the two first individuals, see for example:
            #rs6054250
                #GT: 1|0 .|.
                #AN=2
                #AC=1
                #AC_Hom=0
                #AC_Het=1
                #AF=0.5
                #MAF=0.5
                #NS=1
            #rs6040357
                #GT: 1|1 2|2
                #AN=4
                #AC=2,2
                #AC_Hom=2,2
                #AC_Het=0,0
                #AF=0.5,0.5
                #MAF=0.5
                    #I guess the second most common allele here is the second ALT that also have 0.5.
                #NS=2

#
print_text("remove all previous INFO and FORMAT fields except GT and create the fields you are interested in by using fill-tags but after applying all the filters", header=4)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS AA: %AA AA_upcase: %AA_upcase GTs:[ %GT]\n'")
        #before updating/creating fields, remove all INFO fields (except AA and AA_upcase) and all FORMAT fields (except GT) using annotate and then add the fields we are interested in.
            #--remove LIST
                #List of annotations to remove. 
                    #Use "FILTER" to remove all filters or "FILTER/SomeFilter" to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and "FORMAT" to remove all FORMAT tags. 
                    #To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be abbreviated to "INF" and "FORMAT" to "FMT".
            #https://samtools.github.io/bcftools/bcftools.html#annotate
        #It is ok to have the message "Incorrect number of AC fields". 
            #This happens when we have multiallelic SNPs separated in several lines. You have only 1 ALT in each line, but you still have two AC values separated by comma. Because of this, when you subset by sample, then bcftools tries to update AC, it sees that you have 1 ALT but two AC values, and prints that warning.
            #No problem at all, as AC is correctly updated with the number of ALT alleles we have in the subset genotypes. The same goes for bcftools +fill-tags.
            #In addition, 1KGDP data is already split and the AC/AN fields seem to be updated, with just one value per line. So we should not see this warning.

#
print_text("compare GT between applying or not the +fill-tags commands and the removal of fields", header=4)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --no-header")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools view \
        --no-header")
        #I have checked that the genotypes remain the same despite removing all previous INFO fields and all FORMAT fields (except GT) and then adding new INFO fields, so we are good here.


print_text("extract the position (index) of several columns in the VCF file, so we can be sure we are selecting these columns in later steps", header=3)
print_text("obtain the position of these columns using awk", header=4)
indexes_chrom_pos = run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --header | \
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

print_text("the chromosome name is correct? We do this check here because in the previous line we obtained the number of the column of CHROM", header=4)
print("calculate the number of times each chromosome appears")
chrom_count = run_bash(" \
    bcftools view \
        --drop-genotypes \
        --no-header \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{ \
            for(i=1;i<=NF;i++){ \
                if(i==" + index_chrom + "){ \
                    a[$i]++ \
                } \
            } \
        }END{for(i in a){print i,a[i]\n}}'", return_value=True).strip()
    #extract all rows of the VCF file without header and genotypes
    #load into awk using tabs as separator. In this way, we get separated the main columns, i.e., CHROM, POS...
        #iterate over i, from 1 to the number of columns, adding 1 to i in each iteration
            #if the number of the column is that of CHROM (previously calculated)
                #add an entry to array "a" using the chromosome name ($i) as index and adding 1 to the previous value
                #if the same chromosome appears 2 times, it will have a value of 2...
        #at the END, print each index of "a" (i; chromosome name) and its value (a[$i]; the count). Both values are separated by tab as OFS=\t, and each new pair is a new line (\n)
print("convert the output of awk into a pandas DF")
chrom_count_df = pd.DataFrame([chrom.split("\t") for chrom in chrom_count.split("\n")], columns=["chrom", "count"])
    #split each pair of values using "\n" and for each one
        #split each pair into the two value using \t
        #this is the input for pandas
    #select the column names
    #https://stackoverflow.com/a/54103026/12772630
print("check that we only have 1 chromosome and that is the selected chromosome")
if (chrom_count_df.shape[0]==1) & (chrom_count_df["chrom"].to_numpy()=="chr20"):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WITH THE CHROMOSOMES INCLUDED IN THE VCF FILE")

print_text("check that no SNP has filter different from 'PASS' or 'q10', which are the cases we have in our dummy example. In the case of the real data, FILTER should be '.'", header=4)
    #info for filter in 1KGP
        #according to the specific readme of the dataset we are using (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf), prior phasing they applied several filters, being one of them that all variants has to be PASS for FILTER. I understand that, because of this, all variants in the chromosome have now ".", being this the unique character.
print("calculate the number of SNPs for which FILTER is not 'PASS' or 'q10'")
problematic_fiter = run_bash(" \
    bcftools view \
        --no-header \
        --drop-genotypes \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{ \
            for(i=1;i<=NF;i++){ \
                if(i==" + index_filter + "){ \
                    if($i!=\"PASS\" && $i!=\"q10\"){count++} \
                } \
            } \
        }END{print count}'", return_value=True).strip()
    #load the VCF file after VEP without header and genotypes
    #open in awk using tabs as separator, so we separate the main columns, i.e., CHROM, POS, ID, REF, ALT....FILTER
        #iterate over i from 1 to the number of columns, adding 1 to "i" in each iteration
        #if the number of the column is that of FILTER (previously calculated)
            #if the value of that column ($i) is NOT ".", then add 1 to the array called count
        #END by printing the array count.
print("check that the number of SNPs with FILTER different from 'PASS' or 'q10' is 0")
if problematic_fiter=="":
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! WE HAVE SNPS FOR WHICH FILTER IS NOT '.'!!!")


print_text("switch REF/ALT columns for which REF is not AA", header=3)
print_text("see first cases where REF nor ALT are AA. You can see how we get a case with REF=G and AA_upcase=C. These could be cases where a multiallelic SNP has lost one of the ALT alleles in the subset population and that ALT allele is the ancestral. Therefore, there is no more ancestral in the population. Note, however, that I have found 200K cases like this across all chromosomes without subsetting when running 01b_vep_ancestral.py. Therefore, these could be caused by problems between SNPs according to VEP and out VCF files, so it should not be a high number", header=4)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --include 'REF!=AA_upcase && ALT!=AA_upcase' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA:%AA AA_upcase:%AA_upcase GTs:[ %GT]\n'")

print_text("check now what happens if we add a dummy SNP with AA=g and REF=G. As you can see, this SNP is considered to have AA different from REF and ALT because bcftools in case sensitive. Therefore, we need to use AA_upcase if we want to consider both high and low-confidence SNPs", header=4)
run_bash(" \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}END{ \
            print \"chr20\t1110691\trsdummy\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=g;AA_upcase=G\tGT\t1|0\t1|1\t0|0\" \
        }'\
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --include 'REF!=AA && ALT!=AA' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA:%AA AA_upcase:%AA_upcase GTs:[ %GT]\n'")
    #print all the rows of the dummy VCF file and, at the END, add a new row with a dummy snp having AA=g and REF=G. It has position inside the third interval of the dummy_mask so it passes all filters. The new row is printed using \t as separator between the VCF columns, i.e., REF, ALT, INFO....

print_text("count these cases but in the original dummy vcf file without filtering and just using AWK, so it is faster", header=4)
run_bash(" \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}END{ \
            print \"chr20\t1110691\trsdummy\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=g;AA_upcase=G\tGT\t1|0\t1|1\t0|0\" \
        }'\
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --types snps \
        --no-header \
        --drop-genotypes | \
    awk \
        'BEGIN{FS=\"\t|;|=|,\"}{ \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    aa_row=toupper($(i+1)); \
                    if(index(aa_row, \",\")==0 && length(aa_row)==1){ \
                        if(aa_row!=\"N\" && aa_row!=\"-\" && aa_row!=\".\"){ \
                            if(index($" + index_ref + ", aa_row)==0 && index($" + index_alt + ", aa_row)==0){ \
                               printf \"%s, REF=%s, ALT=%s, AA=%s\\n\", $3, $4, $5, $(i+1); \
                               count++ \
                            } \
                        } \
                    } else {exit 1} \
                } \
            } \
        }END{printf \"The number of these cases is %s\", count}'")
    #add again a line with the dummy SNP (REF=G and AA=g)
    #get all fields removing genotypes and header using bcftools
    #then in awk
        #use as separators "\t", ";", "=" and ",". In this way we separate not only the main fields (REF, ALT, INFO), but we also separate the INFO fields (AF, AA...) and their values (AA=g). The last delimiter is to get the first copy of the ancestral allele if the SNP have multiple consequences (AA=g,g,g). In this way, we get just one value that can be compared with REF, as REF has only one value.
        #run loop across each field
            #if we are in the field corresponding to AA
                #save the next field, i.e., the ancestral allele (first copy is multiple consequence) but in upper case, to consider both low and high-confidence alleles. This will be "aa_row"
                #if aa_row does not include ",", i.e., we do not have several ancestral alleles separated by comma and its number of characters is 1, i.e., we do not have more than 1 base (i.e., not TTG...)
                    #if aa_row is not missing (N, -, .). We can just do aa_row!= because we should not have more than 1 ancestral allele
                        #if the field 4 (REF) and the field 5 (ALT) do not include aa_row, 
                            #then print the row and add 1 to the count
                            #we use index() because if ALT has two alleles, we want to know if one of them is aa_row, i.e., whether aa_row is included in ALT.
                            #we do the same for REF, but it is not necessary because we should have only 1 REF allele always.
                #else, then we failed to extract just 1 ancestral allele per SNP from AA=G,G,G,G.... and we may have more than 1 ancestral allele (e.g., microsatellite TT...) So we need to stop and get an exist status of "1", which is considered as error by run_bash. If we do not add "1", the awk will end without raising an error, and then the python script will continue. We do not want that.
                    #https://unix.stackexchange.com/a/16567

print_text("exclude now cases where REF nor ALT are AA_upcase. This filter in cases where AA_upcase is the REF (rs6040358) or the ALT (rs6040356).", header=4)
print("Remember that one of the reasons for these cases to exist is that a multiallelic SNP lose one ALT when subsetting, and that very ALT is the ancestral. Therefore, we have to do the multiallelic-monomorphic filter before. Indeed, we are doing all filters before. As this filter acts only on REF, ALT and AA_upcase, not considering genotypes or samples, other filters should not affect. If one SNP is not phased or biallelic, then it will be removed before, because we do not want to do the switch for it")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --exclude 'REF!=AA_upcase && ALT!=AA_upcase' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA:%AA AA_upcase:%AA_upcase GTs:[ %GT]\n'")

print_text("create a new vcf file with a new case with AA=a and ALT=A, to check behavior with lower-case. Also add two cases with AA equals to ',' and '' to check behaviour", header=4)
run_bash(" \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}END{ \
            print \"chr20\t1110691\trsdummy\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=a;AA_upcase=A\tGT\t1|0\t1|1\t0|0\"; \
            print \"chr20\t1110692\trsdummy2\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=,;AA_upcase=,\tGT\t1|0\t1|1\t0|0\" \
        }'\
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf; \
        bcftools view \
            --no-header \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf")
    #print all lines of the VCF file
    #add at the end a line with the main columns (ALT, REF, CHROM, genotypes...) separated by tabs and indicate that with FS

print_text("Before the next step, check that REF and ALT are always ACGT, because our awk script relies on the fact that we discard those rows for which the unique ancestral allele in upper case is NOT equal to REF or ALT, removing in that way variants with ancestral allele equal to '.', 'N', '-', etc.... Also check that REF and ALT are NOT the same", header=4)
problematic_cases_ref_alt = run_bash(" \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}END{ \
            print \"chr20\t1110693\trsdummy3\t.\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=a;AA_upcase=A\tGT\t1|0\t1|1\t0|0\"; \
            print \"chr20\t1110694\trsdummy4\tG\t.\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=a;AA_upcase=A\tGT\t1|0\t1|1\t0|0\"; \
            print \"chr20\t1110704\trsdummy5\tG\t \t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=,;AA_upcase=A\tGT\t1|0\t1|1\t0|0\"; \
            print \"chr20\t1110705\trsdummy6\tG\tG\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=,;AA_upcase=A\tGT\t1|0\t1|1\t0|0\" \
        }'\
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools view \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t\"; \
            OFS=\"\t\"; \
            index_ref=" + index_ref + "; \
            index_alt=" + index_alt + "; \
        }{ \
            if($index_ref !~ /A|C|T|G/ || $index_alt !~ /A|C|T|G/ || $index_ref==$index_alt){ \
                count++ \
            } \
        }END{print count}'", return_value=True).strip()
if(problematic_cases_ref_alt == "4"):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM, REF OR ALT ARE NOT ALWAYS ACGT OR REF IS EQUAL TO ALT!!")
    #on the fly, modify the VCF file adding new dummy SNPs for which REF or ALT are not ACGT
    #then split multiallelics and remove the header
    #awk
        #begin with tabs as delimiter and create variables with the index of REF and ALT columns
        #if REF or ALT does NOT include ACGT, count
            #use regex expression to have multiple conditions using "|"
            #negate with "!"
            #https://stackoverflow.com/a/8481180/12772630
        #print count
    #check that the count is exactly 3, because we have added 3 more dummy SNPs for which REF or ALT does not include ACGT

print_text("create a awk script that switch REF/ALT in those rows where AA==ALT, while not doing anything to rows where REF==AA and discarding those rows where REF nor ALT are equal to AA. First apply it to all variants without filter to check behaviour", header=4)
selected_chrom_dummy="chr20"
run_bash(" \
    last_header_row=$( \
        bcftools norm \
            --multiallelic -snps \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf | \
        awk \
            '{ \
                if($0 ~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/){print NR} \
            }'); \
    if [[ $last_header_row == '' ]]; then \
        exit 1; \
    fi; \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf | \
    awk \
        -v last_header_row=\"$last_header_row\" \
        -v date=\"$(date '+%a %b %d %H:%M:%S %Y')\" \
        'BEGIN{ \
            FS=\"\t\"; \
            OFS=\"\t\"; \
            index_ref=" + index_ref + "; \
            index_alt=" + index_alt + "; \
            index_info=" + index_info + "; \
            selected_chrom_dummy=\"" + selected_chrom_dummy + "\"; \
        }{ \
            if(NR > last_header_row){ \
                if(NR==(last_header_row+1)){ \
                    if($1 != selected_chrom_dummy){exit 1} \
                } \
                for(i=1;i<=length($index_info);i++){ \
                    if(substr($index_info,i,3)==\"AA=\"){ \
                        k=i; \
                        while(substr($index_info, k, 1)!=\";\"){ \
                            end_aa=k; \
                            k=k+1; \
                        }; \
                        if(i+3 != end_aa){ \
                            full_aa = substr($index_info, i+3, end_aa-(i+3)+1) \
                        } else { \
                            full_aa = substr($index_info, i+3, 1) \
                        }; \
                        if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                            full_aa_no_comma = full_aa; \
                            gsub(/,/, \"\", full_aa_no_comma); \
                            count=0; \
                            for(k=1;k<=length(full_aa_no_comma);k++){ \
                                char=substr(full_aa_no_comma, k, 1); \
                                if(!array_alleles[char]++){count++}; \
                            }; \
                            if(count==1){ \
                                check_n_csq=\"false\"; \
                                for(k in array_alleles){ \
                                    if(array_alleles[k]==length(full_aa_no_comma)){ \
                                        check_n_csq=\"true\" \
                                    } \
                                }; \
                                if(check_n_csq==\"true\"){ \
                                    k=1; \
                                    while(substr(full_aa, k, 1)!=\",\"){ \
                                       before_comma=k; \
                                       k=k+1\
                                    }; \
                                    unique_aa = substr(full_aa, 1, before_comma); \
                                } else { \
                                    exit 1 \
                                }; \
                            } else { \
                                exit 1 \
                            }; \
                            delete array_alleles; \
                        } else { \
                            unique_aa=full_aa \
                        }; \
                        unique_aa = toupper(unique_aa); \
                    } else { \
                        if(substr($index_info,i,3)==\"AA;\"){exit 1} \
                    } \
                }; \
                if($index_ref==unique_aa || $index_alt==unique_aa){ \
                    if($index_ref!=unique_aa && $index_alt==unique_aa){ \
                        tmp_ref=$index_ref; \
                        $index_ref=$index_alt; \
                        $index_alt=tmp_ref; \
                        print $0 \
                    } else { \
                        print $0 \
                    }; \
                    if($index_ref!=unique_aa && $index_alt==unique_aa){exit 1} \
                }; \
            } else { \
                if(NR==last_header_row){ \
                    printf \"##awk_script: Switched REF/ALT columns for variants where REF!=AA but ALT==AA. Variants where REF!=AA and ALT!=AA where discarded; Date=%s\\n\", date; \
                    print $0; \
                }else{print $0} \
            } \
        }'; \
    check_status=$(echo $?); \
    if [[ $check_status -ne 0 ]]; then \
        echo 'ERROR! FALSE! WHEN SWITCHING REF/ALT COLUMNS'; \
    fi")
    #load the VCF file with only the header into awk, then extract the number of the row (whole row; $0) starting with "CHROM REF...."
        #https://unix.stackexchange.com/a/72763
    #check that the number of rows of the header is NOT empty, meaning that we have correctly detected the header. If it is empty, exit with non-zero exit status.
    #awk
        #load variables 
            #the number of rows of the header as a variable
            #the date following the format used by bcftools
                #The command is "+FORMAT", where FORMAT controls the output. For example: "%H" gives hours.
                #In order to have spaces between the units, you need to use '' around FORMAT.
                #https://unix.stackexchange.com/a/224976
                #https://stackoverflow.com/a/27337807/12772630
                #https://man7.org/linux/man-pages/man1/date.1.html
        #begin awk
            #loading the data into awk using \t for input and output. In this way, we are going to work with the main columns of the VCF file, REF, ALT, INFO, FORMAT, GT...
            #This means that the whole INFO field is going to be a string for each row, i.e., "AA=g;AA_upcase=G...." because we are not using ";" nor "=" as delimiters.
            #we do this to get as output the original INFO field with the ";",  "=", because if you use these as delimiters, they disappear from the data.
            #create variables with the position of the columns of interest that we will use later
        #if the row is after the header
            #if the row is the first just after the header
                #check if the first column (CHROM) is not "chr20", if that is the case, exit with non-zero exit status
                #https://stackoverflow.com/a/3701095/12772630
            #for each character of the INFO field
                #use substring to extract specific characters from the INFO string. 
                    #substr(string, start [, length ])
                    #Return a length-character-long substring of string, starting at character number start. The first character of a string is character number one. For example, substr("washington", 5, 3) returns "ing".
                        #https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html
                #if a substring starting at i up to two positions after (i.e., length=3) is equals to "AA=" (this avoids "AA_upcase=")
                    #while loop: iterate from the starting position of "AA=", i.e., "A". You save the position of the character "end_aa", overwritting the position of the previous character until we reach ";", which is NOT saved. 
                        #this will give the position of the last character before ";", i.e., before the next INFO field.
                        #https://opensource.com/article/19/11/loops-awk
                    #Use this position to select only the ancestral allele data wihthout "AA=" and ";". If the last character before ";" is the position just after "AA=", i.e., "AA" is equal to a single character (e.g., AA=T or AA=-). In other words, three positions after "AA=" starts (i+3) is equals to the last position before ";" (end_aa)
                        #select just the character in that position
                    #else, means that the we have several characters within "AA="
                        #select a substring starting after "AA=" and ending just before ";"
                        #we obtain that position by substracting the start of "AA=" from "end_aa"
                        #For example: In "AA=G,G,G;", "AA=" starts at 1 (i=1), while the first base is in 1+3=4 (i+3=4), and the last character before ";" is at position 8 (end_aa=8). The length of the substring is 8-4+1=5 (end_aa-(i+3)+1). Therefore, if from the position of the last base (end_aa=8), we subtract the position of the first base just after "AA=" (i+3=4), we get the distance between the two characters without including one of them, so you have to add 1 to include botch extremes of this interval: 8-(1+3)+1=5. We want a substring of length 5 to include all alleles.
                            #Note that, in substr, the start is included so you have to start at 4, including 4, and then add 5 positions starting from 4 to reach the final base. First G is at 4, last G is at 8, so 8-4+1=5. This is 5,6,7 plus the extremes, 4 and 8.
                    #if "full_aa" does not include "," or it is exactly "," but it has no other character, this means that "AA" is just a base or a missing including "," (",", "N", "-", ".").
                        #then you can just save "full_aa" in "unique_aa" because there is only 1 character, so it is already unique.
                        #Note that we use bcftools to extract AA from the CSQ field, and I checked that this makes cases with empty "AA" (AA="") as "AA=.". So we will always have a character after AA, and no space.
                            #We have added other "if" to check we have no empty AA field just in case
                    #if it includes "," but it is not exactly ",", this means that we have several alleles for this SNP, e.g., "A,A,A,A", so we need to get just one. 
                        #first check again (we did before) that we only have one ancestral allele per variant
                            #copy "full_aa" into "full_aa_no_comma"
                            #using gsub, change every "," by empty (""), so now we should have only AA data and no commas, the same allele repeated many times
                                #https://unix.stackexchange.com/a/492502
                            #set a "count" variable as zero
                            #for each character in full_aa_no_comma
                                #extract that character
                                #save it in array "a", 
                                    #if it is previously present then add 1 more to its count. 
                                    #if is not present before, add it as a new entry with value 1.
                                        #Negate the expression so this case is true and then we can add 1 to our "count" variable.
                                        #In other words, count the number of unique characters in full_aa_no_comma
                        #We do operations if count==1, if it is higher than 1 and hence we have more than 1 distinct characters, STOP (non-zero exit status). 
                            #this check will fail for microsatellites with different bases like TAT, but this should be ok because we should have only SNPs. If we have different bases in a variant, I want an error.
                            #we do operations if another check works, if non, stop with non-zero exit status. Take the number of counts we have in array "a", which stores the unique characters seens in "full_aa_no_comma" with its count. This number should be the same than the total length of "full_aa_no_comma", i.e., the number of characters without comma should be the same than the count for the single character we got. Remember that the previous if check if 'count'==1, so we already know we have only 1 different character.
                            #If true:
                                #iterate across characters of "full_aa"
                                    #while the character is NOT ",", save the position of the character in "before_comma" and move to the next position
                                    #if "," it reached, then "before_comma" is not updated more, and hence we have the last position just before ","
                                #now extract a substring from "full_aa" from the start to the position just before the first comma. For example, in TT,TT,TT: would be from to 2, because the first comma is in position 3.
                        #delete "a", so we have it clean for the next iteration of the larger loop (iterated over i)
                    #else, means that we only have 1 character after "AA=" and before ";", therefore just take that character.
                    #convert the resulting unique ancestral allele (unique_aa) into uppercase to consdier both low and high-confidence alleles.
                #else means that we do not have "AA=". This is ok because we can have other INFO fields, BUT
                    #if the INFO field is "AA;" STOP, because we should have always at least "." in case AA is missing.
            #if the REF or the ALT fields are equal to the unique ancestral allele
                #if the REF is NOT the AA, but the ALT it is, we need to switch them
                    #save the REF in a temp variable
                    #overwrite the original REF with the ALT
                    #save as new ALT the REF value stored in the temp variable
                #else just print the row without changes
            #else, we do not have ancestral allele for this variant in the panel, so we can discard the SNP, so no print. Without ancestral allele there is anything we can do with that allele.
        #else means the row belongs to the header
            #if we are in the last row of the header
                #before printing that row, add a comment with a explanation of the awk script and the time, ending with new line (\n)
            #else, meaning we are in the previous header lines, just print these lines
    #save the exist status after awk and check it. If is it not zero, we have a problem.

print_text("save the subset of filtered variants and then apply on it the awk script", header=4)
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools view \
        --targets-file ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned.vcf")

print_text("apply the awk script to the cleaned VCF", header=4)
run_bash(" \
    last_header_row=$( \
        bcftools view \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned.vcf | \
        awk \
            '{ \
                if($0 ~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/){print NR} \
            }'); \
    if [[ $last_header_row == '' ]]; then \
        exit 1; \
    fi; \
    bcftools view \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned.vcf | \
    awk \
        -v last_header_row=\"$last_header_row\" \
        -v date=\"$(date '+%a %b %d %H:%M:%S %Y')\" \
        'BEGIN{ \
            FS=\"\t\"; \
            OFS=\"\t\"; \
            index_ref=" + index_ref + "; \
            index_alt=" + index_alt + "; \
            index_info=" + index_info + "; \
            selected_chrom_dummy=\"" + selected_chrom_dummy + "\"; \
        }{ \
            if(NR > last_header_row){ \
                if(NR==(last_header_row+1)){ \
                    if($1 != selected_chrom_dummy){exit 1} \
                } \
                for(i=1;i<=length($index_info);i++){ \
                    if(substr($index_info,i,3)==\"AA=\"){ \
                        k=i; \
                        while(substr($index_info, k, 1)!=\";\"){ \
                            end_aa=k; \
                            k=k+1; \
                        }; \
                        if(i+3 != end_aa){ \
                            full_aa = substr($index_info, i+3, end_aa-(i+3)+1) \
                        } else { \
                            full_aa = substr($index_info, i+3, 1) \
                        }; \
                        if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                            full_aa_no_comma = full_aa; \
                            gsub(/,/, \"\", full_aa_no_comma); \
                            count=0; \
                            for(k=1;k<=length(full_aa_no_comma);k++){ \
                                char=substr(full_aa_no_comma, k, 1); \
                                if(!array_alleles[char]++){count++}; \
                            }; \
                            if(count==1){ \
                                check_n_csq=\"false\"; \
                                for(k in array_alleles){ \
                                    if(array_alleles[k]==length(full_aa_no_comma)){ \
                                        check_n_csq=\"true\" \
                                    } \
                                }; \
                                if(check_n_csq==\"true\"){ \
                                    k=1; \
                                    while(substr(full_aa, k, 1)!=\",\"){ \
                                       before_comma=k; \
                                       k=k+1\
                                    }; \
                                    unique_aa = substr(full_aa, 1, before_comma); \
                                } else { \
                                    exit 1 \
                                }; \
                            } else { \
                                exit 1 \
                            }; \
                            delete array_alleles; \
                        } else { \
                            unique_aa=full_aa \
                        }; \
                        unique_aa = toupper(unique_aa); \
                    } else { \
                        if(substr($index_info,i,3)==\"AA;\"){exit 1} \
                    } \
                }; \
                if($index_ref==unique_aa || $index_alt==unique_aa){ \
                    if($index_ref!=unique_aa && $index_alt==unique_aa){ \
                        tmp_ref=$index_ref; \
                        $index_ref=$index_alt; \
                        $index_alt=tmp_ref; \
                        print $0 \
                    } else { \
                        print $0 \
                    }; \
                    if($index_ref!=unique_aa && $index_alt==unique_aa){exit 1} \
                }; \
            } else { \
                if(NR==last_header_row){ \
                    printf \"##awk_script: Switched REF/ALT columns for variants where REF!=AA but ALT==AA. Variants where REF!=AA and ALT!=AA where discarded; Date=%s\\n\", date; \
                    print $0; \
                }else{print $0} \
            } \
        }' > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf")

print_text("After switching, REF is always equals to AA_upcase while ALT is not?", header=4)
count_problem_ref_aa = run_bash(" \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"; \
            index_ref=" + index_ref + "; \
            index_alt=" + index_alt + " \
        }{ \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA_upcase\"){ \
                    full_aa=$(i+1); \
                    if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                        k=1; \
                        while(substr(full_aa, k, 1) != \",\"){ \
                            before_comma=k; \
                            k=k+1; \
                        }; \
                        unique_aa=substr(full_aa, 1, before_comma); \
                    } else { \
                        unique_aa=full_aa \
                    }; \
                    if(unique_aa!=$index_ref || unique_aa==$index_alt){ \
                        count++ \
                    } \
                } \
            } \
        }END{print count}' \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf", return_value=True).strip()
if(count_problem_ref_aa == ""):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: AFTER SWITCHING, REF IS NOT ALWAYS AA_upcase OR ALT IS EQUAL TO AA_upcase")
    #this is a very good check because we have previously done the switch without using AA_upcase, and now we check  if REF is always equals to AA_upcase, which was obtained using other approach.

print_text("We get the same rows if we take the original VCF and select only rows where AA is REF or ALT and removing REF/ALT columns (which where the ones switched)?", header=4)
run_bash(" \
    bcftools view \
        --no-header \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned.vcf | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"; \
            index_ref=" + index_ref + "; \
            index_alt=" + index_alt + " \
        }{ \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    full_aa=$(i+1); \
                    if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                        k=1; \
                        while(substr(full_aa, k, 1) != \",\"){ \
                            before_comma=k; \
                            k=k+1; \
                        }; \
                        unique_aa=substr(full_aa, 1, before_comma); \
                    } else { \
                        unique_aa=full_aa \
                    }; \
                    unique_aa = toupper(unique_aa);\
                    if(unique_aa==$index_ref || unique_aa==$index_alt){ \
                        print $0 \
                    } \
                } \
            } \
        }' | \
    cut \
        --delimiter '\t' \
        --fields 1-3,6- > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_1.vcf; \
    bcftools view \
        --no-header \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf | \
    cut \
        --delimiter '\t' \
        --fields 1-3,6- > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_2.vcf; \
    check_status=$( \
        cmp \
            --silent \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_1.vcf \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_2.vcf; \
        echo $?); \
    if [[ $check_status -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi; \
    rm ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_1.vcf; \
    rm ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/test_file_2.vcf")
    #create a file with the VCF file before switching but selecting only those rows for which REF or ALT are AA.
        #you can check previous AWK script for information
        #just note you are assuming here that the position of REF and ALT columns is that of "index_ref" and "index_alt".
            #this is ok because despite using several delimiters besides tabs, because these columns are before the INFO field, where we have ";", "=" and hence more columns are created than just considering tabs.
            #in addition, if this is problematic, it will very likely lead to a different file generating false at the final comparison
        #after awk, select all columns except the 4 and 5 with cut
            #again, we are assuming REF and ALT are 4 and 5, if this is not the case, we will get a different test_1 file respect to test_2 file, and hence false in the final check
            #https://stackoverflow.com/a/13446273/12772630
    #get the switched VCF file but removing REF and ALT columns
    #check byte by byte that both files are the same with cmp, is equal, the exit status should be zero.


print_text("see the new header", header=3)
run_bash(" \
    bcftools view \
        --header \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf")
        #the new header shows 
            #the fields of the previous VCF file that remain.
            #every bcftools command used in order, along with the flags, the version, the date.
            #the new fields generated with +fill-tags
            #version of +fill-tags, flags selected
        #Therefore, the new header includes a history of the changes made to the VCF file inside of the file itself.


print_text("see the variants", header=3)
run_bash(" \
    bcftools view \
        --no-header \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf")


print_text("calculate stats of the VCF file, show them here and then use them to make summary plots", header=3)
run_bash(" \
    bcftools stats \
        --samples - \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf")
        #Produce stats of a VCF file using bcftool stats
            #Parses VCF or BCF and produces stats which can be plotted using plot-vcfstats. 
                #When two files are given, the program generates separate stats for intersection and the complements. 
                #By default only sites are compared, -s/-S must given to include also sample columns.
                    #"-" to include all samples
                    #this gives plots showing data with a datapoint per sample
        #You can make plots of the output with plot-vcfstats
            #just plot-vcfstats --prefix for the dir you want to save the plots and the input file (bcftools stats output with .vchk extension). I am not doing this because I have problems with plot-vcfstats in the container
                #--prefix <dir>: Output directory.
                #there are several plotting options, like title...
        #The final looks can be customized by editing the generated outdir/plot.py' script and re-running manually
            #cd outdir && python3 plot.py && pdflatex summary.tex


print_text("convert VCF file to hap file", header=3)
print("First convert the raw dummy VCF file but adding a filter for SNPs in the same line, i.e. removing the microsatellite. The output shows 0/5/1 no-ALT/non-biallelic/filtered indicating that we lose 5 variants that are not biallelic and 1 due to the filter we applied, i.e., removal of non-snps, the microsatellite. Also, the number of records written is 6, which is the result of subtracting 5+1 from 12, the total number of variants. We have included missing alleles, which are shown as '?', while non-phased alleles are shown as '*'")
run_bash(" \
    bcftools convert \
        --include 'TYPE == \"SNP\"' \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf \
        --hapsample ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_IMPUTE2_raw; \
    gunzip -c ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_IMPUTE2_raw.hap.gz")
print("Then convert the cleaned dummy VCF file. The output shows 0/0/0 because no variant has been removed, they were filtered in previous lines. The number of written records is 3. The input had 3 and no filter was applied, nor variant lose due to lack of alt or non-biallelic, thus we save a hap file with 3 variants")
run_bash(" \
    bcftools convert \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf \
        --hapsample ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_IMPUTE2; \
    gunzip -c ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_IMPUTE2.hap.gz")
        #see hap file calculation of the real data to see details about hap format.

#SUMMARY: 
    #With all these commands, we have recreated the scenario we have in 1KGP data, with multiallelic SNPs separated into different lines, select some samples, we then select snps, exclude those SNPs that have the same allele for all samples (considering alleles in genotypes with missing, e.g., 0|.), remove exact duplicates (this does not touch different lines of the same multiallelic snp because they have different ALT). Then we combine all lines of each multiallelic snp and now they have ALT column with several alleles, so we can filter them using --max-alleles 2. Add filter for selecting phased data only. Select only those variants included in interest regions (mask). We also use bcftools +fill-tags to update important fields for each SNP, so if a SNP was multiallelic, but it is not multiallelic in the subset population (i.e., only REF and 1 ALT), we no longer will have two allele frequencies, two allele counts.... for the remainder biallelic SNP in the subset. Although in 1KGP data, multiallelic SNPs are already separated and have only 1 value for these fields (see below). We have also switched REF/ALT columns for those SNPs whose ancestral allele was not REF, also removed those SNPs for which REF nor ALT are the ancestral allele.

#Note about the update of the INFO fields
    #it is important to be sure that the fields you are using for filtering, are updated after subseting samples. Of course, type="snp" will be always "snp" irrespectively of the samples we select, but this is not the case of the number of alleles, because you can have SNPs with 3 alleles considering all 26 populations, but then in GBR they can have only 2 or 1. We are interested in SNPs that are biallelic within the selected population.
    #The same goes for phasing and genotype ^miss. You have to be sure that these filters only consider data from the filtered samples, not fixed data in fields that are not updated.
    #Because of this we have applied all these filters in order.




################################################################
#### function to clean vcf files and create hap - map files ####
################################################################
print_text("function to clean vcf files and create hap - map files", header=1)
import numpy as np
#chr_pop_combination="IBS_1"; debugging=True; debug_file_size=50000
def master_processor(chr_pop_combination, debugging=False, debug_file_size=None):

    #extract selected population and chromosome
    selected_pop = chr_pop_combination.split("_")[0]
    selected_chromosome = chr_pop_combination.split("_")[1]



    print_text("Initial operations", header=2)
    print_text("create folder for the selected chrom and pop", header=3)
    run_bash(" \
        mkdir \
            -p \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "; \
        mkdir \
            -p \
            ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "; \
        mkdir \
            -p \
            ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "; \
        mkdir \
            -p \
            ./scripts/01_hap_map_calcs_outputs/" + selected_pop)


    print_text("redirect standard output", header=3)
    if (debugging==False):
        import sys
        original_stdout = sys.stdout
            #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
        file_output = open("./scripts/01_hap_map_calcs_outputs/" + selected_pop + "/chr" + selected_chromosome + "_" + selected_pop + ".out", "w")
            #open a file where stdout will be saved
        sys.stdout = file_output
            #redirect stdout to that file
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/
            #https://stackoverflow.com/a/23838153


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see VCF file version", header=3)
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
                #GT : genotype, encoded as allele values separated by either of / or | (UNPHASED AND PHASED RESPECTIVELY). The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc. For haploid calls, e.g. on Y, male nonpseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, ‘.’ should be specified for each missing allele in the GT field (for example ‘./.’ for a diploid genotype and ‘.’ for haploid genotype). The meanings of the separators are as follows (see the PS field below for more details on incorporating phasing information into the genotypes):
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
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA %AA_upcase %INFO/AF\n' \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
            head -1")
        #select the format of the query indicating the columns you want to show per SNP.
            #you can include data from INFO
            #end with \n to have different lines per SNPs


    print_text("chr " + selected_chromosome + ": see genotypes of first samples", header=3)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -2")


    print_text("extract the position (index) of several columns in the VCF file, so we can be sure we are selecting these columns in later steps", header=3)
    print_text("obtain the position of these columns using awk", header=4)
    indexes_chrom_pos = run_bash(" \
        bcftools view \
            --header \
            " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
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


    print_text("extract samples IDs", header=3)
    print_text("select the sample IDs for the selected population", header=4)
    subset_pop = unrelated_samples.loc[unrelated_samples["pop"] == selected_pop, :]
    subset_pop = subset_pop.reset_index(drop=True)
        #reset index
        #drop=True to avoid adding the index as a column

    print_text("chr " + selected_chromosome + " - " + selected_pop + ": check we only have the selected pop", header=4)
    print(subset_pop["pop"].unique() == selected_pop)
    if subset_pop["pop"].unique() != selected_pop:
        raise ValueError("SERIOUS ERROR! WE HAVE NOT CORRECTLY SELECTED THE SAMPLES OF POP '" + selected_pop + "' in chromosome '" + selected_chromosome + "'")

    print_text("select the sample IDs and save as a file", header=4)
    selected_samples = subset_pop["sample"]
    selected_samples.to_csv(
        "results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt", 
        sep="\t", 
        header=None, 
        index=False)
        #it is redundant to do it in each chromosome of the same population (same samples) but I do it anyway to avoid confusion between chromosomes


    print_text("Select the input for the next steps. Use only a subset of the samples. Also select only a subset of SNPs if we are on debugging mode", header=3)
    if debugging==True:
        run_bash(" \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
            head -n " + str(debug_file_size) + " > ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/01_debug_subset_chr" + selected_chromosome + "_" + selected_pop + ".vcf; \
            ls -l")
        input_vcf_file="./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/01_debug_subset_chr" + selected_chromosome + "_" + selected_pop + ".vcf"
    else:
        run_bash(" \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                --output-type z \
                --output ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".vcf.gz \
                " + input_vcfs_path + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz")
        input_vcf_file="./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".vcf.gz"
    print(input_vcf_file)
    #select a few samples and then select biallelic snps
        #IMPORTANT: If filtering (by number of alleles) and subseting (by sample) is in the same line (command), the filtering will be done first. Therefore, you could select SNPs that have 2 alleles when considering the 26 pops, but that are monomorphic (1 allele) for the selected population. Because of this, we have to first subset by sample and then filter by number of alleles whitin the selected samples in separated commands (see dummy example).
    #query multiple fields for each snp


    print_text("check that the samples we have in the VCF subset are the correct ones", header=3)
    print(run_bash(" \
        bcftools query \
            --list-samples \
            " + input_vcf_file, return_value=True).strip().split("\n") == \
    selected_samples.to_list())
        #extract the samples from the new VCF file and then compare them with the list of selected samples previously obtained in python.
            #the bcftools output needs to be split with "\n" in order to separate the IDs of the different samples in a list.


    print_text("check we have only SNPs", header=3)
    print("count the number of rows of the VCF when we exclude SNPs")
    n_rows_no_snps = run_bash(" \
        bcftools view \
            --no-header \
            --drop-genotypes \
            --exclude-types snps \
            " + input_vcf_file + " | \
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


    print_text("chr " + selected_chromosome + ": the chromosome name of is correct?", header=3)
    count_problems_chrom_name = run_bash(" \
        bcftools view \
            --no-header \
            --drop-genotypes \
            " + input_vcf_file + " | \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\"; \
                selected_chromosome=\"chr" + selected_chromosome + "\"; \
                chrom_index=" + index_chrom + "; \
            }{ \
                if($chrom_index!=selected_chromosome){count++} \
            }END{print count}'", return_value=True).strip()
    if(count_problems_chrom_name==""):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! We have a problem with the chromosome names in the VCF file")
        #load the VCF file without header
        #awk
            #load using tabs, so the we have separated CHROM, REF, ALT....
            #save two variables with the name of the chromosome and the position of the column of the chromosome name in the VCF file
            #if the chromosome name is not the selected one, add 1 to count
        #extract count and if it is not empty, then raise error!


    print_text("check that no SNP has filter different from '.' We do this check here because in the previous line we obtained the number of the column of FILTER", header=3)
        #according to the specific readme of the dataset we are using (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf), prior phasing they applied several filters, being one of them that all variants has to be PASS for FILTER. I understand that, because of this, all variants in the chromosome have now ".", being this the unique character.
    print("calculate the number of SNPs for which FILTER is not '.'")
    problematic_fiter = run_bash(" \
        bcftools view \
            --no-header \
            --drop-genotypes \
            " + input_vcf_file + " | \
        awk \
            'BEGIN{FS=OFS=\"\t\"; index_filter=" + index_filter + "}{ \
                if($index_filter!=\".\"){count++} \
            }END{print count}'", return_value=True).strip()
    print("check that the number of SNPs with FILTER different from '.' is 0")
    if problematic_fiter=="":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE SNPS FOR WHICH FILTER IS NOT '.'!!!")



    print_text("calculate hap file within selected pop", header=2)
    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see more SNPs and their allele count to see how this data is presented for multiallelic SNPs", header=3)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA %AA_upcase %INFO/AN %INFO/AC %INFO/AF GTs:[ %GT]\n' \
            " + input_vcf_file + " \
            | \
        head -5")
        #IMPORTANT:
            #It seems that multiallelic variants separated in different lines have already updated the AC field. Therefore, each line does not have two allele counts, but one.
            #For example, 1:10453:A:C and 1:10452:A:C are both in the same position (10452) and have the same reference allele. They seem to be multiallelic, but each one has only one allele count, which is 1. 
            #This count is correctly updated for the subset of samples, but remember that --samples does not remove the second count. As we saw in the dummy example, we have to update the field with +fill-tags.
            #My hypothesis is that they updated this field using +fill-tags because they indeed say in the paper that they used bcftools to split the multiallelic SNPs in different lines.
            #Note that AF has also 1 value but it is not correct because we subset samples with --samples, and this command only updates AC and AN.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see monomorphic. Derive allele count should be 0 or equal to the total number of alleles, while AF should be 0 or 1", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -20")
            #update the INFO/AC and INFO/AN fields so we avoid two allele counts in the different lines of multiallelic snps. --samples update AC but maintaining the connection between the different lines of the same multialllelic SNP.
                #it seems that 1KGP data has already updated the AC field for each line separated line of a multiallelic SNP.
                #Therefore, this line would not be necessary, but we are applying just in case, to ensure we have only 1 allele count per line.
                #updating again the AC fields would not do anything wrong, just adding the same value that is was.
            #we also have to update AF if we want to have the actual frequencies after subseting.
            #select those variants for which the number of ALT alleles (allele count) is equal to 0 (no ALT at all) or equal to the total number of alleles (AN), i.e., all alleles are ALT.
            #See dummy example for further details.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": exclude monomorphic. AF has to be >0 and <1", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -20")


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": check we do not have SNPs with genotype missingness > 0.05. If TRUE, then we do not have to apply further filters about missing genotypes", header=3)
    check_missing = run_bash(" \
        n_snps_above_5=$( \
            bcftools view \
                " + input_vcf_file + " | \
            bcftools +fill-tags \
                -- --tags AN,AC,AF | \
            bcftools view \
                --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
            bcftools view \
                --no-header \
                --drop-genotypes \
                --include 'COUNT(GT=\"mis\")/N_SAMPLES >= 0.05' | \
            wc -l); \
        if [[ $n_snps_above_5 -eq 0 ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi", return_value=True)
            #obtain the ID of all variants after applying or not the filter for genotype missingness < 0.05, then check the numbers are the same
            #return the value to do an explicit check
    #if the check is not TRUE, stop because this is important
    if check_missing.strip() != "TRUE":
        #use strip() to remove "\n" at the end of the string
        raise ValueError("ERROR: FALSE! The filter of < 5% of missingness does not seem to be applied in chr '" + selected_chromosome + "' for pop '" + selected_pop + "'")


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": remove also exact duplicates", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -10")
            #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. See dummy example for behaviour.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": combine multiallelic SNPs in one line and select them", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --min-alleles 3 | \
        bcftools +fill-tags \
            -- --tags AC,AF | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -7")
            #we have SNPs with the same position and chromosome but different ALT alleles. According to the 1KGP paper, they split multiallelic SNPs in different lines, so this makes sense. See dummy for further details.
            #combine the different lines of each multiallelic SNP into one line and update the ALT column to include the different ALT alleles and we can filter with --max-alleles
                #two snps in position 10452 with same REF but different ALT in chromosome 1. They get combined into one line. These are the first multiallelic that appear in the previous command, but they were separated.
            #Also update again AC and AF, so we can see the AC and AF fields with the data for each derived allele


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools +fill-tags \
            -- --tags AC,AF | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -7")


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs that are phased for all samples", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools +fill-tags \
            -- --tags AC,AF | \
        bcftools view \
            --phased | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC AF:%AF GTs:[ %GT]\n' | \
        head -7")


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": update and add some fields using fill-tags", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AN:%AN AC:%AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n' | \
        head -7")
            #we use +fill-tags to update fields that are not updated after subsetting like frequency of alternative allele and create some additional fields.
            #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed sue to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields. In addition, --samples maintains 2 allele counts in each line of a splitted multiallelic SNP. In  the dummy example we had to update these fields with +fill-tags to have only 1 count and be able to use this count to filter out monomorphic SNPs. I have checked that splitted multiallelic SNPs in the 1KGP have only 1 count, so I guess the authors updated with fill-tags, but I update before the monomorphc check just in case. See above.
            #see dummy example for further details.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": remove all previous INFO/Fields, retain only what we are interested in and show the result asking for ALL fields", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --no-header | \
        head -7")
            #remove all INFO fields (except the ancestral allele) and all FORMAT fields (except GT) using annotate and then add the fields we are interested in.
            #see dummy example for further details.
            #use view with option --no-header to see all fields without the header.


    print_text("switch REF/ALT columns for which REF is not AA", header=3)
    print_text("see first cases where REF nor ALT are AA. These could be cases where a multiallelic SNP has lost one of the ALT alleles in the subset population and that ALT allele is the ancestral. Therefore, there is no more ancestral in the population. Note, however, that I have found 200K cases like this across all chromosomes without subsetting when running 01b_vep_ancestral.py. Therefore, these could be caused by problems between SNPs according to VEP and out VCF files, so it should not be a high number", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --include 'REF!=AA_upcase && ALT!=AA_upcase' | \
        bcftools query \
            -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA:%AA AA_upcase:%AA_upcase\n' | \
        head -n 10")


    print_text("count these cases", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --no-header \
            --drop-genotypes | \
        awk \
            'BEGIN{ \
                FS=\"\t|;|=|,\"; \
                index_id=" + index_id + "; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + "; \
            }{ \
                for(i=1;i<=NF;i++){ \
                    if($i==\"AA\"){ \
                        aa_row=toupper($(i+1)); \
                        if(index(aa_row, \",\")==0 && length(aa_row)==1){ \
                            if(aa_row!=\"N\" && aa_row!=\"-\" && aa_row!=\".\"){ \
                                if(index($index_ref, aa_row)==0 && index($index_alt, aa_row)==0){ \
                                    if(NR<10000){ \
                                        printf \"ID=%s, REF=%s, ALT=%s, AA=%s\\n\", $index_id, $index_ref, $index_alt, $(i+1); \
                                    }; \
                                    count++ \
                                } \
                            } \
                        } else {exit 1} \
                    } \
                } \
            }END{printf \"The number of these cases is %s\", count}'")
        #add again a line with the dummy SNP (REF=G and AA=g)
        #get all fields removing genotypes and header using bcftools
        #then in awk
            #use as separators "\t", ";", "=" and ",". In this way we separate not only the main fields (REF, ALT, INFO), but we also separate the INFO fields (AF, AA...) and their values (AA=g). The last delimiter is to get the first copy of the ancestral allele if the SNP have multiple consequences (AA=g,g,g). In this way, we get just one value that can be compared with REF, as REF has only one value.
            #add the position of relevant columns as variables
            #run loop across each field
                #if we are in the field corresponding to AA
                    #save the next field, i.e., the ancestral allele (first copy is multiple consequence) but in upper case, to consider both low and high-confidence alleles. This will be "aa_row"
                    #if aa_row does not include ",", i.e., we do not have several ancestral alleles separated by comma and its number of characters is 1, i.e., we do not have more than 1 base (i.e., not TTG...)
                        #if aa_row is not missing (N, -, .). We can just do aa_row!= because we should not have more than 1 ancestral allele
                            #if the field 4 (REF) and the field 5 (ALT) do not include aa_row,
                                #add 1 to the count
                                #if this is one of the first 10K rows
                                    #then print the row
                                #we use index() because if ALT has two alleles, we want to know if one of them is aa_row, i.e., whether aa_row is included in ALT.
                                #we do the same for REF, but it is not necessary because we should have only 1 REF allele always.
                    #else, then we failed to extract just 1 ancestral allele per SNP from AA=G,G,G,G.... and we may have more than 1 ancestral allele (e.g., microsatellite TT...) So we need to stop and get an exist status of "1", which is considered as error by run_bash. If we do not add "1", the awk will end without raising an error, and then the python script will continue. We do not want that.
                        #https://unix.stackexchange.com/a/16567


    print_text("exclude now cases where REF nor ALT are AA_upcase.", header=3)
    print("Remember that one of the reasons for these cases to exist is that a multiallelic SNP lose one ALT when subsetting, and that very ALT is the ancestral. Therefore, we have to do the multiallelic-monomorphic filter before. Indeed, we are doing all filters before. As this filter acts only on REF, ALT and AA_upcase, not considering genotypes or samples, other filters should not affect. If one SNP is not phased or biallelic, then it will be removed before, because we do not want to do the switch for it")
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --exclude 'REF!=AA_upcase && ALT!=AA_upcase' | \
        bcftools query \
            -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AF AA:%AA AA_upcase:%AA_upcase GTs:[ %GT]\n' | \
        head -n 5")


    print_text("Before the next step, check that REF and ALT are always ACGT, because our awk script relies on the fact that we discard those rows for which the unique ancestral allele in upper case is NOT equal to REF or ALT, removing in that way variants with ancestral allele equal to '.', 'N', '-', etc.... I have checked that regex here is case sensitive, so looking for ACGT does not look for acgt in REF/ALT, meaning we are also checking that REF/ALT are in upper case. Also check that REF and ALT are NOT the same", header=3)
    problematic_cases_ref_alt = run_bash(" \
        bcftools view \
            --no-header \
            " + input_vcf_file + "| \
        awk \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\"\t\"; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + "; \
            }{ \
                if($index_ref !~ /A|C|T|G/ || $index_alt !~ /A|C|T|G/ || $index_ref==$index_alt){ \
                    count++ \
                } \
            }END{print count}'", return_value=True).strip()
    if(problematic_cases_ref_alt == ""):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM, REF OR ALT ARE NOT ALWAYS ACGT OR REF IS EQUAL TO ALT!!")
        #awk
            #begin with tabs as delimiter and create variables with the index of REF and ALT columns
            #if REF or ALT does NOT include ACGT, count
                #use regex expression to have multiple conditions using "|"
                #negate with "!"
                #https://stackoverflow.com/a/8481180/12772630
            #print count
        #check that the count is exactly 3, because we have added 3 more dummy SNPs for which REF or ALT does not include ACGT


    print_text("save the subset of filtered variants and then apply on it the awk script", header=3)
    run_bash(" \
        bcftools view \
            " + input_vcf_file + " | \
        bcftools +fill-tags \
            -- --tags AN,AC,AF | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --output ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz \
            --output-type z \
            --compression-level 1")


    print_text("check AA is never ';' nor ',' because of cases with 'AA=;' or 'AA=,,,;'. bcftools should always put '.' for empty ", header=3)
    problem_comma_semicolon = run_bash(" \
        bcftools view \
            --no-header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz | \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\"; \
                index_info=" + index_info + " \
            }{ \
                for(i=1;i<=length($index_info);i++){ \
                    if(substr($index_info, i, 3)==\"AA=\"){ \
                        k=i+3; \
                        while(substr($index_info, k, 1)!=\",\" && substr($index_info, k, 1)!=\";\"){ \
                            last_base=k; \
                            k=k+1 \
                        }; \
                        unique_aa=toupper(substr($index_info, i+3, last_base-(i+3)+1)) \
                    } \
                }; \
                if(unique_aa==\",\" || unique_aa==\";\"){ \
                    exit 1 \
                } \
            }'; \
        echo $?", return_value=True).strip()
    if problem_comma_semicolon=="0":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM: Some SNPs have as unique AA comma or semi-colon and this should not be the case. At least, they should have '.', as this is the string used for empty in bcftools")
    #load the VCF file before switching and without header
    #in awk
        #extract the unique AA like in previous steps
        #if the allele is "," and ";"
            #stop with non-zero exit status
    #print the exit status
    #if that status is not "0", stop the script


    print_text("apply the awk script to the cleaned VCF", header=3)
    run_bash(" \
        gunzip \
            --stdout \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz | \
        awk \
            -v date=\"$(date '+%a %b %d %H:%M:%S %Y')\" \
            'BEGIN{ \
                FS=\"\t\"; \
                OFS=\"\t\"; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + "; \
                index_info=" + index_info + "; \
                index_format=" + index_format + "; \
                selected_chrom=\"chr" + selected_chromosome + "\"; \
                header=\"yes\"; \
            }{ \
                if(header==\"yes\"){ \
                    if($0 !~ /^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/){ \
                        print $0; \
                    } else { \
                        printf \"##awk_script: Switched REF/ALT columns for variants where REF!=AA but ALT==AA. Genotypes have been also switched (0 becomes 1 and vice versa). Variants where REF!=AA and ALT!=AA where discarded; Date=%s\\n\", date; \
                        print $0; \
                        header=\"no\"; \
                    } \
                } else if(header==\"no\" && $1==selected_chrom){ \
                    for(i=1;i<=length($index_info);i++){ \
                        if(substr($index_info,i,3)==\"AA=\"){ \
                            k=i; \
                            while(substr($index_info, k, 1)!=\";\"){ \
                                end_aa=k; \
                                k=k+1; \
                            }; \
                            if(i+3 != end_aa){ \
                                full_aa = substr($index_info, i+3, end_aa-(i+3)+1) \
                            } else { \
                                full_aa = substr($index_info, i+3, 1) \
                            }; \
                            if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                                full_aa_no_comma = full_aa; \
                                gsub(/,/, \"\", full_aa_no_comma); \
                                count=0; \
                                for(k=1;k<=length(full_aa_no_comma);k++){ \
                                    char=substr(full_aa_no_comma, k, 1); \
                                    if(!array_alleles[char]++){count++}; \
                                }; \
                                if(count==1){ \
                                    check_n_csq=\"false\"; \
                                    for(k in array_alleles){ \
                                        if(array_alleles[k]==length(full_aa_no_comma)){ \
                                            check_n_csq=\"true\" \
                                        } \
                                    }; \
                                    if(check_n_csq==\"true\"){ \
                                        k=1; \
                                        while(substr(full_aa, k, 1)!=\",\"){ \
                                           before_comma=k; \
                                           k=k+1\
                                        }; \
                                        unique_aa = substr(full_aa, 1, before_comma); \
                                    } else { \
                                        exit 1 \
                                    }; \
                                } else { \
                                    exit 1 \
                                }; \
                                delete array_alleles; \
                            } else { \
                                unique_aa=full_aa \
                            }; \
                            unique_aa = toupper(unique_aa); \
                        } else { \
                            if(substr($index_info,i,3)==\"AA;\"){exit 1} \
                        } \
                    }; \
                    if($index_ref==unique_aa || $index_alt==unique_aa){ \
                        if($index_ref!=unique_aa && $index_alt==unique_aa){ \
                            tmp_ref=$index_ref; \
                            $index_ref=$index_alt; \
                            $index_alt=tmp_ref; \
                            for(i=(index_format+1);i<=NF;i++){ \
                                if(index($i, \"X\")!=0){exit 1}; \
                                gsub(/0/, \"X\", $i); \
                                gsub(/1/, \"0\", $i); \
                                gsub(/X/, \"1\", $i) \
                            }; print $0 \
                        } else { \
                            print $0 \
                        }; \
                        if($index_ref!=unique_aa && $index_alt==unique_aa){exit 1} \
                    }; \
                } else { \
                    exit 1 \
                }; \
            }' | \
        bgzip --force > ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_raw.vcf.gz")
        #awk
            #load variables 
                #the date following the format used by bcftools
                    #The command is "+FORMAT", where FORMAT controls the output. For example: "%H" gives hours.
                    #In order to have spaces between the units, you need to use '' around FORMAT.
                    #https://unix.stackexchange.com/a/224976
                    #https://stackoverflow.com/a/27337807/12772630
                    #https://man7.org/linux/man-pages/man1/date.1.html
            #begin awk
                #loading the data into awk using \t for input and output. In this way, we are going to work with the main columns of the VCF file, REF, ALT, INFO, FORMAT, GT...
                #This means that the whole INFO field is going to be a string for each row, i.e., "AA=g;AA_upcase=G...." because we are not using ";" nor "=" as delimiters.
                #we do this to get as output the original INFO field with the ";",  "=", because if you use these as delimiters, they disappear from the data.
                #create variables with the position of the columns of interest that we will use later
                #also set a variable indicating "header=yes", because the file starts with header lines. This variable will change to "no" when the end of the header is reached.
            #if the row belongs to header, i.e., header==yes
                #if the row does NOT starts with CHROM, POS...
                    #just print the whole row
                    #https://unix.stackexchange.com/a/72763
                    #https://stackoverflow.com/a/3701095/12772630
                #else means we are in the last row of the header
                    #before printing that row, add a comment with a explanation of the awk script and the time, ending with new line (\n)
                    #print the whole row
                    #change "header" to "no", so the next line will NOT be considered as header
            #else, if header==no (row is after the header) and the first column is the selected chromosome
                    #https://stackoverflow.com/a/3701095/12772630
                #for each character of the INFO field
                    #use substring to extract specific characters from the INFO string. 
                        #substr(string, start [, length ])
                        #Return a length-character-long substring of string, starting at character number start. The first character of a string is character number one. For example, substr("washington", 5, 3) returns "ing".
                            #https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html
                    #if a substring starting at i up to two positions after (i.e., length=3) is equals to "AA=" (this avoids "AA_upcase=")
                        #while loop: iterate from the starting position of "AA=", i.e., "A". You save the position of the character "end_aa", overwritting the position of the previous character until we reach ";", which is NOT saved. 
                            #this will give the position of the last character before ";", i.e., before the next INFO field.
                            #https://opensource.com/article/19/11/loops-awk
                        #Use this position to select only the ancestral allele data wihthout "AA=" and ";". 
                        #If the last character before ";" is the position just after "AA=", i.e., "AA" is equal to a single character (e.g., AA=T or AA=-). In other words, three positions after "AA=" starts (i+3) is equals to the last position before ";" (end_aa)
                            #select just the character in that position
                        #else, means that the we have several characters within "AA="
                            #select a substring starting after "AA=" and ending just before ";"
                            #we obtain that position by substracting the start of "AA=" from "end_aa"
                            #For example: In "AA=G,G,G;", "AA=" starts at 1 (i=1), while the first base is in 1+3=4 (i+3=4), and the last character before ";" is at position 8 (end_aa=8). The length of the substring is 8-4+1=5 (end_aa-(i+3)+1). Therefore, if from the position of the last base (end_aa=8), we subtract the position of the first base just after "AA=" (i+3=4), we get the distance between the two characters without including one of them, so you have to add 1 to include botch extremes of this interval: 8-(1+3)+1=5. We want a substring of length 5 to include all alleles.
                                #Note that, in substr, the start is included so you have to start at 4, including 4, and then add 5 positions starting from 4 to reach the final base. First G is at 4, last G is at 8, so 8-4+1=5. This is 5,6,7 plus the extremes, 4 and 8.
                        #if "full_aa" does not include "," or it is exactly "," but it has no other character, this means that "AA" is just a base or a missing including "," (",", "N", "-", ".").
                            #then you can just save "full_aa" in "unique_aa" because there is only 1 character, so it is already unique.
                            #Note that we use bcftools to extract AA from the CSQ field, and I checked that this makes cases with empty "AA" (AA="") as "AA=.". So we will always have a character after AA, and no space.
                                #We have added other "if" to check we have no empty AA field just in case
                        #if it includes "," but it is not exactly ",", this means that we have several alleles for this SNP, e.g., "A,A,A,A", so we need to get just one. 
                            #first check again (we did before) that we only have one ancestral allele per variant
                                #copy "full_aa" into "full_aa_no_comma"
                                #using gsub, change every "," by empty (""), so now we should have only AA data and no commas, the same allele repeated many times
                                    #https://unix.stackexchange.com/a/492502
                                #set a "count" variable as zero
                                #for each character in full_aa_no_comma
                                    #extract that character
                                    #save it in array "array_alleles", 
                                        #if it is previously present then add 1 more to its count. 
                                        #if is not present before, add it as a new entry with value 1.
                                            #Negate the expression so this case is true and then we can add 1 to our "count" variable.
                                            #In other words, count the number of unique characters in full_aa_no_comma
                            #We do operations if count==1, if it is higher than 1 and hence we have more than 1 distinct characters, STOP (non-zero exit status). 
                                #this check will fail for microsatellites with different bases like TAT, but this should be ok because we should have only SNPs. If we have different bases in a variant, I want an error.
                                #we do operations if another check works, if not, stop with non-zero exit status. Take the number of counts we have in array "array_alleles", which stores the unique characters seen in "full_aa_no_comma" with its count. This number should be the same than the total length of "full_aa_no_comma", i.e., the number of characters without comma should be the same than the count for the single character we got. Remember that the previous if check if 'count'==1, so we already know we have only 1 different character.
                                #If true:
                                    #iterate across characters of "full_aa"
                                        #while the character is NOT ",", save the position of the character in "before_comma" and move to the next position
                                        #if "," it reached, then "before_comma" is not updated more, and hence we have the last position just before ","
                                    #now extract a substring from "full_aa" from the start to the position just before the first comma. For example, in TT,TT,TT: would be from 1 to 2, because the first comma is in position 3. Remember that the range of awk substring includes both extremes!
                            #delete "array_alleles", so we have it clean for the next iteration of the larger loop (iterated over i)
                        #else, means that we only have 1 character after "AA=" and before ";", therefore just take that character.
                        #convert the resulting unique ancestral allele (unique_aa) into uppercase to consdier both low and high-confidence alleles.
                    #else means that we do not have "AA=". This is ok because we can have other INFO fields, BUT
                        #if the INFO field is "AA;" STOP, because we should have always at least "." in case AA is missing.
                #if the REF or the ALT fields are equal to the unique ancestral allele
                    #if the REF is NOT the AA, but the ALT it is, we need to switch them
                        #save the REF in a temp variable
                        #overwrite the original REF with the ALT
                        #save as new ALT the REF value stored in the temp variable
                        #for each column starting after the first FORMAT column, i.e., after GT, thus targeting each genotype
                            #if X is included in the genotype, then stop with non-zero status, because we need not have any "X" character. We will use that character as a temporal step to exchange 0 for 1 and viceversa.
                            #substitute "0" by "X" in the genotype
                            #substitute "1" by "0" in the genotype
                            #substitute "X" by "1" in the genotype
                        #print the whole row
                    #else just print the row without changes
                #else, we do not have ancestral allele for this variant in the panel, so we can discard the SNP, so no print. Without ancestral allele there is anything we can do with that allele.
            #else means we have an error because rows after the header should all have the selected chromosome in the first column: stop with non-zero exist status        


    print_text("check that the switch of REF/ALT is the same than if we use plink 2", header=3)
    print_text("make the switch between REF and ALT using plink2", header=4)
    run_bash(" \
        plink2 \
            --vcf ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz \
            --ref-allele \
                force \
                <( \
                    bcftools view \
                        --no-header \
                        ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz | \
                    awk \
                        'BEGIN{ \
                            FS=OFS=\"\t\"; \
                            index_id=" + index_id + "; \
                            index_ref=" + index_ref + "; \
                            index_alt=" + index_alt + "; \
                            index_info=" + index_info + " \
                        }{ \
                            for(i=1;i<=length($index_info);i++){ \
                                if(substr($index_info, i, 3)==\"AA=\"){ \
                                    k=i+3; \
                                    while(substr($index_info, k, 1)!=\",\" && substr($index_info, k, 1)!=\";\"){ \
                                        last_base=k; \
                                        k=k+1 \
                                    }; \
                                    unique_aa=toupper(substr($index_info, i+3, last_base-(i+3)+1)) \
                                } \
                            }; \
                            if(unique_aa == $index_ref || unique_aa == $index_alt){ \
                                print unique_aa, $index_id \
                            } \
                        }' \
                ) \
                1 \
                2 \
            --recode vcf bgz \
            --threads 1 \
            --out ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_plink_raw")
        #--ref-allele sets all alleles specified in the file to REF, while --alt1-allele does the same for the first ALT allele. Column and skip parameters work the same way as with --update-chr and friends.
            #--force
                #By default, these error out when asked to change a 'known' reference allele. Add the 'force' modifier to permit that (when e.g. switching to a new reference genome).
                    #we want to change the REF even it is known, because we know the SNPs are not polarized
            #filename
                #using process substitution, we select those SNPs for which the AA is REF or ALT. Then extract the AA and the ID of the SNP.
                #this is the file used for plink to polarize, as the allele indicate in that file will be used as REF
            #REF col. number
                #the column from where the REF allele will be taken
            #variant ID col. number
                #the column from where the ID of the SNP will be taken
        #--recode
            #the output file will be a gzipped VCF file
        #https://www.cog-genomics.org/plink/2.0/data#ref_allele

    print_text("compare the switch with plink2 and with my approach", header=4)
    plink_comparison = run_bash(" \
        cmp \
            --silent \
            <( \
                bcftools view \
                    --no-header \
                    ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_plink_raw.vcf.gz | \
                awk \
                    'BEGIN{ \
                        FS=OFS=\"\t\"; \
                        index_chrom=" + index_chrom + "; \
                        index_ref=" + index_ref + "; \
                        index_alt=" + index_alt + "; \
                        index_info=" + index_info + " \
                    }{ \
                        for(i=1;i<=length($index_info);i++){ \
                            if(substr($index_info, i, 3)==\"AA=\"){ \
                                k=i+3; \
                                while(substr($index_info, k, 1)!=\",\" && substr($index_info, k, 1)!=\";\"){ \
                                    last_base=k; \
                                    k=k+1 \
                                }; \
                                unique_aa=toupper(substr($index_info, i+3, last_base-(i+3)+1)) \
                            } \
                        }; \
                        if(unique_aa == $index_ref){ \
                            $index_chrom = \"chr\"$index_chrom; \
                            print $0 \
                        } \
                    }' \
            ) \
            <( \
                bcftools view \
                    --no-header \
                    ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_raw.vcf.gz \
            ); echo $?", return_value=True).strip()
    if(plink_comparison == "0"):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: Making the switch of REF/ALT with plink 2 does not give the same results")
        #using cmp, compare byte by byte two files with silent output
            #file 1
                #the alleles switched with plink 2, but also
                #select only those rows for which the AA is REF, because plink does not filter out those cases where the AA is not REF nor ALT
                #also add "chr" to the chromosome name because plink tend to remove "chr" from that field
            #file 2
                #the alleles switched with my approach
        #print the exit status

    ##checking code from here

    print_text("Also update the INFO fields like AF, AC... because these are calculated considering the old REF/ALTs, but now the AC should be the count of the old REF not the old ALT", header=3)
    #before this step, Jesús created a tabix index from the VCF file (tabix --preset vcf VCF_FILE). I have checked that the creation of a tabix index does NOT change the VCF file. From what I have read, having a tabix index is useful to make fast queries over a large file without having to read the whole file, but I am not sure if this is useful for our case, as we have to do the update of the INFO fields in each row, i.e., in all rows.
        #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
        #https://www.htslib.org/doc/tabix.html#:~:text=Firstly%2C%20tabix%20directly%20works%20with,most%20SQL%20databases%20do%20not.
    run_bash(" \
        bcftools annotate \
            --remove ^INFO/AA,INFO/AA_upcase,^FORMAT/GT \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_raw.vcf.gz | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --output-type z \
            --compression-level 1 \
            --output ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz")


    print_text("check again that AA_upcase is just AA in uppercase, because we are going to use AA_upcase as a check in the next lines", header=3)
    problematic_aa_aa_upcase = run_bash(" \
        bcftools view \
            --no-header \
            --drop-genotypes \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz | \
        awk \
            'BEGIN{FS=\"\t|;|=\"}{ \
                for(i=1;i<=NF;i++){ \
                    if($i == \"AA\"){ \
                        if(toupper($(i+1)) != $(i+3)){ \
                            count++ \
                        } \
                    } \
                } \
            }END{print count}' \
        ", return_value=True).strip()
    if(problematic_aa_aa_upcase == ""):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: AA_upcase is not just AA in uppercase")


    print_text("After switching, REF is always equals to AA_upcase while ALT is not?", header=3)
    count_problem_ref_aa = run_bash(" \
        gunzip \
            --stdout \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz | \
        awk \
            'BEGIN{ \
                FS=\"\t|;|=\"; \
                OFS=\"\t\"; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + " \
            }{ \
                for(i=1;i<=NF;i++){ \
                    if($i==\"AA_upcase\"){ \
                        full_aa=$(i+1); \
                        if(index(full_aa, \",\")!=0 && full_aa!=\",\"){ \
                            k=1; \
                            while(substr(full_aa, k, 1) != \",\"){ \
                                before_comma=k; \
                                k=k+1; \
                            }; \
                            unique_aa=substr(full_aa, 1, before_comma); \
                        } else { \
                            unique_aa=full_aa \
                        }; \
                        if(unique_aa!=$index_ref || unique_aa==$index_alt){ \
                            count++ \
                        } \
                    } \
                } \
            }END{print count}'", return_value=True).strip()
    if(count_problem_ref_aa == ""):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM: AFTER SWITCHING, REF IS NOT ALWAYS AA_upcase OR ALT IS EQUAL TO AA_upcase")
        #this is a very good check because we have previously done the switch without using AA_upcase, and now we check  if REF is always equals to AA_upcase, which was obtained using other approach.


    print_text("check that the SNPs discarded when switching REF/ALT are indeed SNPs without ancestral data", header=3)
    problematic_ref_alt_aa_count = run_bash(" \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\"; \
                index_id=" + index_id + "; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + "; \
                index_info=" + index_info + " \
            }{ \
                if(FNR==NR){ \
                    a[$index_id]; \
                    next \
                }; \
                if(!($index_id in a)){ \
                    for(i=1;i<=length($index_info);i++){ \
                        if(substr($index_info, i, 3)==\"AA=\"){ \
                            k=i+3; \
                            while(substr($index_info, k, 1)!=\",\" && substr($index_info, k, 1)!=\";\"){ \
                                last_base=k; \
                                k=k+1 \
                            }; \
                            unique_aa=toupper(substr($index_info, i+3, last_base-(i+3)+1)) \
                        } \
                    }; \
                    if($index_ref == unique_aa || $index_alt == unique_aa){ \
                        count++ \
                    } \
                } \
            }END{print count}' \
            <( \
                bcftools view \
                    --no-header \
                    --drop-genotypes \
                    ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz \
            ) \
            <( \
                bcftools view \
                    --no-header \
                    --drop-genotypes \
                    ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz \
            )", return_value=True).strip()
    if problematic_ref_alt_aa_count=="":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM: Some SNPs that have been discarded have indeed ancestral allele equal to REF or ALT, they should be included!!")
    #we are going to select those SNPs from the cleaned VCF that were excluded in the next step as their REF/ALT is not equal to the ancestral allele.
    #we are using awk with two files as input, the first one is the VCF file where SNPs without ancestral allele have been removed. We will use this file to make an array of IDs, and then, in the original VCF (second file), select those SNPs whose ID is NOT included in that array, i.e., they are not in the new VCF file.
        #I am processing the files with bcftools in order to remove the header and then directly using them as input for awk with by process substitution
            #Piping the stdout of a command into the stdin of another is a powerful technique. But, what if you need to pipe the stdout of multiple commands? This is where process substitution comes in.
            #https://tldp.org/LDP/abs/html/process-sub.html
    #script based on
        #https://stackoverflow.com/a/14062579/12772630
    #BEGIN
        #using "\t" as separator. In this way, we will work with the main columns, i.e., CHROM, POS, ID, REF, ALT....
        #create variables with the index of some of these columns. The order should be same than in the original file, because we are not changing the columns...
            #we are checking this anyways at the end of this part of the script (see below)
    #if the row belongs to the first file
        #create an array with the SNP ID as index
        #go to the next row, i.e., do not run anything else from the script
            #we use next because we do not want to run more for the first file, so next force to stop running more code and go to the next row. I usually do not want this as I want to run all the script, but I want it in this case.
            #https://www.gnu.org/software/gawk/manual/html_node/Next-Statement.html
        #we know if a row belongs to the first file using if(FNR=NR).
            #NR starts at 1 with the first row of the first file and ends at the last row of the last file.
                #if you have 2 files with 2 rows each one, NR will end at 4, because we have 2 files x 2 rows/file=4 rows
            #FNR, in contrast, is re-started in every new file.
                #in our previous example, the first row of every file will have FNR=1
            #Therefore, only for the first file FNR will be the same than NR.
            #awk '{print NR, FNR}' <(echo -e "hola\tmundo\nhola\tmundos") <(echo -e "hola\tmundo\nhola\tmundos")
        #see dummy example
            #!awk '{if(FNR==NR && NR>1){a[$1];next}; if(($1 in a) && FNR>1){print $0}}' <(echo -e "id\tname\n2\tManu") <(echo -e "id\tname\n1\tJuan\n2\tManu\n3\tPepe")
    #leave for the second file
        #if the SNP ID of the SNP is NOT included as an index the array, this means this row was discarded and should not have ancestral data
            #use (!()) to negate condition in "if" conditional
                #https://stackoverflow.com/a/55940407/12772630
            #for each character in the info field
                #if the selected character plus the next 2 are "AA="
                    #iterate over the characters after "AA=" and save its position until "," or ";" are found, so we get the position just before the "," or ";", i.e., the first copy of the allele. 
                        #we can have AA=T,T,T; or AA=T; or AA=TT;
                        #in all cases we get the first copy
                    #select the strings between the position just after "AA=" and the position just before "," or ";"
                        #substr($index_info, i+3, last_base-(i+3)+1)
                        #Example:
                            #NS=107;AA=TT,TT,TT;
                            #AA starts at 8, thus the "i" we use is 8
                            #i+3=11, thus the alleles start at 11
                            #last_base=12, i.e., last base before "," or ";" is 12
                            #we want everything between 12 and 11, including both extremes, thus 12-11+1=2, i.e., last_base-(i+3)+1.
                            #We start at 11, and get 11 and 12 (two positions), i.e., TT
                    #convert to upper case
            #REF nor ALT should be NOT equal to the ancestral allele, because if one of them is, then we would be able to do the switch between columns and then, this SNP should not be lost
                #if one of them is, add 1 to count
    #After reading all files (END), print the count


    print_text("Check we get the same rows if we take the original VCF and select only rows where AA is REF or ALT, switching only when ALT==AA", header=3)
    run_bash(" \
        bcftools view \
            --no-header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz | \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\"; \
                index_ref=" + index_ref + "; \
                index_alt=" + index_alt + "; \
                index_info=" + index_info + "; \
                index_format=" + index_format + "; \
            }{ \
                for(i=1;i<=length($index_info);i++){ \
                    if(substr($index_info, i, 3)==\"AA=\"){ \
                        k=i+3; \
                        while(substr($index_info, k, 1)!=\",\" && substr($index_info, k, 1)!=\";\"){ \
                            last_base=k; \
                            k=k+1 \
                        }; \
                        unique_aa=toupper(substr($index_info, i+3, last_base-(i+3)+1)) \
                    } \
                }; \
                if(unique_aa==$index_ref || unique_aa==$index_alt){ \
                    if(unique_aa!=$index_ref && unique_aa==$index_alt){ \
                        tmp_ref=$index_ref; \
                        $index_ref=$index_alt; \
                        $index_alt=tmp_ref; \
                        for(i=(index_format+1);i<=NF;i++){ \
                            if(index($i, \"X\")!=0){exit 1}; \
                            gsub(/0/, \"X\", $i); \
                            gsub(/1/, \"0\", $i); \
                            gsub(/X/, \"1\", $i) \
                        }; \
                    }; print $0 \
                } \
            }' | \
        cut \
            --delimiter '\t' \
            --fields 1-7,9- > ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_1.vcf; \
        bcftools view \
            --no-header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz | \
        cut \
            --delimiter '\t' \
            --fields 1-7,9- > ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_2.vcf; \
        check_status=$( \
            cmp \
                --silent \
                ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_1.vcf \
                ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_2.vcf; \
            echo $?); \
        if [[ $check_status -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi; \
        rm ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_1.vcf; \
        rm ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/test_file_2.vcf")
        #create a file with the VCF file before switching but selecting only those rows for which REF or ALT are AA and then switch REF/ALT if ALT==AA. Also substitute 1 by 0 and viceversa in these cases for the genotypes to have as 0 the new REF. We do this from the first column after GT (field 8, FORMAT), because that is where genotypes start. We also check that "X" is not present in the genotype because we set the "0" genotypes as "X" in order to change "1" to "0" and then "X" to "1", avoiding converting all genotypes to "0" in the first step. Select all fields except the INFO field (8), which is not updated by awk. Therefore, this fields is not going to be the same in the new cleaned file where INFO has been updated using bcftools.
            #you can check previous AWK script for information
        #get the switched VCF file and select all fields except the INFO fields (8)
        #NOTE that, in cut, we are assuming that the position of the INFO field is 8, so we are selecting all fields except this one. If INFO is not 8, the script will fail at the beginning anyways, so it should be ok.
        #check byte by byte that both files are the same with cmp, is equal, the exit status should be zero.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see the header of the final VCF file", header=3)
    run_bash(" \
        bcftools head \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz")
            #you can see that after the list of contigs, the only field shown before the list of my commands is just GT (phased genotypes) because I have removed all INFO and FORMAT fields with the exception of FORMAT/GT.
            #then we have all the commands I have run in order to subset individuals and filter SNPs.
            #see dummy example for further details.


    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see the genotypes of a few individuals from the final VCF file", header=3)
    run_bash(" \
        bcftools view \
            --no-header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz | \
        head -5")
            #As expected, we have 
                #CHROM
                #POS
                #ID
                #REF
                #ALT
                #QUAL (empty)
                #FILTER (empty)
                #INFO (with the fields I specifically created)
                #FORMAT (only GT)


    print_text("check we still have the correct column names in order. Both the fixed fields (CHROM, POS...) and the genotype columns", header=3)
    problematic_vcf_columns_header = run_bash(" \
        bcftools view \
            --header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz | \
        awk \
            'BEGIN{ \
                FS=OFS=\"\t\" \
            }END{ \
                if($0 ~ /#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(selected_samples) + "/){ \
                    print \"TRUE\" \
                } else { \
                    print \"FALSE\" \
                }; \
            }'", return_value=True).strip()
    if problematic_vcf_columns_header=="TRUE":
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM: The header of the final VCF file does not have the correct columns. The problem can be the FIXED fields like CHROM, POS... or the columns for the sample genotypes")
        #get the header of the last VCF file
        #open in awk
            #using "\t" as delimiter
            #in the last row, which is the column names
                #if the row has exactly the columns names we expect, i.e., the 8 fixed fields and the FORMAT genotype columns
                    #print "TRUE"
                #else
                    #print "FALSE"
        #stop the script if we do NOT have TRUE


    #summary of the filters applied with vcftools in vcf files from slim simulations
        #--max-alleles 2 --min-alleles 2 --max-missing 1 --phased
            #https://vcftools.sourceforge.net/man_latest.html
        #--max-alleles / --min-alleles
            #Include only sites with a number of alleles greater than or equal to the "--min-alleles" value and less than or equal to the "--max-alleles" value. One of these options may be used without the other.
                #For example, to include only bi-allelic sites, one could use: vcftools --vcf file1.vcf --min-alleles 2 --max-alleles 2
        #--phased
            #Excludes all sites that contain unphased genotypes
        #--max-missing 1
            #Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
            #we used 1 because we were working with simulations, we have for sure data for all samples and snps.

    #my questions to David about SNP filters
        #Inbreeding in 1KGP
            #The new dataset has 3202 samples instead of 2504. It is not obvious if all the new 698 samples (3202-2504=698) are related. Only 600 are implicated in duos/trios, but both Jesús and I think that the other 90 have some degree of inbreeding. They repeat several times in the paper that they added 698 samples related to the original 2504 samples.
            #Within the original set of 2504 samples, I have detected something confusing. According to the new pedigree, 4 samples have their mother/father within the same 2504 samples. In other words, we have 4 duos within the original set of 2504 unrelated samples. This is not new, as Jesús told me, there is a 2015 scientific report ([link](https://www.nature.com/articles/srep17453 "‌")) showing cryptic relatedness in phase 3. The strange thing is that they generate a new pedigree showing duos within the original dataset, but they continue referring to that original dataset as the 2504 unrelated samples.
            #I think we should avoid all the new samples and focus on the original 2504 dataset. Within that original dataset, we could use it in full or remove the 4 samples that are related according to the new pedigree leaving only the parent of each duo. What do you think?
        #Filtering of variants within populations
            #I am using bcftools to select **biallelic** **SNPs that are not monorphic, no duplicates** (i.e., exact same position and REF-ALT alleles) and phased. Importantly, I am doing this within each population, so if a SNP has 3 alleles considering all populations of the panel but only 2 alleles within the selected population, then that SNP is retained for the selected population. This is what we want, right?
            #I have also filtered by percentage of **missing**
                #1KGP authors selected only variants with missingness < 5%.
                #I have additionally removed variants with **any missing genotype**. Is that ok or should I retain all SNPs with missing<5%?
            #Filter by accessibility
                #I am also using the accessibility masks. Jesus told me that it could be relevant to use the mask if we are going to use data from the whole genome. This is the case at least for the climate project because we were thinking of using non-overlapping across the whole genome, not only centered around genes.
                #It is important to note that these masks select regions that are accessible based on the alignment of whole genome **low coverage** data of 2691 phase3 samples to hg38. Do you think it is ok if the masks are based on low instead of high coverage data?
                #Also, one of the filters they apply to consider a region accessible or not is if "_base is an N in the reference genome GRCh37_".
                #Finally, I have checked how many SNPs are lost in chromosome 1 as an example. The VCF file without any filter has 5,013,617 SNPs, with the less stringent mask this number decreases to 4,616,062, while with the more stringent mask it gets reduced to 3,576,231. Therefore, we lose 1.5 million SNPs.
                #I am not sure whether we should use these masks given they are based in the low-coverage genomes and they seem to reduce the number of SNPs a lot. What do you think?
            #filter by MAF?
                #At least for iHS, we usually remove SNPs with MAF<0.05 directly in hapbin.
                #I guess I should **not** apply a MAF filter right now in the vcf files, but then Elise will apply the required MAF filters when calculating the different summary statistics, right?
            #filter by HWE within the population?
                #1KGP authors selected SNPs that pass HWE filter (p-value>1e-10) in at least one superpopulation, but it is possible that some of these SNPs violate HWE within a specific population.
                #Should I also filter by HWE within each specific population?
    #answer of David
        #You should apply the filter with less than 5% missing just as the 1KGP authors.
            #Not add any filter about missing, because the data is already filtered for less 5% of missing.
        #For HWE you should also only apply what the 1KGP authors did and not filter further for specific populations.
            #Not apply any filter about HWE, because the authors already did at the superpopulation level.
        #Biallelic SNPs per specific population are ok.
            #so I can discard SNPs with less than 2 alleles and more than 2 alleles within each specific population.
            #I understand this includes remove monomorphic. For what reason do you want a SNP that is fixed within your population?
            #Also I guess I can remove SNPs with unphased genotypes within the selected population.
        #Masks based on low coverage are ok, you should use the less stringent masks.
            #The less stringent is the pilot.
        #You can use the 2,504 individuals and remove the four related individuals yes.
            #We left 2500 samples.
        #And you are right for MAF, the filtering at MAF>5% is done by the scripts for summary statistics. It is only for the focal SNPs methods like iHS still use the SNPs with lowe MAF around the focal SNPs so it would be an error to remove all SNPs with MAF<5%.
            #I understand that not all summary statistics use the SNP with lower MAF, so we would lose information for these statistics.
            #It is better to specifically filter for MAF when calculating iHS.
        #my follow-up questions
            #When you said "You should apply the filter with less than 5% missing just as the 1KGP authors", you mean that I should not add any further filter about missingness within each population and just use what they already did, i.e., < 5% missing considering all pops, right? Because this filter was already applied in the data I downloaded. 
            #The same would go for HWE. The filter across superpopulations is already applied in the data, so I will just leave it as it is in that regard.
            #Just to double check, I remove a SNP that is monomorphic within a given population, even if it is not monomorphic considering the whole panel. The same goes for phasing, I select SNPs with all genotypes being phased for the specific population, irrespectively of the phasing in other pops.
            #david answered YES to all this.
        #reduction in the number of SNPs after applying filters
            #For example
                #chromosome 1 has 5,013,617 SNPs
                #removing monomorphic SNPs within the british reduces the number to 4,067,698 SNPs.
                #after applying the rest of filters I get 852,072 SNPs for the british, i.e., 4 million of SNPs lost.
            #Is this ok? or it is too much?
                #ok 
        #as I am filtering within pop, each pop will have different snps
            #is this ok?
                #ok



    ####HERE IMPORTANT
    ## USE THE GENETIC MAP OF THE SELECTED CHROMSOOME PREVIOUSLY CREATED, USE IT TO
    #FILTER HAP 
    #FILTER THEN MAP ALSO!!

    #the input VCF file has was changed in the script, check it and put back the cleaned_ref_alt_switeched VCF

    print_text("Create the final HAP and MAP files by selecting only SNPs with genetic position for both files", header=2)
    print_text("chr " + selected_chromosome + ": filter the already cleaned VCF with bcftools", header=3)
    #this file is cleaned regarding biallelic snps, duplicates... but need to retain only SNPs with genetic position
    #We could do this by just creating before the hap file, extract snp positions from there, calculate genetic position and then remove from the hap those rows of SNPs without genetic position. The problem is that we would do that by row index instead of SNP ID, at least if we use the final hap file, so we are going for this option better. In addition, we would have snps that cannot be used in the VCF file because they do not have genetic position. With the other approach we would have a VCF with all SNPs filtered and another one with only snps with genetic position.
    
    print("filter")
    run_bash("\
        bcftools view \
            --include ID==@./results/00b_map_files/chr" + selected_chromosome + "/list_snps_with_gen_pos_" + selected_chromosome + ".tsv.gz \
            --output ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz \
            --output-type z \
            --compression-level 1 \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz")
            #include those SNPs for which ID is included in the list of SNPs with genetic position and save the resulting VCF file
                #https://www.biostars.org/p/373852/

    print("see header of the fully filtered VCF file and some genotypes")
    run_bash(" \
        bcftools head \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz")
    run_bash(" \
        bcftools view \
            --no-header \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz | \
        head -5")

    print_text("check that IDs in the filtered VCF file are the same than the ones in the list of IDs used as input to filter", header=4)
    run_bash(" \
        STATUS=$( \
            cmp \
                --silent \
                <( \
                    bcftools view \
                        --no-header \
                        ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz | \
                    awk \
                        'BEGIN{ \
                            FS=OFS=\"\t\"; \
                            index_id=" + index_id + "} \
                        { \
                            print $index_id \
                        }' \
                ) \
                <( \
                    gunzip \
                        --stdout \
                        ./results/00b_map_files/chr" + selected_chromosome + "/list_snps_with_gen_pos_" + selected_chromosome + ".tsv.gz \
                ); \
            echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #compare two files with cmp in silent mode, as we will use the exit code to make the check
            #both files are directly processed and used as input with process substitution
                #Piping the stdout of a command into the stdin of another is a powerful technique. But, what if you need to pipe the stdout of multiple commands? This is where process substitution comes in.
                #https://tldp.org/LDP/abs/html/process-sub.html
            #file 1
                #load the VCF file with only the SNPs having genetic position using bcftools, remove the header
                #awk: for each row, print the ID
            #file 2:
                #just ungzip the file with the list of SNPs having genetic position
        #if the exist status is 0, we are good. Both files are identical, byte by byte.
            #An exit status of 0 means no differences were found, 1 means some differences were found, and 2 means trouble.
            #https://www.gnu.org/software/diffutils/manual/diffutils.html#Invoking-cmp


    ##load here the map file of the selected chromosome and select only snps included in the VCF
    ##FILTER THE MAP!!
    ##I guess we no longer need a raw map file...
    #remove the old id from the map, use only the new one!!!

    print_text("create the hap file", header=3)
    print("chr " + selected_chromosome + " - " + selected_pop + ": convert to hap file")
    run_bash(" \
        bcftools convert \
            --hapsample ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw \
            ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz")
            #load the cleaned vcf file
            #--vcf-ids 
                #when this option is given, the second column is set to match the ID column of the VCF.
                #Without the option, the format follows https://www.cog-genomics.org/plink/2.0/formats#haps with ids (the second column) of the form "CHR:POS_REF_ALT[_END], with the _END being optional for defining the INFO/END tag when ALT is a symbolic allele.
                #we are not using --vcf-ids so we can have "CHR:POS_REF_ALT". This can avoid strand swaps (see below).
            #--hapsample
                #convert from VCF to hap/sample format used by IMPUTE2 and SHAPEIT. The columns of .hap file begin with ID,RSID,POS,REF,ALT. In order to prevent strand swaps, the program uses IDs of the form "CHROM:POS_REF_ALT".
                #save the results and using the name of the chromosome and the pop
            #https://www.htslib.org/doc/bcftools.html#convert
        #IMPUTE2 hap format
            #https://www.cog-genomics.org/plink/2.0/formats#haps
            #we are going to use the reference panel haplotype file format for IMPUTE2. 
            #This is a text file with no header line, and either 2N+5 or 2N fields where N is the number of samples. In other words, you can have 5 initial columns with data about the variants, or just the haplotype columns. In the former case, the first five columns are:
                #Chromosome code
                #Variant ID
                    #This is in the format CHROM:POS_REF_ALT to prevent strand swaps. I guess you avoid strand swaps by having information about the REF and ALT alleles, so you can be sure the strand you are using.
                #Base-pair coordinate (POS)
                #REF allele
                #ALT allele
            #This is followed by a pair of 0/1-valued haplotype columns for the first sample, then a pair of haplotype columns for the second sample, etc. (For male samples on chrX, the second column may contain dummy '-' entries; otherwise, missing genotype calls are not permitted.)
            #Previous hap format used by us
                #hap IMPUTE format is the one required by hapbin
                    #https://github.com/evotools/hapbin#input-file-formats
                #The format I used for hap files used as input to hapbin was the second option, i.e., ONLY THE HAPLOTYPE COLUMNS.
                    #/xdisk/denard/dftortosa/genomic_determinants_recent_selection/hapbin_inputs/
                #In addition, hap files in yoruba folder of flexsweep are just our hap files for iHS, i.e., ONLY THE HAPLOTYPE COLUMNS.
                    #/xdisk/denard/lauterbur/yoruba/hapmap_YRI_hg19_decode2019/

    print_text("as the IDs are created during the hap file creation (see docs above), we should have the same IDs in hap and map files. Also chromosome and position should be the same between both files", header=4)
    run_bash(" \
        STATUS=$( \
            cmp \
                --silent \
                <( \
                    gunzip \
                        --stdout \
                        ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
                    awk \
                        'BEGIN{ \
                            FS=\" \"; \
                            OFS=\"\t\"} \
                        { print $1, $2, $3 }' \
                ) \
                <( \
                    gunzip \
                        --stdout \
                        ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
                    awk \
                        'BEGIN{ \
                            FS=\" \"; \
                            OFS=\"\t\"} \
                        { print $1, $2, $4 }' \
                ); \
            echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #see previous cmp comparison for explanations about the code

    #note about the output of hap conversion
        #I get an output like this
            #Hap file: ./results/02_hap_map_files_raw/chr1_GBR_IMPUTE2_raw.hap.gz
            #Sample file: ./results/02_hap_map_files_raw/chr1_GBR_IMPUTE2_raw.samples
            #852072 records written, 0 skipped: 0/0/0 no-ALT/non-biallelic/filtered

        #I guess non-biallelic and non-alt should be removed to meet impute requirements, but in our case we already selected those SNPs with 2 alleles only. This explains why we get 0/0/0. This is the number of SNPs with no-ALT, no-biallelic and filtered. All the filters were previously applied.
            #I have checked this in the dummy example.
        
        #if you take the VCF file, filter and then count lines, you get 852072 records, but then say that total is 945919, as when writting the file. "Lines   total/split/realigned/skipped" is produced for some commands like bcftools norm and I have seen in the dummy example that total is the total number of SNPs, while split are those splitting due to multiallelic is the multiallelic flag is used.
            #852072
            #Lines   total/split/realigned/skipped:  945919/0/0/0
            #Lines   total/split/realigned/skipped:  945919/0/0/0

        #the total number of snps in the raw vcf file of chr1 is 5759173, instead of 945919. What it can be happening here is that bcftools norm (the command generating the line output) is run after some (but not all) filters have been applied, thus the input number of SNPs is smaller but the final number (852072). Indeed, if you apply these filters applied before norm, you get 945919 variants, so it makes sense that bcftools norm, which is run after these filters, gives 945919 as total number of snps.

        #The smallest number (852072) is the number of SNPs after we have completely cleaned the vcf file. I have checked that for GBR.

    print_text("chr " + selected_chromosome + " - " + selected_pop + ": see first variants of the hap file", header=4)
    run_bash(
        "gunzip \
            --stdout \
            ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        head -3")
            #decompress the hap file and show in stdout 
    
    print_text("chr " + selected_chromosome + " - " + selected_pop + ": remove first columns of hap file to leave only haplotype columns", header=4)
    run_bash(
        "gunzip \
            --stdout \
            ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        cut \
            --complement \
            --delimiter ' ' \
            --fields 1-5 \
        > ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap; \
        gzip \
            --force \
            ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap")
            #extract the content of compressed hap file
            #remove the first 5 columns
                #--complement keeps the columns other than the ones specified in -f
                #--delimiter specifies the delimiter; in this case, a space
                #--fields specifies the columns to cut (rather than the columns to keep, since --complement is being used);
                    #https://unix.stackexchange.com/questions/222121/how-to-remove-a-column-or-multiple-columns-from-file-using-shell-command
            #save the result as a hap file
            #compress that file
                #-f option : Sometimes a file cannot be compressed. Perhaps you are trying to compress a file called “myfile1” but there is already a file called “myfile1.gz”. In this instance, the “gzip” command won’t ordinarily work. To force the “gzip” command to do its stuff simply use -f option
                #https://www.geeksforgeeks.org/gzip-command-linux/

    print_text("chr " + selected_chromosome + " - " + selected_pop + ": the clean hap file is just the raw hap file but without the first 5 columns?", header=4)
    run_bash(" \
        STATUS=$( \
            cmp \
                --silent \
                <( \
                    gunzip \
                        --stdout \
                        ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
                    awk \
                        'BEGIN{FS=OFS=\" \"}{ \
                            $1=$2=$3=$4=$5=\"\"; \
                            sub(\"     \", \"\"); \
                            print $0 \
                        }' \
                ) \
                <( \
                    gunzip \
                        --stdout \
                        ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz \
                ); \
            echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #check two files are identical using cmp, we are directly processing them as inputs on the fly with process substitution
            #file 1
                #decompress the raw hap file
                #awk
                    #use space as delimiter for input and outputs (this is the delimiter used by hap files)
                    #in each row
                        #set as "" the first 5 columns, which are the ones removed when cleaning the hap file
                        #this creates five empty spaces at the beginning of the row, so we substitute these spaces by nothing "". This is the equivalent of what we did when cleaning the hap file with cut, because we directly removed fields (columns), not just their content.
                        #print the resulting row
                        #https://stackoverflow.com/a/13446273/12772630
            #file 2
                #just decompress the cleaned hap file

    print_text("chr " + selected_chromosome + " - " + selected_pop + ": check we have the same number of rows (variants) in the cleaned vcf, hap and map files", header=4)
    run_bash(" \
        n_snps_vcf=$( \
            bcftools view \
                --no-header \
                ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz | \
            awk 'END{print NR}'); \
        n_snps_map=$( \
            gunzip \
                --stdout \
                ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            awk 'END{print NR}'); \
        n_snps_hap=$( \
            gunzip \
                --stdout \
                ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            awk 'END{print NR}'); \
        if [[ $n_snps_vcf -eq $n_snps_hap && $n_snps_vcf -eq $n_snps_map ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #load the cleaned vcf file with bcftools and show only the snps, with no header, then count the number of lines and save the result
            #decompress the hap and map files and count the number of line
            #check the counts are the same

    
    print_text("chr " + selected_chromosome + " - " + selected_pop + ": do we have the correct number of samples in the hap file?", header=4)
    run_bash(" \
        n_fields_hap=$( \
            gunzip \
                --stdout \
                ./results/03_hap_map_files/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            awk \
                'BEGIN{ \
                    FS=\" \"; \
                    n_snps=" + str(final_genetic_pos_map_file_pruned.shape[0]) + "} \
                { \
                    n_fields[NF]++ \
                }END{ \
                    if(length(n_fields)==1){ \
                        for(i in n_fields){ \
                            if(n_fields[i] == n_snps){ \
                                print i \
                            } \
                        } \
                    } \
                }' \
            ); \
        n_samples_hap=$(echo 'scale=0;' $n_fields_hap '/2'| bc); \
        n_samples_original=$( \
            awk \
                'END{print NR}' \
                ./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt \
            ); \
        if [[ $n_samples_hap -eq $n_samples_original ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #compare the number of samples between the final hap file and the txt file saved with the sample IDs
            #hap file
                #decompress the hap file
                #awk
                    #begin with " " as delimiter and create a variable with the number of SNPs we have in the final map
                    #create an array ("n_fields") and save there the number of columns adding a new entry for each new number of fields and summing 1 to the current count
                    #at the end
                        #the array should be 1-character long because we should only have 1 distinct number of fields, i.e., all rows have the same number of columns
                        #now we know we have only 1 entry in the array, 1 distinct number of fields, extract it, but only if the number of times it appears is exactly the same than the number of SNPs we have in the final map file
                            #remember that i is the index, i.e., the number of fields n_fields[i] is the count, i.e., the number of times that particular number of fields appeared across rows.
                #the number of fields in the hap file, divided by 2 is the number of samples we have, as we have 2 columns per sample. We use "bc" so we can specify the number of decimals.
            #the txt file saved with the sample IDs 
                #just calculate the number of rows, i.e., samples, using awk

    print_text("chr " + selected_chromosome + " - " + selected_pop + ": check that sample file generated with hap has the correct sample IDs?", header=4)
    #read the sample file
    sample_list_from_hap = pd.read_csv(
        "./results/02_hap_map_files_raw/" + selected_pop + "/chr" + selected_chromosome + "/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.samples",
        header=0, 
        sep=" ", 
        low_memory=False)
    #select those rows whose ID is not zero
    sample_list_from_hap_clean = sample_list_from_hap.loc[(sample_list_from_hap["ID_1"] != "0") | (sample_list_from_hap["ID_2"] != "0"), :]
    #reset the index
    sample_list_from_hap_clean = sample_list_from_hap_clean.reset_index()
    #check that these IDs are identical to those of selected_samples
    print(sample_list_from_hap_clean["ID_1"].equals(selected_samples))
    print(sample_list_from_hap_clean["ID_2"].equals(selected_samples))




    print_text("finishing the script", header=2)
    print_text("calculate the number of SNS removed due to the filters, i.e., the difference in the number of SNPs between the VCF file used as initial input in this script and the final VCF generated from where the hap file was created", header=3)
    run_bash(" \
        n_snps_before=$( \
            bcftools view \
                --no-header \
                " + input_vcf_file + " | \
            awk \
                'END{print NR}' \
        ); \
        n_snps_after=$( \
            bcftools view \
                --no-header \
                ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.only_snps_gen_pos.vcf.gz | \
            awk \
                'END{print NR}' \
        ); \
        n_lost_snps=$(echo $n_snps_before '-' $n_snps_after | bc); \
        echo 'During filtering, we have lost' $n_lost_snps 'out of' $n_snps_before 'SNPs'")
        #calculate the number of rows without header in the input VCF file and in the last VCF file generated in this script. Then calculate the difference and print, this is the number of SNPs lost.


    print_text("remove files not needed anymore", header=3)
    run_bash(" \
        rm --force " + input_vcf_file + "; \
        rm --force ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.vcf.gz; \
        rm --force ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_raw.vcf.gz; \
        rm --force ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched_plink_raw.vcf.gz; \
        rm --force ./results/01_cleaned_vep_vcf_files/" + selected_pop + "/chr" + selected_chromosome + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up." + selected_pop + ".cleaned.ref_alt_switched.vcf.gz")


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
print_text("parallelize", header=1)
print_text("create array with all combinations of pops and chromosomes", header=2)
print_text("get pop and chromosome names", header=3)
pop_names=unrelated_samples["pop"].unique()
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



print_text("run parallel analyses", header=2)
print_text("open the pool", header=3)
import multiprocessing as mp
pool = mp.Pool(60)

print_text("run function across pandas rows", header=3)
pool.map(master_processor, full_combinations_pop_chroms)

print_text("close the pool", header=3)
pool.close()




########################################################
#### Do some checks after analyzing all chromosomes ####
########################################################
print_text("Do some checks after analyzing all chromosomes", header=1)
print_text("check we do NOT have any errors in the output files of all chromosomes*pops combinations, also calculate the percentage of SNPs lost", header=2)
print_text("run loop across chromosomes*pops combinations", header=3)
snps_lost_percentage = []
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

    print_text("split the row, extract the numbers and calculate the percentage", header=4)
    row_results_split = row_results.split(" ")
    snps_lost_percentage.append(int(row_results_split[5])/int(row_results_split[8])*100)

print_text("see the percentiles of percentage of SNPs lost across all chromosomes and populations", header=3)
print("Do we have calculated all combinations")
print(len(snps_lost_percentage)==len(full_combinations_pop_chroms))

print("calculate the percentiles across combinations")
for i in [0.1,0.25,0.4,0.5,0.6,0.75,0.9]:
    print("Percentile " + str(i) + "%: " + str(np.quantile(snps_lost_percentage, i)))




####################
#### Next steps ####
####################
print_text("Next steps", header=1)
#if genetic distance is slow, you could do it in a previous step for each chromosome, and the new maps used as input here within each pop to remove snps without genetic distance from the VCF file and then hap file. Also remove snps from the map file not present in the VCF file, so we have the specific map for each population.
#run again  script in container in HPC and do check of the script in the meantime
#check number of snps lost
#once you are finished here, according to david, you can check whether the REF/ALT alleles match between the old and new hap files, but taking into account we have different coordinated, hg19 vs hg38
    ##i guess you could take the old map files, convert to hg38 coordinates and then see if the REF/ALT columns are the same than in the new map files
