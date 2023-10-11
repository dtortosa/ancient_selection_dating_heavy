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
                if($i ~/^HG/ || $i ~/^NA/){ \
                    n_samples++\
                } \
            }; \
            n_fields=NF; \
            printf \"n_fields=%s,n_samples=%s,chrom=%s,pos=%s,id=%s,ref=%s,alt=%s,filter=%s,info=%s\", n_fields, n_samples, chrom_index, pos_index, id_index, ref_index, alt_index, filter_index, info_index \
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
print("total number of fields: " + n_fields)
print("number of samples: " + n_samples)
print("index of column CHROM: " + index_chrom)
print("index of column POS: " + index_pos)
print("index of column ID: " + index_id)
print("index of column REF: " + index_ref)
print("index of column ALT: " + index_alt)
print("index of column FILTER: " + index_filter)
print("index of column INFO: " + index_info)
print("the total number of fields minus the number of samples should be 9. The number of fixed fields in VCF v4.2 is 8 and then FORMAT, which in our case only has GT, thus we should have 9 fields. Also, the index of CHROM, POS, REF, ALT and FILTER should be 1, 2, 4, 5 and 7, respectively")
if (int(n_fields)-int(n_samples) == 9) & (index_chrom=="1") & (index_pos=="2") & (index_id=="3") & (index_ref=="4") & (index_alt=="5") & (index_filter=="7") & (index_info=="8"):
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
print_text("see first cases where REF nor ALT are AA. You can see how we get a case with REF=G and AA_upcase=C. These could be cases where a multiallelic SNP has lost one of the ALT alleles in the subset population and that ALT allele is the ancestral. Therefore, there is no more ancestral in the population. Note, however, that I have found 200K cases like this across all chromosomes without subsetting when running 01b_vep_ancestral.py. Therefore, these could be also errors and should not be a high number", header=4)
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

print("create a new vcf file with a new case with AA=a and ALT=A, to check behavior with lower-case. Also add two cases with AA equals to ',' and '' to check behaviour", header=4)
run_bash(" \
    awk \
        'BEGIN{FS=OFS=\"\t\"}{print $0}END{ \
            print \"chr20\t1110691\trsdummy\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=a;AA_upcase=A\tGT\t1|0\t1|1\t0|0\"; \
            print \"chr20\t1110692\trsdummy2\tG\tA\t29\tPASS\tNS=3;DP=14;AN=6;AC=3;AF=0.5;DB;H2;AA=,;AA_upcase=A\tGT\t1|0\t1|1\t0|0\" \
        }'\
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf > ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf; \
        bcftools view \
            --no-header \
            ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2.vcf")
    #print all lines of the VCF file
    #add at the end a line with the main columns (ALT, REF, CHROM, genotypes...) separated by tabs and indicate that with FS

print_text("create a awk script that switch REF/ALT in those rows where AA==ALT, while not doing anything to rows where REF==AA and discarding those rows where REF nor ALT are equal to AA. First apply it to all variants without filter to check behaviour", header=4)
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
            index_info=" + index_info + " \
        }{ \
            if(NR > last_header_row){ \
                if(NR==(last_header_row+1)){ \
                    if($1 != \"chr20\"){exit 1} \
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
selected_chrom_dummy="chr20"
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
            index_info=" + index_info +"; \
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
    #this is a very good check because we are calculating problematic cases without AA_upcase, so at the end, we can check if REF is always equals to AA_upcase, which was obtained using other approach.

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

print_text("see the new header", header=4)
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

print_text("see the variants", header=4)
run_bash(" \
    bcftools view \
        --no-header \
        ./data/dummy_vcf_files/01_cleaned_dummy_vcf_files_vep/dummy_example_vep_2_anc_up_2_cleaned_ref_alt_switched.vcf")

print_text("calculate stats of the VCF file, show them here and then use them to make summary plots", header=4)
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

print_text("convert VCF file to hap file", header=4)
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

##por aqui

#SUMMARY: 
    #With all these commands, we have recreated the scenario we have in 1KGP data, with multiallelic SNPs separated into different lines, select some samples, we then select snps, exclude those SNPs that have the same allele for all samples (considering alleles in genotypes with missing, e.g., 0|.), remove those with genotype missingness < 5%, remove exact duplicates (this does not touch different lines of the same multiallelic snp because they have different ALT). Then we combine all lines of each multiallelic snp and now they have ALT column with several alleles, so we can filter them using --max-alleles 2. Add filter for selecting phased data only. Select only those variants included in interest regions (mask). We also use bcftools +fill-tags to update important fields for each SNP, so if a SNP was multiallelic, but it is not multiallelic in the subset population (i.e., only REF and 1 ALT), we no longer will have two allele frequencies, two allele counts.... for the remainder biallelic SNP in the subset. although in 1KGP data, multiallelic SNPs are already separated and have only 1 value for these fields (see below).

#Note about the update of the INFO fields
    #it is important to be sure that the fields you are using for filtering, are updated after subseting samples. Of course, type="snp" will be always "snp" irrespectively of the samples we select, but this is not the case of the number of alleles, because you can have SNPs with 3 alleles considering all 26 populations, but then in GBR they can have only 2 or 1. We are interested in SNPs that are biallelic within the selected population.
    #The same goes for phasing and genotype ^miss. You have to be sure that these filters only consider data from the filtered samples, not fixed data in fields that are not updated.
    #Because of this we have applied all these filters in order.



################################################################
#### function to clean vcf files and create hap - map files ####
################################################################

#chr_pop_combination="GBR_1"
def master_processor(chr_pop_combination):

    #extract selected population and chromosome
    selected_pop = chr_pop_combination.split("_")[0]
    selected_chromosome = chr_pop_combination.split("_")[1]

    #redirect standard output
    import sys
    original_stdout = sys.stdout
        #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
    sys.stdout = open("./scripts/01_hap_map_calcs_outputs/chr" + selected_chromosome + "_" + selected_pop + ".out", "w")
        #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/


    #######################
    # one time operations #
    #######################

    #do first some operations that need to be done just one time per chromosome
    if selected_pop == "GBR":
        
        #see file format
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": see VCF file version")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools head \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            grep -i '^##fileformat'")
                #use bcftools to see the header and then select the row starting with ##fileformat to see the version.
                #https://www.htslib.org/howtos/headers.html

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

        #see samples
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": see first 10 samples")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools query \
                --list-samples \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            head -10")
            #https://samtools.github.io/bcftools/howtos/query.html

        #check
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": do we have as samples as the total number of samples considering unrelated individuals and trios-duos?")
        print("#######################################\n#######################################")
        run_bash(" \
            n_samples=$( \
                bcftools query \
                    --list-samples \
                    " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                    wc -l); \
            if [[ $n_samples -eq " + str(samples_pedigree.shape[0]) + " ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi")
            #list the samples and count them
            #if the number of samples is equal to the number of samples in the pedigree, perfect

        #check
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": the chromosome name of the first variant is correct?")
        print("#######################################\n#######################################")
        run_bash(" \
            chr_vcf=$( \
                bcftools query \
                    " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
                    --format '%CHROM\n' | \
                head -1); \
            if [[ $chr_vcf == 'chr" + selected_chromosome + "' ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi")
            #with bcftools, query for the chromosome of the first variant
            #we should have only case, being equal to the selected chromosome
            #equality for strings using 1 bracket is "="
                #https://stackoverflow.com/questions/18102454/why-am-i-getting-an-unexpected-operator-error-in-bash-string-equality-test

        #inspect
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": show the variant type, ID, chromosome, position, alleles and frequency for the first snps")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools query \
                --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF\n' \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                head -5")
            #select the format of the query indicating the columns you want to show per SNP.
                #you can include data from INFO
                #end with \n to have different lines per SNPs

        #inspect
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": all variants has '.' for filter?")
        print("#######################################\n#######################################")
        run_bash(" \
            uniq_filters=$( \
                bcftools query \
                    --format '%FILTER\n' \
                    " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                uniq); \
            if [[ $uniq_filters == '.' ]]; then \
                echo 'True'; \
            else \
                echo 'False';\
            fi")
            #according to the specific readme of the dataset we are using (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf), prior phasing they applied several filters, being one of them that all variants has to be PASS for FILTER. I understand that, because of this, all variants in the chromosome have now ".", being this the unique character. Let's check this for all chromosomes:
                #make a query asking for the FILTER value of all variants
                #select unique cases
                #the unique cases should be a single string with a dot ("."), if yes, perfect.

        #stats with and without accessibility mask. See dummy example about details of the masks
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": See the stats of the whole VCF file with and without selecting SNPs inside the regions that PASS in the pilot and the strict mask")
        print("#######################################\n#######################################")
        print("\n##### NO MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
        print("\n##### PILOT MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                --targets-file ./data/masks/20160622.allChr.pilot_mask.bed \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
        print("\n##### STRICT MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                --targets-file ./data/masks/20160622.allChr.mask.bed.gz \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
                #load the whole VCf file of the selected chromosome, select only SNPs and in one case select variants inside mask "PASS" regions (pilot and strict), while in the second do nothing more. Then see the stats of the resulting BCF files. See dummy examples for further details.



    #######################
    # extract samples IDs #
    #######################

    #select the sample IDs for the selected population
    subset_pop = unrelated_samples.loc[unrelated_samples["pop"] == selected_pop, :]

    #reset index
    subset_pop = subset_pop.reset_index(
        drop=True)
        #drop=True to avoid adding the index as a column

    #check
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we only have the selected pop")
    print("#######################################\n#######################################")
    print(subset_pop["pop"].unique() == selected_pop)
    #stop if False in this check
    if subset_pop["pop"].unique() != selected_pop:
        raise ValueError("SERIOUS ERROR! WE HAVE NOT CORRECTLY SELECTED THE SAMPLES OF POP '" + selected_pop + "' in chromosome '" + selected_chromosome + "'")

    #select the sample IDs
    selected_samples = subset_pop["sample"]

    #save as txt to use it later
    #it is redundant to do it in each chromosome of the same population (same samples) but I do it anyway to avoid confusion between chromosomes
    selected_samples.to_csv(
        "results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt", 
        sep="\t", 
        header=None, 
        index=False)



    ##########################################
    # calculate hap file within selected pop #
    ##########################################

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see genotypes of selected samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #select a few samples and then select biallelic snps
                #IMPORTANT: If filtering (by number of alleles) and subseting (by sample) is in the same line (command), the filtering will be done first. Therefore, you could select SNPs that have 2 alleles when considering the 26 pops, but that are monomorphic (1 allele) for the selected population. Because of this, we have to first subset by sample and then filter by number of alleles whitin the selected samples in separated commands (see dummy example).
            #query multiple fields for each snp

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": show now only SNPs")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #include only those variants with TYPE=snp

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see more SNPs and their allele count to see how this data is presented for multiallelic SNPs")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AN %INFO/AC %INFO/AF GTs:[ %GT]\n' | \
        head -12")
            #IMPORTANT:
                #It seems that multiallelic variants separated in different lines have already updated the AC field. Therefore, each line does not have two allele counts, but one.
                #For example, 1:10453:A:C and 1:10452:A:C are both in the same position (10452) and have the same reference allele. They seem to be multiallelic, but each one has only one allele count, which is 1. 
                #This count is correctly updated for the subset of samples, but remember that --samples does not remove the second count. As we saw in the dummy example, we have to update the field with +fill-tags.
                #My hypothesis is that they updated this field using +fill-tags because they indeed say in the paper that they used bcftools to split the multiallelic SNPs in different lines.
                #Note that AF has also 1 value but it is not correct because we subset samples with --samples, and this command only updates AC and AN.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see monomorphic")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #update the INFO/AC and INFO/AN fields so we avoid two allele counts in the different lines of multiallelic snps. --samples update AC but maintaining the connection between the different lines of the same multialllelic SNP.
                #it seems that 1KGP data has already updated the AC field for each line separated line of a multiallelic SNP.
                #Therefore, this line would not be necessary, but we are applying just in case, to ensure we have only 1 allele count per line.
                #updating again the AC fields would not do anything wrong, just adding the same value that is was.
            #select those variants for which the number of ALT alleles (allele count) is equal to 0 (no ALT at all) or equal to the total number of alleles (AN), i.e., all alleles are ALT.
            #See dummy example for further details.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": exclude monomorphic")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we do not have SNPs with genotype missingness > 0.05. If TRUE, then we do not have to apply further filters about missing genotypes")
    print("#######################################\n#######################################")
    check_missing = run_bash(" \
        n_snps_with_filter=$( \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools view \
                --types snps | \
            bcftools +fill-tags \
                -- --tags AN,AC | \
            bcftools view \
                --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
            bcftools view \
                --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
            bcftools query \
                --format '%ID\n' | \
            wc -l); \
        n_snps_without_filter=$( \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools view \
                --types snps | \
            bcftools +fill-tags \
                -- --tags AN,AC | \
            bcftools view \
                --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
            bcftools query \
                --format '%ID\n' | \
            wc -l); \
        if [[ $n_snps_with_filter -eq $n_snps_without_filter ]];then \
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove also exact duplicates")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. See dummy example for behaviour.

    #combine multiallelic SNPs in one line and select them
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": combine multiallelic SNPs in one line and select them")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --min-alleles 3 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #we have SNPs with the same position and chromosome but different ALT alleles. According to the 1KGP paper, they split multiallelic SNPs in different lines, so this makes sense. See dummy for further details.
            #combine the different lines of each multiallelic SNP into one line and update the ALT column to include the different ALT alleles and we can filter with --max-alleles
                #two snps in position 10452 with same REF but different ALT in chromosome 1. They get combined into one line. These are the first multiallelic that appear in the previous command, but they were separated.

    #now show only biallelic snps
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs that are phased for all samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")



    ##por aqui checking the code
    #this was done in may, maybe a good idea to start from the begining
    #added pop name to bed file cleaned, check that the mask is correctly called in the next lines
    #error with "list_snps_with_gen_pos.txt"
        #this should have the chromosome and pop name to avoid interefernece
        #check that all files are correctly named

    #be sure to use the vcf files with the AA field


    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": clean BED file before applying the mask selecting only intervals in the selected chromosome")
    print("#######################################\n#######################################")
    run_bash(" \
        awk \
            -F '\t' \
            '{if ($1 == \"chr" + selected_chromosome + "\") print $0}' \
            ./data/masks/20160622.allChr.pilot_mask.bed \
            > ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed; \
        gzip \
            --force \
            ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed; \
        gunzip \
            -c \
            ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        head -5")
            #in the bed file, select those rows for which the chromosome name (first column) is the selected chromosome, printing all fields for these rows. Save as a file and then compress. See the first 5 lines. See dummy examples for further details.

    #
    print("\n#######################################\n#######################################")
    print("check we have selected the correct chromosome and mask type")
    print("#######################################\n#######################################")
    run_bash(" \
        uniq_chrom=$(\
            gunzip \
                -c \
                ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
            awk \
                -F '\t' \
                '{print $1}' | \
            uniq); \
        uniq_mask_type=$(\
            gunzip \
                -c \
                ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
            awk \
                -F '\t' \
                '{print $4}' | \
            uniq); \
        if [[ $uniq_chrom == 'chr" + selected_chromosome + "' && $uniq_mask_type == 'pilot' ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress the bed file of the selected chromosome (previously created), then print the chromosome name (first column) for all intervals, i.e., rows, getting then the unique cases, then save as a variable. Do the same with the 4th column, i.e., the mask type. The first variable should be the selected chromosome and the second the select mask type, i.e., pilot.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": select only regions that are accessible to sequencing according to the accessibility mask")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #Use --targets-file to only select SNPs within the intervals defined by a BED file generated by 1kGDP, which is an accessibility mask. Therefore, we select SNPs that are included in regions accessible to sequencing.
            #see dummy example for further details.
            #I have visually inspected the first 120 SNPs without mask and then check what SNPs are retained after applying the mask and if that makes sense with the first intervals in the BED file for the pilot mask. This works PERFECTLY, selecting only SNPs within the range.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": update and add some fields using fill-tags")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n' | \
        head -7")
            #we sue +fill-tags to update fields that are not updated after subsetting like frequency of alternative allele and create some additional fields.
            #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed sue to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields. In addition, --samples maintains 2 allele counts in each line of a splitted multiallelic SNP. In  the dummy example we had to update these fields with +fill-tags to have only 1 count and be able to use this count to filter out monomorphic SNPs. I have checked that splitted multiallelic SNPs in the 1KGP have only 1 count, so I guess the authors updated with fill-tags, but I update before the monomorphc check just in case. See above.
            #see dummy example for further details.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove all previous INFO/Fields, retain only what we are interested in and show the result asking for ALL fields")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools annotate \
            --remove INFO,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --no-header | \
        head -7")
            #remove all INFO fields and all FORMAT fields (except GT) using annotate and then add the fields we are interested in.
            #see dummy example for further details.
            #use view with option --no-header to see all fields without the header.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see header after applying fill-tags")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools annotate \
            --remove INFO,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools head")
            #you can see that after the list of contigs, the only field shown before the list of my commands is just GT (phased genotypes) because I have removed all INFO and FORMAT fields with the exception of FORMAT/GT.
            #then we have all the commands I have run in order to subset individuals and filter SNPs.
            #see dummy example for further details.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": save the cleaned vcf as a compressed file so we have an historial of the changes made in the vcf file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
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
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools annotate \
            --remove INFO,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --output ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
            --output-type z \
            --compression-level 1")
                #after subseting and filtering, save as a compressed VCF file, selecting the option for best speed (see dummy example for further details)

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
            #we used 1 because we were working with simulations, we have for sure data for all samples and snps. Ask David what to do in this case.

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
            #Also I guess I can remove SNPs that all their genotypes are phased within the selected population.
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see the header of the recently created vcf file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools head \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see the genotypes of a few individuals from the recently created vcf file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
            --no-header | \
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



    ##########################################
    # calculate map file within selected pop #
    ##########################################

    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": STARTING MAP FILE CALCULATION")
    print("#######################################\n#######################################")
        #we are going to calculate the map of each population directly using the SNPs in its vcf file FILTERED WITHIN population. Then, when we know for which of the SNPs we have genetic position, we can further filter the VCF file and convert to hap.
        #an alternative would be just take all the SNPs in the raw VCF file and calculate their genetic position.
            #it should be ok regarding the ID of the SNPs because REF/ALT are split when using multiallelics, and these fields are NOT switched based on the frequency of the SNPs in specific subsets.
        #I am going for the first option just to be completely sure I am using the SNPs (and positions) of the selected population.
            #If it is too slow this option, think about the other one.

    #Instructions david
        #I understand that SNPs with genetic position are NO useful for any summary statistic, right? So I can safely remove these SNPs from the VCF and hap files right?
            #ok
        #format of ID
            #in the map file can I use the format "CHROM:POS_REF_ALT" for the ID? 
            #I think remember that the map files I originally got from you in the previous project (before decode2019 conversion) used as ID just the physical position. Not sure if there is any specific reason for doing that.
            #not asked, by irrelevant question
        #data format decode map
            #in the decode map, they say clearly that the data is aligned to hg38. Also they do not specify if the coordinates are 1 or 0-based so I assume they are 1-based, base 1 in the map is base 1 in the genome, to me this should be the default. 
            #1KGP data is also aligned to hg38 and is 1-based.
            #Therefore I can just use the position of the SNPs in 1KGP to calculate their genetic position in the decode map, right?
                #David said that there is not problem if the map and hap files are in the same format. 
                #that is the case because the map file is calculated with decode map, and hap file with 1KGP VCF file, and as I said, I have no reason to think that the decode map is not 1-based.


    ##extract the SNP positions
    #extract snps from the cleaned VCF file
    run_bash(" \
        bcftools view \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools query \
            --format '%CHROM %ID %POS %REF %ALT\n' | \
        gzip \
            --force \
        > ./results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz")
            #from the cleaned VCF file, extract the chromosome, position, REF/ALT to compare with the positions in the raw hap file (see below), compress and save as a file

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


    ##load the raw_map file
    snp_map_raw = pd.read_csv(\
        "./results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz", \
        sep=" ", \
        header=None)
    print("See raw map file loaded in python")
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the correct number of SNPs in the map file loaded in python")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
            wc -l); \
        if [[ $n_snps -eq " + str(snp_map_raw.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    #rename the columns
    snp_map_raw = snp_map_raw.rename(\
        {0: "chr", 1: "id_old", 2: "pos", 3: "ref", 4: "alt"}, \
        axis=1)
            #we name ID as old because this is the ID coming from the VCF file, which we need for selecting those variants in the VCF file with genetic position. The final ID will be in plink format, see below.
            #use a dict with old and new column names. indicated we are renaming columns (axis=1)
    print("see map file with renamed columns")
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": create a new ID variable following plink 2.0 format and check it was correctly created")
    print("#######################################\n#######################################")
    #the original ID column is in the format "CHR:POS:REF:ALT" but the problem is that some SNPs have their position indicated in the ID is shifted. 1KGP authors separated multiallelic SNPs in different lines with bcftools norm and then shifted their position so they could be phased, combing back to the original position afterwards (see next line). I guess during that process, they updated the IDs using chrom and the shifted position was used in the ID. Therefore, even POS comes back to the original position, the ID remains with the shifted position. I have checked several multiallelic SNPs, and they have all the same issue with the position in the ID.
        #From README ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf"): "SHAPEIT2 does not handle multiallelic variant phasing. To phase both biallelic and multiallelic variants we first split the multiallelics into separate rows while left-aligning and normalizing INDELs using bcftools norm tool (Li, 2011). Next, we shifted the position of multiallelic variants (2nd, 3rd, etc ALT alleles) by 1 or more bp (depending on how many ALT alleles there are at a given position) to ensure a unique start position for all variants, which is required for SHAPEIT2. We shifted the positions back to the original ones after phasing".
    #This is not very convenient because when later we create the hap files, SNPs will have the plink format for ID to avoid avoid strand flips (CHR:POS_REF_ALT), so we are going to have different IDs between hap and map files, making more difficult to do checks.
    #Therefore, we are going to update the ID of each SNP using plink format, and ensuring in this way SNPs will be names the same in hap and map files.
    snp_map_raw["id"] = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
    #check
    check_id = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
        #make a series combining chromosome, pos, ref and alt, and using the corresponding separators
    print(check_id.equals(snp_map_raw["id"]))
        #check it is identical to id
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the ref/alt columns as we have this information already included in the ID")
    print("#######################################\n#######################################")
    snp_map_raw = snp_map_raw.drop(["ref", "alt"], axis=1)
    print(snp_map_raw)


    ##load and explore the decode2019 map

    #I know that the original 2019 decode map is alligned to hg38. Also, I assume that the decode2019 map is 1-based because they do not specify is 0-based. I assume that if you say anything, base 1 in your coordinates is base 1 in the genome. I assume this is the default.
        #Data S3.genetic.map.final.sexavg.gor.gz:
            #average genetic map computed from the paternal and maternal genetic maps, which were in turn computed from the paternal and maternal crossover, respectively. The data columns are as follows: Chr (chromosome), Begin (start point position of interval in GRCh38 coordinates), End (end point position of interval in GRCh38 coordinates), cMperMb (recombination rate in interval), cM (centiMorgan location of END POINT of interval)
            #Page 85 of "aau1043-halldorsson-sm-revision1.pdf"

    #1KGP is aligned to hg38 (see paper) and coordinates are 1-based as VCF format 4.2 has 1-based coordinates.

    #therefore, we have the same position format in both datasets, so we can just use the decode 2019 map to calculate the genetic position of each SNP.

    ##IMPORTANT
    #check if the data is 1 or 0 based
        #DAVID: according to DAvid, you can check that with the assembly file, if the base the map is saying is in position X is indeed at position X, you are good.



    # 
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first lines of the Data S3 of decode paper, which is the sex average map (see above). The file has header")
    print("#######################################\n#######################################")
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        head -20")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the header and save")
    print("#######################################\n#######################################")
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        awk \
            'NR>7' | \
        gzip \
            --force > \
        ./data/decode_2019/aau1043_datas3_no_header.gz; \
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3_no_header.gz | \
        head -5")
            #decompress, send to stdout
            #then select any row whose number is larger than 7, i.e., from 8th row and forward.
                #https://www.baeldung.com/linux/remove-first-line-text-file 
            #compress and save
            #decompress and see first lines

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": load decode 2019 map into python")
    print("#######################################\n#######################################") 
    decode2019_map = pd.read_csv(\
        "./data/decode_2019/aau1043_datas3_no_header.gz", \
        sep="\t", \
        header=0, \
        low_memory=False)
    #rename columns in lower case
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": subset decode map for the selected chromosome")
    print("#######################################\n#######################################") 
    decode2019_map_subset = decode2019_map.loc[decode2019_map["chr"] == "chr"+str(selected_chromosome),:]
    print(decode2019_map_subset)
    print("Do we selected the correct chromosome?")
    print(decode2019_map_subset["chr"].unique() == "chr"+str(selected_chromosome))
    #
    print("remove the full decode map")
    del(decode2019_map)
    import gc
    gc.collect()

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": define function to calculate genetic position per SNP")
    print("#######################################\n#######################################") 
    #selected_snp_id=snp_map_raw.iloc[10000]["id"] #snp with cM value just in its position
    #selected_snp_id=snp_map_raw.iloc[5000]["id"] #snp with cM values at both sides
    #selected_snp_id=snp_map_raw.iloc[0]["id"] #snp without cM values around
    def gen_pos(selected_snp_id):

        #extract the row in the raw map for the selected SNP
        selected_snp_row = snp_map_raw.loc[snp_map_raw["id"] == selected_snp_id,:]

        #check we have the correct chromosome
        check_0 = selected_snp_row["chr"].unique()[0] == "chr"+str(selected_chromosome)

        #extract position of the selected snp
        selected_snp_physical_pos = selected_snp_row["pos"].to_numpy()[0]

        #extract old ID (this follows VFP file format so we can use them to filter it)
        selected_snp_old_id = selected_snp_row["id_old"].to_numpy()[0]

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
            check_1=np.unique(\
                (decode2019_map_subset_around_snp["end"] >= (selected_snp_physical_pos - 10**6)) & \
                (decode2019_map_subset_around_snp["end"] <= (selected_snp_physical_pos + 10**6)))[0]

            #if we dot NOT have an interval with an end coordinate exactly similar to the selected SNP
            if (interval_same_pos.shape[0] == 0):

                #check
                check_2 = np.unique(intervals_lower_end["end"] < selected_snp_physical_pos)[0]
                check_3 = np.unique(intervals_upper_end["end"] > selected_snp_physical_pos)[0]

                #from the intervals below the extreme window, select the biggest and hence closest to the extreme window   
                lowest_interval = intervals_lower_end.loc[intervals_lower_end["end"] == max(intervals_lower_end["end"]),:] 
                    #we cannot have two cases with the same value because the coordinates are in increasing order, the coordinate of an interval is bigger than the previous one.

                #from the intervals above the extreme window, select the smallest and hence closest to the extreme window
                highest_interval = intervals_upper_end.loc[intervals_upper_end["end"] == min(intervals_upper_end["end"]),:] 
                    #we cannot have two cases with the same value because the coordinates are in increasing order, the coordinate of an interval is bigger than the previous one.

                #check that the end coordinate with lowest difference respect the SNP is the selected in the previous step both for the lower and higher intervals
                check_4a = (intervals_lower_end.loc[\
                    np.abs(intervals_lower_end["end"]-selected_snp_physical_pos) == \
                    np.min(np.abs(intervals_lower_end["end"]-selected_snp_physical_pos)), \
                    "end"] == lowest_interval["end"]).to_numpy()[0]
                check_4b = (intervals_upper_end.loc[\
                    np.abs(intervals_upper_end["end"]-selected_snp_physical_pos) == \
                    np.min(np.abs(intervals_upper_end["end"]-selected_snp_physical_pos)), \
                    "end"] == highest_interval["end"]).to_numpy()[0]


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
                    #My explanation: What David is doing is 100 + ((102-100)*20)/(20+80). This gives exactly 100.4. David is using the rule of three (https://en.wikipedia.org/wiki/Cross-multiplication#Rule_of_Three). You have three points, A, B and C. If the physical distance distance A-C is 100 kb (20+80) and the genetic distance between these points is 2 cM (102-100) , what would be the genetic distance between A-B if these points are separated by 20 kb? ((102-100 cM) * 20 kb) / (20+80 kb); ((2 cM) * 20 kb) / (100 kb); ((2 cM) * 20 kb) / (100 kb); (40 cM * kb) / 100 kb; 0.4 cM. 0.4 is the genetic distance between A and B. Now we can sum 0.4 and the genetic position of A, to get the genetic position of B in the genome. 100 cM + 0.4 cM = 100.4 cM.
                    #If you the point for which you calculate the genetic distance is exactly in the middle of the two points with 100 and 102 cM of genetic distance, the resulting genetic distance would be exactly in the middle, i.e., 101: 100 + ((102-100)*50)/(50+50). 50 is the physical distance between cM point and the point of interest. 
                    #This method assumes that relationship between genetic distance and physical distance between two points is lineal and stable, so you can estimate the genetic distance based on the physical distance in the genomic region encompassed by these points. Note that you are using point that are at least 1MB close to the point under study, therefore, we are estimating the genetic distance using the relationship between physical and genetic distance in a specific genomic region, not the whole genome.
                    #see figure 31 for further details.
            else:

                #if not and hence we have an deCODE interval exactly in the SNP position
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": Run the function on just one snp")
    print("#######################################\n#######################################")
    print(gen_pos(snp_map_raw.iloc[5000]["id"]))

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": run function across SNPs")
    print("#######################################\n#######################################")
    #we do not use pool.map because we will only want to use 1 core. The parallelization will be done in the parent function across chromosome*pop combinations
    #map seems to be faster than loop even using just 1 core, although there is not a big difference
        #https://www.linkedin.com/pulse/loops-maps-who-faster-time-space-complexity-we-coming-george-michelon/
    final_genetic_pos = list(map(gen_pos, snp_map_raw["id"]))
    #final_genetic_pos = list(map(gen_pos, snp_map_raw.iloc[5000:5100]["id"]))

    #convert the tuple to DF and add the column names
    final_genetic_pos_df = pd.DataFrame(final_genetic_pos, columns=["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "selected_snp_physical_pos", "check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6", "genetic_distance", "left_cM", "right_cM", "distance_left_end", "distance_right_end"])
    print("see results:")
    print(final_genetic_pos_df)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": all checks of genetic position calculation are True?")
    print("#######################################\n#######################################") 
    print(final_genetic_pos_df[["check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6"]].all())
        #important:
            #all() does not consider nan, so if you have nan and the rest True, the output is True.
            #this is ok for us, because we use nan in some conditions.
    
    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the correct number of SNPs in the calculation of genetic position")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
            wc -l); \
        if [[ $n_snps -eq " + str(final_genetic_pos_df.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that we have the exact same snps than in the raw map")
    print("#######################################\n#######################################")
    print(np.array_equal(
        snp_map_raw["chr"].to_numpy(),
        ("chr" + final_genetic_pos_df["selected_chromosome"]).to_numpy()))
    print(np.array_equal(
        snp_map_raw["id"].to_numpy(),
        final_genetic_pos_df["selected_snp_id"].to_numpy()))
    print(np.array_equal(
        snp_map_raw["id_old"].to_numpy(),
        final_genetic_pos_df["selected_snp_old_id"].to_numpy()))
    print(np.array_equal(
        snp_map_raw["pos"].to_numpy(),
        final_genetic_pos_df["selected_snp_physical_pos"].to_numpy()))
 
    
    ##recalculate genetic distance of each SNP
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": compare the new genetic distance and the distance previously calculated. It is ok to have False here if the next check is ok")
    print("#######################################\n#######################################")
    raw_check_gen_dis = new_genetic_distance == final_genetic_pos_df_check_dist_calc["genetic_distance"]
    print(raw_check_gen_dis.groupby(raw_check_gen_dis).count())
        
    #
    #from results, extract ID of snps with NA for last checks but with data for genetic position
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that the cases with NA for last checks but with genetic distance are the cases of SNPs in a position exactly with deCODe data")
    print("#######################################\n#######################################")
    cases_gen_pos_no_last_checks = final_genetic_pos_df.loc[\
        (final_genetic_pos_df["check_4a"].isna()) & \
        (~final_genetic_pos_df["genetic_distance"].isna()), "selected_snp_id"]
    #extract the ID of snps with a position that have deCODE genetic position
    snps_with_decode_data = final_genetic_pos_df.loc[\
        final_genetic_pos_df["selected_snp_physical_pos"].isin(decode2019_map_subset["end"]), "selected_snp_id"]
    #print the check
    print(cases_gen_pos_no_last_checks.equals(snps_with_decode_data))


    ##final map file
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": prepare final map file")
    print("#######################################\n#######################################")
    #subset only the columns for map files
    final_genetic_pos_map_file = final_genetic_pos_df[["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "genetic_distance", "selected_snp_physical_pos"]]
        
    #add chrom
    final_genetic_pos_map_file["selected_chromosome"] = "chr"+final_genetic_pos_map_file["selected_chromosome"]
        #WARNING HERE
    print(final_genetic_pos_map_file)
        #save the chromosome, ID, genetic position and physical position. This is the format expected by hapbin
            #https://github.com/evotools/hapbin

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the SNPs without genetic position from the VCF file")
    print("#######################################\n#######################################")

    #select old_id from the map that have genetic position
    #this ID is the original retained from the VCF file, so we can use it to subset the VCF file
    snps_id_with_gen_pos = final_genetic_pos_map_file.loc[\
        ~final_genetic_pos_map_file["genetic_distance"].isna(), \
        "selected_snp_old_id"]

    #
    print("see the number of SNPs removed due to the lack of genetic position")
    print(final_genetic_pos_map_file.shape[0] - len(snps_id_with_gen_pos))

    #save the names in a txt file
    with open(r"./results/01_cleaned_vep_vcf_files/list_snps_with_gen_pos.txt", "w") as fp:
        fp.write("\n".join(snps_id_with_gen_pos))
            #each name in a different line so we have to add "\n" to the name
            #https://pynative.com/python-write-list-to-file/
        fp.write("\n")
            #add empty line at the end

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": filter the already cleaned VCF with bcftools")
    print("#######################################\n#######################################")
    #this file is cleaned regarding biallelic snps, duplicates... but need to retain only SNPs with genetic position
    #We could do this by just creating before the hap file, extract snp positions from there, calculate genetic position and then remove from the hap those rows of SNPs without genetic position. The problem is that we would do that by row index instead of SNP ID, at least if we use the final hap file, so we are going for this option better. In addition, we would have snps that cannot be used in the VCF file because they do not have genetic position. With the other approach we would have a VCF with all SNPs filtered and another one with only snps with genetic position.
    
    #filter
    run_bash("\
        bcftools view \
            --include ID==@./results/01_cleaned_vep_vcf_files/list_snps_with_gen_pos.txt\
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools view \
            --output ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --output-type z \
            --compression-level 1")
            #include those SNPs for which ID is included in the list of SNPs with genetic position and save the resulting VCF file
                #https://www.biostars.org/p/373852/

    #
    print("see header of the fully filtered VCF file and some genotypes")
    run_bash(" \
        bcftools head \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz")
    run_bash(" \
        bcftools view \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --no-header | \
        head -5")

    #
    print("check that IDs in the filtered VCF file are the same than the ones in the list of IDs used as input to filter")
    run_bash(" \
        bcftools query \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --format '%ID\n' \
        > ./results/01_cleaned_vep_vcf_files/ids_vcf_after_filter.txt; \
        file1='./results/01_cleaned_vep_vcf_files/list_snps_with_gen_pos.txt'; \
        file2='./results/01_cleaned_vep_vcf_files/ids_vcf_after_filter.txt'; \
        STATUS=$(cmp --silent $file1 $file2; echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi; \
        rm $file2")
        #get the IDs in the finally filtered VCF file and save the file
        #create two variables with the names of this file and also the name of the file with the list of IDs used as input to filter        
        #check byte by byte whether the two files are the same
            #cmp takes two files and compare them until 1 byte is different
            #we make it silent and get the final status
            #remember that "$?" gives the return value of the last run command.
                #For example, 
                    #ls somefile
                    #echo $?
                    #If somefile exists (regardless whether it is a file or directory), you will get the return value thrown by the ls command, which should be 0 (default "success" return value). If it doesn't exist, you should get a number other then 0. The exact number depends on the program.
                #https://stackoverflow.com/a/6834572/12772630
            #the return value of cmp will be 0 if the two files are identical, if not, then we have differences between the files
                #https://stackoverflow.com/a/53529649/12772630
        #remove the file created for this check

    #
    print("remove also the SNPs without genetic position from the map file and check")
    final_genetic_pos_map_file = final_genetic_pos_map_file.loc[\
        ~final_genetic_pos_map_file["genetic_distance"].isna(),:]
    print(final_genetic_pos_map_file["selected_snp_old_id"].equals(snps_id_with_gen_pos))

    #
    print("remove old ID as we have already filtered the VCF file and check")
    final_genetic_pos_map_file = final_genetic_pos_map_file.drop(["selected_snp_old_id"], axis=1)
    print(final_genetic_pos_map_file.columns == ["selected_chromosome", "selected_snp_id", "genetic_distance", "selected_snp_physical_pos"])

    #
    print("see final map and save")
    print(final_genetic_pos_map_file)
        #required format according to hapbin
            #The map files (--map) should be in the same format as used by Selscan with one row per variant and four space-separated columns specifiying 
                #chromosome, 
                #locus ID, 
                #genetic position
                #physical position.
    final_genetic_pos_map_file.to_csv(\
        "./results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz", \
        sep=" ", \
        header=False, \
        index=False)

    #
    run_bash("\
        n_rows=$( \
            gunzip \
                --stdout \
                ./results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NR}'); \
        n_cols=$( \
            gunzip \
                --stdout \
                ./results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NF}'); \
        if [[ $n_cols -eq 4 && $n_rows -eq " + str(final_genetic_pos_map_file.shape[0]) + " ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress map file to stdout and then calculate the number of rows (NR) and fields (NF). Do that only after the whole file has been read (END)
                #https://www.gnu.org/software/gawk/manual/html_node/Using-BEGIN_002fEND.html
            #the number of columns (fields) should be 4 following salescan format, while the number of rows should be equal to the number of snps we have in the map file loaded in python, which was indeed used to write this .map file.


    ##convert to hap 
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": convert to hap file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools convert \
            ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --hapsample ./results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw")
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

        #The smallest number (852072) is the number of SNPs after we have completely cleaned the vcf file, including the accesibility mask. I have checked that.

    #see first variants
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first variants of the hap file")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c ./results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        head -3")
            #decompress the hap file and show in stdout 
    
    #clean
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove first columns of hap file to leave only haplotype columns")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        cut \
            --complement \
            --delimiter ' ' \
            --fields 1-5 \
        > results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap; \
        gzip -f results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap")
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

    #check
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": the clean hap file is just the raw hap file but without the first 5 columns?")
    print("#######################################\n#######################################")
    #create a long string with the index of each column of hap_raw file (without 5 first columns)
    #these indexes will be used by awk to select the corresponding columns so we are creating the cleaned hap file again using another approach
    fields_selected_samples = "".join(
        ["$" + str(i) + "," if i != selected_samples.shape[0]*2+5 else "$" + str(i) for i in range(6, selected_samples.shape[0]*2+5+1, 1)])
        #from index 6 (avoiding non-genotype columns) to the index of the last column, i.e., 2 columns times the number of samples plus 5 (because of the non-genotype columns) and 1 (because index 1 in python is 0, so if you do range from 0 to 10, you get until 9, not 10, you need to add 1 (i.e., 11) to get the last one)
        #if the index is NOT the last one
            #add $ to the index and then comma
        #else
            #it is the last one so we do not need add comma
        #join all strings
    #do comparison
    run_bash(" \
        gunzip \
            -c results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        awk -F ' ' '{print " + fields_selected_samples + "}' > \
        results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap; \
        gunzip \
            -kf results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz; \
        file1='results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap'; \
        file2='results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap'; \
        STATUS=$(cmp --silent $file1 $file2; echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi; \
        rm $file1; \
        rm $file2")
        #create again the hap file cleaned
            #decompress the raw hap file and send it to stdout
            #select only the genotype columns using awk, which needs to know the delimiter is ' '
            #save it as a file
        #decompress the previously cleaned hap file, keeping the compressed file and forcing the decompression in case the decompressed file already exist (-kf)
            #https://linux.die.net/man/1/gunzip
        #create two variables with the names of these new files        
        #check byte by byte whether the two files are the same
            #cmp takes two files and compare them until 1 byte is different
            #we make it silent and get the final status
            #remember that "$?" gives the return value of the last run command.
                #For example, 
                    #ls somefile
                    #echo $?
                    #If somefile exists (regardless whether it is a file or directory), you will get the return value thrown by the ls command, which should be 0 (default "success" return value). If it doesn't exist, you should get a number other then 0. The exact number depends on the program.
                #https://stackoverflow.com/a/6834572/12772630
            #the return value of cmp will be 0 if the two files are identical, if not, then we have differences between the files
                #https://stackoverflow.com/a/53529649/12772630
        #remove the new files

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the same number of rows (variants) in the cleaned vcf, hap and map files")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps_vcf=$( \
            bcftools view \
                ./results/01_cleaned_vep_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
                --no-header | \
            wc -l); \
        n_snps_map=$( \
            gunzip \
                -c ./results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            wc -l); \
        n_snps_hap=$( \
            gunzip \
                -c results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            wc -l); \
        if [[ $n_snps_vcf -eq $n_snps_hap && $n_snps_vcf -eq $n_snps_map ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #load the cleaned vcf file with bcftools and show only the snps, with no header, then count the number of lines and save the result
            #decompress the hap file and count the number of line
            #check both numbers are the same

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": do we have the correct number of samples in the hap file?")
    print("#######################################\n#######################################")
    run_bash(" \
        nfields_hap=$( \
            gunzip \
                -c results/03_hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            awk -F ' ' '{print NF}' | \
            sort | \
            uniq); \
        nfields_hap_nlines=$( \
            echo -n $nfields_hap | \
            grep -c '^'); \
        if [[ $nfields_hap_nlines -eq 1 ]]; then \
            nsamples_hap=$(($nfields_hap/2)); \
            nsamples_bcftools=$( \
                cat results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt | \
                grep -c '^'); \
            if [[ $nsamples_hap -eq $nsamples_bcftools ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi\
        fi")
        #calculate the number of fields/columns in the hap file
            #decompress the hap file and send the output to stdout (flag c of gunzip)
            #use the interpreter of awk language (awk), indicating the delimiter of the columns (-F) and asking to print the number of fields. You could also ask for column 1 typing $1... NF will give the number of columns or fields.
                #https://www.unix.com/shell-programming-and-scripting/110272-no-columns-csv-file.html
                #https://askubuntu.com/questions/342842/what-does-this-command-mean-awk-f-print-4
            #you get the number of fields per row, so you have to sort the output and obtain the uniq cases.
            #save into nfields_hap
        #calculate how many number of unique number of columns we have. 
            #We should have the same number of columns in all rows. Therefore, just 1 line.
            #echo the nfields_hap with -n flag, so you do not output the trailing newline, i.e., there is no last empty line at the end. I guess people use this to avoid counting that last line. If you do echo '' | grep -c '^', without -n, you still get 1, when it should be zero.
                #https://unix.stackexchange.com/questions/579651/what-does-n-flags-stands-for-after-echo
            #get the number of lines that start with any character ("^") using grep -c "^". Remember that ^eso would look for any line starting with "eso"
                #https://stackoverflow.com/questions/6314679/in-bash-how-do-i-count-the-number-of-lines-in-a-variable
                #https://stackoverflow.com/questions/13054227/what-does-grep-mean-in-unix
        #if the number of unique number of columns is 1, we are fine so continue
            #take nfields_hap, divide by 2. Remember that we remove the initial 5 columns, so we only have haplo columns. In these columns, we have 2 genotype columns per sample, so the number of columns/2 is the number of samples.
                #https://phoenixnap.com/kb/bash-math
            #calculate the number of samples from the .txt file with the IDs of all samples.
                #load the file with cat and count the number of lines that start with any character
            #check that the number of samples according to the hap file and the txt are the same.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that sample file generated with hap has the correct sample IDs?")
    print("#######################################\n#######################################")
    #read the sample file
    sample_list_from_hap = pd.read_csv(
        "./results/02_hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.samples",
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

    #restore sys.stdout using the previously saved reference to it
    #This is useful if you intend to use stdout for other things
    sys.stdout = original_stdout
        #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/

    #por aquiii
        #run the script 
        #check the redirection of stdout
        #check the whole script in the meantime



#####################
#### paralellize ####
#####################

##
print("\n#######################################\n#######################################")
print("create array with all combinations of pops and chromosomes")
print("#######################################\n#######################################")
#get pop and chromosome names
pop_names
chromosomes = [i for i in range(1, 23, 1)]
print("we are going to analyze 26 pops and 22 chromosomes?")
print((len(pop_names) == 26) & (len(chromosomes) == 22))
print("See them")
print(pop_names)
print(chromosomes)

#get all the combinations but first make a dummy example
import itertools
#create two dummy lists, one with strings and other with integers
print("dummy example to get all possible combinations of two lists")
dummy_x = ["marbella", "cuzco", "granada"]
dummy_y = [1, 2, 3]
#product get all possible combinations between the two lists
dumm_combinations = [x+"_"+str(y) for x in dummy_x for y in dummy_y]
print(dumm_combinations)
    #first for each each value of X, and then for each value of Y, combine X and Y, so combine X1 with Y1, X1 with Y2, .... X2 with Y1, X2 with Y2 and so on...
    #y has to be converted to string with it is integer
print("Do we have all dummy combinations?")
print(len(dumm_combinations) == len(dummy_x)*len(dummy_y))
#get all combinations from the actual pops and chromosomes
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


##run parallel analyses
#open the pool
import multiprocessing as mp
pool = mp.Pool(len(full_combinations_pop_chroms)/2)

#run function across pandas rows
pool.map(master_processor, full_combinations_pop_chroms)

#close the pool
pool.close()

##por aqui
##mail de jesus
    #https://mail.google.com/mail/u/0/?tab=rm&ogbl#drafts/QgrcJHsHpDRJdfjndBxlCjQHdCNBwJJqNSl
    ## primero termina dummy example
        #extract again AA the first base
            #print substr($0,1,2)
            #https://stackoverflow.com/a/1406061/12772630
        #then check that this base is REF
    ##luego run VEP from jesus container on chrom22 and check we get the same results, same AA, same number of missing...
        #in this way we cna be sure we do not have problems avodiing fasta, smartmach error 
    ##think about mask, jesus says we do not need it



##when checking false/error also check you do not have the following warning in ANY pop-chrom:
    #THIS WARNING SHOULD BE ONLY IN DUMMY DATA NOT IN 1KGP DATA as AC field should have only 1 value per line in the data.


##according to david, you can check whether the REF/ALT alleles match between the old and new hap files, but taking into account we have different coordinated, hg19 vs hg38