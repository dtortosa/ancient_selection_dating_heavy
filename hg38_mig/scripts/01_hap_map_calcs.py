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

#We are going to migrate to hg38 as the new release of 1KGP matches this genome reference version.

#The general repo page for this release (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage) and general readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/README_111822.pdf).

#It says that "20220422_3202_phased_SNV_INDEL_SV" is "the most up-to-date version of the phased panel based on the high-coverage 1kGP WGS data which includes SNV, INDEL, and SV calls across 3,202 samples"

#High coverage and phased data can be found in the general repo, working folder and then 20220422_3202_phased_SNV_INDEL_SV (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/). This phased data includes single nucleotide variants, indels, SV in vcf files per chromosome. See its specific readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf). This is the repo used by Jesus Murga. 

#According to the general readme, the pedigree information is in the file "1kGP.3202_samples.pedigree_info.txt" (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt). This includes the ID of each sample.



##################
#### Starting ####
##################


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
        print("WARNING: THIS WARNING SHOULD BE ONLY IN DUMMY DATA NOT IN 1KGP DATA as AC field should have only 1 value per line in the data. " + complete_process.stderr)

        #return also the value if required
        if return_value==True:
            return complete_process.stdout
    elif ("Warning" in complete_process.stderr) | ("warning" in complete_process.stderr) | ("W:" in complete_process.stderr):

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #print the standard error without stopping
        print("WARNING: " + complete_process.stderr)

        #return also the value if required
        if return_value==True:
            return complete_process.stdout
    elif ("Lines   total/split/realigned/skipped" in complete_process.stderr) | ("no-ALT/non-biallelic/filtered" in complete_process.stderr) | ("Parsing bcftools stats output" in complete_process.stderr):

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #print the standard error, which is not an error
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
input_vcfs_path = "data/vcf_files_hg38"

#create folders to save the results
run_bash(" \
    mkdir \
        -p ./results/hap_map_files_raw; \
    mkdir \
        -p ./results/hap_map_files; \
    mkdir \
        -p ./results/cleaned_vcf_files")
    #-p: no error if the folder already exists, make parent directories as needed



#############
# pops prep #
#############

#load original 2504 unrelated samples from phase 3. This includes sample IDs and pop/superpop codes. This is the definitive ped we will use both for pop codes and sample IDs
#http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
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

#check
print("\n#######################################\n#######################################")
print("check that sample IDs in the original pedigree are all included in peds of high coverage dataset")
print("#######################################\n#######################################")
import numpy as np
print(np.unique(original_unrel_ped["sample"].isin(samples_pedigree["sampleID"])))
print(np.unique(original_unrel_ped["sample"].isin(samples_pedigree_pop["SampleID"])))

#split the new peds between samples that are included in the original 2504 set and not included in it
samples_pedigree_unrel = samples_pedigree.\
    loc[samples_pedigree["sampleID"].isin(original_unrel_ped["sample"]),:]
samples_pedigree_rel = samples_pedigree.\
    loc[~samples_pedigree["sampleID"].isin(original_unrel_ped["sample"]),:]
samples_pedigree_pop_unrel = samples_pedigree_pop.\
    loc[samples_pedigree_pop["SampleID"].isin(original_unrel_ped["sample"]),:]
samples_pedigree_pop_rel = samples_pedigree_pop.\
    loc[~samples_pedigree_pop["SampleID"].isin(original_unrel_ped["sample"]),:]

#
print("\n#######################################\n#######################################")
print("check that we have the same sex values in the three peds")
print("#######################################\n#######################################")
#sex="male"
for sex in original_unrel_ped["gender"].unique():

    #start with selected sex
    print("###### " + sex + " ######")

    #change the format of sex to that of the new peds
    if sex == "male":
        sex_new_peds=1
    elif sex == "female":
        sex_new_peds=2

    #get IDs of the selected sex in the three peds and reset index to let us compare with the other peds
    original_unrel_ped_sex_ids = original_unrel_ped.\
        loc[original_unrel_ped["gender"]==sex, "sample"].\
        reset_index(drop=True)
    samples_pedigree_unrel_sex_ids = samples_pedigree_unrel.\
        loc[samples_pedigree_unrel["sex"]==sex_new_peds, "sampleID"].\
        reset_index(drop=True)
    samples_pedigree_pop_unrel_sex_ids = samples_pedigree_pop_unrel.\
        loc[samples_pedigree_pop_unrel["Sex"]==sex_new_peds, "SampleID"].\
        reset_index(drop=True)

    #compare
    print(np.unique(original_unrel_ped_sex_ids == samples_pedigree_unrel_sex_ids))
    print(np.unique(original_unrel_ped_sex_ids == samples_pedigree_pop_unrel_sex_ids))

#
print("\n#######################################\n#######################################")
print("check that we have the same number of samples per pop")
print("#######################################\n#######################################")
#pop="GBR"
for pop in original_unrel_ped["pop"].unique():

    #start with selected sex
    print("###### " + pop + " ######")

    #get IDs of the sample for the selected pop in peds and reset index to let us compare with the other peds
    samples_old_ped = original_unrel_ped.\
        loc[original_unrel_ped["pop"]==pop, "sample"]\
        .reset_index(drop=True)
    samples_new_ped = samples_pedigree_pop_unrel.\
        loc[samples_pedigree_pop_unrel["Population"]==pop, "SampleID"].\
        reset_index(drop=True) 

    #compare
    print(np.unique(samples_old_ped == samples_new_ped))


##explore these peds
#"As part of this publication, we sequenced 3,202 lymphoblastoid cell line (LCL) samples from the 1kGP collection, including 1,598 males and 1,604 females."
    #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 1598 males in the latest ped?")
print("#######################################\n#######################################")
male_samples = samples_pedigree.loc[samples_pedigree["sex"] == 1, :]
print(male_samples.shape[0] == 1598)
print(male_samples)
print("\n#######################################\n#######################################")
print("Do we have 1604 females in the latest ped?")
print("#######################################\n#######################################")
female_samples = samples_pedigree.loc[samples_pedigree["sex"] == 2, :]
print(female_samples.shape[0] == 1604)
print(female_samples)
    #Sex is the same in the latest ped and the original ped (for shared samples), so we are good with sex in the original ped.

#see the number of samples per superpopulation
#The 3,202 samples were drawn from 26 populations (listed in Table S1) across the following 5 continental ancestry groups: African (AFR, n = 893), European (EUR, n = 633), East Asian (EAS, n = 601), South Asian (SAS, n = 585), and American (AMR, n = 490)
    #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
afr_samples = samples_pedigree_pop.\
    loc[samples_pedigree_pop["Superpopulation"] == "AFR", :]
eur_samples = samples_pedigree_pop.\
    loc[samples_pedigree_pop["Superpopulation"] == "EUR", :]
eas_samples = samples_pedigree_pop.\
    loc[samples_pedigree_pop["Superpopulation"] == "EAS", :]
sas_samples = samples_pedigree_pop.\
    loc[samples_pedigree_pop["Superpopulation"] == "SAS", :]
amr_samples = samples_pedigree_pop.\
    loc[samples_pedigree_pop["Superpopulation"] == "AMR", :]
#check
print("\n#######################################\n#######################################")
print("Do we have the correct number of samples per superpopulation in the high-coverage ped with pop names?")
print("#######################################\n#######################################")
print("## AFR: " + str(afr_samples.shape[0] == 893))
print("## EUR: " + str(eur_samples.shape[0] == 633))
print("## EAS: " + str(eas_samples.shape[0] == 585))
print("## SAS: " + str(sas_samples.shape[0] == 601))
print("## AMR: " + str(amr_samples.shape[0] == 490))
    #the number of samples of populations matches what I have except for the case of EAS and SAS. I think EAS and SAS are flipped, because i got the exact number but in the opposite way.
    #indeed, I have sum the number of samples of EAS and SAS according to Table S1 of the paper and matches what I have
        #EAS: 44+49+46+57+86+77+56+48+60+62=585
        #SAS: 60+71+56+47+61+46+77+69+65+49=601
        #https://www.cell.com/cms/10.1016/j.cell.2022.08.004/attachment/1e0ee9bf-5514-4660-8ff2-6d2c2be97687/mmc1.pdf
    #sample IDs per pop are the same in the high-coverage ped with pop data and the original ped (for shared samples), so we are good with pops in the original ped.

#Among the 3,202 samples, there are 602 father-mother-child trios
trios = samples_pedigree.\
    loc[\
        (samples_pedigree["fatherID"] != "0") &\
        (samples_pedigree["motherID"] != "0"), :]
    #select samples with father and mother
    #This includes 
        #2 trios that are part of a multi-generational family, i.e., one sample is child of another sample but also is parent of another sample. This would be the case of HG00702.
        #10 trios that were split from 5 quads for the purpose of pedigree-based correction applied after haplotype phasing). I GUESS a quad would be a sample that is parent of another sample, that other sample is in turn parent of another sample. This would be the case of HG00656, which is parent of HG00702 that in turn is parent of HG00703.
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 602 trios in the latest ped?")
print("#######################################\n#######################################")
print(trios.shape[0] == 602)
print(trios)

#we also have duos
duos = samples_pedigree.\
    loc[\
        (samples_pedigree["fatherID"] == '0') & (samples_pedigree["motherID"] != '0') | 
        (samples_pedigree["fatherID"] != '0') & (samples_pedigree["motherID"] == '0'), :]
    #select samples that
        #have mother but not father OR
        #have father but not mother
    #there are also 6 parent-child duos. Where a sample only has one parent included here. For example, HG02569 only has mother (HG02568), not father.
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 6 duos in the latest ped?")
print("#######################################\n#######################################")
print(duos.shape[0] == 6)
print(duos)

#the total number of trios and duos is 608, which matches the number of cases with parentID different from zero.
print("\n#######################################\n#######################################")
print("Do we have 608 duos and trios in the latest ped?")
print("#######################################\n#######################################")
print(duos.shape[0] + trios.shape[0] == 608)

#The rest of samples have no parental relationships.
no_trios_duos = samples_pedigree.\
    loc[\
        (samples_pedigree["fatherID"] == '0') & (samples_pedigree["motherID"] == '0'), :]
#check
print("\n#######################################\n#######################################")
print("Do we have 2594 samples with no trios/duos in the latest ped?")
print("#######################################\n#######################################")
print(no_trios_duos.shape[0] == 2594)
print(no_trios_duos)
print("\n#######################################\n#######################################")
print("The 90 samples added on top of the original 2504 are related!! See script for details")
print("#######################################")
    #In the 2022 paper they say 
        #"Using the Illumina NovaSeq 6000 System, we performed WGS of the original 2,504 1kGP unrelated samples and an additional 698 related samples".
        #"Here, we present high-coverage WGS and comprehensive analyses of the original 2,504 1kGP samples, as well as of 698 additional related samples that now complete 602 trios in the 1kGP cohort"
            #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
            #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue#supplementaryMaterial
    #there is an older readme saying that the number of trios is 698 (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/README_2504_plus_additional_698_related_samples.txt).
    #There is also a INFO field in the VCF files that says "2504 unrelated"
    #Jesus and I think that ALL the new samples are related, the 90 samples that are not included as trios/duos maybe are related like cousins or similar. We have to remove these.


##check in detail if we have duos/trios inside the original 2504 sample
print("\n#######################################\n#######################################")
print("check in detail if we have duos/trios inside the original 2504 sample")
print("#######################################")
#we are going to look for samples of the original 2504 set that are sons of samples already included in this dataset.
#these are the problematic cases, we do not care if a sample is in the original set and then, in the high coverage, they add his parents. These parents are not in our original sample, so when we select for 2504 samples, they will be filtered out.

#first merge the original set of samples and the last ped with new duos-trios
merge_parents_original = pd.merge(
    right=original_unrel_ped,
    left=samples_pedigree,
    how="inner",
    right_on="sample",
    left_on="sampleID")
print("\nDo we have 2504 samples when merging the original 2504 sample and the new ped? If True, we have only the original 2504 unrelated samples")
print(merge_parents_original.shape[0] == 2504)
print(np.unique(merge_parents_original["sample"] == merge_parents_original["sampleID"]) == True)
print("\nSee the merged file")
print(merge_parents_original)

#
print("\nNow select those samples of the original dataset whose parents are included that original dataset")
original_unrelated_parents_inside = merge_parents_original.loc[ \
    (\
        (merge_parents_original["motherID"].isin(merge_parents_original["sample"])) | \
        (merge_parents_original["fatherID"].isin(merge_parents_original["sample"]))) & \
    (\
        (merge_parents_original["motherID"] != "0") | \
        (merge_parents_original["fatherID"] != "0")), \
    :]
        #select samples with mother/father ID different from zero and that these parental IDs are included in the sample column, i.e., they are already included as samples in the original 2504 dataset
print(original_unrelated_parents_inside)


##see data about inbreeding in Gazal et al. (2015).
#This paper detect a high level of inbreeding within the phase 1000 genomes project. we propose two different panels of unrelated and outbred individuals to TGP users, such as previously proposed for HapMap III panel8,10. First, for both panels, we removed 7 individuals for quality reasons (Q-score≤50). Then, for the first panel, labeled TGP2457, we removed the 14 individuals involved in first and second degree relationships inferred by RELPAIR and 26 individuals inferred as avuncular offspring (AV) or double first-cousin offspring (2×1C) by FSuite. Finally for the second panel, labeled TGP2261, we removed individuals from 227 relationships up to first-cousins detected by RELPAIR (see Table S3) and the 94 individuals that have been inferred as first-cousin offspring or closer by FSuite. These filters mainly reduced the number of individuals in STU, PJL, ASW and LWK populations, with a sample size decrease of 35%, 31%, 26% and 23%, respectively (Table S5). These 2 lists are provided in Table S4.
    #https://www.nature.com/articles/srep17453
print("\n#######################################\n#######################################")
print("exploring inmbreeding according to Gazal et al. (2015) in phase 3")
print("#######################################")

#I have downloaded Table S4 and convert it to csv in order to make some checks.
gazal_imbreeding = pd.read_csv(
    "./data/pedigrees/gazal_relatedness_phase3_data/41598_2015_BFsrep17453_MOESM3_ESM.csv",
    sep=",",
    header=0,
    low_memory=False)

#we have OUT for non-imbreeding, 2C for second cousins, 1C for first cousins, double first-cousin offspring (2x1C), AV for avuncular offspring (uncle-niece)
print("types of inmbreeding according to Gazal et al. (2015): " + str(gazal_imbreeding.loc[:, "Mating type"].unique()))

#
print("\nSee Gazal data for the 4 samples with parents insides the original unrealted 2504 dataset")
print("\nSee the data of these samples in the new ped")
print(original_unrelated_parents_inside)
print("\nSamples with parents inside the original 2504 unrelated dataset")
print(gazal_imbreeding.loc[ \
    gazal_imbreeding["IID"].isin(original_unrelated_parents_inside["sample"])])
print("\nSee their parents")
print(gazal_imbreeding.loc[ \
    gazal_imbreeding["IID"].isin(original_unrelated_parents_inside["fatherID"]) | \
    gazal_imbreeding["IID"].isin(original_unrelated_parents_inside["motherID"])])
    #Table S3 from Gazal (41598_2015_BFsrep17453_MOESM2_ESM.xls) shows that these 4 relations (PO) are included, correctly showing the corresponding duos.


##summary
#All the new samples added on top of the original 2504 are related in some degree to the previous one, so we should not use them.
#In addition, I have found within the original set of 2504 that some samples are included in trios/duos in the last ped.
#we are going to use unrelated samples only. If we are calculating haplotype homozygosity by counting haplotypes, if a father and child have the same haplotype, this is not caused by selection but just by shared ancestry
#Gazal data does not match duos trios in high coverage, so I do NOT think we can use this data for high coverage.
#Therefore, we should use the original dataset, but also removing the samples included in new trios/duos, reducing the probability to include parents/sons.


##remove trios/duos within the original 2504 set
#select those samples from the original 2504 set that are NOT included in known duos/trios according to last ped
unrelated_samples = original_unrel_ped.\
    loc[\
        ~original_unrel_ped["sample"].isin(original_unrelated_parents_inside["sample"]), :]
print("\n#######################################\n#######################################")
print("see final pedigree data with only unrelated samples")
print("#######################################\n#######################################")
print(unrelated_samples)

#extract the distinct population names
pop_names = unrelated_samples["pop"].unique()

#check
print("\n#######################################\n#######################################")
print("Do we have 26 pops?")
print("#######################################\n#######################################")
print(len(pop_names) == 26)
print(pop_names)



#################
# bcftools prep #
#################

#see version of bcftools
print("\n#######################################\n#######################################")
print("see bcftools version")
print("#######################################\n#######################################")
run_bash("bcftools -v")
    #bcftools cheatsheet
        #https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b#file-bcftools-cheat-sheet



####################################################
# check bcftools's behaviour with a dummy vcf file #
####################################################
print("\n#######################################\n#######################################")
print("We are going to use dummy vcf file to check the behaviour of bcftools: ")
print("#######################################\n#######################################")

#
print("\n#######################################\n#######################################")
print("see the dummy variants")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #We have variants with different characteristics
            #Biallelic SNPs:
                #SNP rs6054248 chr20 14290 C A 5 0 0.2 GTs: 0|0 0|0 1|.
                #SNP rs6054249 chr20 14300 C A 5 1 0.2 GTs: 0|0 0|0 1|.
                #SNP rs6054250 chr20 14310 C A 3 0 0 GTs: 1|0 .|. .|.
                #SNP rs6054251 chr20 14320 C A 5 1 0.2 GTs: 0|0 0|0 1|1
                #SNP rs6054252 chr20 14350 C A 4 2 0.5 GTs: 1|0 1|0 .|.
                #SNP rs6054257 chr20 14370 G A 6 3 0.5 GTs: 1|0 1|1 0|0
                #SNP rs6054255 chr20 14371 G C 6 3 0.667 GTs: 0|1 1|1 1|0
                #SNP rs6040351 chr20 17330 T A 6 2 0.333 GTs: 0|0 1|0 1|0
            #Multiallelic SNPs:
                #SNP rs6040355 chr20 1110696 A G,T 6 2,2 0.333,0.333 GTs: 1|1 2|2 0|0
                #SNP rs6040356 chr20 1110697 A G,T 6 2,1 0.333,0 GTs: 1|1 0|0 0|0
                #SNP rs6040357 chr20 1110698 A G,T 6 3,3 0.5,0.5 GTs: 1|1 2|2 1|2
                #SNP rs6040358 chr20 1110699 A G,T 6 2,0 0.333,0 GTs: 1|1 0|0 .|.
                #SNP rs6040359 chr20 1110700 A G,T 5 0,3 0,0.6 GTs: 2|2 2|2 .|.
            #Exact duplicate SNPs, i.e., pos, chr, REF, alt
                #SNP rs6040360 chr20 1110701 A G 6 3 0.5 GTs: 1|1 0|0 0|1
                #SNP rs6040360_copy chr20 1110701 A G 6 3 0.5 GTs: 1|1 1|0 0|1
            #SNP with unphased data
                #SNP rs6040361 chr20 1110702 A G 6 3 0.5 GTs: 1|1 0/0 0|1
            #microsatellite (indel)
                #INDEL microsat1 chr20 1110703 GTC G,GTCT . . . GTs: 0/1 0/2 1/1

#
print("\n#######################################\n#######################################")
print("split multiallelic SNPs ")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
    #split multiallelic SNPs in different lines
        #--multiallelic -snps
            #split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+).
            #https://samtools.github.io/bcftools/bcftools.html#norm
    #when a multiallelic snp (e.g., REF=A, ALT=G,T) is split in several lines, these lines have the same position and REF allele, but different ALT. This was done in 1000 genomes data. They also add 1 to the position of each new line to do the phasing, but then they put all lines back in the same position.
    #In these new lines, the ALT is now only one, but the AC still has two fields, i.e., count for both ALTs in both lines. the 1000 genomes project SNPs that are multiallelic have their AC field just with one value so I guess they used +fill-tags to update these fields (see below).
            #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf 
        #First scenario: rs6040355
            #the first line shows the genotypes for A and G, so a sample that is G|G would be 1|1, but a T|T sample would be 0|0, so A and T gets 0 in this line. This is ok, because the next line shows genotypes for A and T, and that T|T sample will be 1|1.
                #This case is multiallelic, so we should remove both lines. We can just combine these two lines with norm --multiallelic, as it looks for snps with same pos and similar REF to merge (see when --multiallelic +snps is first used).
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

#note about the output of bcftools norm
    #Lines   total/split/realigned/skipped: 12/5/0/0
    #12 is the total number of snps, while 5 is the number of snps that are multiallelic so have been splitted

#
print("\n#######################################\n#######################################")
print("select only two samples to check whether AC, AN and AF changes")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
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
print("\n#######################################\n#######################################")
print("select SNPs")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools norm \
        --multiallelic -snps | \
    bcftools view \
        --types snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #--types/--exclude-types LIST
            #Select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other

#
print("\n#######################################\n#######################################")
print("see monomorphic snps: rs6054249 (0|0 0|0 1|.) is not considered as monomorphic which is right, given that the last genotype not only has missing but also an ALT allele. rs6054248 (0|0 0|0 1|.) is not included even having AC erronously set as zero, thanks we update AC field with +fill-tags. The cases considered as monomorphic are those with all REF (e.g., rs6040358) or all ALT (e.g., rs6040359) irrespectively if they have missing. For example, rs6040359 has '1|1 1|1 1|.' as GT and it is considered monomorphic despite having a missing allele in the last genotype")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
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
print("\n#######################################\n#######################################")
print("see monomorphic snps after subseting two first samples: rs6054249 is now considered monomorphic because it has 0|0 0|0 1|., being 1|. removed after filtering out the third sample. We get a warning after subsetting because --sampls updated AC, but for splitted multiallelic snps, we have several allele counts in each line, so the first time it sees this (rs6040355), print warning and no more. But AC is correctly updated anyway, leaving just one count, thus breaking the connection between the allele counts of the lines of a given multiallelic snp.")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")

#
print("\n#######################################\n#######################################")
print("exlcude monomorphic snps")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")

#
print("\n#######################################\n#######################################")
print("see what missing definition uses bcftools in stats: Stats considers as missing both '.|.' and '0|.', resulting in 5 missing for the third sample (3 cases with .|. and 2 cases with 1|.). This makes sense because we are looking for the proportion of genotypes with missing data, '0|.' has missing data besides the allele '0'. Therefore, it should be counted as missing. We will consider this using GT='mis' (see below)")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
print("\n#######################################\n#######################################")
print("filter by frequency of missing. This approach counts the number of missing genotypes with GT='mis' ('.|.', '.', '0|.') and divided by the number of samples, so you get the proportion of samples with missing genotypes. The strength of this approach over calculating the tag F_MISSING (+fill-tags) is that you are counting anything with '.', so even '1|.' would be considered missing. To me this is more correct because we are filtering by genotype missingness, i.e., the proportion of genotypes with missing. Therefore, '1|.' has missing even though it has one allele. It is different in previous filters when we want to check if a SNP is monomorphic. In that situation, having 0|0 1|. is not monomorphic if we count the last allele in the genotype with missing. At the same time the genotype has missing (it should be count for missingness) and has an ALT allele (it is not monomorphic), we should consider both aspects")
print("#######################################\n#######################################")
print("#### SNPs with 1/3 or more of missing: we have snps with 1 out of 3 missing genotypes (e.g., rs6040358) and variants with 2 out of 3 missing (rs6054250) ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES >= 1/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
print("#### SNPs with 2/3 of missing: We get a variant with 2 out of 3 missing (rs6054250) ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES = 2/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
print("#### select variants with less than 1/3 of missing: We have lost all variants with missing because we only have variants 2/3 and 1/3 of missing ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 1/3' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #Lines of multiallelic rs6040358 have been removed because the last sample has a missing genotype. Note that rs6040358 has not been excluded as monomorphic as it has 1 and 0, but because it has missing.
print("#### select variants with less than 0.05 of missing: We have lost all variants with missing because we only have variants 2/3 and 1/3 of missing ####")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools view \
        --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #1/2 would be 0.5, i.e., 50%. Therefore, 5% would be 0.05 (i.e., 0.5/10).

#
print("\n#######################################\n#######################################")
print("remove now the exact duplicates: it does not matter if we remove duplicates before or after subsetting samples. Subseting does not update chrom, pos and REF/ALT, which are the fields used by '--rm-dup exact' to remove exact duplicates. In an hypothetically scenario where we have two duplicated SNPs and one of them is monomorphic, we would remove first the monomorphic due to the previous filters and then retain the second. If you use --rm-dup before, you would select the SNP that appears first, that could be the monomorphic or the other one. Removing mono first force to select always the second, while removing duplicates before makes possible to select one or the other. I do not think this is relevant")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. As you can see for rs6040359, it does not matter the ID or the genotypes, this command only targets chr, pos, ref, and alt. 
        #bcftools norm --rm-dup exact selects only the first appearance and remove the next. Consequently, rs6040360_copy is filtered out, retaining rs6040360.
            #https://github.com/samtools/bcftools/issues/1089
        #https://samtools.github.io/bcftools/bcftools.html#norm

#
print("\n#######################################\n#######################################")
print("combine lines of multialleic snps in just one line")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
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
print("\n#######################################\n#######################################")
print("check what happens if we remove the third sample and then rs6040355 have no REF anymore, only two ALTs: Despite losing the REF, both lines of rs6040355 are merged into a multiallelic SNP, thus it can be removed with --max-alleles 2. This makes sense because --multiallelics +snps looks for snps with the same position and same REF or ALT")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #.

#
print("\n#######################################\n#######################################")
print("now we can remove those SNPs with more than 1 allele in ALT. rs6040358 is now included because the second ALT is not present, only REF and ALT, having no missing (it has missing for the last sample that is filtered out)")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #--min-alleles/--max-alleles INT  
            #Minimum/maximum number of alleles listed in REF and ALT
        #This looks at REF/ALT, so we need these columns updated in order to use this. If you just do this when the multiallelic SNPs are splitted across different lines, then you would not remove them, because each one is presented like a biallelic SNP.

#   
print("\n#######################################\n#######################################")
print("select only phased data")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
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
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #rs6040360 has "1|1 0/0 0|1" as GT, thus the second sample does not have phased data.
        #it is removed after applying the phased filter.
        #I separated phased from --max/min-allele filter to be completely sure that first we remove multiallelic snps, and then we remove those snps with unphased data for any sample.
        #I guess it is not mandatory to do it in this order because if a SNP that is multiallelic has also unphased data, will be removed anyways, and biallelic and phased SNPs will not be affected. 
        #The order is more important in the case of samples, as you could remove a SNP that has unphased data for some samples, but any of those samples are included in the subset, so you are removing it while it has complete phased data for the samples you selected.

#   
print("\n#######################################\n#######################################")
print("select only SNPs in regions with high accessibility")
print("#######################################\n#######################################")
#general info
    #The 1000 Genomes Project created what they defined as accessibilty masks for the pilot phase, phase one and phase three of the Project. Some other studies have similar files.
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
#mask data here
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/
    #we will use the bed files with the regions has to be retained
        #"In addition to masked fasta files, a bed file of all passed sites can be found in this directory."
#dummy mask here:
    #/home/dftortosa/singularity/dating_climate_adaptation/hg38_mig/data/dummy_vcf_files/dummy_pilot_mask.bed
#compress the dummy (mask) bed file to match what we will do with the real data. bcftools --targets-file can take a .bed.gz file as input (see below).
run_bash(" \
    gzip \
        --force \
        --keep \
        ./data/dummy_vcf_files/dummy_pilot_mask.bed")
    #--force: force overwrite of output file and compress links
    #--keep: keep (don't delete) input files
#apply the dummy mask
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
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
print("\n#######################################\n#######################################")
print("check you can filter the BED file by chromosome using awk")
print("#######################################\n#######################################")
run_bash(" \
    gunzip -c ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    awk \
        -F '\t' \
        '{if ($1 == \"chr20\") print $0}' \
        > ./data/dummy_vcf_files/dummy_pilot_mask_chr20.bed; \
    gzip --force ./data/dummy_vcf_files/dummy_pilot_mask_chr20.bed; \
    gunzip -c ./data/dummy_vcf_files/dummy_pilot_mask_chr20.bed.gz")
        #decompress the dummy bed file and sent it to stdout
        #process it with awk
            #-F is the delimiter, "\t" in this case
            #for each row, if the first field (column) has a value of "chr20", print all fields and save into a file.
            #compress the file even if a compressed file with the same name exists (--force).
                #https://unix.stackexchange.com/questions/399560/using-awk-to-select-rows-with-specific-value-in-specific-column
                #https://stackoverflow.com/questions/2961635/using-awk-to-print-all-columns-from-the-nth-to-the-last

#
print("\n#######################################\n#######################################")
print("check the cleaning of the BED file with awk")
print("#######################################\n#######################################")
run_bash(" \
    uniq_chrom=$(\
        gunzip -c ./data/dummy_vcf_files/dummy_pilot_mask_chr20.bed.gz | \
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
print("\n#######################################\n#######################################")
print("use +fill-tag to update INFO fields")
print("#######################################\n#######################################")
print("\n#######################################\n#######################################")
print("first, see the sample 1 and 2 to check AF is not updated")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #As expected, after the subset of samples AN and AC are updated considering only the selected samples, while AF is not. For example, rs6040351 has only 1 ALT allele after selecting the two first samples, so AN=4, AC=1 but AF is still 0.333 instead of 0.25 (1/4).

#
print("\n#######################################\n#######################################")
print("update AF with fill-tags")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools +fill-tags \
        -- --tags AN,AC,AF | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #we use +fill-tags to update fields that are not updated after subsetting like frequency of the alternative allele and create some additional fields.
            #+fill-tags has updated AF, so for example, rs6040355 has 0.5 for the frequency of both ALTs, while before it was 0.333. Similarly, the frequency of the ALT is 0.25 in rs6040351, which is correct.
        #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed due to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields.

#
print("\n#######################################\n#######################################")
print("update and create more fields with fill-tags")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        --samples NA00001,NA00002 \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n'")
        #we have create several fields that are included the 1KGDP vcf files or are related:
            #INFO/AN: Total number of alleles in called genotypes
            #INFO/AC: Allele count in genotypes
                #According to vcf specification file v2: 
                    #allele count in genotypes, for each ALT allele, in the same order as listed
            #INFO/AC_Hom: Allele counts in homozygous genotypes
            #INFO/AC_Het: Allele counts in heterozygous genotypes
            #INFO/AF: Allele frequency from FMT/GT or AC,AN if FMT/GT is not present
            #INFO/MAF: Frequency of the second most common allele
            #INFO/ExcHet: Test excess heterozygosity; 1=good, 0=bad
            #INFO/HWE: HWE test (PMID:15789306); 1=good, 0=bad
            #INFO/NS: Number of samples with data
        #This works after selecting the two first individuals, see for example:
            #rs6054257
                #GT: 1|0 1|1
                #AN=4
                #AC=3
                #AC_Hom=2
                #AC_Het=1
                #AF=0.75
                #MAF=0.25
                #NS=2
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
print("\n#######################################\n#######################################")
print("remove all previous INFO and FORMAT fields except GT and create the fields you are interested in by using fill-tags but after applying all the filters")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove INFO,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n'")
        #before updating/creating fields, remove all INFO fields and all FORMAT fields (except GT) using annotate and then add the fields we are interested in.
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
print("\n#######################################\n#######################################")
print("compare GT between applying or not the +fill-tags commands and the removal of fields")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove INFO,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --no-header")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools view \
        --no-header")
        #I have checked that the genotypes remain the same despite removing all previous INFO fields and all FORMAT fields (except GT) and then adding new INFO fields, so we are good here.

#
print("\n#######################################\n#######################################")
print("see the new header")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove INFO,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools head")
        #the new header shows 
            #the fields of the previous VCF file.
            #every bcftools command used in order, along with the flags, the version, the date.
            #the new fields generated with +fill-tags
            #version of +fill-tags, flags selected
        #Therefore, the new header includes a history of the changes made to the VCF file inside of the file itself.

#
print("\n#######################################\n#######################################")
print("save the cleaned vcf file")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/dummy_example.vcf |\
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
        --targets-file ./data/dummy_vcf_files/dummy_pilot_mask.bed.gz | \
    bcftools annotate \
        --remove INFO,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --output ./data/dummy_vcf_files/dummy_example_cleaned.vcf.gz \
        --output-type z \
        --compression-level 1")
        #--output
            #just the path and name of the file
        #--output-type
            #z for compressed VCF file
        #--compression-level
            #Compression level: 0 uncompressed, 1 best speed, 9 best compression
        #you can see in stdout that in total lines we get 12, which is the total number of variants we have before filtering, instead of the final number of variants (3), see below.

#
print("\n#######################################\n#######################################")
print("see the header of the recently created dummy vcf file")
print("#######################################\n#######################################")
run_bash(" \
    bcftools head \
        ./data/dummy_vcf_files/dummy_example_cleaned.vcf.gz")
    #we get first the fields of FILTER and FORMAT that have been not removed (you only removed INFO field and all FORMAT except GT). 
    #Then we get all the commands we have run and the fields we have added.

#
print("\n#######################################\n#######################################")
print("see the genotypes of a few individuals from the recently created dummy vcf file")
print("#######################################\n#######################################")
run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/dummy_example_cleaned.vcf.gz \
        --no-header")

#
print("\n#######################################\n#######################################")
print("calculate stats of the VCF file, show them here and then use them to make summary plots")
print("#######################################\n#######################################")
run_bash(" \
    bcftools stats \
        --samples - \
        ./data/dummy_vcf_files/dummy_example_cleaned.vcf.gz")
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

#
print("\n#######################################\n#######################################")
print("convert VCF file to hap file")
print("#######################################\n#######################################")
print("First convert the raw dummy VCF file but adding a filter for SNPs in the same line, i.e. removing the microsatellite. The output shows 0/5/1 no-ALT/non-biallelic/filtered indicating that we lose 5 variants that are not biallelic and 1 due to the filter we applied, i.e., removal of non-snps, the microsatellite. Also, the number of records written is 6, which is the result of subtracting 5+1 from 12, the total number of variants")
run_bash(" \
    bcftools convert \
        --include 'TYPE == \"SNP\"' \
        ./data/dummy_vcf_files/dummy_example.vcf \
        --hapsample ./data/dummy_vcf_files/dummy_example_IMPUTE2_raw; \
    gunzip -c ./data/dummy_vcf_files/dummy_example_IMPUTE2_raw.hap.gz")
print("Then convert the cleaned dummy VCF file. The output shows 0/0/0 because no variant has been removed, they were filtered in previous lines. The number of written records is 3. The input had 3 and no filter was applied, nor variant lose due to lack of alt or non-biallelic, thus we save a hap file with 3 variants")
run_bash(" \
    bcftools convert \
        ./data/dummy_vcf_files/dummy_example_cleaned.vcf.gz \
        --hapsample ./data/dummy_vcf_files/dummy_example_IMPUTE2_raw; \
    gunzip -c ./data/dummy_vcf_files/dummy_example_IMPUTE2_raw.hap.gz")
        #see hap file calculation of the real data to see details about hap format.

#SUMMARY: 
    #With all these commands, we have recreated the scenario we have in 1KGP data, with multiallelic SNPs separated into different lines, select some samples, we then select snps, remove those with genotype missingness < 5%, remove exact duplicates (this does not touch different lines of the same multiallelic snp because they have different ALT), exclude those SNPs that have the same allele for all samples (considering alleles in genotypes with missing, e.g., 0|.). Then we combine all lines of each multiallelic snp and now they have ALT column with several alleles, so we can filter them using --max-alleles 2. Add filter for selecting phased data only. Select only those variants included in interest regions. We can also use bcftools +fill-tags to update important fields for each SNP, so if a SNP was multiallelic, but it is not multiallelic in the subset population (i.e., only REF and 1 ALT), we no longer will have two allele frequencies, two allele counts.... for the remainder biallelic SNP in the subset.

#Note about the update of the INFO fields
    #it is important to be sure that the fields you are using for filtering, are updated after subseting samples. Of course, type="snp" will be always "snp" irrespectively of the samples we select, but this is not the case of the number of alleles, because you can have SNPs with 3 alleles considering all 26 populations, but then in GBR they can have only 2 or 1. We are interested in SNPs that are biallelic within the selected population.
    #The same goes for phasing and genotype ^miss. You have to be sure that these filters only consider data from the filtered samples, not fixed data in fields that are not updated.
    #Because of this we have applied all these filters in order.



################################################################
#### function to clean vcf files and create hap - map files ####
################################################################

#selected_chromosome="1"; selected_pop = "GBR"
def master_processor(selected_chromosome, selected_pop):



    #######################
    # one time operations #
    #######################

    #do first some operations that need to be done just one time per chromosome
    if selected_pop == pop_names[0]:
        
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
                #load the whole VCf file of the selected chromosome, select only SNPs and in one case select variants inside mask "PASS" regions, while in the second do nothing more. Then see the stats of the resulting BCF files. See dummy examples for further details.



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

    #select the sample IDs
    selected_samples = subset_pop["sample"]

    #save as txt to use it later
    #it is redundant to do it in each chromosome of the same population (same samples) but I do it anyway to avoid confusion between chromosomes
    selected_samples.to_csv(
        "results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt", 
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
    run_bash(" \
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
        fi")
            #obtain the ID of all variants after applying or not the filter for genotype missingness < 0.05, then check the numbers are the same

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


    ##por aqui
    #check genotype missingness
        #run :  check we do not have SNPs with genotype missingness > 0.05.
    #use pilot mask
        #check all the lines where pilot has been added


    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": clean BED file before applying the mask selecting only intervals in the selected chromosome")
    print("#######################################\n#######################################")
    run_bash(" \
        awk \
            -F '\t' \
            '{if ($1 == \"chr" + selected_chromosome + "\") print $0}' \
            ./data/masks/20160622.allChr.pilot_mask.bed \
            > ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed; \
        gzip --force ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed; \
        gunzip -c ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
        head -5")
            #in the bed file, select those rows for which the chromosome name (first column) is the selected chromosome, printing all fields for these rows. Save as a file and then compress. See the first 5 lines. See dummy examples for further details.

    #
    print("\n#######################################\n#######################################")
    print("check the cleaning of the BED file with awk")
    print("#######################################\n#######################################")
    run_bash(" \
        uniq_chrom=$(\
            gunzip -c ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
            awk \
                -F '\t' \
                '{print $1}' | \
            uniq); \
        if [[ $uniq_chrom == 'chr" + selected_chromosome + "' ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress the bed file of the selected chromosome (previously created), then print the chromosome name (first column) for all intervals, i.e., rows, getting then the unique cases, then save as a variable. The variable should have only the name of the selected chromosome.

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
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #Use --targets-file to only select SNPs within the intervals defined by a BED file generated by 1kGDP, which is an accessibility mask. Therefore, we select SNPs that are included in regions accessible to sequencing.
            #see dummy example for further details.
    
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
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n' | \
        head -7")
            #we sue +fill-tags to update fields that are not updated after subsetting like frequency of alternative allele and create some additional fields.
            #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed sue to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields.
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
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
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
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
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
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + ".pilot_mask.bed.gz | \
        bcftools annotate \
            --remove INFO,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --output ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
            --output-type z \
            --compression-level 1")
            #after subseting and filtering, save as a compressed VCF file, selecting the option for best speed (see dummy example for further details)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see the header of the recently created vcf file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools head \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see the genotypes of a few individuals from the recently created vcf file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
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

    #    
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": convert to hap file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools convert \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
            --hapsample ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw")
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
            #Hap file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.hap.gz
            #Sample file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.samples
            #938126 records written, 0 skipped: 0/0/0 no-ALT/non-biallelic/filtered
        
        #I guess non-biallelic and non-alt should be removed to meet impute requirements, but in our case we already selected those SNPs with 2 alleles only. This explains why we get 0/0/0. This is the number of SNPs with no-ALT, no-biallelic and filtered. All the filters were previously applied.
            #I have checked this in the dummy example.
        
        #if you take the VCF file, filter and then count lines, you get 938126 records, but then say that total is 945919, as when writting the file. "Lines   total/split/realigned/skipped" is produced for some commands like bcftools norm and I have seen in the dummy example that total is the total number of SNPs, while split are those splitting due to multiallelic is the multiallelic flag is used.
            #938126
            #Lines   total/split/realigned/skipped:  945919/0/0/0
            #Lines   total/split/realigned/skipped:  945919/0/0/0

        #the total number of snps in the raw vcf file of chr1 is 5759173, instead of 945919. What it can be happening here is that bcftools norm (the command generating the line output) is run after some filters have been applied, thus the input number of SNPs is smaller. Indeed, if you apply these previous filters, you get 945919 variants, so it makes sense that bcftools norm, which is run after these filters, gives 945919 as total number of snps.

    #see first variants
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first variants of the hap file")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        head -3")
            #decompress the hap file and show in stdout 
    
    #clean
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove first columns of hap file to leave only haplotype columns")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        cut \
            --complement \
            --delimiter ' ' \
            --fields 1-5 \
        > results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap; \
        gzip -f results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap")
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
            -c results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        awk -F ' ' '{print " + fields_selected_samples + "}' > \
        results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap; \
        gunzip \
            -kf results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz; \
        file1='results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap'; \
        file2='results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap'; \
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
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the same number of rows (variants) in the cleaned vcf and hap files")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps_vcf=$( \
            bcftools view \
                ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
                --no-header | \
            wc -l); \
        n_snps_hap=$( \
            gunzip \
                -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            wc -l); \
        if [[ $n_snps_vcf -eq $n_snps_hap ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #load the cleaned vcf file with bcftools and show only the snps, with no header, then count the number of lines and save the result
            #decompress the hap file and count the number of line
            #check both numbers are the same

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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": do we have the correct number of samples in the hap file?")
    print("#######################################\n#######################################")
    run_bash(" \
        nfields_hap=$( \
            gunzip \
                -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            awk -F ' ' '{print NF}' | \
            sort | \
            uniq); \
        nfields_hap_nlines=$( \
            echo -n $nfields_hap | \
            grep -c '^'); \
        if [[ $nfields_hap_nlines -eq 1 ]]; then \
            nsamples_hap=$(($nfields_hap/2)); \
            nsamples_bcftools=$( \
                cat results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt | \
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
        "./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.samples",
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



    ##########################################
    # calculate map file within selected pop #
    ##########################################

    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": STARTING MAP FILE CALCULATION")
    print("#######################################\n#######################################")

    ##extract the SNP positions
    #first extract from the raw hap file, the clean hap file only have genotypes
    run_bash(" \
        gunzip \
            --stdout \
            ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        awk \
            -F ' ' \
            '{print $1,$2,$3,$4,$5}' | \
        gzip \
            --force \
        > ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz")
            #decompress raw hap file, select first columns with position data, IDs, and allele names, then compress and save as a file

    #then from the cleaned VCF file
    run_bash(" \
        bcftools view \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools query \
            --format '%CHROM %POS %REF %ALT\n' | \
        gzip \
            --force \
        > ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw_check.map.gz")
            #from the cleaned VCF file, extract the chromosome, position, REF/ALT to compare with the positions in the raw hap file (see below), compress and save as a file

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that the two raw map files obtained from hap_raw and vcf_clean are the same")
    print("#######################################\n#######################################")
    run_bash("\
        cd ./results/hap_map_files_raw/; \
        gunzip \
            --stdout \
            ./chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz | \
        awk \
            -F ' ' \
            '{print $1,$3,$4,$5}' > \
        ./chr" + selected_chromosome + "_" + selected_pop + "_raw.map; \
        gunzip \
            --force \
            ./chr" + selected_chromosome + "_" + selected_pop + "_raw_check.map.gz; \
        status=$(\
            cmp \
                --silent \
                ./chr" + selected_chromosome + "_" + selected_pop + "_raw.map \
                ./chr" + selected_chromosome + "_" + selected_pop + "_raw_check.map;\
            echo $?); \
        if [[ $status -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi; \
        rm ./chr" + selected_chromosome + "_" + selected_pop + "_raw.map; \
        rm ./chr" + selected_chromosome + "_" + selected_pop + "_raw_check.map")
            #get chromosome, position and REF/ALT for each SNP from the raw hap file 
                #We are not using ID because it is complicated to match the format of VCF and hap files. Importantly, we do not need that because the ID of the hap file follows the format CHROM:POS_REF_ALT, i.e., it is based in chromosome, position and REF/ALT, which is the data we are comparing.
                #decompress the map file and send to standard output, select the corresponding columns and save as a file
            #decompress the map file obtained from the VCF file, and do not keep the compressed file because that is only for the check
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
            #if status is zero (both files are the same), perfect
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": This check also tell us that the coordinates in hap file are 1 based. Hap and VCF file have the same positions, and the pos in VCF files v4.2 is 1-based according to the specification file (this is the format of 1KGP data)")
    print("#######################################\n#######################################")
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


    ##calculate the genetic position
    #load the raw_map file
    snp_map_raw = pd.read_csv(\
        "./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz", \
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
                ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
            wc -l); \
        if [[ $n_snps -eq " + str(snp_map_raw.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    #rename the columns
    snp_map_raw = snp_map_raw.rename(\
        {0: "chr", 1: "id", 2: "pos", 3: "ref", 4: "alt"}, \
        axis=1)
            #use a dict with old and new column names. indicated we are renaming columns (axis=1)
    print("see map file with renamed columns")
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that the ID column in the raw map file is exactly the combination of chromosome, pos, ref, alt")
    print("#######################################\n#######################################")
    check_id = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
        #make a series combining chromosome, pos, ref and alt, and using the corresponding separators
    print(check_id.equals(snp_map_raw["id"]))
        #check it is identical to id

    #
    print("remove the ref/alt columns as we have this information already included in the ID")
    snp_map_raw = snp_map_raw.drop(["ref", "alt"], axis=1)
    print(snp_map_raw)


    ##load and explore the decode2019 map

    #I know that the original 2019 decode map is alligned to hg38. Also, I assume that the decode2019 map is 1-based because they do not specify is 0-based
        #Data S3.genetic.map.final.sexavg.gor.gz:
            #average genetic map computed from the paternal and maternal genetic maps, which were in turn computed from the paternal and maternal crossover, respectively. The data columns are as follows: Chr (chromosome), Begin (start point position of interval in GRCh38 coordinates), End (end point position of interval in GRCh38 coordinates), cMperMb (recombination rate in interval), cM (centiMorgan location of end point of interval)
            #Page 85 of "aau1043-halldorsson-sm-revision1.pdf"

    #1KGP is aligned to hg38 (see paper) and coordinates are 1-based as VCF format 4.2 has 1-based coordinates.

    #therefore, we have the same position format in both datasets, so we can just use the decode 2019 map to calculate the genetic position of each SNP.

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


    #selected_snp=snp_map_raw.iloc[0,:]
    def gen_pos(selected_snp):

        #select snp id
        selected_snp_id = selected_snp["id"]

        #selected chromosome
        selected_chrom = selected_snp["chr"]

        #extract position of the selected snp
        selected_snp_physical_pos = selected_snp["pos"]

        decode2019_map_subset = decode2019_map.loc[decode2019_map["chr"] == selected_chrom,:]

        np.unique(decode2019_map_subset["chr"]) == selected_chrom
        
        decode2019_map_subset


        #it makes sense to repeat this for each populations? maybe do it per chromosome and then select those snps included in each pop?
            #maybe you can calculate the maps in a different script run before this and then run this, selecting snps of the corresponid chromosome map that are included in the vcf file of the pop,
            #remember that remove filter snps wihitn pop, so each pop will have different snps
            #the ID would be different because REF/ALT are included in the ID in order to avoid strand flips
            #if it is not very slow, we could just do it per popuatilion-chrom is we did for iHS

        #cM (centiMorgan location of END POINT of interval)


        ##por aquii,
            #use new monomorphic filter
            #select missing < 5%
            #use pilot mask
            #use 2504 samples less 4 related. 

            #You should apply the filter with less than 5% missing just as the 1KGP authors. For HWE you should also only apply what the 1KGP authors did and not filter further for specific populations. Biallelic SNPs per specific population are ok. Masks based on low coverage are ok, you should use the less stringent masks. You can use the 2,504 individuals and remove the four related individuals yes. And you are right for MAF, the filtering at MAF>5% is done by the scripts for summary statistics. It is only for the focal SNPs methods like iHS still use the SNPs with lowe MAF around the focal SNPs so it would be an error to remove all SNPs with MAF<5%.





#############
# ask enard #
#############

#as i filter within pop, each pop can have different snps.
#in the map file we can use the format "CHROM:POS_REF_ALT" for the ID? I think remember that the map files I originally got from you in the previous project (before decode2019 conversion) used as ID just the physical position. Not sure if there is any specific reason for doing that.
