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



##############################################################
######## SELECT THE PEDEGREE AND SAMPLES FOR HG38 MIG ########
##############################################################

#We are going to migrate to hg38 as the new release of 1KGP matches this genome reference version (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

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

#set working dir
working_dir = "/home/dftortosa/singularity/dating_climate_adaptation/hg38_mig"
import os
os.chdir(working_dir)



#############
# pops prep #
#############

#load original 2504 unrelated samples from phase 3. This includes sample IDs and pop/superpop codes. This is the definitive ped we will use both for pop codes and sample IDs
#Data:
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#Readme
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220
import pandas as pd
original_unrel_ped = pd.read_csv(
    "./data/pedigrees/integrated_call_samples_v3.20130502.ALL.panel.txt", 
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
#these are the problematic cases, we do not care if a sample is in the original set and then, in the high coverage, they add his parents. These parents are not in our original sample, so when we select for 2504 samples, these parents will be filtered out.

#first merge the original set of samples and the last ped with new duos-trios selecting only samples included in both datasets
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
print("\nNow select those samples of the original dataset whose parents are included that original dataset. One of the cases (NA20318) is mentioned in known-issues of phase 3")
    #known issues: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_known_issues_20200731
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
    #These 4 samples are not included in the cleaned datasets of Gazal, suggesting that they were detected as inbreeding cases. The parents, in contrast, were retained in the dataset.


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

#
print("\n#######################################\n#######################################")
print("check that final pedigree does NOT have any previously known related individual according to phase 3 data. These should be out from the original 2504 dataset")
print("#######################################\n#######################################")
#load the data
#these are samples known to be related in phase 3. These should be out from the 2504 samples. Indeed, they say in the readme file of phase 3 that "We have removed the genotypes for 31 individuals who have blood relationship with the 2504 samples in the main release. This was done to ensure we do not over estimate allele frequency. The 31 related samples are listed in 20140625_related_individuals.txt. This has resulted in a small number of AC=0 sites from rare alleles only present in one or more of these 31 individuals."
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220
#I downloaded from the same folder than the original pedigree (2504 samples): http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt
known_related_in_original_dataset = pd.read_csv(\
    "./data/pedigrees/20140625_related_individuals.txt", \
    sep="\t", \
    header=0, \
    low_memory=False)
#solve a problem in Sample column name
known_related_in_original_dataset = known_related_in_original_dataset.rename({"Sample ": "Sample"}, axis=1)
#
print("known related samples in the original dataset")
print(known_related_in_original_dataset)
#see if no sample in the original dataset is included in the dataset with known related samples 
print("Do the check:")
print(unrelated_samples.loc[unrelated_samples["sample"].isin(known_related_in_original_dataset["Sample"]), :].shape[0] == 0)

#
print("\n#######################################\n#######################################")
print("Do we have 26 pops in the final pedigree?")
print("#######################################\n#######################################")
#extract the distinct population names
pop_names = unrelated_samples["pop"].unique()
#check
print(len(pop_names) == 26)
print(pop_names)


##Note about phase 3 pedigree of 2020
#In the folder of phase 3, there is a pedegree more recent than "integrated_call_samples_v3.20130502.ALL.panel", which is the ped I have used to obtain the original 2504 samples and the one that Jesus recommended me to use.
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped
#This more recent ped has strange stuff
    #there are 3691 samples in total!
    #samples label as unrelated are only 1181
    #samples with parental and maternal ID equal to zero are 2927
#I am not going to use this pedigree for anything.


#
print("\n#######################################\n#######################################")
print("save the final pedigree")
print("#######################################\n#######################################")
unrelated_samples.to_csv( \
    "./data/pedigrees/unrelated_samples.tsv", \
    sep='\t', \
    header=True, \
    index=False)
run_bash("ls ./data/pedigrees")
