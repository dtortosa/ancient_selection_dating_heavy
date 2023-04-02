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
    elif ("Lines   total/split/realigned/skipped" in complete_process.stderr):

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #print the standard error, which is not an error
        print(complete_process.stderr)

        #return also the value if required
        if return_value==True:
            return complete_process.stdout 
    elif ("no-ALT/non-biallelic/filtered" in complete_process.stderr):
        
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

#load pedigree of the latest version of the phased data that has sample IDs, sex and parents but no pop names. 
#Downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
import pandas as pd
samples_pedigree = pd.read_csv(
    "data/pedigrees/1kGP.3202_samples.pedigree_info.txt", 
    sep=" ", 
    header=0, 
    low_memory=False)

#load also a pedigree present in the main directory of the high coverage data. This has sample and pop IDs, but parents and sex are different with respect to the pedigree of the new sample
#downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
samples_pedigree_pop = pd.read_csv(
    "data/pedigrees/20130606_g1k_3202_samples_ped_population.txt", 
    sep=" ", 
    header=0, 
    low_memory=False)

#check
print("\n#######################################\n#######################################")
print("check that sample IDs in the pedigree of latest data are all included in the ped_sample_pop file of the general high coverage dataset")
print("#######################################\n#######################################")
print(all(samples_pedigree["sampleID"].isin(samples_pedigree_pop["SampleID"])))

#check
print("\n#######################################\n#######################################")
print("explicitly test whether the columns of the two datasets are the same")
print("#######################################\n#######################################")
print("sample ID: " + str(samples_pedigree["sampleID"].equals(samples_pedigree_pop["SampleID"])))
print("father ID: " + str(samples_pedigree["fatherID"].equals(samples_pedigree_pop["FatherID"])))
print("mother ID: " + str(samples_pedigree["motherID"].equals(samples_pedigree_pop["MotherID"])))
print("sex ID: " + str(samples_pedigree["sex"].equals(samples_pedigree_pop["Sex"])))
    #differneces in parents and sex, so we have to use the ped of the latest version

#merge both datasets using the sample ID
ped_merged = samples_pedigree.merge(
    samples_pedigree_pop[["SampleID", "Population", "Superpopulation"]],
    left_on="sampleID", #use the column with IDs in sample_map
    right_on="SampleID", #use the column with IDs in pheno_data
    suffixes=["_sample_ped", "_sample_ped_pop"], #set the suffix for repeated columns
    how="inner") #only IDs included in both datasets

#checks
print("\n#######################################\n#######################################")
print("the number of rows in the merged dataset using inner, i.e., shared rows, is equals to 3202?")
print("#######################################\n#######################################")
print(ped_merged.shape[0] == 3202)
print("\n#######################################\n#######################################")
print("the sampleID columns of both datasets are identical in the merged dataset?")
print("#######################################\n#######################################")
print(ped_merged["sampleID"].equals(ped_merged["SampleID"]))

#remove the SampleID column from one of the merged datasets
ped_merged = ped_merged.drop("SampleID", axis=1)

#change some column names
ped_merged = ped_merged.rename(columns={"Population": "population", "Superpopulation": "superpopulation"})


##explore the pedigree in more detail
#As part of this publication, we sequenced 3,202 lymphoblastoid cell line (LCL) samples from the 1kGP collection, including 1,598 males and 1,604 females.
    #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 1598 males?")
print("#######################################\n#######################################")
male_samples = ped_merged.loc[ped_merged["sex"] == 1, :]
print(male_samples.shape[0] == 1598)
print(male_samples)
print("\n#######################################\n#######################################")
print("Do we have 1604 females?")
print("#######################################\n#######################################")
female_samples = ped_merged.loc[ped_merged["sex"] == 2, :]
print(female_samples.shape[0] == 1604)
print(female_samples)

#see the number of samples per superpopulation
#The 3,202 samples were drawn from 26 populations (listed in Table S1) across the following 5 continental ancestry groups: African (AFR, n = 893), European (EUR, n = 633), East Asian (EAS, n = 601), South Asian (SAS, n = 585), and American (AMR, n = 490)
    #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
afr_samples = ped_merged.loc[ped_merged["superpopulation"] == "AFR", :]
eur_samples = ped_merged.loc[ped_merged["superpopulation"] == "EUR", :]
eas_samples = ped_merged.loc[ped_merged["superpopulation"] == "EAS", :]
sas_samples = ped_merged.loc[ped_merged["superpopulation"] == "SAS", :]
amr_samples = ped_merged.loc[ped_merged["superpopulation"] == "AMR", :]
#check
print("\n#######################################\n#######################################")
print("Do we have the correct number of samples per superpopulation?")
print("#######################################\n#######################################")
print(afr_samples.shape[0] == 893)
print(afr_samples)
print(eur_samples.shape[0] == 633)
print(eur_samples)
print(eas_samples.shape[0] == 585)
print(eas_samples)
print(sas_samples.shape[0] == 601)
print(sas_samples)
print(amr_samples.shape[0] == 490)
print(amr_samples)
        #the number of samples of populations matches what I have except for the case of EAS and SAS. I think EAS and SAS are flipped, because i got the exact number but in the opposite way.
        #indeed, I have sum the number of samples of EAS and SAS according to Table S1 of the paper and matches what I have
            #EAS: 44+49+46+57+86+77+56+48+60+62=585
            #SAS: 60+71+56+47+61+46+77+69+65+49=601
        #https://www.cell.com/cms/10.1016/j.cell.2022.08.004/attachment/1e0ee9bf-5514-4660-8ff2-6d2c2be97687/mmc1.pdf
#Among the 3,202 samples, there are 602 father-mother-child trios
trios = ped_merged.loc[(ped_merged["fatherID"] != '0') & (ped_merged["motherID"] != '0'), :]
    #select samples with father and mother
    #This includes 
        #2 trios that are part of a multi-generational family, i.e., one sample is child of another sample but also is parent of another sample. This would be the case of HG00702.
        #10 trios that were split from 5 quads for the purpose of pedigree-based correction applied after haplotype phasing). I GUESS a quad would be a sample that is parent of another sample, that other sample is in turn parent of another sample. This would be the case of HG00656, which is parent of HG00702 that in turn is parent of HG00703.
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 602 trios?")
print("#######################################\n#######################################")
print(trios.shape[0] == 602)
print(trios)

#we also have duos
duos = ped_merged.loc[
    (ped_merged["fatherID"] == '0') & (ped_merged["motherID"] != '0') | 
    (ped_merged["fatherID"] != '0') & (ped_merged["motherID"] == '0'), :]
    #select samples that
        #have mother but not father OR
        #have father but not mother
    #there are also 6 parent-child duos. Where a sample only has one parent included here. For example, HG02569 only has mother (HG02568), not father.
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue
#check
print("\n#######################################\n#######################################")
print("Do we have 6 duos?")
print("#######################################\n#######################################")
print(duos.shape[0] == 6)
print(duos)

#the total number of trios and duos is 608, which matches the number of cases with parentID different from zero.
print("\n#######################################\n#######################################")
print("Do we have 608 duos and trios?")
print("#######################################\n#######################################")
print(duos.shape[0] + trios.shape[0] == 608)

#The rest of samples are unrelated.
unrelated_samples = ped_merged.loc[(ped_merged["fatherID"] == '0') & (ped_merged["motherID"] == '0'), :]
#check
print("\n#######################################\n#######################################")
print("Do we have 2594 unrelated samples?")
print("#######################################\n#######################################")
print(unrelated_samples.shape[0] == 2594)
print(unrelated_samples)
    #STRANGE THING HERE: They say "Using the Illumina NovaSeq 6000 System, we performed WGS of the original 2,504 1kGP unrelated samples and an additional 698 related samples". But we have 3202-608=2594 unrelated, instead of 2504. If they have 602+6 trios and duos, what other samples are related (Table S1 says there are 602 trios)? 2504 is the sample size of phase 3, so maybe they have added more unrelated samples?
        #ASKING JESUS
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue#supplementaryMaterial

#there is an older readme saying that the number of trios is 698 (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/README_2504_plus_additional_698_related_samples.txt), but this comes from older data, a dataset released in 2020. But they say these are samples mainly completing trios and they say there are only 602+8 trios+duos! so not sure what is happening here
    #ASKING JESUS

#we are going to use unrelated samples only. If we are calculating haplotype homozygosity by counting haplotypes, if a father and child have the same haplotype, this is not caused by selection but just by shared ancestry
print("\n#######################################\n#######################################")
print("see final pedigree data with only unrelated samples")
print("#######################################\n#######################################")
print(unrelated_samples)

#extract the distinct population names
pop_names = unrelated_samples["population"].unique()

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
                #SNP rs6054257 20 14370 G A 6 3 0.5 GTs: 1|0 1|1 0|0
                #SNP rs6054255 20 14371 G C 6 3 0.667 GTs: 0|1 1|1 1|0
                #SNP rs6040351 20 17330 T A 6 2 0.333 GTs: 0|0 1|0 1|0
            #Multiallelic SNPs:
                #SNP rs6040355 20 1110696 A G,T 6 2,2 0.333,0.333 GTs: 1|1 2|2 0|0
                #SNP rs6040356 20 1110697 A G,T 6 2,1 0.333,0 GTs: 1|1 0|0 0|0
                #SNP rs6040357 20 1110698 A G,T 6 3,3 0.5,0.5 GTs: 1|1 2|2 1|2
                #SNP rs6040358 20 1110699 A G,T 6 2,0 0.333,0 GTs: 1|1 0|0 0|.
                #SNP rs6040359 20 1110700 A G,T 5 0,3 0,0.6 GTs: 2|2 2|2 .|2
            #Exact duplicate SNPs, i.e., pos, chr, REF, alt
                #SNP rs6040360 20 1110701 A G 6 3 0.5 GTs: 1|1 0|0 0|1
                #SNP rs6040360_copy 20 1110701 A G 6 3 0.5 GTs: 1|1 1|0 0|1
            #SNP with unphased data
                #SNP rs6040361 20 1110702 A G 6 3 0.5 GTs: 1|1 0/0 0|1
            #microsatellite (indel)
                #INDEL microsat1 20 1110703 GTC G,GTCT . . . GTs: 0/1 0/2 1/1

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
    #In these new lines, the ALT is now only one, but the AC still has two fields, i.e., count for both ALTs in both lines. the 1000 genomes project SNPs that are multiallelic have their AC field just with one value so I guess they used +flags to update these fields (see below).
            #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf 
        #First scenario: rs6040355
            #the first line shows the genotypes for A and G, so a sample that is G|G would be 1|1, but a T|T sample would be 0|0, so A and T gets 0 in this line. This is ok, because the next line shows genotypes for A and T, and that T|T sample will be 1|1.
                #This case is multiallelic, so we should remove both lines by removing duplicates. We can just combine these two lines with norm --multiallelic, as it looks for snps with same pos and similar REF or ALT to merge (see when --multiallelic +snps is first used).
        #Second scenario: rs6040356
            #we could also have a case where the second ALT does not exist in our subpop, then the first line would be 0s (for A) and 1s (for G), while the second one would be all zeros, because T (1) is not present.
                #We need to retain the first row but not the second. This can be solved removing monomorphic so the second row with all 0 is removed, then pass duplicates filter, which should not affect the first line, as its "sister" row has been removed. When then we combine lines with multiallelic +snps, this remaining line should not have any other line with the same position and REF/ALT.
        #Third scenario: rs6040357
            #we could also have a case where the first line is all 1 (G) except for the cases that are T (0), zero in this case, not A (REF). In the second row, all 1 (T), except those that are G (0) in this case, not A (REF). Therefore, we do not have the REF in the subset, and only the two ALT alleles.
                #We have several ALT alleles, and no REF, so we should remove these SNPs, all lines. This case would be solved just by merging these lines with --multiallelic and then filter out those with number of ALT>2.
        #Fourth scenario: rs6040358
            #we could also have a scenario similar where the second ALT is not present in the subset, and we also have missing data.

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
            #Similarly, rs6040351 has AN=6, 2 copies of the allele and frequency of 2/6=0.33, but after removing the third sample (1|0), AN becomes 4, AC is 1, BUT AF remains 0.333, when it should be 0.25.
        #Indeed, there is a command (--no-update) to avoid (re)calculating INFO fields for the subset, currently including only INFO/AC and INFO/AN. Therefore, they are clearly saying that AN and AC are the only INFO fields updated for the subset.
            #You can update other fields using bcftools +fill
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
print("see monomorphic snps")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --include 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #select those variants for which the number of ALT|ALT or REF|REF is equal to the number of samples (3 in this case), so there is no variability.
            #we use count() to count the number of homozygous AA and RR
                #count is a function that gives you the number of cases, it can be used on FORMAT tags (over samples) and INFO tags (over vector fields).
                #we are counting number of genotypes that satisfy the condition, given that GT is a format field, we are using count over samples (see below).
            #We look for genotypes including both copies of reference (GT="RR") or alternative (GT="AA") alleles
                #sample genotype: reference (haploid or diploid), alternate (hom or het, haploid or diploid), missing genotype, homozygous, heterozygous, haploid, ref-ref hom, alt-alt hom, ref-alt het, alt-alt het, haploid ref, haploid alt (case-insensitive)
                #It seems this consider the genotype columns, not AC field, as we are using GT. Indeed, I have purposely set AC field for rs6040356 as 2,1 indicating that there is 1 copy of the second ALT allele, while this is not true in the genotype. The expression still excludes this line with all 0|0 in the genotype, so we are good.
                #Also, if you have missing (.), this is not included in GT=AA or GT=RR, because not all genotypes are AA or RR. Because of this, we removed first any SNP with missing genotypes.
                    #https://www.biostars.org/p/360620/
            #We also look for genotypes that are missing "."
                #missing genotypes can be matched regardless of phase and ploidy (".|.", "./.", ".", "0|.", "1|.") using these expressions
                    #GT="mis", GT~"\.", GT!~"\."
            #This is second condition that has to be met, we are interested in selecting also those cases where everything is the same allele except missing cases. 
                #We have to use "|" instead of "||". 
                    #You look to a given sample and check in the same sample whether the genotype is A|A, A|. or .|A. 
                    #If it satisfies any of these conditions, then we can count it
                    #Therefore, the conditions are checked within the same sample, so we have to use "|" instead of "||".
                    #In case it is not clear see and example from the manual:
                        #Say our VCF contains the per-sample depth and genotype quality annotations and we want to include only sites where one or more samples have big enough coverage (DP>10) and genotype quality (GQ>20). The expression -i 'FMT/DP>10 & FMT/GQ>20' selects sites where the conditions are satisfied WITHIN THE SAME SAMPLE: 'FMT/DP>10 & FMT/GQ>20'
                        #On the other hand, if we need to include sites where both conditions met but not necessarily in the same sample, we use the && operator rather than &: 'FMT/DP>10 && FMT/GQ>20'
                            #http://samtools.github.io/bcftools/howtos/filtering.html
            #check for equality with N_SAMPLES. This is the number of samples and it is calculated on the fly by bcftools
                #variables calculated on the fly if not present: number of alternate alleles; number of samples; count of alternate alleles; minor allele count (similar to AC but is always smaller than 0.5); frequency of alternate alleles (AF=AC/AN); frequency of minor alleles (MAF=MAC/AN); number of alleles in called genotypes; number of samples with missing genotype; fraction of samples with missing genotype; indel length (deletions negative, insertions positive)
            #details about the use of expressions can be found here
                #https://samtools.github.io/bcftools/bcftools.html#expressions
        #we get all SNPs that are monomorphic. This includes 
            #all 0
                #SNP rs6040356 20 1110697 A T 6 2,1 0 GTs: 0|0 0|0 0|0
            #all 0 or all 1 except missing in the first or second position. This includes cases with multiallelic SNPs, where only the lines with the same alelle (+ missing) are considered monomorphic.
                #SNP rs6040358 20 1110699 A T 6 2,0 0 GTs: 0|0 0|0 0|.
                #SNP rs6040359 20 1110700 A G 5 0,3 0 GTs: 0|0 0|0 .|0
                #SNP rs6040359 20 1110700 A T 5 0,3 0.6 GTs: 1|1 1|1 .|1

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
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")

#
print("\n#######################################\n#######################################")
print("remove SNPs with missing")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
    #remove missing genotypes
        #--genotype
            #Require one or more hom/het/missing genotype or, if prefixed with "^", exclude such sites
    #Lines of multiallelic rs6040358 have been removed because the last sample has a missing genotype. Note that rs6040358 has not been excluded as monomorphic as it has 1 and 0, but because it has missing.

#
print("\n#######################################\n#######################################")
print("remove now the exact duplicates")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
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
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
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
            #the normal scenario is to have same REF and different ALTs, so it just add the different ALTs as different field in ALT column
            #I have checked that it also merge when REF is different but ALT is the same, in that case it takes the REF of the first SNP.
            #It merges with REF and ALT are the same.
            #when REF and ALT are different, no merging is done and get and error.
            #This should not be a problem for us. 
                #If we have SNPs with the same REF and ALT and position, these will be removed with dup (see below).
                #If we have snps with same position and same ALT, they will get merged and then removed by our filters. We do not have want these SNPs.
                #if REF and ALT different, we get error and we will see.
        #some genotypes appear as not phased within the new fused multiallelic SNPs, but that is not a problem because these will be removed anyway.

print("\n#######################################\n#######################################")
print("check what happens if we remove the third sample and then rs6040355 have no REF anymore, only two ALTs, thus it is not monomorphic")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #Despite losing the REF, both lines of rs6040355 are merged into a multiallelic SNP, thus it can be removed with --max-alleles 2. This makes sense because --multiallelics +snps looks for snps with the same position and same REF or ALT.

#
print("\n#######################################\n#######################################")
print("now we can remove those SNPs with more than 1 allele in ALT")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
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
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
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
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
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
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
    bcftools annotate \
        --remove INFO,^FORMAT/GT | \
    bcftools +fill-tags \
        -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
    bcftools view \
        --no-header")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
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
        --no-header")
    #I have checked that the genotypes remain the same despite removing all previous INFO fields and all FORMAT fields (except GT) and then adding new INFO fields, so we are good here.

#
print("\n#######################################\n#######################################")
print("see the new header")
print("#######################################\n#######################################")
run_bash(" \
    bcftools norm \
        --multiallelic -snps \
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
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
        data/dummy_vcf_files/dummy_example.vcf | \
    bcftools view \
        --samples NA00001,NA00002 | \
    bcftools view \
        --types snps | \
    bcftools view \
        --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
    bcftools view \
        --genotype ^miss |\
    bcftools norm \
        --rm-dup exact | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools view \
        --phased | \
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

#SUMMARY: 
    #With all these commands, we have recreated the scenario we have in 1KGP data, with multiallelic SNPs separated into different lines, select some samples, we then select snps, remove those with missing genotypes, remove exact duplicates (this does not touch different lines of the same multiallelic snp because they have different REF/ALT), exclude those SNPs that have the same allele for all samples, and we can do that using RR (REF|REF) and AA (A|A), because there are no missing (.) genotypes, as they have been removed. Then we combine all lines of each multiallleic snp and now they have ALT column with several alleles, so we can filter them using --max-alleles 2. Add last filter for selecting phased data only. We can also use bcftools +fill-tags to update important fields for each SNP, so if a SNP was multiallelic, but it is not multiallelic in the subset population (i.e., only REF and 1 ALT), we no longer will have two allele frequencies, two counts.... for the remainder biallelic SNP in the subset.

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
            if [[ $n_samples -eq " + str(ped_merged.shape[0]) + " ]]; then \
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



    #######################
    # extract samples IDs #
    #######################

    #select the sample IDs for the selected population
    subset_pop = unrelated_samples.loc[unrelated_samples["population"] == selected_pop, :]

    #reset index
    subset_pop = subset_pop.reset_index(
        drop=True)
        #drop=True to avoid adding the index as a column

    #check
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we only have the selected pop and unrelated samples")
    print("#######################################\n#######################################")
    print(subset_pop["population"].unique()[0] == selected_pop)
    print(all((subset_pop["fatherID"] == "0") & (subset_pop["motherID"] == "0")))

    #select the sample IDs
    selected_samples = subset_pop["sampleID"]

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
    print("chr " + selected_chromosome + " - " + selected_pop + ": see monomorphic")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --include 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #select those variants for which the number of ALT|ALT or REF|REF is equal to the number of samples, so there is no variability. All ALT + missing or all REF + missing are also considered. See dummy example for further details.

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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove SNPs with missing")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #exclude (^) those genotypes with missing. See dummy example to check behaviour.
                #there are multiple ways to filter by missing
                    #https://www.biostars.org/p/362060/

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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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

    #show now only biallelic SNPs that are phased for all samples
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs that are phased for all samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        bcftools view \
            --exclude 'COUNT(GT=\"AA\" | GT=\"mis\")=N_SAMPLES || COUNT(GT=\"RR\" | GT=\"mis\")=N_SAMPLES' |\
        bcftools view \
            --genotype ^miss | \
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
        if [[ $n_snps_vcf == $n_snps_hap ]];then \
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
        if [[ $nfields_hap_nlines == 1 ]]; then \
            nsamples_hap=$(($nfields_hap/2)); \
            nsamples_bcftools=$( \
                cat results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt | \
                grep -c '^'); \
            if [[ $nsamples_hap == $nsamples_bcftools ]]; then \
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

    #WHEN PREPARING MAP FILE, REMEMBER THAT VCF IS 1 BASED!

    #


    #loading the position and chromosome of each snp, VCf file or better impute file

    "./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw"

    run_bash(" \
        bcftools view \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools query \
            --format '%ID %CHROM %POS\n' |\
        head -5 > eso.txt; \
        cat eso.txt")

    snp_map_raw = pd.read_csv("eso.txt", sep=" ", header=None)
    snp_map_raw

    snp_map_raw.rename()


    #last conversion hap, I guess non-billalte and non-alt shoudl be removed to meet impute requeriments
    #Hap file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.hap.gz
    #Sample file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.samples
    #938126 records written, 0 skipped: 0/0/0 no-ALT/non-biallelic/filtered





    #if you take the VCF file, filter and then count lines, you get 938126 records, but then say that total is 945919, as when writting the file. Maybe this is the total originally, CHECK
    #938126
    #Lines   total/split/realigned/skipped:  945919/0/0/0
    #Lines   total/split/realigned/skipped:  945919/0/0/0

    #
    #In [4]: run_bash(" \
    #...:     bcftools view \
    #...:         " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
    #...:     wc -l")
    #5759173


    #In [5]: run_bash(" \
    #...:     bcftools view \
    #...:             --samples " + ",".join(selected_samples) + " \
    #...:         " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
    #...:     wc -l")
    #...: 
    #5759173