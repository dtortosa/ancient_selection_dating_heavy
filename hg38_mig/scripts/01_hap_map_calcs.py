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
    elif ("Warning" in complete_process.stderr) | ("warning" in complete_process.stderr):

        #print the standard output without "\n" and other characters
        print(complete_process.stdout)

        #print the standard error without stopping
        print("WARNING! FALSE!: " + complete_process.stderr)

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
    else:
        #print the standard error and stop
        raise ValueError("ERROR! FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

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
run_bash(
    "mkdir \
        -p ./results/hap_map_files")
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
    #They do not give an exact number of unrelated samples. I think it is ok if 3202-608 gives 2594 instead 2504. 2504 is the sample size of phase 3, but there are more samples in this new phase of the project.
        #https://www.cell.com/cell/fulltext/S0092-8674(22)00991-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422009916%3Fshowall%3Dtrue#supplementaryMaterial

#there is an older readme saying that the number of trios is 698 (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/README_2504_plus_additional_698_related_samples.txt), but this comes from older data, a dataset released in 2020. In contrast, our data here matches with the paper of the latest release

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
        print("chr " + selected_chromosome + ": the chromosome name in the vcf file is correct?")
        print("#######################################\n#######################################")
        run_bash(" \
            chr_vcf=$( \
                bcftools query \
                    " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
                    --format '%CHROM\n' | \
                uniq); \
            if [[ $chr_vcf == 'chr" + selected_chromosome + "' ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi")
            #with bcftools, query for the chromosome of each variant and get the unique cases
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
            #according to the specific readme of the dataset we are using (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf), prior phasing they applied several filters, being one of them that all variants has to be PASS for FILTER. I understand that, because of this, all variants in the chromosome 1 have now ".", being this the unique character. Let's check this for all chromosomes:
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
    selected_samples.to_csv(
        "results/hap_map_files/samples_to_bcftools_" + selected_pop + ".txt", 
        sep="\t", 
        header=None, 
        index=False)



    ##################################
    # operations within selected pop #
    ##################################

    #see genotypes of few samples
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": see genotypes of few samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #select a few samples and then select biallelic snps
                #IMPORTANT: If filtering (by number of alleles) and subseting (by sample) is in the same line (command), the filtering will be done first. Therefore, you could select SNPs that have 2 alleles when considering the 26 pops, but that are monomorphic (1 allele) for the selected population. Because of this, we have to first subset by sample and then filter by number of alleles whitin the selected samples in separated commands.
            #query the genotype of each variant

    #now show only snps
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": show now only SNPs")
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
            #query multiple fields for each snp

    #remove SNPs with missing
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": remove SNPs with missing")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #exclude (^) those genotypes with missing. See dummy example to check behaviour.
                #there are multiple ways to filter by missing
                    #https://www.biostars.org/p/362060/

    #remove also exact duplicates
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": remove also exact duplicates")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools norm \
            --rm-dup exact | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. See dummy example for behaviour.

    #exclude monomorphic
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": exclude monomorphic")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools norm \
            --rm-dup exact | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\")=" + str(len(selected_samples)) + " || COUNT(GT=\"RR\")=" + str(len(selected_samples)) + "' |\
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #exclude those variants for which the number of ALT|ALT or REF|REF is equal to the number of samples, so there is no variability. See dummy example for further details.

    #combine multiallelic SNPs in one line and select them
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": combine multiallelic SNPs in one line and select them")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools norm \
            --rm-dup exact | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\")=" + str(len(selected_samples)) + " || COUNT(GT=\"RR\")=" + str(len(selected_samples)) + "' | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --min-alleles 3 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
        #we have SNPs with the same position and chromosome but different ALT alleles. According to the 1KGP paper, they split multiallelic SNPs in different lines, so this makes sense. See below for a way to deal with this
        #combine the different lines of each multiallelic SNP into one line and update the ALT column to include the different ALT alleles.
            #--multiallelic +snps
                #Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]
            #I have checked that this command does NOT combine different lines that have the same chrom, pos, REF and ALT. These are considered exact duplicates and removed (see dummy example).
            #by combining the different lines of each multiallelic SNP, we update the ALT column, which now includes several ALT alleles. Therefore, we can select/filter out these SNPs with --max-alleles. In other words, --multiallelic +snps lets us target those lines that have the same chromosome and position, but different ALT
        #select those SNPs with 3 or more alleles in ALT

    #now show only biallelic snps
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": show now only biallelic SNPs")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools norm \
            --rm-dup exact | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\")=" + str(len(selected_samples)) + " || COUNT(GT=\"RR\")=" + str(len(selected_samples)) + "' | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #we can filter out those SNPs with more than 2 alleles in ALT because --multiallelic +snps combined lines of the same multiallelic SNP into one single line, updating ALT column. Consequently, in chr1, 1:10453:A:C and 1:10452:A:G are filtered out. 
                #--min-alleles/--max-alleles INT  
                    #Minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)

    #show now only biallelic SNPs that are phased for all samples
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + ": show now only biallelic SNPs that are phased for all samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools view \
            --genotype ^miss | \
        bcftools norm \
            --rm-dup exact | \
        bcftools view \
            --exclude 'COUNT(GT=\"AA\")=" + str(len(selected_samples)) + " || COUNT(GT=\"RR\")=" + str(len(selected_samples)) + "' | \
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


        ###POR AQUIIII
        #check COUNT(GT=\"AA\")=" + str(len(selected_samples)





        #IMPORTANT NOTES AND CHECKS ABOUT FILTERING
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": We are going to use dummy vcf files in order to check that the subset and filtering methods we are using are correct: ")
        print("#######################################\n#######################################")
        
        #see dummy snps
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": see dummy SNPs: ")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools view \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #we have 4 SNPs, with different number of alleles, frequencies, etc...
                #SNP rs6054257 20 14370 G A 6 3 0.5 GTs: 1|0 1|1 0|0
                #SNP rs6040351 20 17330 T A 6 2 0.333 GTs: 0|0 1|0 1|0
                #SNP rs6040355 20 1110696 A G,T 6 2,2 0.333,0.333 GTs: 1|1 2|2 0|0
                #SNP rs6054255 20 14371 G C 6 3 0.667 GTs: 0|1 1|1 1|0

        #primero 
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": select only the first and second sample, to reduce the number of alleles of the multiallelic snp to 2, also check how AN and AC changes but AF, ALT... ")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools view \
                --samples NA00001,NA00002 \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --include 'TYPE=\"snp\"' | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #now, the multiallelic snp only have allele 1 and 2 but no 0, in the subset
            #AN AND AC ARE UPDATED: you can see in the last snp how the total number of alleles is 4 now, and the count for ALT is 2, while in the whole set is 6 and 3.
            #AF IS NOT UPDATED: The frequency of ALT should be now 0.75, because is 3/4, but the old frequency is retained.

        #now check if the --max-allele filter works, it should not because ALT is not updated after subsetting
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": now check if the --max-allele filter works, it should not because ALT is not updated after subsetting: ")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools view \
                --samples NA00001,NA00002 \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --include 'TYPE=\"snp\"' \
                --max-alleles 2 \
                --min-alleles 2 | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #we have removed the third sample, so the multiallelic snp is no longer multiallelic because there are only two alleles in the subset, but --max-alleles 2 still removes this snp because ALT remains the same despite the subsetting.
            #therefore, we cannot just use this to filter a subset of samples

        #try --trim-alt-alleles
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": try --trim-alt-alleles")
        print("#######################################\n#######################################")
        run_bash(" \
            bcftools view \
                --samples NA00001,NA00003\
                --trim-alt-alleles \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --include 'TYPE=\"snp\"' \
                --max-alleles 2 \
                --min-alleles 2 | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #--trim-alt-alleles
                #remove alleles not seen in the genotype fields from the ALT column. Note that if no alternate allele remains after trimming, the record itself is not removed but ALT is set to ".". If the option -s or -S is given, removes alleles not seen in the subset. INFO and FORMAT tags declared as Type=A, G or R will be trimmed as well.
                    #https://samtools.github.io/bcftools/bcftools.html
                #therefore, we have to subset the samples and apply this option in the same step, if not we do not get rid of snps with the lacking allele?
                    #CHECK THIS

            #PROBLEM, if the ancestral dissapear, the REF remains the same but now you have 0 for one of the alt..
                #not sure if this is very frequent...

            #ALSO, multiallelic snps in 1000 genomes are split in several lines, so this will not work







        #POR AQUII

        #it is important to be sure that the fields you are using for filtering, are updated after subseting samples. Of course, type="snp" will be always "snp" irrespectively of the samples we select, but this is not the case of the number of alleles, because you can have SNPs with 3 alleles considering all 26 populations, but then in GBR they can have only 2 or 1. We are interested in SNPs that are biallelic within the selected populaion, CHECK DAVID
        #The same goes for phasing and genotype ^miss. You have to be sure that these filters only consider data from the filtered samples, not fixed data in fields that are not updated.






        #CHECK ABOUT PHASING AND NUMBER OF ALLELES?
            #PS is hpasing state and is NOT a field

        #https://www.biostars.org/p/216114/



        #por aquiiiii
        #next lines do the most updated filters
        #clean previus lines and update the lines with selected subset samples


        #CHECK THAT WE HAVE MULTIALLEIC SPLIT IN 1000 GENOMES



        ###PON DUMMY AL PRINCIPIO QUE SE CORRA UNA SOLA VEZ, O INCLUSO, SACALLO DE LA FUNCION!



        #see dummy vcf file including multiallelic snps, but split those in different lines
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
            #when a multiallelic snp (e.g., REF=A, ALT=G,T) is split in several lines, having all the same ID and position, thus being considered duplicates. This was done in 1000 genomes data. They also add 1 to the position of each new line to do the phasing, but then they put all lines in the same position.
                    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf 
                #First scenario: rs6040355
                    #the first line shows the genotypes for A and G, so a sample that is G|G would be 1|1, but a T|T sample would be 0|0, although 0 should be only A. Indeed, A|A would also be 0|0. This is ok, because the next line shows genotypes for A and T, and that T|T sample will be 1|1.
                        #This case is multiallelic, so we should remove both lines by removing duplicates.
                #Second scenario: rs6040356
                    #we could also have a case where the second ALT does not exist in our subpop, then the first line would be 0s (for A) and 1s (for G), while the second one would be all zeros, because T (1) is not present.
                        #We need to retain the first row but not the second. This can be solved removing monomorphic so the second row with all 0 is removed, then pass duplicates filter, which should not affect the first line, as its "sister" row has been removed.
                #Third scenario: rs6040357
                    #we could also have a case where the first line is all 1 (G) except for the cases that are T (0), zero in this case, not A (REF). In the second row, all 1 (T), except those that are G (0) in this case, not A (REF). Therefore, we do not have the REF in the subset, and only the two ALT alleles.
                        #We have several ALT alleles, and no REF, so we should remove these SNPs, all lines. This case would be solved just by removing duplicated SNPs, i.e., similar ID and position.
                #Fourth scenario: rs6040358
                    #we could also have a scenario similar to the previous one but with a missing for the second chromosome of the last sample, so the second line is no longer 0|0, but it has also 0|.
            #we also have a exact duplicate, i.e., same pos, chr, REF, alt (rs6040359)
                #we can select only the appearance and remove the next using bcftools norm --rm-dup exact
            #we also have a microsatellite

        #important note about AC here.
            #the output of the previous line shows that the number of alleles and the counts of ALT are not updated in INFO because the last multiallelic SNP (rs6040358) lacks one allele for one sample but thte AN is still 6, not 5.

        #select only SNPs
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

        #remove first the SNPs with missing
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss |\
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #remove missing genotypes
                #--genotype
                    #Require one or more hom/het/missing genotype or, if prefixed with "^", exclude such sites
            #The last two lines corresponding with the last multiallelic SNP have been removed because the last sample has a missing genotype.

        #remove now the exact duplicates
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss |\
            bcftools norm \
                --rm-dup exact | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #remove those snps that are exact duplicates, meaning identical chr, pos, ref, and alt. As you can see for rs6040359, it does not matter the ID or the genotypes, this command only targets chr, pos, ref, and alt.
                #https://github.com/samtools/bcftools/issues/1089
                #https://samtools.github.io/bcftools/bcftools.html#norm
            #when using --multiallelic +snps, you also get rid of these duplicated SNPs because they are not considered different alleles of a multiallelic snps given that they have the same ALT. We add this step anyways just in case to be sue we remove first exact duplicates.

        #exclude also those monomorphic snps
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss |\
            bcftools norm \
                --rm-dup exact | \
            bcftools view \
                --exclude 'COUNT(GT=\"AA\")=" + str(3) + " || COUNT(GT=\"RR\")=" + str(3) + "' |\
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
            #exclude those variants for which the number of ALT|ALT or REF|REF is equal to the number of samples, so there is no variability.
                #we use count() to count the number of homozygous and then exclude those
                #It seems this consider the genotype columns, not AC field. Indeed, I have purposely set AC field for rs6040356 as 2,1 indicating that there is 1 copy of the second ALT allele, while this is not true in the genotype. The expression still excludes this line with all 0|0 in the genotype, so we are good.
                #Also, if you have missing (.), this is not included in GT=AA or GT=RR, because not all genotypes are AA or RR. Because of this, we removed first any SNP with missing genotypes.
                    #https://www.biostars.org/p/360620/

        #combine lines of multialleic snps in just one line
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss |\
            bcftools norm \
                --rm-dup exact | \
            bcftools view \
                --exclude 'COUNT(GT=\"AA\")=" + str(3) + " || COUNT(GT=\"RR\")=" + str(3) + "' |\
            bcftools norm \
                --multiallelic +snps | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
                #by combining the different lines of each multiallelic SNP, we update the ALT column, which now includes several ALT alleles. Therefore, we can filter out these SNPs with --max-alleles.
                #--multiallelic +snps
                    #split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps
                        #https://samtools.github.io/bcftools/bcftools.html#norm
                #some genotypes appear as not phased within the new fused multiallelic SNPs, but that is not a problem because these will be removed anyway.

        #now we can remove those SNPs with more than 1 allele in ALT
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss | \
            bcftools norm \
                --rm-dup exact | \
            bcftools view \
                --exclude 'COUNT(GT=\"AA\")=" + str(3) + " || COUNT(GT=\"RR\")=" + str(3) + "' |\
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
            
        #select only phased data
        run_bash(" \
            bcftools norm \
                --multiallelic -snps \
                data/dummy_vcf_files/dummy_example.vcf | \
            bcftools view \
                --types snps | \
            bcftools view \
                --genotype ^miss | \
            bcftools norm \
                --rm-dup exact | \
            bcftools view \
                --exclude 'COUNT(GT=\"AA\")=" + str(3) + " || COUNT(GT=\"RR\")=" + str(3) + "' |\
            bcftools norm \
                --multiallelic +snps | \
            bcftools view \
                --max-alleles 2 \
                --min-alleles 2 \
                --phased | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #SUMMARY: With all these commands, we have recreated the scenario we have in 1KGP data, with multiallelic SNPs separated into different lines, we then select snps, remove those with missing genotypes, remove exact duplicates (this does not touch different lines of the same multiallelic snp because they have different REF/ALT), exclude those SNPs that have the same allele for all samples, and we can do that using RR (REF|REF) and AA (A|A), because there are no missing (.) genotypes, as they have been removed. Then we combine all lines of each multiallleic snp and now they have ALT column with several alleles, so we can filter them usin --max-alleles 2. Add last filter for selecting phased data only and finished.





    #remove next lines and combine with the previous



    ##################################
    # operations within selected pop #
    ##################################

    #see genotypes
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first genotypes of biallelic snps for selected samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --include 'TYPE=\"snp\"' \
            --max-alleles 2 \
            --min-alleles 2 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -5")
            #select of the samples of the selected population using a list of them separated by ","

    #see complete data for some of the selected samples
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see complete data for some of the selected samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --include 'TYPE=\"snp\"' \
            --max-alleles 2 \
            --min-alleles 2 \
            --no-header | \
        head -3")
        #view complete row, i.e., all INFO fields, genotypes...
            #for selected samples
            #for each variant with type=SNP and 2 alleles in the selected samples
            #and avoid header
        #show only 10 first

    #in case you want to save the filtered dataset, but it is slow the file is very big, so it is not worth it, we can directly go to hap
    #run_bash(
    #    "bcftools view \
    #        " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    #        -i 'TYPE=\"snp\"' \
    #        -H \
    #        -s " + ",".join(selected_samples) + " \
    #        -o data/" + selected_pop + ".bcf.gz")

    #convert VCF file to hap file    
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": convert to hap file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --include 'TYPE=\"snp\"' \
            --max-alleles 2 \
            --min-alleles 2 \
            --phased \
            --genotype ^miss | \
        bcftools convert \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
            --hapsample ./results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw")
        #for those samples of the selected pop
        #select those variants that are 
            #biallelic 
            #snps
            #that are phased for all the selected samples
            #that do not have missing data
        #convert to hap format saving in results and using the name of the chromosome and the pop
            #--vcf-ids 
                #when this option is given, the second column is set to match the ID column of the VCF.
                #Without the option, the format follows https://www.cog-genomics.org/plink/2.0/formats#haps with ids (the second column) of the form "CHR:POS_REF_ALT[_END]"
                #https://www.htslib.org/doc/bcftools.html#convert
        #IMPUTE2 hap format
            #https://www.cog-genomics.org/plink/2.0/formats#haps
            #we are going to use the reference panel haplotype file format for IMPUTE2. 
            #This is a text file with no header line, and either 2N+5 or 2N fields where N is the number of samples. In other words, you can have 5 initial columns with data about the variants, or just the happlotype columns. In the former case, the first five columns are:
                #Chromosome code
                #Variant ID
                    #This is in the format CHROM:POS_REF_ALT to prevent strand swaps.
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
        "gunzip -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        head -2")
    
    #clean
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove first columns of hap file to leave only haplotype columns")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
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
    #create a long string with the index of each column of hap_raw file
    #these indexes will be used by awk to select the corresponding columns
    fields_selected_samples = "".join(
        ["$" + str(i) + "," if i != selected_samples.shape[0]*2+5 else "$" + str(i) for i in range(6, selected_samples.shape[0]*2+5+1, 1)])
        #from index 6 (avoiding non-genotype columns) to the index of the last column, i.e., 2 columns times the number of samples plus 5 (because of the non-genotype columns) and 1 (because index 0 is 1)
        #if the index is the last one
            #add $ to the index and then comma
        #else
            #it is the last one so we do not need add comma
        #join all strings
    #do comparison
    run_bash(
        "gunzip \
            -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
        awk -F ' ' '{print " + fields_selected_samples + "}' > \
        results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap; \
        gunzip \
            -kf results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz; \
        file1='results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap'; \
        file2='results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw_clean_check.hap'; \
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
            #remember that $? gives the return value of the last run command.
                #For example, 
                    #ls somefile
                    #echo $?
                    #If somefile exists (regardless whether it is a file or directory), you will get the return value thrown by the ls command, which should be 0 (default "success" return value). If it doesn't exist, you should get a number other then 0. The exact number depends on the program.
                #https://stackoverflow.com/a/6834572/12772630
            #the return value of cmp will be 0 if the two files are identical, if not, then we have differences between the files
                #https://stackoverflow.com/a/53529649/12772630
        #remove the new files




    ##POR AQUII

    #filters applied in vcftools
        #--max-alleles 2 --min-alleles 2 --max-missing 1 --phased
            #https://vcftools.sourceforge.net/man_latest.html
        #--max-alleles / --min-alleles
            #Include only sites with a number of alleles greater than or equal to the "--min-alleles" value and less than or equal to the "--max-alleles" value. One of these options may be used without the other.
                #For example, to include only bi-allelic sites, one could use: vcftools --vcf file1.vcf --min-alleles 2 --max-alleles 2
            #in bcftools
                #-m/M, --min-alleles/--max-alleles INT  
                    #Minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
        #--phased
            #Excludes all sites that contain unphased genotypes
            #in bcftools
                #-p/P, --phased/--exclude-phased
                    #Select/exclude sites where all samples are phased
        #--max-missing 1
            #Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
            #in bcftools
                #

        #bcftools filter --include 'AN=2*N_SAMPLES' [VCF/BCF]
            #https://www.biostars.org/p/362060/

        #add this as a second condition in --include
            #https://github.com/samtools/bcftools/issues/822

        ##CHECK WHY YOU SELECTED THESE FITLERS OF VCFTOOLS WHEN CREATING HAP FILES


    #YOU NEED TO SELECT SNPS WITH A MXIMUMG OF 2 ALLELES IN ORDER TO AVOID MULTIALLELIC, WE DID THAT WITH IHS
        #this is done WITHIN EACH POPUALTION, NOT CONSIDERING ALL OF THEM !!! 
         #--max-alleles 2
         #ASK JESUS AND DAIVD  
         #https://www.biostars.org/p/141156/

    #CHECK MAF AND FILTERS IN SLIM
        #check filters applied in 1000 genomes project


        #-p/P, --phased/--exclude-phased        Select/exclude sites where all samples are phased

    #you can add mask with bcftools filter
        #    -M, --mask-file [^]FILE        Soft filter regions listed in a file, "^" to negate



    #CHCK THIS!!!
        #5013617 records written, 745443 skipped: 0/0/745443 no-ALT/non-biallelic/filtered




    #check we are only selecting SNPs
    os.system(
        "n_no_snps=$( \
            bcftools query \
                --include 'TYPE=\"snp\"' \
                --samples " + ",".join(selected_samples) + " \
                --format '%TYPE\n' \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            grep -v 'SNP' | \
            wc -l) ; \
        if [ $n_no_snps -eq 0 ]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
        #make a query that
            #includes only snps and the selected samples
            #returns the variant type, which should be only snp
        #from the list of variant types
            #select those that do not include SNP (-v flag)
            #https://stackoverflow.com/questions/29489050/unix-command-to-get-lines-not-containing-certain-text
        #the total number of these lines should be zero


    #check that generated hap file has the correct number of SNPs, i.e., rows
    os.system(
        'n_unique_snps=$( \
            bcftools query \
                --include "TYPE=\'snp\'" \
                --samples ' + ','.join(selected_samples) + ' \
                --format "%ID\n" \
                ' + input_vcfs_path + '/1kGP_high_coverage_Illumina.chr' + selected_chromosome + '.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            uniq | \
            wc -l); \
        echo "We have" "${n_unique_snps}" "unique SNPs in the query of VCF file"; \
        n_lines_hap=$( \
            gunzip \
                -c \
                results/hap_map_files/chr' + selected_chromosome + '_' + selected_pop + '_IMPUTE2.hap.gz | \
            wc -l); \
        echo "We have" "${n_lines_hap}" "lines in the hap file"; \
        echo "Both numbers are the same?"; \
        if [ $n_unique_snps -eq $n_lines_hap ]; then \
            echo "TRUE"; \
        else \
            echo "FALSE"; \
        fi')
        #calculate the number of unique snps
            #querying variants 
                #with type=snp
                #for selected samples
                #get the ID of the variants
            #select unique cases
            #count them
        #calculate the number of lines in the previously generate hap file
            #gunzip with "c" flag to simply write the output stream to stdout. This will leave the compressed file untouched. In other words, you do not create a file, but just show the content of the compressed file into the terminal. We will pipe that to wc
                #https://superuser.com/questions/135329/count-lines-in-a-compressed-file
            #count lines of the file
        #check that both numbers are the same

    #do we have the correct number of samples in the hap file?
    os.system(
        'nfields_hap=$( \
            gunzip \
                -c results/hap_map_files/chr' + selected_chromosome + '_' + selected_pop + '_IMPUTE2.hap.gz | \
            awk -F " "  "{print NF}" | \
            sort | \
            uniq); \
        nfields_hap_nlines=$( \
            echo -n "$nfields_hap" | \
            grep -c "^"); \
        if [ "${nfields_hap_nlines}" -eq 1 ]; then \
            nsamples_hap=$(($nfields_hap/2)); \
            nsamples_bcftools=$( \
                cat results/hap_map_files/samples_to_bcftools_' + selected_pop + '.txt | \
                grep -c "^"); \
            if [ "${nsamples_hap}" -eq "${nsamples_bcftools}" ]; then \
                echo "TRUE"; \
            else \
                echo "FALSE"; \
            fi\
        fi')
        #calculate the number of fields/columns in the hap file
            #decompress the hap file and send the output to stdout (flag c of gunzip)
            #use the interpreter of awk language (awk), indicating the delimiter of the columns (-F) and asking to print the number of field. You could also ask for column 1 typing $1... NF will give the number of columns or fields.
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
            #calculate the number of samples from the .txt file generated.
                #load the file with cat and count the number of lines that start with any character
            #check that the number of samples according to the hap file and the txt are the same.



    sample_list_from_hap = pd.read_csv(
        './results/hap_map_files/IBS.samples',
        header=0, 
        sep=' ', 
        low_memory=False)

    sample_list_from_hap_clean = sample_list_from_hap.loc[(sample_list_from_hap["ID_1"] != '0') | (sample_list_from_hap["ID_2"] != '0'), :]

    sample_list_from_hap_clean = sample_list_from_hap_clean.reset_index()

    sample_list_from_hap_clean["ID_1"].equals(selected_samples)
    sample_list_from_hap_clean["ID_2"].equals(selected_samples)



    #check for duplicates in snp ID and position
        #https://www.biostars.org/p/274824/


    #more CHECKS??
        #https://www.biostars.org/p/270381/


    #ask jesus about the filters he applied besides the ones applied by 1000 genomes project according to readme file.

    #ON FILTERING SNPS WITH BCFTOOLS
        #https://www.htslib.org/workflow/filter.html

    #BCFTOOLS CHEATSHEET
        #https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b#file-bcftools-cheat-sheet


    #filter by accesibility mask?



##POR AQUIII

#ASK JESUS ABOUT 
    #THE MASKS
    #the pedigree used
#ASK ENARD ABOUT REMOVING RELATED INDIVIDUALS
    #I should remove indels and structural variations, right?
    #¿Estás usando los individuos emparentados? Hasta ahora yo he trabajado con unrelated, pero en el último dataset han metido también 698 que son parientes de los 2504 originales. No sé si lo has hablado con David, pero imagino que habrá que eliminar estos individuos ya que lógicamente sus genotipos van a estar muy correlacionados con individuos ya incluidos en la muestra. Si, por ejemplo, estás calculando la homozigosidad, estarías metiendo secuencias muy correlacionadas que podrían inflar la homozigosidad por parentesco, no porque haya habido selección en la población y ha hecho esa haplotipo más frecuente.



#WHEN PREPARING MAP FILE, REMEMBER THAT VCF IS 1 BASED!

for index, variant in enumerate(VCF("data/1KGP_vcf_files/1kGP_high_coverage_Illumina."+selected_chromosome+".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")): # or VCF('some.bcf')
    if(index == 0):
        print(len(variant)) # e.g. REF='A', ALT=['C', 'T']

#each column is the genotype of a sample in the same order than in the pedigree?
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/README_111822.pdf
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf

    import vcf
    vcf_reader = vcf.Reader(open("data/1KGP_vcf_files/1kGP_high_coverage_Illumina."+selected_chromosome+".filtered.SNV_INDEL_SV_phased_panel.vcf.gz", 'r'))
    for record in vcf_reader:
        print record




    #create temporary folder to save
    import tempfile
    temp_dir = tempfile.TemporaryDirectory()
    #print(temp_dir.name)
    
    #to remove the dir
    #temp_dir.cleanup()
        #https://stackoverflow.com/questions/3223604/how-to-create-a-temporary-directory-and-get-its-path-file-name
    
    #read only the zip and get list of files in it
    import zipfile
    zipdata = zipfile.ZipFile("data/1KGP_vcf_files/1kGP_high_coverage_Illumina."+selected_chromosome+".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    zipinfos = zipdata.infolist()

#get the name of each zip file
names_files = [zipinfo.filename for zipinfo in zipinfos]

#iterate across files and get only final reports
print("\n#######################################\n#######################################")
print("Unzipping data: ")
print("#######################################\n#######################################")
#we are using a loop because it is a bit faster than zipdata.extractall. Parallelizing does not seems to be a good solution for unzipping and I got problems with pool. This takes a few minutes anyways.
    #parallel slower
        #https://stackoverflow.com/questions/43313666/python-parallel-processing-to-unzip-files
    #count time
        #https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
import numpy as np
#zipinfo=zipinfos[1]
#zipinfo=zipinfos[np.where([element == zip_name + "/SNP_Map.txt" for element in names_files])[0][0]]
#zipinfo=zipinfos[np.where([element == zip_name + "/Sample_Map.txt" for element in names_files])[0][0]]
import os
for zipinfo in zipinfos:

    #fir the file name starts with the adequate zip (batch) name and FinalReport
    if (zipinfo.filename.startswith(zip_name + "/" + zip_name + "_FinalReport")) | (zipinfo.filename.startswith(zip_name + "/SNP_Map.txt")) | (zipinfo.filename.startswith(zip_name + "/Sample_Map.txt")):

        #rename by removing the name of the parent folder, so we only get the file not the parent folder
        zipinfo.filename = zipinfo.filename.split(zip_name+"/")[1]
        print(zipinfo.filename)

        #extract the file in the temp dict
        zipdata.extract(zipinfo, temp_dir.name)
            #https://stackoverflow.com/questions/44079913/renaming-the-extracted-file-from-zipfile/56362289#56362289
            #https://stackoverflow.com/questions/13765486/how-to-unzip-specific-folder-from-a-zip-with-python

        #if the selected file is a FinalReport
        if zipinfo.filename.startswith(zip_name + "_FinalReport"):

            #remove the first 10 lines of the file, which is the head 
            os.system("cd " + temp_dir.name + "; tail -n +11 " + zipinfo.filename + " > tmp.txt && mv tmp.txt " + zipinfo.filename)
                #if you use tail with "+number_line" you can get all lines starting from the selected number of line
                    #tail is faster than sed
                #then save the result in a temporal file and overwrite the original file. The && will make sure that the file doesn't get overwritten when there is a problem.
                    #https://stackoverflow.com/questions/339483/how-can-i-remove-the-first-line-of-a-text-file-using-bash-sed-script
                    #https://www.baeldung.com/linux/remove-first-line-text-file



#parse large vcf file?
#https://www.biostars.org/p/101444/


#ASK DAVID
    #remove unrelated samples
        #we are going to use unrelated samples only. If we are calculating haplotype homozygosity by counting haplotypes, if a father and child have the same haplotype, this is not caused by selection but just by shared ancestry




### 
#finish hap and map files as soon as possible and send to elise
#while she calculates flex sweep prob, you calculate predictors in non-overlapping windwos for hg38, ask david about non-overlapping. It can be good to explore the whole genome, have more data, reduce correlation between windows, and we can just add the number of interest genes in each window, but then who we do the enrichment test? maybe taking the closest non-overlapping window to the center of each gene?

#this should be done in April, in may you could spend time doing modeling and get first results in hg38, ready for ambizione and conferences.

#in these two months, maybe you could have ready the cleaning of combat genes data, and you can show just a preliminary manhantan plot for the change in fitness? 