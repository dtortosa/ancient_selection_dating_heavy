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




######################################################
#### function to calculate map files for ALL SNPs ####
######################################################
print_text("function to calculate map files for ALL SNPs", header=1)
#selected_chromosome="1"; debugging=True
def master_processor(selected_chromosome, debugging=False):

    #redirect standard output ONLY when running for production
    if debugging == False:
        import sys
        original_stdout = sys.stdout
            #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
        sys.stdout = open("./scripts/00b_map_calcs_outputs/chr" + selected_chromosome + ".out", "w")
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/



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
            head -n 5000 > ./results/00b_map_files/chr" + selected_chromosome + "/00b_map_calc_debug_subset_chr" + selected_chromosome + ".vcf; \
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




    print_text("remove files not required anymore", header=3)
    run_bash(" \
        rm " + input_vcf_file_map_calc + "; \
        ls -l ./results/00b_map_files/chr" + selected_chromosome)


    print_text("FINISH", header=3)


    print_text("restore sys.stdout using the previously saved reference to it. This is useful if you intend to use stdout for other things only required if we are in production, as we changed stdout only in that case", header=3)
    if debugging==False:
        sys.stdout = original_stdout
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/



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
import multiprocessing as mp
pool = mp.Pool(len(chromosomes))


print_text("run function across pandas rows", header=3)
pool.map(master_processor, chromosomes)

print_text("close the pool", header=3)
pool.close()




########################################################
#### Do some checks after analyzing all chromosomes ####
########################################################
print_text("Do some checks after analyzing all chromosomes", header=1)
print_text("Number of SNPs with low- and high-confidence ancestral alleles across chromosomes", header=2)
print_text("empty lists to save the counts of high and low confidence ancestral alleles across chromosomes", header=3)
count_snps_acgt_list = []
count_snps_acgt_upper_anc_list = []
count_snps_acgt_lower_anc_list = []

print_text("run loop across chromosomes", header=3)
#chrom=1
for chrom in chromosomes:
    print("Doing chromosome " + str(chrom))

    print("from the output file of the selected chromosome, extract the row with the counts of high and low confidence ancestral alleles")
    row_results = run_bash(" \
        grep \
            'count_snps_acgt' \
            ./scripts/00b_map_calcs_outputs/chr" + str(chrom) + ".out", return_value=True).strip()

    print("see the row")
    print(row_results)

    print("split the row")
    row_results_split = row_results.split(",")
    print(row_results_split)

    print("check that we only have one row, and it can be split in three parts with comma")
    if ("\n" not in row_results) & (len(row_results_split)==3):
        print("YES! GOOD TO GO!")
    else: 
        raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs with low and high confidence ancestral alleles")

    print("append to each list the corresponding count")
    count_snps_acgt_list.append(row_results_split[0].replace("count_snps_acgt=", ""))
    count_snps_acgt_upper_anc_list.append(row_results_split[1].replace("count_snps_acgt_upper_anc=", ""))
    count_snps_acgt_lower_anc_list.append(row_results_split[2].replace("count_snps_acgt_lower_anc=", ""))


print_text("process the results", header=3)
print_text("convert to int each count if the count is NOt zero", header=4)
count_snps_acgt_list = [int(x) if x!="" else 0 for x in count_snps_acgt_list]
count_snps_acgt_upper_anc_list = [int(x) if x!="" else 0 for x in count_snps_acgt_upper_anc_list]
count_snps_acgt_lower_anc_list = [int(x) if x!="" else 0 for x in count_snps_acgt_lower_anc_list]

print_text("check we have all chromosomes", header=4)
#selected_list=[count_snps_acgt_list, count_snps_acgt_upper_anc_list, count_snps_acgt_lower_anc_list][0]
check_count = [len(selected_list) == 22 for selected_list in [count_snps_acgt_list, count_snps_acgt_upper_anc_list, count_snps_acgt_lower_anc_list]]
if sum(check_count) == len(check_count):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs with low and high confidence ancestral alleles")

print_text("make the sum across chromosomes", header=4)
count_snps_acgt_total = sum(count_snps_acgt_list)
count_snps_acgt_upper_anc_total = sum(count_snps_acgt_upper_anc_list)
count_snps_acgt_lower_anc_total = sum(count_snps_acgt_lower_anc_list)

print_text("check the total is equal to the sum of high and low-confidence", header=4)
if count_snps_acgt_total == count_snps_acgt_upper_anc_total+count_snps_acgt_lower_anc_total:
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs with low and high confidence ancestral alleles")

print_text("see the sums", header=4)
print("Total number of SNPs with ancestral allele: " + str(count_snps_acgt_total))
print("Total number of SNPs with high-confidence ancestral allele: " + str(count_snps_acgt_upper_anc_total))
print("Total number of SNPs with low-confidence ancestral allele: " + str(count_snps_acgt_lower_anc_total))

print_text("calculate the percentage of SNPs with low-confidence ancestral allele across all chromosomes", header=4)
print(count_snps_acgt_lower_anc_total/count_snps_acgt_total*100)

print_text("IMPORTANT: if you need the number of SNPs without ancestral allele data, you can find it per chromosome in the line starting with 'The number of SNPs without ancestral allele data is'", header=4)



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
#for 01d_hap_map_calcs.py
    #You will do the polarization within each populatin. We need only biallelic snps (to easily exchange ref by alt if needed), and we need to do this within pop, porque an allele can be biallelic for one pop but not for other.
    #list of SNPs to for which we have to switch REF/ALT
        #WE HAVE TO EXCLUDE CASES WHERE REF NOR ALT ARE THE ANCESTRAL ALLELE
            #these can be multiallelic SNPs for which one of the ALTs is not present in the selected population and that very ALT is the ancestral. We need these SNPs OUT, if no ancestral, we cannot estimate selection.
            #you hav eto do it after removing truly multiallelci snps in the population
            #COUNT THESE CASES TO SEE THE IMPACT
            #see line 812 of 01b_vep_ancestral.py
        #then create the list of SNPs for which REF is not AA
        #alternatively, you can directly export chrom, pos, REF, ALT, for all SNPs, then select rows with different REF than Ancestral using awk and generate a BED file, which is also accpeted by annotate
            #http://www.htslib.org/doc/bcftools.html#annotate
    #make the switch
        #then you can go to -fixref and use it to switch these snps in the original VCF file
            #bcftools +fixref file.bcf -Ob -o out.bcf -- -i List_of_1.vcf.gz 
                #https://www.biostars.org/p/411202/
                #https://samtools.github.io/bcftools/howtos/plugin.fixref.html
            #fixref ES PELIGROSO!!
        #I think you could also use --derived to set the ancestral as reference in bcftools. If you need to change the name of AA_upcase, I think you can use something like "bcftools annotate -c INFO/1kg_v2a_AF:=INFO/AF"
            #https://www.biostars.org/p/304979/
            #https://github.com/samtools/bcftools/issues/1695
    #check and then go to the real data
