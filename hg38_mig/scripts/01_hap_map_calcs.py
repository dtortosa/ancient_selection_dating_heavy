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
import os
os.system("mkdir -p ./results/clean_vcf")
os.system("mkdir -p ./results/hap_map_files")
    #-p: no error if the folder already exists, make parent directories as needed



#############
# pops prep #
#############

#load pedigree of the latest version of the phased data that has sample IDs, sex and parents but no pop names
import pandas as pd
samples_pedigree = pd.read_csv("data/pedigrees/1kGP.3202_samples.pedigree_info.txt", 
    sep=" ", 
    header=0, 
    low_memory=False)

#load also a pedigree present in the main directory of the high coverage data. This has sample and pop IDs, but parents and sex are different with respect to the pedigree of the new sample
samples_pedigree_pop = pd.read_csv("data/pedigrees/20130606_g1k_3202_samples_ped_population.txt", 
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

#look
print("\n#######################################\n#######################################")
print("final pedigree data")
print("#######################################\n#######################################")
print(ped_merged)

#extract the distinct population names
pop_names = ped_merged["population"].unique()

#check
print("\n#######################################\n#######################################")
print("see distinct population names and check they are 26")
print("#######################################\n#######################################")
print(pop_names)
print(len(pop_names) == 26)



#################
# bcftools prep #
#################

#see version of bcftools
print("\n#######################################\n#######################################")
print("see bcftools version")
print("#######################################\n#######################################")
print(os.system("bcftools -v"))



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
        os.system("bcftools head " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | grep -i '^##fileformat'")
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
        os.system(
            "bcftools query \
                -l \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                head -10")
            #https://samtools.github.io/bcftools/howtos/query.html

        #check
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": do we have the correct number of total samples?")
        print("#######################################\n#######################################")
        os.system(
            "n_samples=$( \
                bcftools query \
                    -l \
                    " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                    wc -l); \
            if [ $n_samples -eq " + str(ped_merged.shape[0]) + " ]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi")
            #list the samples and count them
            #if the number of samples is equal to the number of samples in the pedigree, perfect

        #inspect
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": show the variant type, ID, chromosome, position, alleles and frequency for the first snps")
        print("#######################################\n#######################################")
        os.system(
            "bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF\n' \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
                head -5")
            #select the format of the query indicating the columns you want to show per SNP.
                #you can include data from INFO
                #end with \n to have different lines per SNPs

        #now show only snps
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": show now only SNPs")
        print("#######################################\n#######################################")
        os.system(
            "bcftools filter \
                -i 'TYPE=\"snp\"' \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF\n' | \
            head -5")
            #you can directly filter in query
            #https://davetang.github.io/learning_vcf_file/#filtering-variant-types

        #see snps of only three samples
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": see genotypes of these snps for a few samples: ")
        print("#######################################\n#######################################")
        os.system(" \
            bcftools filter \
                -i 'TYPE=\"snp\"' \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools query \
                -s HG00096,HG00097,HG00099 \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
            head -5")
            #select only snps
                #you can directly filter in query
            #make a query making a subset of individuals and then ask in each SNP for the type, chrom... and genotypes. As we are only selecting three samples, we will only see three genotypes
            #show only first 5 snps



    #######################
    # extract samples IDs #
    #######################

    #select the sample IDs for the selected population
    subset_pop = ped_merged.loc[ped_merged["population"] == selected_pop, :]

    #reset index
    subset_pop = subset_pop.reset_index(
        drop=True)
        #drop=True to avoid adding the index as a column

    #check
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we only have the selected pop")
    print("#######################################\n#######################################")
    print(all(subset_pop["population"] == selected_pop))

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

    #see genotypes
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first genotypes for selected samples")
    print("#######################################\n#######################################")
    os.system(" \
        bcftools query \
            -i 'TYPE=\"snp\"' \
            -s " + ",".join(selected_samples) + " \
            -f '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        head -5")
            #make bcftools query
                #include snps
                #subset the samples of the selected population
                #format the query asking for variant type, ID, chromosome, genotype of each variant
            #show only the first 5 snps

    #POR AQUIII

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see complete data for some of the selected samples")
    print("#######################################\n#######################################")
    os.system("bcftools filter -i 'TYPE=\"snp\"' " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | bcftools view -H -s " + ",".join(selected_samples) + " | head -10")
        #-H to avoid header,




    os.system("bcftools filter -i 'TYPE=\"snp\"' " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | bcftools view -H -s " + ",".join(selected_samples) + " -o data/eso.bcf.gz")
        #https://www.htslib.org/workflow/filter.html


    #THIS SHOULD GO IN THE GENERAL SECTION THAT IT IS ONLY RUN ONE TIME
    #check that the chromosome name in the vcf file is correct
    os.system(
        'chr_vcf=$(bcftools query \
            ' + input_vcfs_path + '/1kGP_high_coverage_Illumina.chr' + selected_chromosome + '.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
            --format "%CHROM\n" | \
            head -1); \
        if [ $chr_vcf = "chr' + selected_chromosome + '" ]; then \
            echo "TRUE"; \
        else \
            echo "FALSE"; \
        fi')
        #with bcftools, query for the chromosome of each variant but get only the first line
        #if the chromosome number is equal to selected_chromosome, perfect
        #equality for strings using 1 bracket is "="
            #https://stackoverflow.com/questions/18102454/why-am-i-getting-an-unexpected-operator-error-in-bash-string-equality-test

    #more CHECKS??

    #convert VCF file to hap file
    os.system(
        "bcftools convert " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        --include 'TYPE=\"snp\"' \
        --samples " + ",".join(selected_samples) + " \
        --hapsample ./results/hap_map_files/" + selected_pop)
        #convert 
            #those variants that are SNPs
            #for those samples of the selected pop
            #to hap format saving in results and using the name of the pop
        #to hap file (IMPUTE2 format)
            #https://www.htslib.org/doc/bcftools.html#convert

        

        #CHECK THIS IS THE CORRECT FORMAT FOR FLEXSEEP
            #look at the files and compare with your hap files for iHS

        #CHECK CONVERT FUNCTION
            #https://www.htslib.org/doc/bcftools.html#convert

        #https://www.biostars.org/p/270381/
    

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
                results/hap_map_files/' + selected_pop + '.hap.gz | \
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
                -c results/hap_map_files/' + selected_pop + '.hap.gz | \
            awk -F " "  "{print NF}" | \
            sort | \
            uniq); \
        nfields_hap_nlines=$( \
            echo -n "$nfields_hap" | \
            grep -c "^"); \
        if [ "${nfields_hap_nlines}" -eq 1 ]; then \
            nsamples_hap=$((($nfields_hap - 5)/2)); \
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
            #take nfields_hap, substract 5 because the first 5 columns of the hap file are not genotypes and then divide by 2. Remember that after the alleles names, we have 2 genotype columns per sample, so the number of columns/2 is the number of samples.
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


    #check chromosome number correspond with the name of the file

    #check for duplicates in snp ID and position
        #https://www.biostars.org/p/274824/

    #ask jesus about the filters he applied besides the ones applied by 1000 genomes project according to readme file.

    #ON FILTERING SNPS WITH BCFTOOLS
        #https://www.htslib.org/workflow/filter.html

    #BCFTOOLS CHEATSHEET
        #https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b#file-bcftools-cheat-sheet


    #filter by accesibility mask?


    for index, variant in enumerate(VCF("data/vcf_files_hg38/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")): # or VCF('some.bcf')
        if(index == 0):
            print(variant) # e.g. REF='A', ALT=['C', 'T']


#run in paralle across combinations of chromosomes and populations

##YOU CAN USEE SPARK IN A SINGLE DF!!! THE FILE WILL NOT BE FULLY LOADED IN MEMORY

##clean vcf, then impute to create the hap file as in the slim simulations

#they say this package is faster than the original becuase uses cython (like C but in python)
#https://stackoverflow.com/questions/72574865/how-to-read-a-vcf-gz-file-in-python
#https://brentp.github.io/cyvcf2/
#https://github.com/brentp/cyvcf2
from cyvcf2 import VCF

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





#