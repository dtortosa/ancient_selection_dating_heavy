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
ped_merged.shape[0] == 3202

ped_merged["sampleID"].equals(ped_merged["SampleID"])

ped_merged = ped_merged.drop("SampleID", axis=1)

ped_merged = ped_merged.rename(columns={"Population": "population", "Superpopulation": "superpopulation"})

ped_merged


ped_merged.loc[ped_merged["population"] == "IBS", :]


#see version of bcftools
print(os.system("bcftools -v"))



################################################################
#### function to clean vcf files and create hap - map files ####
################################################################

#selected_chromosome="22"
def master_processor(selected_chromosome):


    ped_merged.loc[ped_merged["population"] == "IBS", :]

    #extract info from columns
    #https://davetang.github.io/learning_vcf_file/#extracting-info-from-columns

    #subset by sample
    #https://davetang.github.io/learning_vcf_file/#subset-samples-from-a-multi-sample-vcf-file

    #filter by variant type, i.e., remove indels, structurla variations...
    #https://davetang.github.io/learning_vcf_file/#filtering-variant-types



    #get list of samples per pop  to make the subset

    "vcf-subset -e -c list_of_Punjabi_persons.txt 1000_genomes_chr3.vcf > PJL_chr3.vcf"
        #https://www.biostars.org/p/241810/
        #https://manpages.ubuntu.com/manpages/bionic/man1/vcf-subset.1.html
        #if vcftools is too slow, you can use bcftools
            #https://www.biostars.org/p/184950/



    for index, variant in enumerate(VCF("data/vcf_files_hg38/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz")): # or VCF('some.bcf')
        if(index == 0):
            print(variant) # e.g. REF='A', ALT=['C', 'T']



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





