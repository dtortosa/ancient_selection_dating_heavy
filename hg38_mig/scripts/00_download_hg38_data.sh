#!/bin/bash 
	#to run this script: chmod +x script.sh; ./script.sh
	#!/bin/sh does not work with my terminal en msi of David.
	#if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
	#you can save the output and the errors
		#./script.sh > script.out #only output
		#./script.sh 2> error.out #only error
		#./script.sh > script.out 2> error.out #both in different files
		#./script.sh > script.out 2>&1 #both in the same file
		#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/



############################################################################
############################ DOWNLOAD HG38 DATA ############################
############################################################################

#Bash script dedicated to download hg38 data of 1KGP high coverage

#We are going to migrate to hg38 as the new release of 1KGP matches this genome reference version.

#The general repo page for this release (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage) and general readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/README_111822.pdf). 

#It says that "20220422_3202_phased_SNV_INDEL_SV" is "the most up-to-date version of the phased panel based on the high-coverage 1kGP WGS data which includes SNV, INDEL, and SV calls across 3,202 samples"

#High coverage data (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38). There you can get the phased data, including single nucleotide variants, indels, SV in vcf files per chromosome (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/). See readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf). This is the repo used by Jesus Murga. 



##################
#### Starting ####
##################


######
# wd #
######

#save the working directory
path_working_dir="/home/dftortosa/singularity/dating_climate_adaptation/hg38_mig"

#set the working directory
cd $path_working_dir



#######################
#### Download data ####
#######################

#define the URL where download the data
path_hg38_data="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"

#make a loop to download the data of each chromosome
#i=22
for i in $(seq 1 22); do

	#show the chromosome
	echo "#############################"
	echo "Downloading chromosome ${i}"
	echo "#############################"
	echo "see the path:"
	echo "${path_hg38_data}1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

	#download the vcf file
	wget "${path_hg38_data}1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" --directory-prefix=./data/vcf_files_hg38/
		#https://linuxize.com/post/bash-concatenate-strings/
done
	#there are other options for looping over a sequence that are better for memory
		#https://stackoverflow.com/questions/169511/how-do-i-iterate-over-a-range-of-numbers-defined-by-variables-in-bash


#CHECK WE HAVE THE CORRECT FILES
	#IF THE FILE ALREAYD EXIST YOU GET THE SAME NAME WITH .1



#ASK JESUS ABOUT THE MASKS