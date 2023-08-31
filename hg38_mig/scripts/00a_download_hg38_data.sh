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

#Bash script dedicated to download hg38 data of 1KGP high coverage (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

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

#save the working directory
path_working_dir="./data/vcf_files_hg38/"

#set the working directory
cd $path_working_dir



#######################
#### Download data ####
#######################



############################################
# loop across chromosomes to get VCF files #
############################################

#define the URL where download the data
path_hg38_data="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"

#set the total number of chromosomes you want
n_chrom=22

#make a loop to download the data of each chromosome
#i=22
for i in $(seq 1 $n_chrom); do

	#show the chromosome
	echo "#############################"
	echo "Downloading chromosome ${i}"
	echo "#############################"
	
	#set the file name
	vcf_file="1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

	#see the path
	echo "see the path:"
	echo "${path_hg38_data}${vcf_file}"

	#if the file we are trying to download exists
		#https://linuxize.com/post/bash-check-if-file-exists/
	if [[ -f "$vcf_file" ]]; then

		#remove it
		rm "./${vcf_file}"
	fi

	#download the vcf file
	wget -nv "${path_hg38_data}${vcf_file}" --directory-prefix=./
		#-nv: turn off verboseness, without being quiet, so you get only a bit of output
		#https://linuxize.com/post/bash-concatenate-strings/
done
	#there are other options for looping over a sequence that are better for memory
		#https://stackoverflow.com/questions/169511/how-do-i-iterate-over-a-range-of-numbers-defined-by-variables-in-bash



###############################
# check we have all the files #
###############################

#list files with the name for vcf across all chromosomes and count them
n_files=$(ls "1kGP_high_coverage_Illumina.chr"*".filtered.SNV_INDEL_SV_phased_panel.vcf.gz" | wc -l)

#check
echo "#############################"
echo "DO WE HAVE VCF FILES FOR ALL CHROMOSOMES?"
echo "#############################"
if [[ $n_files -eq $n_chrom ]]; then
	echo "TRUE, WE HAVE ${n_files} FILES"
else
	echo "FALSE, WE HAVE ${n_files} FILES"
fi



########################
# get strict mask file #
########################

#download the strict accessibility mask to select variant in regions that are accessible to short-reads sequencing (see 01_hap_map_calcs.py for further details)
wget -nv http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed --directory-prefix=../masks/
	#-nv: turn off verboseness, without being quiet, so you get only a bit of output
	#https://linuxize.com/post/bash-concatenate-strings/

#compress the file
gzip --force --keep ../masks/20160622.allChr.mask.bed
	#--force: force overwrite of output file and compress links
	#--keep: keep (don't delete) input files
