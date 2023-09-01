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



#####################################################################
#################### DOWNLOAD FASTA ANCESTRAl #######################
#####################################################################

#This bash script will download the fasta file with the estimations of ancestral alleles. This fasta file is used by the AncestraAllele plugin in order to polarize SNPs that are provided in a VCF file.
	#https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele



##################
#### Starting ####
##################


######
# wd #
######

#go to the general folder
cd /home/dftortosa/singularity/dating_climate_adaptation/hg38_mig

#save the working directory
path_working_dir="./data/fasta_ancestral/fasta_ancestral"

#set the working directory
mkdir -p $path_working_dir
cd $path_working_dir



#######################
#### Download data ####
#######################
wget https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz
