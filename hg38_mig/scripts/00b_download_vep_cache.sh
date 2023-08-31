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
############################ DOWNLOAD VEP CACHE ############################
############################################################################

#This bash script will download the cache required for VEP.
	#https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

#VEP can use a variety of annotation sources to retrieve the transcript models used to predict consequence types.
	#Cache - a downloadable file containing all transcript models, regulatory features and variant data for a species
	#GFF or GTF - use transcript models defined in a tabix-indexed GFF or GTF file
	#Database - connect to a MySQL database server hosting Ensembl databases

#Using a cache is the most efficient way to use VEP; we would encourage you to use a cache wherever possible. Caches are easy to download and set up using the installer. 
	#Using a cache (--cache) is the fastest and most efficient way to use VEP, as in most cases only a single initial network connection is made and most data is read from local disk. Use offline mode to eliminate all network connections for speed and/or privacy.

#we are going to download the cache independently, i.e., not with INSTALL.PL because we want run INSTALL.PL within the singularity container, so If we ask it for download the cache, it will be downloaded every time we build the container!! the whole 26GB! 

#we need to downloaded one time, selecting the current version of the cache, i.e., if the current version of VEP is 110, you need cache 110. Then, you can upload it to the HPC and have it there ready to use it with vep. If you use a folder different from $HOME/.vep, then you have to use --dir-cache when running vep to indicate where the cache is located.

#CAUTION:
	#We strongly recommend that you download/use the VEP Cache version which corresponds to your Ensembl VEP installation, i.e. the VEP Cache version 110 should be used with the Ensembl VEP tool version 110.
	#This is mainly due to the fact that the VEP Cache (data content and structure) is generated every Ensembl release, regarding the data and API updates for this release, therefore the cache data format might differ between versions (and be incompatible with a newer version of the Ensembl VEP tool).

#Automatic download of caches
	#Ensembl creates cache files for every species for each Ensembl release. They can be automatically downloaded and configured using INSTALL.pl.
	#If interested in RefSeq transcripts you may download an alternate cache file (e.g. homo_sapiens_refseq), or a merged file of RefSeq and Ensembl transcripts (eg homo_sapiens_merged); remember to specify --refseq or --merged when running VEP to use the relevant cache. See documentation for full details.

#manual download
	#It is also simple to download and set up caches without using the installer. By default, VEP searches for caches in $HOME/.vep; to use a different directory when running VEP, use --dir_cache.
	#Indexed cache
		#Essential for human and other species with large sets of variant data
		#(https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/) - requires Bio::DB::HTS (setup by INSTALL.pl) or tabix.



##################
#### Starting ####
##################


######
# wd #
######

#go to the general folder
cd /home/dftortosa/singularity/dating_climate_adaptation/hg38_mig

#set the cache version
cache_version="110"

#save the working directory
path_working_dir="./data/vep_cache/cache_${cache_version}"

#set the working directory
mkdir -p $path_working_dir
cd $path_working_dir



#######################
#### Download data ####
#######################
curl -O https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache/homo_sapiens_vep_${cache_version}_GRCh38.tar.gz
	#-O: Use HTTP 1.0



####################
#### decompress ####
####################
tar xzf homo_sapiens_vep_${cache_version}_GRCh38.tar.gz
	#x: extract files from an archive
	#z: extract gzip
	#f: we decompress a file indicating its name
