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



##########################################################
######## CALCULATE ANCESTRAL ALLELES FOR EACH SNP ########
##########################################################

#We are going to migrate to hg38 as the new release of 1KGP matches this genome reference version (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38). 

#This script will estimate which allele is the ancestral for each SNP included in the new panel. REF/ALT in 1000 genomes project is not referring to Ancestral/Derived, so we need to polarize the alleles according to Jesus. I have checked the README and the paper and there is no information about ancestral, so it seems to be the case that they did not polarize the alleles. According to Jesus, I have to use a combination of VEP and bcftools to update the VCF files.


#The general repo page for this release (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage) and general readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/README_111822.pdf).

#It says that "20220422_3202_phased_SNV_INDEL_SV" is "the most up-to-date version of the phased panel based on the high-coverage 1kGP WGS data which includes SNV, INDEL, and SV calls across 3,202 samples"

#High coverage and phased data can be found in the general repo, working folder and then 20220422_3202_phased_SNV_INDEL_SV (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/). This phased data includes single nucleotide variants, indels, SV in vcf files per chromosome. See its specific readme (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf). This is the repo used by Jesus Murga. 

#For information about the download of the pedigrees, see 01a_selecting_pedegree.py.



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
print_text("Preparate pedigree data", header=1)

#save name of path including input vcf files
input_vcfs_path = "data/vcf_files_hg38"

#create folders to save the results
run_bash(" \
    mkdir \
        -p ./data/dummy_vcf_files/00_dummy_vcf_files_vep \
    mkdir \
        -p ./results/00_vep_vcf_files")
    #-p: no error if the folder already exists, make parent directories as needed

#add readme file
run_bash(" \
    cd ./results/00_vep_vcf_files; \
    echo 'These VCF files include an INFO field with the Ancestral Allele of each variant, BUT this has not been updated in the REF/ALT fields, meaning that the REF alleles in these files can be derived. We will polarize variants in the next script where VCF cleaning will be done within each population.' > README.txt")



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



####################################################
# check bcftools's behaviour with a dummy vcf file #
####################################################
print_text("Check bcftools's behaviour with a dummy vcf file", header=1)
print_text("We are going to use dummy vcf file to check the behaviour of bcftools", header=2)
print_text("see the dummy variants", header=3)
run_bash(" \
    bcftools view \
        data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example.vcf | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %INFO/AF GTs:[ %GT]\n'")
        #We have variants with different characteristics. In some cases, the AN, AC fields are not correct, but we use this to check whether the different commands we use can correctly update these fields.
            #Biallelic SNPs:
                #SNP rs6054247 chr20 14280 T A 5 0 0.2 GTs: 0|0 0|0 1|1
                #SNP rs6054248 chr20 14290 C A 5 0 0.2 GTs: 0|0 0|0 1|.
                #SNP rs6054249 chr20 14300 C A 5 1 0.2 GTs: 0|0 0|0 1|.
                #SNP rs6054250 chr20 14310 C A 3 0 0 GTs: 1|0 .|. .|.
                #SNP rs6054251 chr20 14320 C A 6 2 0.333 GTs: 0|0 0|0 1|1
                #SNP rs6054252 chr20 14350 C A 4 2 0.5 GTs: 1|0 1|0 .|.
                #SNP rs6054257 chr20 14370 G A 6 3 0.5 GTs: 1|0 1|1 0|0
                #SNP rs6054255 chr20 14371 G C 6 3 0.667 GTs: 0|1 1|1 1|0
                #SNP rs6040351 chr20 17330 T A 6 2 0.333 GTs: 0|0 1|0 1|0
            #Multiallelic SNPs:
                #SNP rs6040355 chr20 1110696 A G,T 6 2,2 0.333,0.333 GTs: 1|1 2|2 0|0
                #SNP rs6040356 chr20 1110697 A C,T 6 2,1 0.333,0 GTs: 1|1 0|0 0|0
                #SNP rs6040357 chr20 1110698 A G,T 6 3,3 0.5,0.5 GTs: 1|1 2|2 1|2
                #SNP rs6040358 chr20 1110699 G A,T 6 2,0 0.333,0 GTs: 1|1 0|0 .|.
                #SNP rs6040359 chr20 1110700 A G,T 4 0,3 0,0.6 GTs: 2|2 2|2 2|.
            #Exact duplicate SNPs, i.e., pos, chr, REF, alt
                #SNP rs6040360 chr20 1110701 A G 6 3 0.5 GTs: 1|1 0|0 0|1
                #SNP rs6040360_copy chr20 1110701 A G 6 3 0.5 GTs: 1|1 1|0 0|1
            #SNP with unphased data
                #SNP rs6040361 chr20 1110702 A G 6 3 0.5 GTs: 1|1 0/0 0|1
            #microsatellite (indel)
                #INDEL microsat1 chr20 1110703 GTC G,GTCT . . . GTs: 0/1 0/2 1/1


print_text("compress the dummy file", header=3)
run_bash("\
    cd ./data/dummy_vcf_files/00_dummy_vcf_files_vep/; \
    bgzip \
        --keep \
        --force \
        --threads 1 \
        dummy_example.vcf; \
    ls -l")



#
print_text("Instructions for installing and using VEP in order to polarize alleles", header=2)
print("REF/ALT in 1000 genomes project is not referring to Ancestral/Derived, so we need to polarize the alleles according to Jesus. I have checked the README and the paper and there is no information about ancestral, so it seems to be the case that they did not polarize the alleles. According to Jesus, I have to use a combination of VEP and bcftools to update the VCF files.")
#VEP (Variant Effect Predictor) and AncestralAllele plugin (https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele):
    #VEP plugin that retrieves ancestral allele sequences from a FASTA file.
    #Ensembl produces FASTA file dumps of the ancestral sequences of key species, including humans
        #Data files for GRCh37 are available from https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/
        #Data files for GRCh38 are available from https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/
    #In the VEP paper
        #they say "In Ensembl, we infer genome-wide ancestral sequences (Paten et al., 2008) for different groups of species. Select the “Ancestral Allele” option (Figure 5) to obtain the ancestral ALLELE PREDICTED FROM THE ALIGNMENT OF 12 PRIMATE SPECIES, INCLUDING HOMO SAPIENS."
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7613081/
    #Indeed, in the supplementary file of the paper for 1KGP phase 3 (previous hg19 version; section 8), they say they polarized alleles using ensembl data from "ftp://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2", which is an old version of the data indicated in VEP webpage.
        #https://static-content.springer.com/esm/art%3A10.1038%2Fnature15393/MediaObjects/41586_2015_BFnature15393_MOESM86_ESM.pdf
    #In both, the supplementary and the README of hg19 (https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.README), they say these ancestral alleles are obtained using the EPO pipeline:
        #In the EPO (Enredo-Pecan-Ortheus) pipeline, Ortheus infers ancestral states from the Pecan alignments. The confidence in the ancestral call is determined by comparing the call to the ancestor of the ancestral sequence as well as the 'sister' sequence of the query species. For instance, using a human-chimp- macaque alignment to get the ancestral state of human, the human-chimp ancestor sequence is compared to the chimp and to the human-chimp-macaque ancestor. A high-confidence call is made when all three sequences agree. If the ancestral sequence agrees with one of the other two sequences only, we tag the call as a low-confidence call. If there is more disagreement, the call is not made.
        #The convention for the sequence is:
            #ACTG: high-confidence call, ancestral state supported by the other two sequences
            #actg: low-confidence call, ancestral state supported by one sequence only
            #N: failure, the ancestral state is not supported by any other sequence
            #-: the extant species contains an insertion at this postion
            #.: no coverage in the alignment

#install VEP
    #the next steps show the commands to install VEP in David's laptop. If you need to install it in a singularity container, then go an check the .def file "01_ubuntu_20_04_hg38_mig_hap_map.def"
        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
    #install HTSLib (from SAMtools)
        #download the last version in March 2023
            #wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
        #decompress
            #tar -xf htslib-1.17.tar.bz2
                #https://linuxize.com/post/how-to-extract-unzip-tar-bz2-file/
        #go to the new folder
            #cd htslib-1.17
        #indicate where to install the executable program
            #./configure --prefix=/usr/local
                #if not /usr/bin, add the path to bashrc in laptop or %environment in the singularity
                #export PATH=/opt/htslib-1.17:$PATH
        #compile
            #make
            #sudo make install
    #install required perl modules
        #tool for easily installing modules
            #sudo cpan App::cpanminus
                #you have to answer some questions here, so maybe this is problematic for container
                #http://www.cpan.org/modules/INSTALL.html
        #create folder to save modules
            #mkdir -p $HOME/cpanm
            #echo -e "export PERL5LIB=$PERL5LIB:$HOME/cpanm/lib/perl5\n" >> $HOME/.bashrc
                #indicate path where to install the perl modules
                #https://superuser.com/questions/154936/echo-text-with-new-line-in-bash
        #modules required for VEP
            #cpanm --local-lib $HOME/cpanm DBI
            #cpanm --local-lib $HOME/cpanm DBD::mysql
            #cpanm --local-lib $HOME/cpanm Archive::Zip
        #modules maybe useful for AncestralAllele plugin
            #sudo cpanm --local-lib $HOME/cpanm Bio::DB::HTS
            #sudo cpanm --local-lib $HOME/cpanm Bio::DB::Fasta
            #sudo cpanm --local-lib $HOME/cpanm Bio::DB::HTS::Faidx
            #sudo cpanm --local-lib $HOME/cpanm PerlIO::gzip
                #indicate the local lib where to install
    #install VEP (https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer):
        #git clone https://github.com/Ensembl/ensembl-vep.git
        #cd ensembl-vep
        #perl INSTALL.pl --AUTO ap --PLUGINS AncestralAllele --ASSEMBLY GRCh38 --SPECIES "homo_sapiens"
            #https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html
            #--ASSEMBLY GRCh38
                #If using the --AUTO functionality to install without prompts, remember to add the assembly version required using e.g. "--ASSEMBLY GRCh38"
                        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#assembly
            #--AUTO ap
                #Run installer without prompts
                #install also
                    #a: API + perl module Bio::DB::HTS/htslib
                    #c: download cache, not needed here because I have downloaded the last cache version in the a separated bash script in my laptop. This cache has been then sent to the HPC, just one time. YOU HAVE TO USE THE SAME VERSION FOR VEP AND CACHE
                        #Using a cache (--cache) is the fastest and most efficient way to use VEP, as in most cases only a single initial network connection is made and most data is read from local disk. Use offline mode to eliminate all network connections for speed and/or privacy.
                        #By default the installer will download the latest version of VEP caches and FASTA files.
                        #We strongly recommend that you download/use the VEP Cache version which corresponds to your Ensembl VEP installation, i.e. the VEP Cache version 109 should be used with the Ensembl VEP tool version 109.
                        #SO INSTALL CACHE EACH TIME YOU INSTALL VEP
                            #https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
                    #f: download fasta files
                        #I do not think I have to download all fasta files as we are only going to use ancestral alleles, which are in a specific fasta file.
                        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta
                    #p: install plugins indicated with --PLUGINS
            #--PLUGINS AncestralAllele
                #Comma-separated list of plugins to install when using --AUTO. To install all available plugins, use --PLUGINS all.
            #--PLUGINSDIR
                #By default the script will install the plugins files in the "Plugins" subdirectory of the --CACHEDIR directory. This option configures where the plugins files are installed.
                #The --dir_plugins flag must be passed when running the VEP if a non-default plugins directory is given.
            #--SPECIES "homo_sapiens"
                #Comma-separated list of species to install when using --AUTO. To install the RefSeq cache, add "_refseq" to the species name, e.g. "homo_sapiens_refseq", or "_merged" to install the merged Ensembl/RefSeq cache. Remember to use --refseq or --merged when running the VEP with the relevant cache!
                #Use all to install data for all available species.    #IMPORTANT information about the cache
        #correspondence between versions
            #We strongly recommend that you download/use the VEP Cache version which corresponds to your Ensembl VEP installation,i.e. the VEP Cache version 110 should be used with the Ensembl VEP tool version 110.
            #This is mainly due to the fact that the VEP Cache (data content and structure) is generated every Ensembl release, regarding the data and API updates for this release, therefore the cache data format might differ between versions (and be incompatible with a newer version of the Ensembl VEP tool).
            #Therefore, the cache should be installed in "/home/dftortosa/.vep/homo_sapiens/" and should be named as "110_GRCh38" at least by August 30th.
        #download cache
            #automatically
                #Ensembl creates cache files for every species for each Ensembl release. They can be automatically downloaded and configured using INSTALL.pl.
            #manually
                #It is also simple to download and set up caches without using the installer. By default, VEP searches for caches in $HOME/.vep; to use a different directory when running VEP, use --dir_cache.
                #Essential for human and other species with large sets of variant data. These lines download cache 110 (August 2023)! check the version you need to download if update VEP!!
                    #cd $HOME/.vep
                    #curl -O https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
                    #tar xzf homo_sapiens_vep_110_GRCh38.tar.gz
                    #CHANGE VERSION CACHE TO THE NEW VEP VERSION IF YOU UPDATE!
        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache



print_text("set cache version", header=2)
print("IMPORTANT: We are setting the cache version manually, using this time version 110, which is ok for VEP 110 (August 2023). If you update VEP, you will need to download the new cache version")
cache_version="110"
    #set the cache version

print_text("check we have the same VEP and cache versions", header=2)
vep_version_raw = run_bash(" \
    vep | \
    grep ensembl-vep", return_value=True)
    #run VEP and get the line with the VEP version
vep_version=vep_version_raw.strip().split("ensembl-vep          : ")[1].split(".")[0]
    #extract the version removing unnecessary text
if vep_version == cache_version:
    print("GOOD TO GO! Both VEP and the cache have the same version")
else:
    raise ValueError("SERIOUS ERROR! The version of VEP and the cache are not the same. You HAVE TO download the new version of the cache. Check script ('install VEP' section) and go to 'https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache'")
    #do the check



#
print_text("download the fasta files with the ancestral alleles", header=2)
#https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
#we have done it in a separate bash script in David's laptop ("00c_download_fasta_ancestral.sh")


print_text("list the files in the compressed file", header=3)
run_bash("\
    tar \
        --gunzip \
        --list \
        --file \
        ./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38.tar.gz")


print_text("see fasta file of the second chromosome", header=3)
#you can extract one specific file but typing its path and name within the compressed file
run_bash("\
    tar \
        -zxOvf \
        ./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38.tar.gz \
        homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_2.fa | \
    head -500")
    #tar
        #z for gzip files, x to extract files, O for sending to stdout, v for verbose and f for specifying the compressed file name
            #https://unix.stackexchange.com/questions/61461/how-to-extract-specific-files-from-tar-gz
            #https://www.gnu.org/software/tar/manual/html_node/Writing-to-Standard-Output.html
    #fasta file format (https://en.wikipedia.org/wiki/FASTA_format):
        #A sequence begins with a greater-than character (">") followed by a description of the sequence (all in a single line).
            #in our case is ">ANCESTOR_for_chromosome:GRCh38:2:1:242193529:1" 
        #The next lines immediately following the description line are the sequence representation, with one letter per amino acid or nucleic acid, and are typically no more than 80 characters in length.
            #In our case, shows support for considering a base as ancestral, based on the comparison with the human-chimp-macaque ancestor (see above).
                #ACTG: high-confidence call, ancestral state supported by the other two sequences
                #actg: low-confidence call, ancestral state supported by one sequence only
                #N: failure, the ancestral state is not supported by any other sequence
                #-: the extant species contains an insertion at this postion
                #.: no coverage in the alignment


print_text("create one single, bgzipped file with all fasta files", header=3)
#For optimal retrieval speed, you should pre-process the FASTA files into a single bgzipped file that can be accessed via the perl module Bio::DB::HTS::Faidx (installed by VEP's INSTALL.pl). Therefore, we need to use bgzip (from HTSlib) to compress the fast file. We are using the uncompressed file because id does not work if we compress with bgzip
    #https://www.htslib.org/doc/bgzip.html
print_text("create one single file with all fasta files", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    tar \
        --extract \
        --gunzip \
        --file ./homo_sapiens_ancestor_GRCh38.tar.gz")

print_text("check we have fasta files for all chromosomes", header=4)
listing_fastas_raw = run_bash("\
    cd ./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38; \
    ls *.fa", return_value=True)
listing_fastas = listing_fastas_raw.split("\n")[0:-1]
    #get a list of chromosomes with fasta
from natsort import natsorted
listing_fastas_ordered = natsorted(listing_fastas)
    #sort
expected_chr_raw = [str(i) for i in range(1, 23)]
[expected_chr_raw.append(non_auto) for non_auto in ["MT", "X", "Y", "scaffold"]]
expected_chr = ["homo_sapiens_ancestor_" + chromosome + ".fa" for chromosome in expected_chr_raw] 
    #get a list with the expected chromosomes
expected_chr == listing_fastas_ordered
    #check

print_text("combine all fasta files and remove the unzipped folder", header=4)
run_bash(" \
    cd ./data/fasta_ancestral/; \
    cat ./homo_sapiens_ancestor_GRCh38/*.fa \
    > ./homo_sapiens_ancestor_GRCh38_final.fa; \
    rm \
        --recursive \
        --force \
        ./homo_sapiens_ancestor_GRCh38/; \
    ls -l")
        #use tar to unzip a tar.gz file generating a folder with all the fasta files
        #print all fasta files and bgzipped them sending to standard output, save as a file
        #remove the folder with all the fasta files
        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele

print_text("see the number of lines of the combined fasta file", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    awk \
        'END{print NR}' \
        homo_sapiens_ancestor_GRCh38_final.fa")

print_text("see first 500 lines of the combined fasta file", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    awk \
        '{if(NR>=1 && NR <=500){print $0}}' \
        homo_sapiens_ancestor_GRCh38_final.fa")

print_text("see the number of the line of the second chromosome so we can then check the first lines for that chromosome in the combined file", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    awk \
        '{if(/>ANCESTOR_for_chromosome:GRCh38:2:/){print NR}}' \
        homo_sapiens_ancestor_GRCh38_final.fa")
        #https://stackoverflow.com/questions/17908555/printing-with-sed-or-awk-a-line-following-a-matching-pattern

print_text("print the 500 rows after the second chromosome starts in the combined fasta file", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    awk \
        '{if(NR>=12871036 && NR<=12871536){print $0}}' \
        homo_sapiens_ancestor_GRCh38_final.fa")

print_text("bgzip the fasta file with ancestral alleles because this should make things faster when running VEP (see documentation below)", hader=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    bgzip \
        --keep \
        --force \
        --threads 1 \
        homo_sapiens_ancestor_GRCh38_final.fa; \
    ls -l")
        #--force
            #overwrite files without asking
        #--keep
            #don't delete input files during operation
        #--threads
            #number of compression threads to use

print_text("check integrity of the bgzipped file", header=4)
fasta_ancestra_integrity_test = run_bash("\
    cd ./data/fasta_ancestral/; \
    bgzip \
        --test \
        homo_sapiens_ancestor_GRCh38_final.fa.gz", return_value=True)
        #--test
            #test integrity of compressed file
if fasta_ancestra_integrity_test == "":
    print("The integrity of the bgzipped fasta file with ancestral allele estimations is OK")
else:
    raise ValueError("SERIOUS ERROR! There is a problem with the integrity of the bgzipped fasta file with ancestral allele estimations")

print_text("remove previous indexes of the fasta file if they are present. These are created the first time VEP is run on the files, so just to be sure we are using the correct index, we remove those existing previously", header=4)
run_bash(" \
    cd ./data/fasta_ancestral/; \
    n_prev_fasta_index=$(ls -l | grep 'fai\|gzi' | wc -l); \
    if [[ $n_prev_fasta_index -gt 0 ]];then \
        echo 'We have fasta index previously generated, removing related files...'; \
        rm *fai; \
        rm *gzi; \
    fi; \
    ls -l")


print_text("run VEP's plugin ancestral allele on the cleaned VCF file", header=3)
print_text("make a run of VEP", header=4)
run_bash("\
    vep \
        --offline \
        --verbose \
        --species 'homo_sapiens' \
        --assembly GRCh38 \
        --input_file ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example.vcf.gz \
        --format vcf \
        --output_file ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --vcf \
        --compress_output gzip \
        --force_overwrite \
        --fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,STRAND,AA \
        --cache \
        --dir_cache ./data/vep_cache/cache_" + cache_version + " \
        --plugin AncestralAllele,./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38_final.fa.gz \
        --dir_plugins /opt/ensembl-vep/vep_plugins")
        #https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html
        #--offline
            #Enable offline mode. No database connections will be made, and a cache file or GFF/GTF file is required for annotation. Add --refseq to use the refseq cache (if installed). Not used by default
        #--verbose
            #Print out a bit more information while running. Not used by default
        #--species
            #Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can speed up initial database connection as the registry does not have to load all available database aliases on the server. Default = "homo_sapiens"
        #assembly
            #Select the assembly version to use if more than one available. If using the cache, you must have the appropriate assembly's cache file installed. If not specified and you have only 1 assembly version installed, this will be chosen by default. Default = use found assembly version.
        #--input_file
            #Input file name. If not specified, VEP will attempt to read from STDIN. Can use compressed file (gzipped).
        #--format
            #Input file format - one of "ensembl", "vcf", "hgvs", "id", "region", "spdi". By default, VEP auto-detects the input file format. Using this option you can specify the input file is Ensembl, VCF, IDs, HGVS, SPDI or region format. Can use compressed version (gzipped) of any file format listed above. Auto-detects format by default
        #--output-file
            #Output file name. Results can write to STDOUT by specifying 'STDOUT' as the output file name - this will force quiet mode. Default = "variant_effect_output.txt"
        #--vcf
            #Writes output in VCF format. Consequences are added in the INFO field of the VCF file, using the key "CSQ". Data fields are encoded separated by "|"; the order of fields is written in the VCF header. Output fields in the "CSQ" INFO field can be selected by using --fields.
            #IF THE INPUT FORMAT WAS VCF, THE FILE WILL REMAIN UNCHANGED SAVE FOR THE ADDITION OF THE CSQ FIELD (unless using any filtering).
            #Custom data added with --custom are added as separate fields, using the key specified for each data file.
            #Commas in fields are replaced with ampersands (&) to preserve VCF format.
        #--compress_output
            #Writes output compressed using either gzip or bgzip. Not used by default
        #--force_overwrite
            #By default, VEP will fail with an error if the output file already exists. You can force the overwrite of the existing file by using this flag.
        #--fields
            #Configure the output format using a comma separated list of fields. Can only be used with tab (--tab) or VCF format (--vcf) output. For the tab format output, the selected fields may be those present in the default output columns, or any of those that appear in the Extra column (including those added by plugins or custom annotations) if the appropriate output is available (e.g. use --show_ref_allele to access 'REF_ALLELE'). Output remains tab-delimited. For the VCF format output, the selected fields are those present within the "CSQ" INFO field.
                #https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html
        #cache
            #Enables use of the cache. Add --refseq or --merged to use the refseq or merged cache, (if installed).
        #--dir_cache
            #Specify the cache directory to use. Default = "$HOME/.vep/"
            #The --dir_cache flag must be passed when running the VEP if a non-default cache directory is given
        #--plugin
            #Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory (defaults to $HOME/.vep/) or where indicated in --dir-plugins. Multiple plugins can be used by supplying the --plugin flag multiple times. See plugin documentation. Not used by default
        #--dir_plugins
            #Specify the plugin directory to use. Default = "$HOME/.vep/"
            #The --dir_plugins flag must be passed when running the VEP if a non-default plugins directory is given
        #AncestralAllele plugin
            #https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
            #A VEP plugin that retrieves ancestral allele sequences from a FASTA file.
            #Ensembl produces FASTA file dumps of the ancestral sequences of key species.
                #Data files for GRCh37 are available from https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/
                #Data files for GRCh38 are available from https://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/
            #For optimal retrieval speed, you should pre-process the FASTA files into a single bgzipped file that can be accessed via Bio::DB::HTS::Faidx (installed by VEP's INSTALL.pl):
            #Data file is only available for GRCh38.
            #The plugin is also compatible with Bio::DB::Fasta and an uncompressed FASTA file.
            #Note the first time you run the plugin with a newly generated FASTA file it will spend some time indexing the file. DO NOT INTERRUPT THIS PROCESS, particularly if you do not have Bio::DB::HTS installed.
                #the indexing creates two files in the folder of the fasta file: homo_sapiens_ancestor_GRCh38_final.fa.gz.fai and homo_sapiens_ancestor_GRCh38_final.fa.gz.gzi
                #thanks to these files, the second time you run it, it is much faster.
            #Special cases:
                #"-" represents an insertion
                #"?" indicates the chromosome could not be looked up in the FASTA
#we get different transcripts for each SNP
    #VEP me da para cada SNP información de la strand, alelo ancestral, e impacto para diferentes transcritos de un mismo gen en los que "cae" dicho SNP. Así, algunas filas pone protein coding, nonsense... y tienen diferentes strands (1/-1). Para los casos que he mirado, todas las filas del mismo SNP tienen el mismo alelo Ancestral, así que entiendo que podría coger cualquier

print_text("see first lines of the generated VCF file", header=4)
run_bash(" \
    gunzip \
        --stdout \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz | \
    head -50")

print_text("do some checks about about whether we have used the correct data using the row added to the header in the VCF file by VEP", header=4)
line_checks_after = run_bash(" \
    gunzip \
        --stdout \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz | \
    grep '##VEP='", return_value=True)
    #extract the row of the header with that information
line_checks_after_split = line_checks_after.split(" ")
    #split the different fields
print("Do we have the correct VEP version in the VCF file?")
print('##VEP="v'+vep_version+'\"' in line_checks_after_split)
print("Do we have the correct cache and assembly version in the VCF file?")
import numpy as np
print(line_checks_after_split[np.where(['ensembl=' in field for field in line_checks_after_split])[0][0]].split(".")[0] == 'ensembl='+cache_version)
print(line_checks_after_split[np.where(['assembly=' in field for field in line_checks_after_split])[0][0]].split(".")[0] == 'assembly="GRCh38')
    #select the field including "ensemble" or "assembly", then split by "." to have only the general number version and then check that number is the correct one.

print_text("see what happens with multiallelic snps", header=4)
run_bash(" \
    bcftools view \
        --max-alleles 10  \
        --min-alleles 3 \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz | \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %CSQ\n'")
print("The ancestral allele is selected even if it is one of the multiple derived alleles. Therefore, no problem.")

print_text("According to the VEP manual, vep does not change the vcf file, only add the CSQ field. Therefore, we could process the original vcf files and then process with our cleaning. Let's check if the VCF is exactly the same if we avoid the CSQ field, which is the one added by AncestralAllele plugin of VEP", header=4)
run_bash(" \
    cd ./data/dummy_vcf_files/00_dummy_vcf_files_vep; \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO/AN %INFO/AC %INFO/AF %INFO/NS %INFO/DP %INFO/AA %INFO/DB %INFO/H2 GTs:[ %GT]\n' \
        ./dummy_example.vcf.gz > file_check_1.txt; \
    bcftools query \
        -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO/AN %INFO/AC %INFO/AF %INFO/NS %INFO/DP %INFO/AA %INFO/DB %INFO/H2 GTs:[ %GT]\n' \
        ./dummy_example_vep.vcf.gz > file_check_2.txt; \
    echo '##See head first check file:'; head file_check_1.txt; \
    echo '##See head second check file:'; head file_check_2.txt; \
    check_status=$(cmp --silent file_check_1.txt file_check_2.txt; echo $?); \
    echo '##Do the check'; \
    if [[ $check_status -eq 0 ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")


print_text("Create a new AA field using the ancestral allele stored in CSQ/AA", header=3)
print_text("see the header of the VCF file after VEP processing", header=4)
run_bash("\
    bcftools head \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz")

print_text("see all the tags in the CSQ field", header=4)
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --list")

print_text("make tags inside CSQ available", header=4)
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --annotation CSQ \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
        #--annotation
            #INFO annotation to parse, being CSQ the default
        #you now can directly ask for AA and other tags added by VEP in the CSQ field using --columns (see below)

print_text("check we can call AA from CSQ and use it to filter", header=4)
print("exclude those SNPs for which the REF allele IS NOT the ancestral")
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --annotation CSQ \
        --exclude REF==AA \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
print("exclude those SNPs for which the REF allele IS the ancestral")
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --annotation CSQ \
        --include 'REF=AA'\
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
print("As you can see, we can just filter by the ancestral allele using AA and it does not consider C, C, C.... but only C, see below")

print_text("extract AA tag from CSQ and then remove CSQ", header=4)
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ | \
    bcftools view")
    #bcftools +split-vep
        #--columns
            #Extract the fields listed either as indexes or names. The default type of the new annotation is String but can be also Integer/Int or Float/Real.
    #bcftools annotate
        #--remove
            #List of annotations to remove. Use "FILTER" to remove all filters or "FILTER/SomeFilter" to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and "FORMAT" to remove all FORMAT tags except GT. To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be abbreviated to "INF" and "FORMAT" to "FMT".

print_text("select SNPs for which the ancestral allele is ACGT or acgt, avoiding cases where AA='.' or '-'", header=4)
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep.vcf.gz \
        --annotation CSQ \
        --include 'AA=\"A,C,G,T,a,c,g,t\"'\
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
        #we get rid of the SNP with AA="."
        #Comma in strings is interpreted as a separator and when multiple values are compared, the OR logic is used. Consequently, the following two expressions are equivalent but not the third:
            #-i 'TAG="hello,world"'
            #-i 'TAG="hello" || TAG="world"'
            #-i 'TAG="hello" && TAG="world"'
                #https://samtools.github.io/bcftools/bcftools.html#expressions
        #This is exactly what we want, check that our TAG ("AA") has ACGT both in upper (ACGT) and lower (acgt) case. See next line to understand why we need lower case.
        #Remember that the convention for the sequence in ancestral allele determination is:
            #ACTG: high-confidence call, ancestral state supported by the other two sequences
            #actg: low-confidence call, ancestral state supported by one sequence only
            #N: failure, the ancestral state is not supported by any other sequence
            #-: the extant species contains an insertion at this postion
            #.: no coverage in the alignment
                #supplementary file of the paper for 1KGP phase 3 (section 8)
                    #https://static-content.springer.com/esm/art%3A10.1038%2Fnature15393/MediaObjects/41586_2015_BFnature15393_MOESM86_ESM.pdf
                #README ancestral fasta file for hg19
                    #https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.README
        #In addition, the plugin AncestralAllele has two special conditions:
            #"-" represents an insertion
            #"?" indicates the chromosome could not be looked up in the FASTA
                #https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
        #Therefore, we only want to consider cases for which the ancestral allele is inferred.

print_text("now manually change ancestral of the first SNP from '.' to 'c,c,c' to check whether our expression catch it: We can see how the first SNP now with AA=c is included by the expression AA='ACTGactg', so we are targeting ancestral alleles with high and low confidence", header=4)


#find and replace
run_bash(" \
    sed \
        --in-place \
        --expression 's/chr20:14370|A||||intergenic_variant||./chr20:14370|A||||intergenic_variant||c/g' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf")
        #change "." by "c,c,c" in a specific SNP, modify in place
        #https://stackoverflow.com/questions/525592/find-and-replace-inside-a-text-file-from-a-bash-command
#filter
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf \
        --annotation CSQ \
        --include 'AA=\"A,C,G,T,a,c,g,t\"'\
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
#reverse the modification of the VCF file
run_bash(" \
    sed \
        --in-place \
        --expression 's/chr20:14370|A||||intergenic_variant||c/chr20:14370|A||||intergenic_variant||./g' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf")

#
print("\n#######################################\n#######################################")
print("exclude SNPs whose ancestral allele is ACGT or acgt: We get a SNP with ancestral = '.'")
print("#######################################\n#######################################")
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf \
        --annotation CSQ \
        --exclude 'AA=\"A,C,G,T,a,c,g,t\"'\
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")

#
print("\n#######################################\n#######################################")
print("after selecting SNPs with ACGT-acgt ancestral allele, select those whose REF is not equal to ancestral")
print("#######################################\n#######################################")
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf \
        --annotation CSQ \
        --include 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT=AA' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
    #select those SNPs for which 
        #we have actual ancestral allele infered ("N", "."...) AND 
        #the REF is different from that ancestral allele AND 
        #the ALT is the ancestral. These are the SNPs we have to switch.
    #I use && instead of & because only one & looks for variants satisfying the condition within sample.
        #For example, say our VCF contains the per-sample depth and genotype quality annotations per SNP and we want to include only sites where ONE OR MORE samples have big enough coverage (DP>10) and genotype quality (GQ>20). The expression -i 'FMT/DP>10 & FMT/GQ>20'. This would select variants for which at least one sample has good coverage and genotype quality.  
        #but we do not want this. We want to select SNPs that fullfill certain requirements across all samples, not just within at least one sample. We want to match the whole record, using features that are similar across samples, like REF, ALT and AA.
            #http://samtools.github.io/bcftools/howtos/filtering.html




#ANSWER JESUS WHEN FINISHED THIS PART


##por aqui
    #we have a problem with case sensitive: if ALT=C and AA=c, ALT is NOT equal to AA
    #so I am creating a new AA field with all alleles as upper case, so we avoid this problem.
    #check what to do with cases like "." because they are not updated in the new AA field, but this should not be a problem, right? They would be out after selecting for ACTG
    #check the new field is the same with case insensitive
    #then you should update the previous code with the new options like extracting AA from CSQ and then remove CSQ
    #then create the list of SNPs for which REF is not AA
    #then you can go to -fixref and use it to switch these snps in the original VCF file
    #check and then go to the real data



#create a tab separated file with the position info and ancestral allele in upper case of each SNP. Then create an index for this tab-delimited file using tabix.
run_bash(" \
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ | \
    bcftools view \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    $(i+1) = toupper($(i+1)); \
                    print $1, $2, $(i+1); \
                    next \
                } \
            } \
        }' \
    > ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv; \
    sed  \
        --in-place \
        --expression '1i #CHROM\tPOS\tAA_upcase' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv; \
    bgzip \
        --force \
        --keep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv; \
    echo 'See the tab file with upper case ancestral alleles'; \
    gunzip \
        --stdout \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv.gz; \
    tabix \
        --sequence 1 \
        --begin 2 \
        --end 2 \
        --force \
        --comment '#' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv.gz; \
    echo 'See the tab file with upper case ancestral alleles after creating the index'; \
    gunzip \
        --keep \
        --stdout \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv.gz")
        #convert to upper case the ancestral alleles using awk
            #Based on "#https://www.biostars.org/p/304979/"
            #start and set 
                #the field separator
                    #"\t", ";" and "=" so we separate the fields like POS, CHROM... and the tags we have inside INFO that are in the format AC=3;AF=1... so we have the tag and its value as a different field
                        #https://www.gnu.org/software/gawk/manual/html_node/Field-Separators.html
                    #you can have sevaral delimiters separated with "|"
                        #https://stackoverflow.com/questions/12204192/using-multiple-delimiters-in-awk
                #the output field separator
                    #"\t" because we want a tsv
            #for i to the number of delimited fields
                #if the field the ancestral allele (AA)
                    #update the next field (value of AA) with the same string that was present but in upper case.
                    #then print the fields 1 (chrom), 2 (pos), i+1 (AA value)
                        #this is the format required by bcftools annotate to create a new field based on a tab delimited file
                            #--annotations: VCF file or tabix-indexed FILE with annotations: CHR\tPOS[\tVALUE]
                    #go to next row
        #save the result as a file
        #add a header to that file (tab separated names) but annotating that new line with "#" to avoid problems with tabix
            #https://unix.stackexchange.com/a/401673
        #compress using bgzip (from samtools) so the file is recognised by tabix, see below
            #https://github.com/samtools/bcftools/issues/668
        #create index of the tab delimited file with SNP positions and upper ancestral alleles
            #Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
            #After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)
            #tabix
                #--sequence: column number for sequence names
                    #chromosome in our case (first column)
                #--begin: column number for region start
                    #position of the SNP in our case (second column)
                #--end: column number for region end (if no end, set INT to -b)
                    #again the SNP position in our case (second column)
                #--force: overwrite existing index without asking
                #--comment: skip comment lines starting with CHAR
                    #in our case we use "#" to comment the first line with the header
                #http://www.htslib.org/doc/tabix.html

#take the indexed and tab-delimited file with the ancestral alleles in upper case and the position to create a new field with upper ancestral alleles, save the new VCF file
run_bash(" \
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep.vcf \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ \
        --annotations ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv.gz \
        --columns CHROM,POS,AA_upcase \
        --header-line '##INFO=<ID=AA_upcase,Number=.,Type=String,Description=\"The AA field from INFO/CSQ after converting alleles to uppercase\">' > \
    ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep_anc_up.vcf")
        #From the CSQ field added by VEP, extract the tag "AA" as a string, which is the ancestral state.
        #remove the CSQ field
        #add a new INFO/Tag using the tab delimited file previously created
        #select the columns from the tab file in which we are interested
        #add the header line for this new tag
        #save as a new file

#
print("\n#######################################\n#######################################")
print("check that the new INFO/TAG with ancestral alleles in upper case is exactly the same than the original AA tag but in uppercase always")
print("#######################################\n#######################################")
count_aa_no_upper = run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep_anc_up.vcf \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    if(toupper($(i+1))== \"A\" || toupper($(i+1))== \"C\" || toupper($(i+1))== \"T\" || toupper($(i+1))== \"G\" && toupper($(i+1))!=$(i+3)){ \
                        count++; \
                        next \
                    } \
                } \
            } \
        }; \
        END{ \
            print count\
            }'", return_value=True).strip()

if (count_aa_no_upper == ""):
    print("TRUE")
else:
    raise ValueError("SERIOUS ERROR! We have not correctly converted to upper case all ancestral alleles")

#
print("\n#######################################\n#######################################")
print("check that REF and ALT are always in upper case")
print("#######################################\n#######################################")
count_ref_alt_no_upper = run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep_anc_up.vcf \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            if(toupper($4)!=$4 || toupper($5)!=$5){ \
                count++; \
                next \
            } \
        }; \
        END{ \
            print count \
        }'", return_value=True).strip()
        #remove genotypes and header
        #process with awk
            #begin by defining the delimiter of the input and the output. In the input we have a the different fields separated with "\t". The INFO field with tag separated by ";" and each tag name followed by "=" and its value.
            #if in a given row, the fourth OR fifth fields (REF and ALT, respectively), do not have the value in upper case, add 1 to the count
                #I understand that 4 and 5 columns should be always REF and ALT
                #https://stackoverflow.com/questions/13067532/awk-and-operator
                #https://stackoverflow.com/questions/12809909/efficient-way-to-count-the-amount-lines-obeying-some-condition
            #go to the next row because we have already checked the fields we are interested in
            #when done, print the count
        #save the count as an object in python without "\n" and the end of the line (using strip for that)

#if the count is not empty, then we have a problem
if (count_ref_alt_no_upper == ""):
    print("TRUE")
else:
    raise ValueError("SERIOUS ERROR! We do not have the all REF and ALT alleles in upper case so we cannot correctly select those SNPs whose REF is not AA, because our ancestral allele data is always upper case")













        #then extract AA tag from the CSQ field 
        #using annotate
            #remove the CSQ field
            #add the new AA_upcase field using the tab delimited file including also a header to be appended to the VCF file about the new field
                #--annotations: VCF file or tabix-indexed FILE with annotations: CHR\tPOS[\tVALUE]
                    #this command can be used to transfer values from a tab-delimited file into a new INFO/TAG annotation. Note that if the TAG is not defined in the VCF header, a header fragment with the definition must be provided via the -h option.
                #--header-line: Header line which should be appended to the VCF header, can be given multiple times.
                #--columns: List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details

                #https://samtools.github.io/bcftools/howtos/annotate.html


            #https://stackoverflow.com/a/44932064/12772630


    #the first snp with "." does not get AA_upcase

    #check more tabix
        #maybe add more examples to check behaviour

    #replace value with awk?
        #https://stackoverflow.com/questions/51258235/replace-particular-column-value-using-awk-if-found

    #the original script with awk used next. It is used to go to the next line once you ahve fullfill the condition. It makes things faster if you have another if after, because if you already satisifed the condition, you do not need to do more stuff there.
        #https://www.tecmint.com/use-next-command-with-awk-in-linux/
        #https://www.biostars.org/p/304979/


    #mete en los codigos previos y en los siguientes lo de sacar AA y quitar CSQ


    #check AA and AA_upcase are the same with case insensitive




#filter
run_bash("\
    bcftools query \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_cleaned_vep_anc_up.vcf \
        --include 'AA_upcase=\"A,C,G,T\" && REF!=AA_upcase && ALT=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n'")

    #check different scenarios
        #just change dummy_example_cleaned_vep_anc_up.vcf
        #check "count_aa_no_upper" fails when setting "." to "a" or "c"
            #it seems awk does not cosnider "A"=="A"




    #check options this plugin, you can select specific data from CSQ...
        #https://samtools.github.io/bcftools/howtos/plugin.split-vep.html



#POR AQUII
#atcg is not considered using ACTG, so we I had to add actg
#do checks about the use of the plugin on VEP data
    #maybe change AA for ancestral? it can be confusing with AA (ALT | ALT)
#generate a list of SNPs with different REF - ancestral
#then use fixref to switch REF/ALT using that list




#with bcftools
    #--include those SNPs for which R != INFO/AA
        #you have to find a way to create INFO/AA because that information is right now, and ask for chr, pos, REF ALT, and export as vcf
    #alternatively, you can directly export chrom, pos, REF, ALT, for all SNPs, then select rows with different REF than Ancestral using awk and generate a BED file, which is also accpeted by annotate
        #http://www.htslib.org/doc/bcftools.html#annotate
    #use +fixref to switch REF/ALT for these SNPs
        #bcftools +fixref file.bcf -Ob -o out.bcf -- -i List_of_1.vcf.gz 
            #https://www.biostars.org/p/411202/
            #https://samtools.github.io/bcftools/howtos/plugin.fixref.html
        #fixref ES PELIGROSO!!



    #WHERE THE CACHE IS BEING INSTALLED?
        #.vep folder
    #WE HAVE TO USE CACHE
        #the question is if we do this step at the end, so we need to upload the whole cache to the HPC, with the risk that you have to be sure you are using the correct version of the cache respect to the version of VEP
        #the alternative is doing locally with the raw VCF files, and then clean the resulting VCF files.


##STRAND, SOME CASES 1 AND OTHER -1!!!!
    #check if you get the same strand in the output non-VCF than in the VCF


#check usage of vep
    #https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html
    #there is no output




#finish checking options of vep installing
    #https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html



#dbi is installed in the container but I get 
    #ERROR: DBI module not found. VEP requires the DBI perl module to function
    #maybe environmental variables is the problem? check your bashrc to check if folder of modules is indicated


#check if we really need Bio::DB::HTS, it just mentioned in the page of AncestralAllele
    #https://useast.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#ancestralallele
    #trying to stop the container so I get the log file with details about the error
        #cpanm --local-lib /opt/cpanm Bio::DB::HTS &
        #sleep 60 
        #cat /root/.cpanm/work/1682462531.14977/build.log
        #thisdoes not work becausr the folder number changes
            #cat: /root/.cpanm/work/1682462531.14977/build.log: No such file or directory
            #baybe force to open with widlcard ?


    #https://github.com/Ensembl/Bio-DB-HTS/issues/91


#several months ago I started the pollarization of the alleles. I installed VEP and used one of its plugins to estimate the ancestral allele of each SNP within a new VCF file. Then, I started working with AWK to convert to upper case the new column with the ancestral allele, so we have high and low-confidence ancestral alleles and they can be used in conditionals. In other words, I am selecting those snps for which REF is not Ancestral in bcftools, and for that I need the same case, as the conditions of bcftools are case sensitive.

#see email from Jesus about VEP installation, and answer once you have solved the upper case problem
    #https://mail.google.com/mail/u/0/?tab=rm&ogbl#drafts/QgrcJHsHpDRJdfjndBxlCjQHdCNBwJJqNSl 

#then you should go to the actual data


#you should also check the number of variants for which polarization could not be done and hence are lost due to this.







#Note about the update of the INFO fields
    #it is important to be sure that the fields you are using for filtering, are updated after subseting samples. Of course, type="snp" will be always "snp" irrespectively of the samples we select, but this is not the case of the number of alleles, because you can have SNPs with 3 alleles considering all 26 populations, but then in GBR they can have only 2 or 1. We are interested in SNPs that are biallelic within the selected population.
    #The same goes for phasing and genotype ^miss. You have to be sure that these filters only consider data from the filtered samples, not fixed data in fields that are not updated.
    #Because of this we have applied all these filters in order.



################################################################
#### function to clean vcf files and create hap - map files ####
################################################################

#chr_pop_combination="GBR_1"
def master_processor(chr_pop_combination):

    #extract selected population and chromosome
    selected_pop = chr_pop_combination.split("_")[0]
    selected_chromosome = chr_pop_combination.split("_")[1]

    #redirect standard output
    import sys
    original_stdout = sys.stdout
        #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
    sys.stdout = open("./scripts/01_hap_map_calcs_outputs/chr" + selected_chromosome + "_" + selected_pop + ".out", "w")
        #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/


    #######################
    # one time operations #
    #######################

    #do first some operations that need to be done just one time per chromosome
    if selected_pop == "GBR":
        
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
            if [[ $n_samples -eq " + str(samples_pedigree.shape[0]) + " ]]; then \
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

        #stats with and without accessibility mask. See dummy example about details of the masks
        print("\n#######################################\n#######################################")
        print("chr " + selected_chromosome + ": See the stats of the whole VCF file with and without selecting SNPs inside the regions that PASS in the pilot and the strict mask")
        print("#######################################\n#######################################")
        print("\n##### NO MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
        print("\n##### PILOT MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                --targets-file ./data/masks/20160622.allChr.pilot_mask.bed \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
        print("\n##### STRICT MASK ####")
        run_bash(" \
            bcftools view \
                --types snps \
                --targets-file ./data/masks/20160622.allChr.mask.bed.gz \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools stats")
                #load the whole VCf file of the selected chromosome, select only SNPs and in one case select variants inside mask "PASS" regions (pilot and strict), while in the second do nothing more. Then see the stats of the resulting BCF files. See dummy examples for further details.



    #######################
    # extract samples IDs #
    #######################

    #select the sample IDs for the selected population
    subset_pop = unrelated_samples.loc[unrelated_samples["pop"] == selected_pop, :]

    #reset index
    subset_pop = subset_pop.reset_index(
        drop=True)
        #drop=True to avoid adding the index as a column

    #check
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we only have the selected pop")
    print("#######################################\n#######################################")
    print(subset_pop["pop"].unique() == selected_pop)
    #stop if False in this check
    if subset_pop["pop"].unique() != selected_pop:
        raise ValueError("SERIOUS ERROR! WE HAVE NOT CORRECTLY SELECTED THE SAMPLES OF POP '" + selected_pop + "' in chromosome '" + selected_chromosome + "'")

    #select the sample IDs
    selected_samples = subset_pop["sample"]

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
    print("chr " + selected_chromosome + " - " + selected_pop + ": see more SNPs and their allele count to see how this data is presented for multiallelic SNPs")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AN %INFO/AC %INFO/AF GTs:[ %GT]\n' | \
        head -12")
            #IMPORTANT:
                #It seems that multiallelic variants separated in different lines have already updated the AC field. Therefore, each line does not have two allele counts, but one.
                #For example, 1:10453:A:C and 1:10452:A:C are both in the same position (10452) and have the same reference allele. They seem to be multiallelic, but each one has only one allele count, which is 1. 
                #This count is correctly updated for the subset of samples, but remember that --samples does not remove the second count. As we saw in the dummy example, we have to update the field with +fill-tags.
                #My hypothesis is that they updated this field using +fill-tags because they indeed say in the paper that they used bcftools to split the multiallelic SNPs in different lines.
                #Note that AF has also 1 value but it is not correct because we subset samples with --samples, and this command only updates AC and AN.

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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --include 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #update the INFO/AC and INFO/AN fields so we avoid two allele counts in the different lines of multiallelic snps. --samples update AC but maintaining the connection between the different lines of the same multialllelic SNP.
                #it seems that 1KGP data has already updated the AC field for each line separated line of a multiallelic SNP.
                #Therefore, this line would not be necessary, but we are applying just in case, to ensure we have only 1 allele count per line.
                #updating again the AC fields would not do anything wrong, just adding the same value that is was.
            #select those variants for which the number of ALT alleles (allele count) is equal to 0 (no ALT at all) or equal to the total number of alleles (AN), i.e., all alleles are ALT.
            #See dummy example for further details.

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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we do not have SNPs with genotype missingness > 0.05. If TRUE, then we do not have to apply further filters about missing genotypes")
    print("#######################################\n#######################################")
    check_missing = run_bash(" \
        n_snps_with_filter=$( \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools view \
                --types snps | \
            bcftools +fill-tags \
                -- --tags AN,AC | \
            bcftools view \
                --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
            bcftools view \
                --include 'COUNT(GT=\"mis\")/N_SAMPLES < 0.05' | \
            bcftools query \
                --format '%ID\n' | \
            wc -l); \
        n_snps_without_filter=$( \
            bcftools view \
                --samples " + ",".join(selected_samples) + " \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            bcftools view \
                --types snps | \
            bcftools +fill-tags \
                -- --tags AN,AC | \
            bcftools view \
                --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
            bcftools query \
                --format '%ID\n' | \
            wc -l); \
        if [[ $n_snps_with_filter -eq $n_snps_without_filter ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi", return_value=True)
            #obtain the ID of all variants after applying or not the filter for genotype missingness < 0.05, then check the numbers are the same
            #return the value to do an explicit check
    #if the check is not TRUE, stop because this is important
    if check_missing.strip() != "TRUE":
        #use strip() to remove "\n" at the end of the string
        raise ValueError("ERROR: FALSE! The filter of < 5% of missingness does not seem to be applied in chr '" + selected_chromosome + "' for pop '" + selected_pop + "'")

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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
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
                #two snps in position 10452 with same REF but different ALT in chromosome 1. They get combined into one line. These are the first multiallelic that appear in the previous command, but they were separated.

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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
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

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": show now only biallelic SNPs that are phased for all samples")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
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

    ##por aqui checking the code
    #added pop name to bed file cleaned, check that the mask is correctly called in the next lines
    #error with "list_snps_with_gen_pos.txt"
        #this should have the chromosome and pop name to avoid interefernece
        #check that all files are correctly named

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": clean BED file before applying the mask selecting only intervals in the selected chromosome")
    print("#######################################\n#######################################")
    run_bash(" \
        awk \
            -F '\t' \
            '{if ($1 == \"chr" + selected_chromosome + "\") print $0}' \
            ./data/masks/20160622.allChr.pilot_mask.bed \
            > ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed; \
        gzip \
            --force \
            ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed; \
        gunzip \
            -c \
            ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        head -5")
            #in the bed file, select those rows for which the chromosome name (first column) is the selected chromosome, printing all fields for these rows. Save as a file and then compress. See the first 5 lines. See dummy examples for further details.

    #
    print("\n#######################################\n#######################################")
    print("check we have selected the correct chromosome and mask type")
    print("#######################################\n#######################################")
    run_bash(" \
        uniq_chrom=$(\
            gunzip \
                -c \
                ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
            awk \
                -F '\t' \
                '{print $1}' | \
            uniq); \
        uniq_mask_type=$(\
            gunzip \
                -c \
                ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
            awk \
                -F '\t' \
                '{print $4}' | \
            uniq); \
        if [[ $uniq_chrom == 'chr" + selected_chromosome + "' && $uniq_mask_type == 'pilot' ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress the bed file of the selected chromosome (previously created), then print the chromosome name (first column) for all intervals, i.e., rows, getting then the unique cases, then save as a variable. Do the same with the 4th column, i.e., the mask type. The first variable should be the selected chromosome and the second the select mask type, i.e., pilot.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": select only regions that are accessible to sequencing according to the accessibility mask")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools view \
            --samples " + ",".join(selected_samples) + " \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        bcftools view \
            --types snps | \
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' | \
        head -7")
            #Use --targets-file to only select SNPs within the intervals defined by a BED file generated by 1kGDP, which is an accessibility mask. Therefore, we select SNPs that are included in regions accessible to sequencing.
            #see dummy example for further details.
            #I have visually inspected the first 120 SNPs without mask and then check what SNPs are retained after applying the mask and if that makes sense with the first intervals in the BED file for the pilot mask. This works PERFECTLY, selecting only SNPs within the range.

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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AN %AC %AC_Hom %AC_Het %AF %MAF %ExcHet %HWE %NS GTs:[ %GT]\n' | \
        head -7")
            #we sue +fill-tags to update fields that are not updated after subsetting like frequency of alternative allele and create some additional fields.
            #I understand that when using --multiallelic + or -, there is no update because the genotypes should not change, you are just spliting or merging the different ALT alleles. If AC/AN has changed sue to the subset, this is updated in the AC/AN fields and these are used to do the combine/split AC/AN fields. The problem is that only AC/AN are updated, not the rest of fields. In addition, --samples maintains 2 allele counts in each line of a splitted multiallelic SNP. In  the dummy example we had to update these fields with +fill-tags to have only 1 count and be able to use this count to filter out monomorphic SNPs. I have checked that splitted multiallelic SNPs in the 1KGP have only 1 count, so I guess the authors updated with fill-tags, but I update before the monomorphc check just in case. See above.
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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
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
        bcftools +fill-tags \
            -- --tags AN,AC | \
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
        bcftools norm \
            --rm-dup exact | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2  \
            --min-alleles 2 | \
        bcftools view \
            --phased | \
        bcftools view \
            --targets-file ./data/masks/20160622.chr" + selected_chromosome + "_" + selected_pop + ".pilot_mask.bed.gz | \
        bcftools annotate \
            --remove INFO,^FORMAT/GT | \
        bcftools +fill-tags \
            -- --tags AN,AC,AC_Hom,AC_Het,AF,MAF,ExcHet,HWE,NS | \
        bcftools view \
            --output ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz \
            --output-type z \
            --compression-level 1")
                #after subseting and filtering, save as a compressed VCF file, selecting the option for best speed (see dummy example for further details)

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

    #my questions to David about SNP filters
        #Inbreeding in 1KGP
            #The new dataset has 3202 samples instead of 2504. It is not obvious if all the new 698 samples (3202-2504=698) are related. Only 600 are implicated in duos/trios, but both Jesús and I think that the other 90 have some degree of inbreeding. They repeat several times in the paper that they added 698 samples related to the original 2504 samples.
            #Within the original set of 2504 samples, I have detected something confusing. According to the new pedigree, 4 samples have their mother/father within the same 2504 samples. In other words, we have 4 duos within the original set of 2504 unrelated samples. This is not new, as Jesús told me, there is a 2015 scientific report ([link](https://www.nature.com/articles/srep17453 "‌")) showing cryptic relatedness in phase 3. The strange thing is that they generate a new pedigree showing duos within the original dataset, but they continue referring to that original dataset as the 2504 unrelated samples.
            #I think we should avoid all the new samples and focus on the original 2504 dataset. Within that original dataset, we could use it in full or remove the 4 samples that are related according to the new pedigree leaving only the parent of each duo. What do you think?
        #Filtering of variants within populations
            #I am using bcftools to select **biallelic** **SNPs that are not monorphic, no duplicates** (i.e., exact same position and REF-ALT alleles) and phased. Importantly, I am doing this within each population, so if a SNP has 3 alleles considering all populations of the panel but only 2 alleles within the selected population, then that SNP is retained for the selected population. This is what we want, right?
            #I have also filtered by percentage of **missing**
                #1KGP authors selected only variants with missingness < 5%.
                #I have additionally removed variants with **any missing genotype**. Is that ok or should I retain all SNPs with missing<5%?
            #Filter by accessibility
                #I am also using the accessibility masks. Jesus told me that it could be relevant to use the mask if we are going to use data from the whole genome. This is the case at least for the climate project because we were thinking of using non-overlapping across the whole genome, not only centered around genes.
                #It is important to note that these masks select regions that are accessible based on the alignment of whole genome **low coverage** data of 2691 phase3 samples to hg38. Do you think it is ok if the masks are based on low instead of high coverage data?
                #Also, one of the filters they apply to consider a region accessible or not is if "_base is an N in the reference genome GRCh37_".
                #Finally, I have checked how many SNPs are lost in chromosome 1 as an example. The VCF file without any filter has 5,013,617 SNPs, with the less stringent mask this number decreases to 4,616,062, while with the more stringent mask it gets reduced to 3,576,231. Therefore, we lose 1.5 million SNPs.
                #I am not sure whether we should use these masks given they are based in the low-coverage genomes and they seem to reduce the number of SNPs a lot. What do you think?
            #filter by MAF?
                #At least for iHS, we usually remove SNPs with MAF<0.05 directly in hapbin.
                #I guess I should **not** apply a MAF filter right now in the vcf files, but then Elise will apply the required MAF filters when calculating the different summary statistics, right?
            #filter by HWE within the population?
                #1KGP authors selected SNPs that pass HWE filter (p-value>1e-10) in at least one superpopulation, but it is possible that some of these SNPs violate HWE within a specific population.
                #Should I also filter by HWE within each specific population?
    #answer of David
        #You should apply the filter with less than 5% missing just as the 1KGP authors.
            #Not add any filter about missing, because the data is already filtered for less 5% of missing.
        #For HWE you should also only apply what the 1KGP authors did and not filter further for specific populations.
            #Not apply any filter about HWE, because the authors already did at the superpopulation level.
        #Biallelic SNPs per specific population are ok.
            #so I can discard SNPs with less than 2 alleles and more than 2 alleles within each specific population.
            #I understand this includes remove monomorphic. For what reason do you want a SNP that is fixed within your population?
            #Also I guess I can remove SNPs that all their genotypes are phased within the selected population.
        #Masks based on low coverage are ok, you should use the less stringent masks.
            #The less stringent is the pilot.
        #You can use the 2,504 individuals and remove the four related individuals yes.
            #We left 2500 samples.
        #And you are right for MAF, the filtering at MAF>5% is done by the scripts for summary statistics. It is only for the focal SNPs methods like iHS still use the SNPs with lowe MAF around the focal SNPs so it would be an error to remove all SNPs with MAF<5%.
            #I understand that not all summary statistics use the SNP with lower MAF, so we would lose information for these statistics.
            #It is better to specifically filter for MAF when calculating iHS.
        #my follow-up questions
            #When you said "You should apply the filter with less than 5% missing just as the 1KGP authors", you mean that I should not add any further filter about missingness within each population and just use what they already did, i.e., < 5% missing considering all pops, right? Because this filter was already applied in the data I downloaded. 
            #The same would go for HWE. The filter across superpopulations is already applied in the data, so I will just leave it as it is in that regard.
            #Just to double check, I remove a SNP that is monomorphic within a given population, even if it is not monomorphic considering the whole panel. The same goes for phasing, I select SNPs with all genotypes being phased for the specific population, irrespectively of the phasing in other pops.
            #david answered YES to all this.
        #reduction in the number of SNPs after applying filters
            #For example
                #chromosome 1 has 5,013,617 SNPs
                #removing monomorphic SNPs within the british reduces the number to 4,067,698 SNPs.
                #after applying the rest of filters I get 852,072 SNPs for the british, i.e., 4 million of SNPs lost.
            #Is this ok? or it is too much?
                #ok 
        #as I am filtering within pop, each pop will have different snps
            #is this ok?
                #ok

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



    ##ANCESTRAL/DERIVED...

    #por aquiii





    ##########################################
    # calculate map file within selected pop #
    ##########################################

    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": STARTING MAP FILE CALCULATION")
    print("#######################################\n#######################################")
        #we are going to calculate the map of each population directly using the SNPs in its vcf file FILTERED WITHIN population. Then, when we know for which of the SNPs we have genetic position, we can further filter the VCF file and convert to hap.
        #an alternative would be just take all the SNPs in the raw VCF file and calculate their genetic position.
            #it should be ok regarding the ID of the SNPs because REF/ALT are split when using multiallelics, and these fields are NOT switched based on the frequency of the SNPs in specific subsets.
        #I am going for the first option just to be completely sure I am using the SNPs (and positions) of the selected population.
            #If it is too slow this option, think about the other one.

    #Instructions david
        #I understand that SNPs with genetic position are NO useful for any summary statistic, right? So I can safely remove these SNPs from the VCF and hap files right?
            #ok
        #format of ID
            #in the map file can I use the format "CHROM:POS_REF_ALT" for the ID? 
            #I think remember that the map files I originally got from you in the previous project (before decode2019 conversion) used as ID just the physical position. Not sure if there is any specific reason for doing that.
            #not asked, by irrelevant question
        #data format decode map
            #in the decode map, they say clearly that the data is aligned to hg38. Also they do not specify if the coordinates are 1 or 0-based so I assume they are 1-based, base 1 in the map is base 1 in the genome, to me this should be the default. 
            #1KGP data is also aligned to hg38 and is 1-based.
            #Therefore I can just use the position of the SNPs in 1KGP to calculate their genetic position in the decode map, right?
                #David said that there is not problem if the map and hap files are in the same format. 
                #that is the case because the map file is calculated with decode map, and hap file with 1KGP VCF file, and as I said, I have no reason to think that the decode map is not 1-based.


    ##extract the SNP positions
    #extract snps from the cleaned VCF file
    run_bash(" \
        bcftools view \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools query \
            --format '%CHROM %ID %POS %REF %ALT\n' | \
        gzip \
            --force \
        > ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz")
            #from the cleaned VCF file, extract the chromosome, position, REF/ALT to compare with the positions in the raw hap file (see below), compress and save as a file

    #Note about the format of the positions
    #pos in VCF files v4.2 is 1-based according to the specification file (this is the format of 1KGP data). Therefore, we have here 1-based coordinates.
        #POS - position: The reference position, with the 1st base having position 1. Positions are sorted numerically, in increasing order, within each reference sequence CHROM. It is permitted to have multiple records with the same POS. Telomeres are indicated by using positions 0 or N+1, where N is the length of the corresponding chromosome or contig. (Integer, Required)
            #https://samtools.github.io/hts-specs/VCFv4.2.pdf

    #required format according to hapbin
        #The map files (--map) should be in the same format as used by Selscan with one row per variant and four space-separated columns specifiying 
            #chromosome, 
            #locus ID, 
            #genetic position
            #physical position.
        #therefore, we still need to 
            #include the genetic position
            #remove the allele names


    ##load the raw_map file
    snp_map_raw = pd.read_csv(\
        "./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_raw.map.gz", \
        sep=" ", \
        header=None)
    print("See raw map file loaded in python")
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the correct number of SNPs in the map file loaded in python")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
            wc -l); \
        if [[ $n_snps -eq " + str(snp_map_raw.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    #rename the columns
    snp_map_raw = snp_map_raw.rename(\
        {0: "chr", 1: "id_old", 2: "pos", 3: "ref", 4: "alt"}, \
        axis=1)
            #we name ID as old because this is the ID coming from the VCF file, which we need for selecting those variants in the VCF file with genetic position. The final ID will be in plink format, see below.
            #use a dict with old and new column names. indicated we are renaming columns (axis=1)
    print("see map file with renamed columns")
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": create a new ID variable following plink 2.0 format and check it was correctly created")
    print("#######################################\n#######################################")
    #the original ID column is in the format "CHR:POS:REF:ALT" but the problem is that some SNPs have their position indicated in the ID is shifted. 1KGP authors separated multiallelic SNPs in different lines with bcftools norm and then shifted their position so they could be phased, combing back to the original position afterwards (see next line). I guess during that process, they updated the IDs using chrom and the shifted position was used in the ID. Therefore, even POS comes back to the original position, the ID remains with the shifted position. I have checked several multiallelic SNPs, and they have all the same issue with the position in the ID.
        #From README ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf"): "SHAPEIT2 does not handle multiallelic variant phasing. To phase both biallelic and multiallelic variants we first split the multiallelics into separate rows while left-aligning and normalizing INDELs using bcftools norm tool (Li, 2011). Next, we shifted the position of multiallelic variants (2nd, 3rd, etc ALT alleles) by 1 or more bp (depending on how many ALT alleles there are at a given position) to ensure a unique start position for all variants, which is required for SHAPEIT2. We shifted the positions back to the original ones after phasing".
    #This is not very convenient because when later we create the hap files, SNPs will have the plink format for ID to avoid avoid strand flips (CHR:POS_REF_ALT), so we are going to have different IDs between hap and map files, making more difficult to do checks.
    #Therefore, we are going to update the ID of each SNP using plink format, and ensuring in this way SNPs will be names the same in hap and map files.
    snp_map_raw["id"] = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
    #check
    check_id = snp_map_raw["chr"] + ":" + snp_map_raw["pos"].astype("str") + "_" + snp_map_raw["ref"] + "_" + snp_map_raw["alt"]
        #make a series combining chromosome, pos, ref and alt, and using the corresponding separators
    print(check_id.equals(snp_map_raw["id"]))
        #check it is identical to id
    print(snp_map_raw)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the ref/alt columns as we have this information already included in the ID")
    print("#######################################\n#######################################")
    snp_map_raw = snp_map_raw.drop(["ref", "alt"], axis=1)
    print(snp_map_raw)


    ##load and explore the decode2019 map

    #I know that the original 2019 decode map is alligned to hg38. Also, I assume that the decode2019 map is 1-based because they do not specify is 0-based. I assume that if you say anything, base 1 in your coordinates is base 1 in the genome. I assume this is the default.
        #Data S3.genetic.map.final.sexavg.gor.gz:
            #average genetic map computed from the paternal and maternal genetic maps, which were in turn computed from the paternal and maternal crossover, respectively. The data columns are as follows: Chr (chromosome), Begin (start point position of interval in GRCh38 coordinates), End (end point position of interval in GRCh38 coordinates), cMperMb (recombination rate in interval), cM (centiMorgan location of END POINT of interval)
            #Page 85 of "aau1043-halldorsson-sm-revision1.pdf"

    #1KGP is aligned to hg38 (see paper) and coordinates are 1-based as VCF format 4.2 has 1-based coordinates.

    #therefore, we have the same position format in both datasets, so we can just use the decode 2019 map to calculate the genetic position of each SNP.

    # 
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first lines of the Data S3 of decode paper, which is the sex average map (see above). The file has header")
    print("#######################################\n#######################################")
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        head -20")

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the header and save")
    print("#######################################\n#######################################")
    run_bash("\
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3.gz | \
        awk \
            'NR>7' | \
        gzip \
            --force > \
        ./data/decode_2019/aau1043_datas3_no_header.gz; \
        gunzip \
            --stdout \
            ./data/decode_2019/aau1043_datas3_no_header.gz | \
        head -5")
            #decompress, send to stdout
            #then select any row whose number is larger than 7, i.e., from 8th row and forward.
                #https://www.baeldung.com/linux/remove-first-line-text-file 
            #compress and save
            #decompress and see first lines

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": load decode 2019 map into python")
    print("#######################################\n#######################################") 
    decode2019_map = pd.read_csv(\
        "./data/decode_2019/aau1043_datas3_no_header.gz", \
        sep="\t", \
        header=0, \
        low_memory=False)
    #rename columns in lower case
    decode2019_map=decode2019_map.rename({"Chr":"chr", "Begin":"begin", "End":"end", "cMperMb":"cM_Mb", "cM":"cM"}, axis=1)
    print(decode2019_map)
        #Average genetic map computed from the paternal and maternal genetic maps.
        #The data columns are as follows:
        #Chr (chromosome)
        #Begin (start point position of interval in GRCh38 coordinates)
        #End (end point position of interval in GRCh38 coordinates)
        #cMperMb (recombination rate in interval)
        #cM (centiMorgan location of END POINT of interval)

        #I did a lot of checks on this map regarding overlapping of the intervals etc, check recomb_v3.R in method_deep paper for further details.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": subset decode map for the selected chromosome")
    print("#######################################\n#######################################") 
    decode2019_map_subset = decode2019_map.loc[decode2019_map["chr"] == "chr"+str(selected_chromosome),:]
    print(decode2019_map_subset)
    print("Do we selected the correct chromosome?")
    print(decode2019_map_subset["chr"].unique() == "chr"+str(selected_chromosome))
    #
    print("remove the full decode map")
    del(decode2019_map)
    import gc
    gc.collect()

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": define function to calculate genetic position per SNP")
    print("#######################################\n#######################################") 
    #selected_snp_id=snp_map_raw.iloc[10000]["id"] #snp with cM value just in its position
    #selected_snp_id=snp_map_raw.iloc[5000]["id"] #snp with cM values at both sides
    #selected_snp_id=snp_map_raw.iloc[0]["id"] #snp without cM values around
    def gen_pos(selected_snp_id):

        #extract the row in the raw map for the selected SNP
        selected_snp_row = snp_map_raw.loc[snp_map_raw["id"] == selected_snp_id,:]

        #check we have the correct chromosome
        check_0 = selected_snp_row["chr"].unique()[0] == "chr"+str(selected_chromosome)

        #extract position of the selected snp
        selected_snp_physical_pos = selected_snp_row["pos"].to_numpy()[0]

        #extract old ID (this follows VFP file format so we can use them to filter it)
        selected_snp_old_id = selected_snp_row["id_old"].to_numpy()[0]

        #select those deCODE intervals that are at least 1 MB close to the selected SNP: I have search for genetic position data around each SNP of each population, 1MB at each side. I think remember you told me that if we do not find data points at both side of the SNP, we can safely remove it. In that way, we include areas with low recombination (possible haplotypes) but not areas with a lot of missing data.
        decode2019_map_subset_around_snp = decode2019_map_subset.loc[\
            (decode2019_map_subset["end"] >= (selected_snp_physical_pos - 1000000)) & \
            (decode2019_map_subset["end"] <= (selected_snp_physical_pos + 1000000)), :]
                #we are only interested in the END coordinate because the cM data of each interval came from the end of the interval. Indeed, the start coordinate of the next interval is the same of the end of the previous one. Therefore, we focus on end coordinate.
                #Note that we are using 1000000 directly. If a SNP is at 1000001, then 1000001-1000000=1, length(1:1000001) is equal to 1000001, which is not exactly 1MB, but this is only 1 base of difference. This is not important.
                #Also note that for SNPs between base 0 and 1000kb, the difference between SNP position and 1000kb will be negative, but this is OK:
                    #If a SNP is before base 1000kb, then there is less than 1000kb bases to look for decode intervals before the SNP, reaching base 0.
                    #therefore, having a negative value would mean the same than just look up to zero (there are not decode intervals below zero).
                        #SNP at position 500kb
                            #500kb+1000kb=1500kb
                            #500kb-1000kb=-500kb
                            #there is no enough space at the left of the SNP to look for decode intervals up to 1000kb, so we have to reach 0, which is the same than looking for values equal or higher than a negative given that no decode interval has a negative position.
                    #this will be the case until a SNP in position 1001kb, as 1001kb-1000kb would be 1, so we do not look for decode intervals below 1.
                    #As we move foward from 1000kb, the lower limit of the window starts moving away from 1.
                    #indeed, using the absolute value would not work
                        #SNP at position 500kb
                            #1000kb+500kb=1500kb
                            #1000kb-500kb=500kb
                            #the lower limit cannot be 500kb, when the SNP is at 500kb. We would automatically lose this SNP.
                            #you would need -500 to 1500kb.
                    #one concern about this is that for some SNPs we are looking for decode intervals in a smaller region, but SNPs below 1MB are not frequent. 
                        #For example, in chromosome 1, only 1339 out of 800K are at a coordinate below 1000kb. Therefore, this does not seem to be a problem. I have not tested it, but I guess the same would go for SNPs close to the end of the chromosome, this would be a small proportion of the total number of SNPs.
                        #More important, there are NO decode interval below base 500kb, so we will discard any SNP before that base because no cM value will be available to the left in order to interpolate. Therefore, the importance of this issue is very limited.       

        #select those intervals before and after the selected SNP
        intervals_lower_end = decode2019_map_subset_around_snp.loc[(decode2019_map_subset_around_snp["end"] < selected_snp_physical_pos), :]
        intervals_upper_end = decode2019_map_subset_around_snp.loc[(decode2019_map_subset_around_snp["end"] > selected_snp_physical_pos), :]

        #select those decode intervals with the same position than the selected SNP
        interval_same_pos = decode2019_map_subset_around_snp.loc[decode2019_map_subset_around_snp["end"] == selected_snp_physical_pos, :]

        #if we have deCODE intervals 1MB around the selected SNP, i.e., we have intervals at both sides, intervals ending before and after the selected SNP OR we have a deCODE interval ending exactly at the SNP. In the second case if you have cM value in the exact position of the selected SNP, then you do not need intervals at both sides.
        if (intervals_lower_end.shape[0]>0) & (intervals_upper_end.shape[0]>0) | (interval_same_pos.shape[0]>0): 
            #the first condition do not need equal because in the next condition (after "|") we consider the option of equal coordinate between window extreme and deCODE end interval.

            #checks
            check_1=np.unique(\
                (decode2019_map_subset_around_snp["end"] >= (selected_snp_physical_pos - 10**6)) & \
                (decode2019_map_subset_around_snp["end"] <= (selected_snp_physical_pos + 10**6)))[0]

            #if we dot NOT have an interval with an end coordinate exactly similar to the selected SNP
            if (interval_same_pos.shape[0] == 0):

                #check
                check_2 = np.unique(intervals_lower_end["end"] < selected_snp_physical_pos)[0]
                check_3 = np.unique(intervals_upper_end["end"] > selected_snp_physical_pos)[0]

                #from the intervals below the extreme window, select the biggest and hence closest to the extreme window   
                lowest_interval = intervals_lower_end.loc[intervals_lower_end["end"] == max(intervals_lower_end["end"]),:] 
                    #we cannot have two cases with the same value because the coordinates are in increasing order, the coordinate of an interval is bigger than the previous one.

                #from the intervals above the extreme window, select the smallest and hence closest to the extreme window
                highest_interval = intervals_upper_end.loc[intervals_upper_end["end"] == min(intervals_upper_end["end"]),:] 
                    #we cannot have two cases with the same value because the coordinates are in increasing order, the coordinate of an interval is bigger than the previous one.

                #check that the end coordinate with lowest difference respect the SNP is the selected in the previous step both for the lower and higher intervals
                check_4a = (intervals_lower_end.loc[\
                    np.abs(intervals_lower_end["end"]-selected_snp_physical_pos) == \
                    np.min(np.abs(intervals_lower_end["end"]-selected_snp_physical_pos)), \
                    "end"] == lowest_interval["end"]).to_numpy()[0]
                check_4b = (intervals_upper_end.loc[\
                    np.abs(intervals_upper_end["end"]-selected_snp_physical_pos) == \
                    np.min(np.abs(intervals_upper_end["end"]-selected_snp_physical_pos)), \
                    "end"] == highest_interval["end"]).to_numpy()[0]


                ##calculate cM value of the snp
                #extract the centimorgan of each the closest deCODE intervals to the SNP
                left_cM = lowest_interval["cM"].to_numpy()[0]
                right_cM = highest_interval["cM"].to_numpy()[0]

                #calculate the distance from each interval to the SNP
                distance_left_end = (selected_snp_physical_pos - lowest_interval["end"]).to_numpy()[0]
                distance_right_end = (highest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]
                    #We do not need to include both extremes, we want the distance from one point to another. Imagine the window begins at 1 and ends at 3. Including both extremes, the size of the window is 3, you have 3 bases. However, the distance from the point 1 to 3 is 2 (3-1=2). We want the distance between two points with centiMorgan values.

                #check that calculating the distance with abs and changing order gives the same result
                check_5a = distance_left_end == abs(lowest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]
                check_5b = distance_right_end == abs(highest_interval["end"] - selected_snp_physical_pos).to_numpy()[0]

                #check that the sum of the physical distance of each deCODE end point to the SNP is the same than the total distance between the deCODE end points
                check_6 = distance_left_end + distance_right_end == np.abs(lowest_interval["end"].to_numpy()[0] - highest_interval["end"].to_numpy()[0])

                #calculate the genetic distance using the formula of David
                genetic_distance = left_cM + (right_cM - left_cM) * distance_left_end / (distance_left_end + distance_right_end)
                    #Explanation of David: In that case, you have to find the genetic map position of a single SNP (or in general a position of the genome). To assign a position to each SNP, you can use the two genetic map positions (end points as you described) left and right form the SNP. You can then consider that the genetic position increases linearly between the two left and right positions. For example, if a SNP is between a genetic position on the left at 100 cM, and a genetic position on the right at 102 cM, and the SNP is located 20 kb from the left genetic position but 80kb from the right position, then the SNP will be located at genetic position: 100 cM + (102 cM-100 cM) * 20kb / (20kb+80kb) = 100.4 cM. Of course, if the SNP is right on the coordinate of an end point, then just use the genetic map position directly for that SNP.
                    #My explanation: What David is doing is 100 + ((102-100)*20)/(20+80). This gives exactly 100.4. David is using the rule of three (https://en.wikipedia.org/wiki/Cross-multiplication#Rule_of_Three). You have three points, A, B and C. If the physical distance distance A-C is 100 kb (20+80) and the genetic distance between these points is 2 cM (102-100) , what would be the genetic distance between A-B if these points are separated by 20 kb? ((102-100 cM) * 20 kb) / (20+80 kb); ((2 cM) * 20 kb) / (100 kb); ((2 cM) * 20 kb) / (100 kb); (40 cM * kb) / 100 kb; 0.4 cM. 0.4 is the genetic distance between A and B. Now we can sum 0.4 and the genetic position of A, to get the genetic position of B in the genome. 100 cM + 0.4 cM = 100.4 cM.
                    #If you the point for which you calculate the genetic distance is exactly in the middle of the two points with 100 and 102 cM of genetic distance, the resulting genetic distance would be exactly in the middle, i.e., 101: 100 + ((102-100)*50)/(50+50). 50 is the physical distance between cM point and the point of interest. 
                    #This method assumes that relationship between genetic distance and physical distance between two points is lineal and stable, so you can estimate the genetic distance based on the physical distance in the genomic region encompassed by these points. Note that you are using point that are at least 1MB close to the point under study, therefore, we are estimating the genetic distance using the relationship between physical and genetic distance in a specific genomic region, not the whole genome.
                    #see figure 31 for further details.
            else:

                #if not and hence we have an deCODE interval exactly in the SNP position
                genetic_distance = interval_same_pos["cM"].to_numpy()[0]

                #set NA for the rest of results. They are not needed.
                check_2 = np.nan
                check_3 = np.nan                
                check_4a = np.nan
                check_4b = np.nan
                check_5a = np.nan
                check_5b = np.nan
                check_6 = np.nan
                left_cM = np.nan
                right_cM = np.nan
                distance_left_end = np.nan
                distance_right_end = np.nan
        else:

            #if not, and hence we cannot calculate the cM of SNP
            genetic_distance = np.nan

            #set NA for the rest of results. They are not needed.
            check_1 = np.nan
            check_2 = np.nan
            check_3 = np.nan
            check_4a = np.nan
            check_4b = np.nan
            check_5a = np.nan
            check_5b = np.nan
            check_6 = np.nan
            left_cM = np.nan
            right_cM = np.nan
            distance_left_end = np.nan
            distance_right_end = np.nan

        #save results
        return(tuple([selected_chromosome, selected_snp_id, selected_snp_old_id, selected_snp_physical_pos, check_0, check_1, check_2, check_3, check_4a, check_4b, check_5a, check_5b, check_6, genetic_distance, left_cM, right_cM, distance_left_end, distance_right_end]))

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": Run the function on just one snp")
    print("#######################################\n#######################################")
    print(gen_pos(snp_map_raw.iloc[5000]["id"]))

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": run function across SNPs")
    print("#######################################\n#######################################")
    #we do not use pool.map because we will only want to use 1 core. The parallelization will be done in the parent function across chromosome*pop combinations
    #map seems to be faster than loop even using just 1 core, although there is not a big difference
        #https://www.linkedin.com/pulse/loops-maps-who-faster-time-space-complexity-we-coming-george-michelon/
    final_genetic_pos = list(map(gen_pos, snp_map_raw["id"]))
    #final_genetic_pos = list(map(gen_pos, snp_map_raw.iloc[5000:5100]["id"]))

    #convert the tuple to DF and add the column names
    final_genetic_pos_df = pd.DataFrame(final_genetic_pos, columns=["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "selected_snp_physical_pos", "check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6", "genetic_distance", "left_cM", "right_cM", "distance_left_end", "distance_right_end"])
    print("see results:")
    print(final_genetic_pos_df)

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": all checks of genetic position calculation are True?")
    print("#######################################\n#######################################") 
    print(final_genetic_pos_df[["check_0", "check_1", "check_2", "check_3", "check_4a", "check_4b", "check_5a", "check_5b", "check_6"]].all())
        #important:
            #all() does not consider nan, so if you have nan and the rest True, the output is True.
            #this is ok for us, because we use nan in some conditions.
    
    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the correct number of SNPs in the calculation of genetic position")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps=$(\
            bcftools view \
                --no-header \
                ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
            wc -l); \
        if [[ $n_snps -eq " + str(final_genetic_pos_df.shape[0]) + " ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #count the number of lines in the cleaned VCF file without the header, and check that number is equal to the number of SNPs we have in the map file loaded in python 

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that we have the exact same snps than in the raw map")
    print("#######################################\n#######################################")
    print(np.array_equal(
        snp_map_raw["chr"].to_numpy(),
        ("chr" + final_genetic_pos_df["selected_chromosome"]).to_numpy()))
    print(np.array_equal(
        snp_map_raw["id"].to_numpy(),
        final_genetic_pos_df["selected_snp_id"].to_numpy()))
    print(np.array_equal(
        snp_map_raw["id_old"].to_numpy(),
        final_genetic_pos_df["selected_snp_old_id"].to_numpy()))
    print(np.array_equal(
        snp_map_raw["pos"].to_numpy(),
        final_genetic_pos_df["selected_snp_physical_pos"].to_numpy()))
 
    
    ##recalculate genetic distance of each SNP
    #but we exclude SNPs that have cM exactly in their position in the decode map or SNPs without decode data around
    final_genetic_pos_df_check_dist_calc = final_genetic_pos_df.loc[\
        (~final_genetic_pos_df["genetic_distance"].isna()) & \
        (~final_genetic_pos_df["left_cM"].isna()), :]
    #calculate the physical distance between the end of the two deCODE ranges, that is, the distance from the closest interval to the SNP of the left PLUS the distance from the SNP to the closest deCODE interval from the right
    phys_distance_decode_intervals = final_genetic_pos_df_check_dist_calc["distance_left_end"] + final_genetic_pos_df_check_dist_calc["distance_right_end"]

    #calculate the genetic distance between the end of the two deCODE interval closest to the SNP (on both sides, right and left)
    gen_distance_decode_intervals = final_genetic_pos_df_check_dist_calc["right_cM"] - final_genetic_pos_df_check_dist_calc["left_cM"]

    #extract the physical distance from the SNP to the closest deCODE interval to the left
    phy_distance_left_decode = final_genetic_pos_df_check_dist_calc["distance_left_end"]

    #calculate the genetic distance: If the physical distance between the end of deCODE intervals (phys_distance_decode_intervals) corresponds with a known genetic distance (gen_distance_decode_intervals), the physical distance from the closest deCODE intervals from the left to the selected SNP (phy_distance_left_decode) would correspond with X; thus X = (gen_distance_decode_intervals*phy_distance_left_decode)/phys_distance_decode_intervals X is the genetic distance from the closest deCODE interval from the left to the SNP If you sum this to the genetic position of that closest deCODE interval from the left (final_genetic_pos$left_cM), you would have the genetic distance of the selected SNP. The genetic position of that interval gives the cM value until that point, and you just calculated the rest of cM increase until the SNP 
    new_genetic_distance = ((gen_distance_decode_intervals * phy_distance_left_decode) / phys_distance_decode_intervals) + final_genetic_pos_df_check_dist_calc["left_cM"]
        #see figure 31 and the calculation of genetic distance in "recomb_calc" function for the full explanation using the words of David and also my interpretation.

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": compare the new genetic distance and the distance previously calculated. It is ok to have False here if the next check is ok")
    print("#######################################\n#######################################")
    raw_check_gen_dis = new_genetic_distance == final_genetic_pos_df_check_dist_calc["genetic_distance"]
    print(raw_check_gen_dis.groupby(raw_check_gen_dis).count())
        
    #
    #from results, extract ID of snps with NA for last checks but with data for genetic position
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": check that the cases with NA for last checks but with genetic distance are the cases of SNPs in a position exactly with deCODe data")
    print("#######################################\n#######################################")
    cases_gen_pos_no_last_checks = final_genetic_pos_df.loc[\
        (final_genetic_pos_df["check_4a"].isna()) & \
        (~final_genetic_pos_df["genetic_distance"].isna()), "selected_snp_id"]
    #extract the ID of snps with a position that have deCODE genetic position
    snps_with_decode_data = final_genetic_pos_df.loc[\
        final_genetic_pos_df["selected_snp_physical_pos"].isin(decode2019_map_subset["end"]), "selected_snp_id"]
    #print the check
    print(cases_gen_pos_no_last_checks.equals(snps_with_decode_data))


    ##final map file
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": prepare final map file")
    print("#######################################\n#######################################")
    #subset only the columns for map files
    final_genetic_pos_map_file = final_genetic_pos_df[["selected_chromosome", "selected_snp_id", "selected_snp_old_id", "genetic_distance", "selected_snp_physical_pos"]]
        
    #add chrom
    final_genetic_pos_map_file["selected_chromosome"] = "chr"+final_genetic_pos_map_file["selected_chromosome"]
        #WARNING HERE
    print(final_genetic_pos_map_file)
        #save the chromosome, ID, genetic position and physical position. This is the format expected by hapbin
            #https://github.com/evotools/hapbin

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": remove the SNPs without genetic position from the VCF file")
    print("#######################################\n#######################################")

    #select old_id from the map that have genetic position
    #this ID is the original retained from the VCF file, so we can use it to subset the VCF file
    snps_id_with_gen_pos = final_genetic_pos_map_file.loc[\
        ~final_genetic_pos_map_file["genetic_distance"].isna(), \
        "selected_snp_old_id"]

    #
    print("see the number of SNPs removed due to the lack of genetic position")
    print(final_genetic_pos_map_file.shape[0] - len(snps_id_with_gen_pos))

    #save the names in a txt file
    with open(r"./results/cleaned_vcf_files/list_snps_with_gen_pos.txt", "w") as fp:
        fp.write("\n".join(snps_id_with_gen_pos))
            #each name in a different line so we have to add "\n" to the name
            #https://pynative.com/python-write-list-to-file/
        fp.write("\n")
            #add empty line at the end

    #
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": filter the already cleaned VCF with bcftools")
    print("#######################################\n#######################################")
    #this file is cleaned regarding biallelic snps, duplicates... but need to retain only SNPs with genetic position
    #We could do this by just creating before the hap file, extract snp positions from there, calculate genetic position and then remove from the hap those rows of SNPs without genetic position. The problem is that we would do that by row index instead of SNP ID, at least if we use the final hap file, so we are going for this option better. In addition, we would have snps that cannot be used in the VCF file because they do not have genetic position. With the other approach we would have a VCF with all SNPs filtered and another one with only snps with genetic position.
    
    #filter
    run_bash("\
        bcftools view \
            --include ID==@./results/cleaned_vcf_files/list_snps_with_gen_pos.txt\
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + ".vcf.gz | \
        bcftools view \
            --output ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --output-type z \
            --compression-level 1")
            #include those SNPs for which ID is included in the list of SNPs with genetic position and save the resulting VCF file
                #https://www.biostars.org/p/373852/

    #
    print("see header of the fully filtered VCF file and some genotypes")
    run_bash(" \
        bcftools head \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz")
    run_bash(" \
        bcftools view \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --no-header | \
        head -5")

    #
    print("check that IDs in the filtered VCF file are the same than the ones in the list of IDs used as input to filter")
    run_bash(" \
        bcftools query \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
            --format '%ID\n' \
        > ./results/cleaned_vcf_files/ids_vcf_after_filter.txt; \
        file1='./results/cleaned_vcf_files/list_snps_with_gen_pos.txt'; \
        file2='./results/cleaned_vcf_files/ids_vcf_after_filter.txt'; \
        STATUS=$(cmp --silent $file1 $file2; echo $?); \
        if [[ $STATUS -eq 0 ]]; then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi; \
        rm $file2")
        #get the IDs in the finally filtered VCF file and save the file
        #create two variables with the names of this file and also the name of the file with the list of IDs used as input to filter        
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
        #remove the file created for this check

    #
    print("remove also the SNPs without genetic position from the map file and check")
    final_genetic_pos_map_file = final_genetic_pos_map_file.loc[\
        ~final_genetic_pos_map_file["genetic_distance"].isna(),:]
    print(final_genetic_pos_map_file["selected_snp_old_id"].equals(snps_id_with_gen_pos))

    #
    print("remove old ID as we have already filtered the VCF file and check")
    final_genetic_pos_map_file = final_genetic_pos_map_file.drop(["selected_snp_old_id"], axis=1)
    print(final_genetic_pos_map_file.columns == ["selected_chromosome", "selected_snp_id", "genetic_distance", "selected_snp_physical_pos"])

    #
    print("see final map and save")
    print(final_genetic_pos_map_file)
        #required format according to hapbin
            #The map files (--map) should be in the same format as used by Selscan with one row per variant and four space-separated columns specifiying 
                #chromosome, 
                #locus ID, 
                #genetic position
                #physical position.
    final_genetic_pos_map_file.to_csv(\
        "./results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz", \
        sep=" ", \
        header=False, \
        index=False)

    #
    run_bash("\
        n_rows=$( \
            gunzip \
                --stdout \
                ./results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NR}'); \
        n_cols=$( \
            gunzip \
                --stdout \
                ./results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            awk \
                -F ' ' \
                'END {print NF}'); \
        if [[ $n_cols -eq 4 && $n_rows -eq " + str(final_genetic_pos_map_file.shape[0]) + " ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #decompress map file to stdout and then calculate the number of rows (NR) and fields (NF). Do that only after the whole file has been read (END)
                #https://www.gnu.org/software/gawk/manual/html_node/Using-BEGIN_002fEND.html
            #the number of columns (fields) should be 4 following salescan format, while the number of rows should be equal to the number of snps we have in the map file loaded in python, which was indeed used to write this .map file.


    ##convert to hap 
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": convert to hap file")
    print("#######################################\n#######################################")
    run_bash(" \
        bcftools convert \
            ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
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

    #note about the output of hap conversion
        #I get an output like this
            #Hap file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.hap.gz
            #Sample file: ./results/hap_map_files_raw/chr1_GBR_IMPUTE2_raw.samples
            #852072 records written, 0 skipped: 0/0/0 no-ALT/non-biallelic/filtered

        #I guess non-biallelic and non-alt should be removed to meet impute requirements, but in our case we already selected those SNPs with 2 alleles only. This explains why we get 0/0/0. This is the number of SNPs with no-ALT, no-biallelic and filtered. All the filters were previously applied.
            #I have checked this in the dummy example.
        
        #if you take the VCF file, filter and then count lines, you get 852072 records, but then say that total is 945919, as when writting the file. "Lines   total/split/realigned/skipped" is produced for some commands like bcftools norm and I have seen in the dummy example that total is the total number of SNPs, while split are those splitting due to multiallelic is the multiallelic flag is used.
            #852072
            #Lines   total/split/realigned/skipped:  945919/0/0/0
            #Lines   total/split/realigned/skipped:  945919/0/0/0

        #the total number of snps in the raw vcf file of chr1 is 5759173, instead of 945919. What it can be happening here is that bcftools norm (the command generating the line output) is run after some (but not all) filters have been applied, thus the input number of SNPs is smaller but the final number (852072). Indeed, if you apply these filters applied before norm, you get 945919 variants, so it makes sense that bcftools norm, which is run after these filters, gives 945919 as total number of snps.

        #The smallest number (852072) is the number of SNPs after we have completely cleaned the vcf file, including the accesibility mask. I have checked that.

    #see first variants
    print("\n#######################################\n#######################################")
    print("chr " + selected_chromosome + " - " + selected_pop + ": see first variants of the hap file")
    print("#######################################\n#######################################")
    run_bash(
        "gunzip -c ./results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2_raw.hap.gz | \
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
    print("chr " + selected_chromosome + " - " + selected_pop + ": check we have the same number of rows (variants) in the cleaned vcf, hap and map files")
    print("#######################################\n#######################################")
    run_bash(" \
        n_snps_vcf=$( \
            bcftools view \
                ./results/cleaned_vcf_files/chr" + selected_chromosome + "_" + selected_pop + "_only_snps_gen_pos.vcf.gz \
                --no-header | \
            wc -l); \
        n_snps_map=$( \
            gunzip \
                -c ./results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_selscan.map.gz | \
            wc -l); \
        n_snps_hap=$( \
            gunzip \
                -c results/hap_map_files/chr" + selected_chromosome + "_" + selected_pop + "_IMPUTE2.hap.gz | \
            wc -l); \
        if [[ $n_snps_vcf -eq $n_snps_hap && $n_snps_vcf -eq $n_snps_map ]];then \
            echo 'TRUE'; \
        else \
            echo 'FALSE'; \
        fi")
            #load the cleaned vcf file with bcftools and show only the snps, with no header, then count the number of lines and save the result
            #decompress the hap file and count the number of line
            #check both numbers are the same

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
        if [[ $nfields_hap_nlines -eq 1 ]]; then \
            nsamples_hap=$(($nfields_hap/2)); \
            nsamples_bcftools=$( \
                cat results/hap_map_files_raw/chr" + selected_chromosome + "_" + selected_pop + "_samples_to_bcftools.txt | \
                grep -c '^'); \
            if [[ $nsamples_hap -eq $nsamples_bcftools ]]; then \
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

    #restore sys.stdout using the previously saved reference to it
    #This is useful if you intend to use stdout for other things
    sys.stdout = original_stdout
        #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/



#####################
#### paralellize ####
#####################

##
print("\n#######################################\n#######################################")
print("create array with all combinations of pops and chromosomes")
print("#######################################\n#######################################")
#get pop and chromosome names
pop_names
chromosomes = [i for i in range(1, 23, 1)]
print("we are going to analyze 26 pops and 22 chromosomes?")
print((len(pop_names) == 26) & (len(chromosomes) == 22))
print("See them")
print(pop_names)
print(chromosomes)

#get all the combinations but first make a dummy example
import itertools
#create two dummy lists, one with strings and other with integers
print("dummy example to get all possible combinations of two lists")
dummy_x = ["marbella", "cuzco", "granada"]
dummy_y = [1, 2, 3]
#product get all possible combinations between the two lists
dumm_combinations = [x+"_"+str(y) for x in dummy_x for y in dummy_y]
print(dumm_combinations)
    #first for each each value of X, and then for each value of Y, combine X and Y, so combine X1 with Y1, X1 with Y2, .... X2 with Y1, X2 with Y2 and so on...
    #y has to be converted to string with it is integer
print("Do we have all dummy combinations?")
print(len(dumm_combinations) == len(dummy_x)*len(dummy_y))
#get all combinations from the actual pops and chromosomes
full_combinations_pop_chroms = [pop+"_"+str(chrom) for pop in pop_names for chrom in chromosomes]
print("Do we have all combinations of chromosomes and populations?")
print(len(full_combinations_pop_chroms) == len(pop_names) * len(chromosomes))
print("is this equivalent to itertools.product?")
print(\
    full_combinations_pop_chroms == \
    [i[0]+"_"+str(i[1]) for i in list(itertools.product(pop_names, chromosomes))])
    #itertools.product gives a tuple for each combination, so you can extract both elements and bind them with join.
        #https://stackoverflow.com/a/34032549/12772630
        #https://docs.python.org/3/library/itertools.html#itertools.product


##run parallel analyses
#open the pool
import multiprocessing as mp
pool = mp.Pool(len(full_combinations_pop_chroms)/2)

#run function across pandas rows
pool.map(master_processor, full_combinations_pop_chroms)

#close the pool
pool.close()

#por aquiii
    #run the script 
    #check the redirection of stdout
    #check the whole script in the meantime


#you need to save in upper case the ancestral allele and save the VCF files in a new folder, indicate in that folder that the REF is not yet ancestral. You will do the polarization within each populatin. We need only biallelic snps (to easily exchange ref by alt if needed), and we need to do this within pop, porque an allele can be biallelic for one pop but not for other.
