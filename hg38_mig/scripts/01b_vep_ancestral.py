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
        -p ./data/dummy_vcf_files/00_dummy_vcf_files_vep; \
    mkdir \
        -p ./scripts/00_ancestral_calcs_outputs; \
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

print_text("calculate the total number of rows across all fasta files to check later with the combined file. We do it using a loop with awk and wc -l to check both give the same", header=4)
print("calculate with awk")
n_rows_across_fastas_awk = run_bash("\
    cd ./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38; \
    list_fasta=$(ls *.fa); \
    for fasta in \"${list_fasta[@]}\"; do \
        awk 'END{print NR}' $fasta; \
    done", return_value=True).strip()
    #it seems that awk consider all files in the loop like continuous!!
print("calculate with wc -l")
n_rows_across_fastas_wc_raw = run_bash("\
    cd ./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38; \
    list_fasta=$(ls *.fa); \
    for fasta in \"${list_fasta[@]}\"; do \
        wc -l $fasta; \
    done", return_value=True)
n_rows_across_fastas_wc = n_rows_across_fastas_wc_raw.split("\n")[-2].replace(" total", "").strip()
    #indeed wc -l gives the number of rows of each file and at the end gives the total number!
        #get the raw list of counts and the select only the total count and 
#I have no checked the behavior in detail, I have just compared two approaches to have more confidence the total number of rows is correct
print("both awk and wc -l give the same number of rows across all fastas?")
if n_rows_across_fastas_awk == n_rows_across_fastas_wc:
    print("YES! GOOD TO GO!")
else:
    print("ERROR! FALSE! WE HAVE A PROBLEM")

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

print_text("check the number of lines in the combined file is the same than the total number of lines across all individual fasta files", header=4)
run_bash("\
    cd ./data/fasta_ancestral/; \
    total_n_rows=$( \
        awk \
            'END{print NR}' \
            homo_sapiens_ancestor_GRCh38_final.fa); \
    if [[ $total_n_rows -eq " + n_rows_across_fastas_awk + " ]];then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi")

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

print_text("bgzip the fasta file with ancestral alleles because this should make things faster when running VEP (see documentation below)", header=4)
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

print_text("remove previous indexes of the fasta file if they are present. These are created the first time VEP is run on the files, so just to be sure we are using the correct index, we remove those existing previously. these indexes will be used estimating the ancestral allele of all chromosomes", header=4)
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
    #I have checked in the information obtained from VEP for this dummy file that the same SNP can be in different transcripts with different strand, but the ancestral allele is always the same.

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
    fi; \
    rm file_check_1.txt; rm file_check_2.txt; \
    ls -l")


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
    #--list: 
        #Parse the VCF header and list the annotation fields

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

print_text("now manually change ancestral of rs6054257 and rs6054252 from '.' to 'c' to check whether our expression catch it: We can see how the first SNP now with AA=c is included by the expression AA='ACTGactg', so we are targeting ancestral alleles with high and low confidence. We are also adding a few cases of '-' and 'N' to see how we deal with that (see below)", header=4)
run_bash(" \
    cd ./data/dummy_vcf_files/00_dummy_vcf_files_vep/; \
    gunzip \
        --stdout \
        dummy_example_vep.vcf.gz | \
    sed \
        --expression 's/chr20:14310|A||||intergenic_variant||./chr20:14310|A||||intergenic_variant||N/g' | \
    sed \
        --expression 's/chr20:14320|A||||intergenic_variant||./chr20:14320|A||||intergenic_variant||-/g' | \
    sed \
        --expression 's/chr20:14350|A||||intergenic_variant||./chr20:14350|A||||intergenic_variant||c/g' | \
    sed \
        --expression 's/chr20:14370|A||||intergenic_variant||./chr20:14370|A||||intergenic_variant||c/g' > dummy_example_vep_2.vcf")
        #decompress the VCF file to avoid problems with sed
        #change "." by "c" in two SNPs and also change in other two to N and -
            #ACTG: high-confidence call, ancestral state supported by the other two sequences
            #actg: low-confidence call, ancestral state supported by one sequence only
            #N: failure, the ancestral state is not supported by any other sequence
            #-: the extant species contains an insertion at this postion
            #.: no coverage in the alignment
            #https://stackoverflow.com/questions/525592/find-and-replace-inside-a-text-file-from-a-bash-command
        #then save as a new file
run_bash("\
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2.vcf \
        --annotation CSQ \
        --include 'AA=\"A,C,G,T,a,c,g,t\"'\
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
print("see how rs6054252 and rs6054257 are now included between those with ancestral allele. We will maintain these changes to have lower case ancestral allele and learn to deal with that (see below)")
run_bash(" \
    cd ./data/dummy_vcf_files/00_dummy_vcf_files_vep/; \
    bgzip \
        --force \
        --keep \
        dummy_example_vep_2.vcf; \
    ls -l")

print_text("see the new VCF file", header=4)
run_bash(" \
    bcftools +split-vep \
        --annotation CSQ \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AF AA=%AA\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2.vcf.gz")

print_text("WE HAVE TO EXCLUDE CASES WHERE REF NOR ALT ARE THE ANCESTRAL ALLELE. This can be the case for a multiallelic SNP that only has two alleles in the selected populations, being the missing one the ancestral. We will just remove actual multiallelic SNPs in this dummy population, i.e., more than 2 alleles are present in the genotype data. Then, we will exclude those SNPs where AA is not REF nor ALT.", header=4)
run_bash("\
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2.vcf.gz |\
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools +split-vep \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ | \
    bcftools view \
        --exclude 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT!=AA' | \
    bcftools view \
        --include 'AA=\"A,C,G,T,a,c,g,t\"' | \
    bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")
    #split multiallelic SNPs
    #update the AN and AC fields to be sure we have the new allele count considering only each pair of REF/ALT
    #remove those monomorphic SNPs. This will remove REF/ALT pairs of a multiallelic SNP if one of the ALTs is not present in the population. For example, the pair A/C in the example will be removed because C is not present.
        #snp_1 REF=A ALT=C A|. A|. A|A AA=C
        #snp_1 REF=A ALT=G A|G G|G A|A AA=C
    #bind again the multiallelic SNPs and remove those that still have more than 1 ALT, because these are truly multiallelic SNPs in the dummy population.
    #Extract AA data from CSQ and create a new, independent field with that.
    #exclude SNPs where AA is not REF nor ALT
        #in the previous example, A/C has been already removed, thus only A/G is present. Therefore, this SNP now does not have the ancestral allele. AA is not present in this dummy population, thus REF!=AA (A!=C) and ALT!=AA (G!=C) and this SNP is excluded.
        #snp_1 REF=A ALT=G A|G G|G A|A AA=C
    #then select only for SNPs with actual ancestral allele estimated, instead of "." or "N"
    #I use && instead of & because only one & looks for variants satisfying the condition within sample.
        #For example, say our VCF contains the per-sample depth and genotype quality annotations per SNP and we want to include only sites where ONE OR MORE samples have big enough coverage (DP>10) and genotype quality (GQ>20). The expression -i 'FMT/DP>10 & FMT/GQ>20'. This would select variants for which at least one sample has good coverage and genotype quality.  
        #but we do not want this. We want to select SNPs that fulfill certain requirements across all samples, not just within at least one sample. We want to match the whole record, using features that are similar across samples, like REF, ALT and AA.
            #http://samtools.github.io/bcftools/howtos/filtering.html
print("IMPORTANT: We can see how rs6054252 is not included in this list even having AA=c and REF=C, thus at least one of the REF/ALT columns includes the AA and then should not be removed by the filter REF!=AA && ALT!=AA. There is a problem with the case because bcftools filters are case sensitive: C is not the same than c")

print_text("select now those cases where REF!=AA but ALT=AA, because these are the cases that can be switched", header=4)
run_bash("\
    bcftools norm \
        --multiallelic -snps \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2.vcf.gz |\
    bcftools +fill-tags \
        -- --tags AN,AC | \
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools +split-vep \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ | \
    bcftools view \
        --exclude 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT!=AA' | \
    bcftools view \
        --include 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT=AA' | \
    bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n'")

print_text("Problem upper vs. lower case ancestral allele", header=4)
print("We have a problem with case sensitive: if ALT=C and AA=c, ALT is NOT equal to AA. So I am creating a new AA field with all alleles as upper case, so we avoid this problem.")


print_text("working on the upper vs. lower case problem", header=3)
print_text("create a tab separated file with the position info and ancestral allele in upper case of each SNP. Then create an index for this tab-delimited file using tabix.", header=4)
run_bash(" \
    bcftools +split-vep \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2.vcf.gz \
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
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/anc_alleles_uppercase.tsv.gz; \
    ls -l ./data/dummy_vcf_files/00_dummy_vcf_files_vep/")
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
                        #next is used to go to the next line once you ahve fullfill the condition. It makes things faster if you have another if after, because if you already satisifed the condition, you do not need to do more stuff there.
                            #https://www.tecmint.com/use-next-command-with-awk-in-linux/
                            #https://www.biostars.org/p/304979/
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
                    #CHECK THAT CHROMSOME AND POS ARE THE FIRST AND SECOND COLUMNS WHEN USING THIS ON REAL DATA
                #--begin: column number for region start
                    #position of the SNP in our case (second column)
                    #I am using the coordinates from the VCF file, so they are 1-based, so I have not to use "--zero-based".
                #--end: column number for region end (if no end, set INT to -b). The end column can be the same as the start column.
                    #again the SNP position in our case (second column)
                #--force: overwrite existing index without asking
                #--comment: skip comment lines starting with CHAR
                    #in our case we use "#" to comment the first line with the header
                #http://www.htslib.org/doc/tabix.html

print_text("take the indexed and tab-delimited file with the ancestral alleles in upper case and the position to create a new field with upper ancestral alleles, save the new VCF file", header=4)
run_bash(" \
    cd ./data/dummy_vcf_files/00_dummy_vcf_files_vep/; \
    bcftools +split-vep \
        ./dummy_example_vep_2.vcf.gz \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ \
        --annotations ./anc_alleles_uppercase.tsv.gz \
        --columns CHROM,POS,.AA_upcase \
        --header-line '##INFO=<ID=AA_upcase,Number=.,Type=String,Description=\"The AA field from INFO/CSQ after converting alleles to uppercase\">' > ./dummy_example_vep_2_anc_up.vcf; \
    cat ./dummy_example_vep_2_anc_up.vcf")
        #From the CSQ field added by VEP, extract the tag "AA" as a string, which is the ancestral state.
        #remove the CSQ field
        #add a new INFO/Tag using the tab delimited file previously created
        #select the columns from the tab file in which we are interested
            #We do ".AA" because we want to include also missing values (i.e., 'AA=.')
                #.TAG 
                    #Add TAG even if the source value is missing. This can overwrite non-missing values with a missing value and can create empty VCF fields (TAG=.)
                #TAG
                    #Add TAG if the source value is not missing (“.”). If TAG exists in the target file, it will be overwritten
                #+TAG
                    #Add TAG if the source value is not missing and TAG is not present in the target file.
                #.+TAG
                    #Add TAG even if the source value is missing but only if TAG does not exist in the target file; existing tags will not be overwritten.
                #https://samtools.github.io/bcftools/howtos/annotate.html
        #add the header line for this new tag
        #save as a new file
print("We can how the new tag AA_upcase has ACTG in uppercase, while the rest of characters ('.', '-', 'N') remain the same")

print_text("check that the new INFO/TAG with ancestral alleles in upper case is exactly the same than the original AA tag but in uppercase always", header=4)
print("calculate the number of SNPs for which AA is just AA_upcase but in lowercase")
count_aa = run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    if(toupper($(i+1))==$(i+3)){ \
                        count++; \
                        next \
                    } \
                } \
            } \
        }; \
        END{ \
            print count\
            }'", return_value=True).strip()
print("then check that this number is equal to the total number of SNPs")
check_aa = run_bash(" \
    n_variants=$( \
        bcftools view \
            ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
            --drop-genotypes \
            --no-header | \
        awk \
            'END{print NR}'); \
    if [[ " + count_aa +  " -eq $n_variants ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi", return_value=True).strip()
    #load the VCF file without the header and the genotypes, having one row per SNP
    #count the number of rows at the end of the file with awk
    #if the total number of rows (i.e., SNPs) is equal to the number of SNPs for which AA is just AA_upcase in lowercase, then we are good, all SNPs are ok.
if (check_aa == "TRUE"):
    print("GOOD TO GO! The new tag with ancestral alleles is exactly the same than the original AA field but in upper case")
else:
    raise ValueError("SERIOUS ERROR! We have not correctly converted to upper case all ancestral alleles")

print_text("check that REF and ALT are always in upper case", header=4)
print("calculate the number of SNPs for which the REF and ALT alleles are in uppercase")
count_ref_alt_upper = run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            if(toupper($4)==$4 && toupper($5)==$5){ \
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
print("then check that this number is equal to the total number of SNPs")
check_ref_alt_upper = run_bash(" \
    n_variants=$( \
        bcftools view \
            ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
            --drop-genotypes \
            --no-header | \
        awk \
            'END{print NR}'); \
    if [[ " + count_ref_alt_upper +  " -eq $n_variants ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi", return_value=True).strip()
    #load the VCF file without the header and the genotypes, having one row per SNP
    #count the number of rows at the end of the file with awk
    #if the total number of rows (i.e., SNPs) is equal to the number of SNPs for which AA is just AA_upcase in lowercase, then we are good, all SNPs are ok.
if (check_ref_alt_upper == "TRUE"):
    print("GOOD TO GO! The REF and ALT alleles are always in uppercase")
else:
    raise ValueError("SERIOUS ERROR! We do not have the all REF and ALT alleles in upper case so we cannot correctly select those SNPs whose REF is not AA, because our ancestral allele data is always upper case and bcftools filters are case sensitive")

print_text("Use the new AA_upcase tag to filter and check the behavior", header=4)
print("SNPs where REF is ancestral")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("SNPs where REF is NOT ancestral but ALT is")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF!=AA_upcase && ALT=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("SNPs where REF nor ALT is ancestral")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF!=AA_upcase && ALT!=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("We get as cases with AA not being the ALT those multiallelics where 1 of the ALT alleles are indeed the ancestral. For example, REF=A;ALT=T,G;AA=G is considered when filtering by ALT!=AA_upcase, when indeed G is the ancestral and it is one of the alternative alleles. Oddly enough, these SNPs are again filtered in when doing ALT=AA_upcase. In both cases makes sense, because we have ALT alleles that are the AA but other are not. For things like this, we should do ancestral filtering after dealing with multiallelic snps.")

print_text("Check the number of SNPs without ancestral allele", header=4)
run_bash(" \
    n_snps_missing_ancestral=$( \
        bcftools view \
            --drop-genotypes \
            --no-header \
            --include 'AA_upcase=\".,-,N\"' \
            ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
        awk 'END{print NR}'); \
    n_snps_total=$( \
        bcftools view \
            --drop-genotypes \
            --no-header \
            ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf | \
        awk 'END{print NR}'); \
    printf 'We have a total of %s SNPs\n' \"$n_snps_total\"; \
    printf 'We have %s SNPs WITHOUT ancestral alleles\n' \"$n_snps_missing_ancestral\"; \
    loss_percent=$( \
        awk \
            -v x=$n_snps_missing_ancestral \
            -v y=$n_snps_total \
            'BEGIN{print (x/y)*100}'); \
    if [[ $loss_percent < 40 ]]; then \
        printf 'GOOD TO GO! The percentage of SNPs without ancestral allele is: %s' \"$loss_percent\"; \
    else \
        echo 'ERROR: FALSE! There are more than 40% of SNPs without ancestral allele';\
    fi")
    #Count the number of rows without header (i.e., SNPs) for which the ancestral allele is ".", "-" or "N". We do not use --exclude "A,C,T,G" because this lead to include indels with ancestral state like "TT" because "TT" is not "T". Use awk to count with print NR at the end of the file.
        #from the fasta file info
            #ACTG: high-confidence call, ancestral state supported by the other two sequences
            #actg: low-confidence call, ancestral state supported by one sequence only
                #these are no longer present in AA_upcase, they have been converted to uppercase
            #N: failure, the ancestral state is not supported by any other sequence
            #-: the extant species contains an insertion at this postion
            #.: no coverage in the alignment
    #Count the total number of SNPs, i.e., no filtering.
    #load both numbers into AWK and divide SNPs with missing ancestrals by the total number of SNPs, then multiply by 100 to get percentage
    #check this number is not very high
        #we can use < inside [[ ]] as it will not be considered redirection
        #-gt/-lt does not work for floats



print_text("See the dummy VCF file after all operations with the Ancestral Alleles as a new field with all bases in uppercase to match that of REF and ALT columns, so we will be able to apply filters in next steps", header=2)
run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
        --drop-genotypes")




################################################################
#### function to estimate the ancestral allele for ALL SNPs ####
################################################################
print_text("function to estimate the ancestral allele for ALL SNPs", header=1)
#selected_chromosome="1"; debugging=True
def master_processor(selected_chromosome, debugging=False):

    #redirect standard output ONLY when running for production
    if debugging == False:
        import sys
        original_stdout = sys.stdout
            #save off a reference to sys.stdout so we can restore it at the end of the function with "sys.stdout = original_stdout"
        sys.stdout = open("./scripts/00_ancestral_calcs_outputs/chr" + selected_chromosome + ".out", "w")
            #https://www.blog.pythonlibrary.org/2016/06/16/python-101-redirecting-stdout/



    print_text("Initial operations", header=2)        
    print_text("chr" + selected_chromosome + ": see VCF file version", header=3)
    vcf_version = run_bash(" \
        bcftools head \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
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
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        head -10")
        #https://samtools.github.io/bcftools/howtos/query.html


    print_text("chr " + selected_chromosome + ": the number of samples is equal to the number of samples considering unrelated individuals and trios-duos?", header=3)
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


    print_text("chr " + selected_chromosome + ": the chromosome name of the first variant is correct?", header=3)
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


    print_text("chr " + selected_chromosome + ": show the variant type, ID, chromosome, position, alleles and frequency for the first snps")
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF\n' \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            head -5")
        #select the format of the query indicating the columns you want to show per SNP.
            #you can include data from INFO
            #end with \n to have different lines per SNPs

    print_text("chr " + selected_chromosome + ": check that all variants has '.' for filter, but only if you are NOT in debugging mode because it takes a few minutes", header=3)
    if debugging==False:
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


    print_text("chr " + selected_chromosome + ": see genotypes of first samples", header=3)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF GTs:[ %GT]\n' \
            " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
        head -2")



    print_text("run VEP's plugin ancestral allele on the VCF file", header=2)
    print_text("Select the input VCF for VEP. Use only a subset of the data if we are on debugging mode. Also select only SNPs in both cases, as we do not need indels and we can remove 1/6 of variants", header=3)
        #number of variants
            #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf
    print_text("apply the SNP filter always, but in the case of debugging=True reduce the sample size more", header=4)
    if debugging==True:
        run_bash(" \
            bcftools view \
                --types snps \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
            head -n 5000 > 00_ancestral_debug_subset_chr" + selected_chromosome + ".vcf; \
            ls -l")
        input_vcf_file_vep="00_ancestral_debug_subset_chr" + selected_chromosome + ".vcf"
    else:
        run_bash(" \
            bcftools view \
                --types snps \
                --output-type z \
                --output ./1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vcf.gz \
                " + input_vcfs_path + "/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_INDEL_SV_phased_panel.vcf.gz; \
            ls -l")
        input_vcf_file_vep="./1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vcf.gz"
    print(input_vcf_file_vep)
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
            " + input_vcf_file_vep + " | \
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


    print_text("make a run of VEP", header=3)
    run_bash("\
        vep \
            --offline \
            --verbose \
            --species 'homo_sapiens' \
            --assembly GRCh38 \
            --input_file " + input_vcf_file_vep + " \
            --format vcf \
            --output_file ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz\
            --vcf \
            --compress_output gzip \
            --force_overwrite \
            --fields Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,STRAND,AA \
            --cache \
            --dir_cache ./data/vep_cache/cache_" + cache_version + " \
            --plugin AncestralAllele,./data/fasta_ancestral/homo_sapiens_ancestor_GRCh38_final.fa.gz \
            --dir_plugins /opt/ensembl-vep/vep_plugins")
        #see dummy example for details about the arguments


    print_text("explore the generated VCF file and do checks", header=3)
    print_text("view the file", header=4)
    run_bash(" \
        bcftools view \
            --header \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz")

    print_text("do some checks about about whether we have used the VEP/ensemble/cache/assembly versions using the row added to the header in the VCF file by VEP", header=4)
    line_checks_after = run_bash(" \
        gunzip \
            --stdout \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        grep '##VEP='", return_value=True)
    #extract the row of the header with that information
    line_checks_after_split = line_checks_after.split(" ")
        #split the different fields
    import numpy as np
    #extract the different versions
    vep_version_vcf = line_checks_after_split[np.where(['##VEP=' in field for field in line_checks_after_split])[0][0]].replace('##VEP="v', '').replace('"', '')
    ensembl_version_vcf = line_checks_after_split[np.where(['ensembl=' in field for field in line_checks_after_split])[0][0]].split(".")[0].replace("ensembl=", "")
    assembly_version_vcf = line_checks_after_split[np.where(['assembly=' in field for field in line_checks_after_split])[0][0]].split(".")[0].replace('assembly="', '')
        #select the field of interest, then split and/or replace to correctly extract the version number
    import re
    cache_version_vcf = re.split("_|/", line_checks_after_split[np.where(['cache=' in field for field in line_checks_after_split])[0][0]])[-2]
    cache_version_assembly_vcf = re.split("_|/", line_checks_after_split[np.where(['cache=' in field for field in line_checks_after_split])[0][0]])[-1].replace('"', "")
        #in the case of the cache line, we need to split using two separators to obtain both the cache and assembly version
    print("The vep version in the VCF file is the same than the one I selected in the script?")
    print(vep_version_vcf == vep_version)
    print("The cache version in the VCF (according to cache and ensemble lines) matches the version I selected in the script?")
    print(cache_version_vcf == cache_version)
    print(ensembl_version_vcf == cache_version)
    print("The cache version of the VCF file is the same than the VEP version used according to the VCF file?")
    print(cache_version_vcf == vep_version_vcf)
    print("The cache version of the VCF file is the same than the VEP version used according to my script?")
    print(cache_version_vcf == vep_version)
    print("The assembly version according to the cache and assembly lines in the VCF matches?")
    print(cache_version_assembly_vcf == assembly_version_vcf)
    print("The assembly is GRCh38?")
    print(assembly_version_vcf == 'GRCh38')
    print("The cache version of the VCF file according to the ensemble line is equal to the vep version according to the VCf file and according to my script?")
    print(ensembl_version_vcf == vep_version_vcf)
    print(ensembl_version_vcf == vep_version)

    print_text("see CSQ field for the first variants", header=4)
    run_bash(" \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AF CSQ: %CSQ\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        head -n 1000")

    print_text("According to the VEP manual, vep does not change the vcf file, only add the CSQ field. Therefore, we could process the original vcf files and then process with our cleaning. Let's check if the vep VCF is exactly the same after we removed the CSQ field, which is the one added by AncestralAllele plugin of VEP. Only do it if NOT in debugging mode because we need the whole dataset to this check to work", header=4)
    if debugging==False:
        run_bash(" \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO %FORMAT\n' \
                " + input_vcf_file_vep + " > chr" + selected_chromosome + "_file_check_1.txt; \
            bcftools annotate \
                --remove INFO/CSQ \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO %FORMAT\n' > chr" + selected_chromosome + "_file_check_2.txt; \
            echo '##See head first check file:'; head -n 5 chr" + selected_chromosome + "_file_check_1.txt; \
            echo '##See head second check file:'; head -n 5 chr" + selected_chromosome + "_file_check_2.txt; \
            check_status=$(cmp --silent chr" + selected_chromosome + "_file_check_1.txt chr" + selected_chromosome + "_file_check_2.txt; echo $?); \
            echo '##Do the check'; \
            if [[ $check_status -eq 0 ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi; \
            rm chr" + selected_chromosome + "_file_check_1.txt; rm chr" + selected_chromosome + "_file_check_2.txt; \
            ls -l")
    
    
    print_text("Create a new AA field using the ancestral allele stored in CSQ/AA", header=3)
    print_text("see again the header of the VCF file after VEP processing", header=4)
    run_bash("\
        bcftools head \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz")
    
    print_text("see all the tags in the CSQ field", header=4)
    run_bash("\
        bcftools +split-vep \
            --list \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz")

    print_text("make tags inside CSQ available and then extract AA as a new column outside CSQ. You can then make a query and call AA without +split-vep", header=4)
    run_bash("\
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 1000")


    print_text("why I am directly calling AA using +split-vep?", header=4)
    print("If you query %AA using  +split-vep and --annotation CSQ, you would not get variants where there is empty AA, i.e., '', you would get those with A='.', but not with AA=''")
    print("In the script, if you run the next lines for chromosome 1, you will see how we have more variants if we call AA as an independent INFO/AA field with bcftools query than if we just call it inside CSQ with +split-vep. 1:10398:C:CCCCTAA has AA='', i.e., there is no data, no space for AA. This is the variant lost when directly making the query with +split-vep. If we use +split-vep but do not call AA in the query, we DO get the missing variant. It seems that bcftools add '.' to these empty cases when querying outside of CSQ as shown for 1:10398:C:CCCCTAA. I prefer to get all variants and then deal with those without data, so we will use the first approach")
    print("We are not going to run this because we do not longer have INDELS in the data, we have previously filtered out. I leave the code as a reminder of why I am using AA from INFO.")
    if False:
        run_bash(" \
            bcftools +split-vep \
                --annotation CSQ \
                --columns AA:String \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
            bcftools query \
                --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
            head -n 4")
        run_bash(" \
            bcftools +split-vep \
                --annotation CSQ \
                --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
            head -n 3")
        run_bash(" \
            bcftools +split-vep \
                --annotation CSQ \
                --format '%TYPE %ID %CHROM %POS %REF %ALT\n' \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
            head -n 3")

    print_text("use AA to filter", header=4)
    print("exclude those SNPs for which the REF allele IS NOT the ancestral. Note that 'A' and 'a' are considered different, this is something we will deal later")
    run_bash("\
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools query \
            --exclude 'REF=AA' \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 70")
    print("exclude those SNPs for which the REF allele IS the ancestral")
    run_bash("\
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools query \
            --include 'REF=AA' \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 70")
    print("As you can see, we can just filter by the ancestral allele using AA and it does not consider C, C, C.... but only C. Also note that C is not considered equal to 'c', see below")

    print_text("extract AA tag from CSQ and then remove CSQ", header=4)
    run_bash("\
        bcftools +split-vep \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz \
            --annotation CSQ \
            --columns AA:String | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 50")


    print_text("IMPORTANT CHECK: Check that all the consequences the same variant has across transcripts have ALWAYS the same ancestral allele", header=3)
    print_text("First, extract the AA field, where we have 1 row per SNP and each consequence across transcripts is separated by ',', these are going to by our columns", header=4)
    run_bash(" \
        bcftools +split-vep \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz \
            --annotation CSQ \
            --columns AA:String | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools query \
            --format '%AA\n' > check_ancestral_transcripts_chr" + selected_chromosome + ".csv; \
        head -n 50 check_ancestral_transcripts_chr" + selected_chromosome + ".csv")
    
    print_text("get the number of rows of the file", header=4)
    total_ancestral_check_before = run_bash(" \
        awk \
            'END{print NR}'\
            check_ancestral_transcripts_chr" + selected_chromosome + ".csv", return_value=True).strip()
    
    print_text("Then add a new lines at the beginning with two different alleles to check our approach is able to discard these cases. Also add two lines with the same alleles, i.e., two spaces and two empty fields, which should be considered as equal", header=4)
    run_bash(" \
        sed \
            --in-place \
            --expression '1i A,a' \
            --expression '1i A,1' \
            --expression '1i A,!' \
            --expression '1i \\ ,A' \
            --expression '1i ,A' \
            --expression '1i \\ ,' \
            --expression '1i \\ , ' \
            --expression '1i ,' \
            check_ancestral_transcripts_chr" + selected_chromosome + ".csv; \
        head -n 50 check_ancestral_transcripts_chr" + selected_chromosome + ".csv")
        #Change directly in the file (in place), and add 
            #"A,a": 
                #So we can be sure awk differentiate even the same letter in different case
            #"A,1"
                #to check awk can differentiate letters from numbers
            #"A,!"
                #to check awk can differentiate letters from symbols
            #" ,A": 
                #Space and then A (we have to add "\" in sed in order to add spaces in the expression) to check letters are differentiated from spaces
            #",A"
                #empty field and then A to check letters are differentiated from blank fields
            #" ,"
                #space and empty to check that both are clearly differentiated
            #" , "
                #two spaces to check whether two spaces are considered equal
            #","
                #two empty fields to check that two different fields that are empty are considered as identical
        #https://unix.stackexchange.com/a/99351
        #https://stackoverflow.com/a/43420225/12772630
        #https://stackoverflow.com/a/39788362/12772630
    
    print_text("Calculate the number of SNPs for which there is only 1 unique ancestral allele. In some rows, we will have more consequences (columns) and other will have less consequences and hence less columns. We will deal with that using awk", header=4)
    n_pass_ancestral_check = run_bash("\
        awk \
            'BEGIN{ \
                FS=\",\" \
            }{ \
                if(NF != 1){ \
                    for(i=2; i<=NF; i++){ \
                        if($(i-1)==$i){count_1++} \
                    };\
                    if((NF-1)==count_1){count_2++}; \
                    count_1=\"\" \
                }else if(NF == 1){ \
                    count_2++ \
                }; \
            }END{ \
                print count_2++ \
            }' \
            check_ancestral_transcripts_chr" + selected_chromosome + ".csv", return_value=True).strip()
        #load the file using the commas as field separator
        #if the row has more than 1 column:
            #remember that each row can have a different number of columns because each SNP can have more or less consequences. Therefore our loop across columns only work only for those rows with more than 1 column
            #from the field 2 to the last one
                #if the previous field (i-1) and the current field are identical, then add 1 to the count_1
                #this is able to detect as different an space " " in the field 1 vs an blank value "" in the next field or upper vs lower case (see new lines added to the csv).
                #in previous versions,
                    #I considered as OK those cases where $i or $(i-1) were blank, for example, "A, ,".
                    #I did that but adding an additional if else, so if both fields were NOT blank, I do the if previously explained, but if one of them is blank, no check was done and I gave it for good. 
                        #if($(i-1)!=\"\" && $i!=\"\"){ \
                            #if($(i-1)==$i){count_1++} \
                        #}else{count_1++} \
                    #but this is not ok because "A" and "" are not the same. I would expect that all consequences have EXACTLY the same ancestral allele. So we should get TRUE for this check without this line of code
            #then check whether the count_1 across fields is equal to the number of fields minus 1, e.g., we have 3 columns, you compare 1 with 2 and 2 with 3. If all are equal, we would have 2 counts (3 fields - 1). If True, then add 1 to count_2
            #then set the count_1 to "", so it is empty when we repeat the process with the next row
        #if not, and the number of fields is only 1
            #just add 1 to count_2
                #if there is only one consequence, i.e., one ancestral allele, there is really no need to compare. We know for sure there is only 1 unique ancestral allele, which is the main goal of this script.
        #return count_2 at the END, which is the total number of rows, i.e., SNPs, for which there is only 1 unique ancestral allele.
        #https://stackoverflow.com/a/57984015/12772630
    
    print_text("now check that the number of SNPs with 1 unique ancestral allele is exactly the total number of variantes before adding the new lines PLUS 2. Remember we added several lines with different ancestral alleles in each one. These new lines should NOT be counted by our approach, except the last 2, which included two spaces and two empty fields, respectively, and hence are identical. Therefore, we should have the original number of variants plus 2", header=4)
    if (int(total_ancestral_check_before)+2) == int(n_pass_ancestral_check):
        print("YES! GOOD TO GO! All variants have 1 unique ancestral allele, so we can just select 1 per variant")
    else:
        raise ValueError("SERIOUS ERROR! THE CHECK FOR 1 UNIQUE ANCESTRAL ALLELE PER VARIANT FAILED. WE NEED TO CAREFULLY CHECK THIS BECAUSE WE CAN ONLY SELECT ONE ANCESTRAL ALLELE PER VARIANT")
    print("remove the csv files used to do the check")
    run_bash(" \
        rm check_ancestral_transcripts_chr" + selected_chromosome + ".csv; \
        ls -l")


    print_text("make more checks on the vep generated file", header=3)
    print_text("there is any problem with multiallelic variants?", header=4)
    run_bash(" \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --min-alleles 3 | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %INFO/AA\n' | \
        head -n 50")
    print("As you can see in these few cases, the different alleles of the same variant get the same ancestral allele. The same goes for the rest of CSQ fields. This makes sense because the different lines of a multiallelic variant have the same position, thus they get the same ancestral allele.")

    print_text("select SNPs for which the ancestral allele is ACGT or acgt, avoiding cases where AA='.' or '-'", header=4)
    run_bash("\
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools query\
            --include 'AA=\"A,C,G,T,a,c,g,t\"'\
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 70")
            #we get rid of the SNP with AA=".", but also of AA="", i.e., empty. The variant 1:10398:C:CCCCTAA in chromosome 1 has no string for AA and it is removed from here.
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

    print_text("WE HAVE TO EXCLUDE CASES WHERE REF NOR ALT ARE THE ANCESTRAL ALLELE. This can be the case for a multiallelic SNP that only has two alleles in the selected populations, being the missing one the ancestral. After multiallelic snps have been removed, we have to exclude those for which AA is not REF nor ALT. In this example, we are going to specifically select these cases to see how bcftools deals with letters with different case. Note that the VCF files of the 1KGP already have multiallelix snps separated, so we not need to do --multiallelic -snps", header=4)
    run_bash("\
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2 \
            --min-alleles 2 | \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools view \
            --include 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT!=AA' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 50")
            #remove those monomorphic SNPs. This will remove REF/ALT pairs of a multiallelic SNP if one of the ALTs is not present in the population.
            #bind the multiallelic SNPs and remove those that still have more than 1 ALT, because these are truly multiallelic SNPs in the dummy population.
            #Extract AA data from CSQ and create a new, independent field with that.
        #select SNPs where AA is not REF nor ALT. These can be cases where one line of the multiallelic has been removed (e.g., because it was all REF, i.e., monomorphic) and the ALT allele of that line was actually the AA. After the removal the remaining line has REF and ALT, being none of them the AA.
            #this filter removes problematic cases where the AA was a allele considered as ALT, maybe because was in low frequency as it was surpassed by other alleles.
        #then select only for SNPs with actual ancestral allele estimated, instead of "." or "N"
        #I use && instead of & because only one & looks for variants satisfying the condition within sample.
            #For example, say our VCF contains the per-sample depth and genotype quality annotations per SNP and we want to include only sites where ONE OR MORE samples have big enough coverage (DP>10) and genotype quality (GQ>20). The expression -i 'FMT/DP>10 & FMT/GQ>20'. This would select variants for which at least one sample has good coverage and genotype quality.  
            #but we do not want this. We want to select SNPs that fulfill certain requirements across all samples, not just within at least one sample. We want to match the whole record, using features that are similar across samples, like REF, ALT and AA.
                #http://samtools.github.io/bcftools/howtos/filtering.html
    print("IMPORTANT: You can see many cases for which the AA is like the REF but in lower case. These cases should not be included here because REF does not meet REF!=AA. There is a problem with the case here in bcftools.")

    print_text("select now those cases where REF!=AA but ALT=AA, because these are the cases where REF and ALT can be switched", header=4)
    run_bash("\
        bcftools view \
            --exclude 'INFO/AC=INFO/AN || INFO/AC=0' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools norm \
            --multiallelic +snps | \
        bcftools view \
            --max-alleles 2 \
            --min-alleles 2 | \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools view \
            --exclude 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT!=AA' | \
        bcftools view \
            --include 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT=AA' | \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
        head -n 20")


    print_text("Problem upper vs. lower case ancestral allele", header=3)
    print("We have a problem with case sensitive. We have previously seen that if REF=C and AA=c, REF is NOT considered to be equal to AA. Therefore, if ALT=C and AA=c, ALT is NOT considered to be equal to AA. So I am creating a new AA field with all alleles as upper case, so we avoid this problem and we can consider ancestral alleles with high and low confidence (upper and lower case, respectively)")
    
    print_text("extract the position (index) of the columns CHROM and POS in the VCF file, so we can be sure we are selecting these two columns. We need these columns to create a file with chrom, pos and upper ancestral allele that will be included as a new field in the VCF file using tabix (see below)", header=4)
    print("Important: From this step, we are going to remove non-SNPs. We can have a SNP and a INDEL in the same position. This generates errors when the index (position and ancestral allele) tries to add the Ancestral allele of a given variant. For example, it can give the AA to the SNP and leave the INDEL without AA or viceversa, even though both the SNP and the INDEL have AA. The problem is that they are in the same position.")
    print("obtain the position of these columns using awk")
    indexes_chrom_pos = run_bash(" \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools view \
            --types snps | \
        bcftools view \
            --header | \
        awk \
            'BEGIN{FS=\"\t\"}; \
            END{ \
                for(i=1;i<=NF;i++){ \
                    if($i == \"#CHROM\"){ \
                        chrom_index=i \
                    }; \
                    if($i == \"POS\"){ \
                        pos_index=i \
                    }; \
                    if($i ~/^HG/ || $i ~/^NA/){ \
                        n_samples++\
                    } \
                }; \
                n_fields=NF; \
                printf \"n_fields=%s,n_samples=%s,chrom=%s,pos=%s\", n_fields, n_samples, chrom_index, pos_index \
            }'", return_value=True).strip()
        #get the header of the VCF file after extracting AA from CSQ, removing CSQ and select SNPs, just like we are going to do when we replace lower for upper case in the next line
        #open the header with AWK
            #when you reach the last line, which includes the headers
                #run loop across fields, i.e., headers
                    #if the header is CHROM or POS then save the index of the field as a new variable, chrom_index and pos_index, respectively.
                        #"i" is the index, like 1, 2, 3... because of this, when you do $i is like you are doing $1.
                        #if we save "i", we are saving the index, the number
                        #https://unix.stackexchange.com/a/616495
                    #if the header starts with HG or NA, add 1 to the count of n_samples, because this is a GT column for a given sample
                        #"~" let you use regular expression
                        #"/.../" is a regular expression to match text that meet condition
                        #"^" text that starts with...
                        #https://unix.stackexchange.com/a/72763
                #save the number of fields
            #then print the number of fields, the number of sample and the index of both CHROM and POS
                #https://www.gnu.org/software/gawk/manual/html_node/Printf-Examples.html
    print("extract the numbers")
    n_fields = indexes_chrom_pos.split(",")[0].replace("n_fields=", "")
    n_samples = indexes_chrom_pos.split(",")[1].replace("n_samples=", "")
    index_chrom = indexes_chrom_pos.split(",")[2].replace("chrom=", "")
    index_pos = indexes_chrom_pos.split(",")[3].replace("pos=", "")
    print("total number of fields: " + n_fields)
    print("number of samples: " + n_samples)
    print("index of column CHROM: " + index_chrom)
    print("index of column POS: " + index_pos)
    print("the total number of fields minus the number of samples should be 9. The number of fixed fields in VCF v4.2 is 8 and then FORMAT, which in our case only has GT, thus we should have 9 fields. Also, the index of CHROM and POS should be 1 and 2, respectively")
    if (int(n_fields)-int(n_samples) == 9) & (index_chrom=="1") & (index_pos=="2"):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WITH THE FIELDS OF THE VCF FILE BEFORE CONVERTING TO UPPER CASE ANCESTRAL ALLELES")

    print_text("create a tab separated file with the position info and ancestral allele in upper case of each SNP", header=4)
    run_bash(" \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools view \
            --types snps | \
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
                        print $" + index_chrom + ", $" + index_pos + ", $(i+1); \
                        next \
                    } \
                } \
            }' \
        > ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv; \
        sed  \
            --in-place \
            --expression '1i #CHROM\tPOS\tAA_upcase' \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv; \
        bgzip \
            --force \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv; \
        echo 'See first lines of the tab file with upper case ancestral alleles'; \
        gunzip \
            --stdout \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz | head -n 70")
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
                    #if the field is the ancestral allele (AA)
                        #update the next field (value of AA) with the same string that was present but in upper case.
                        #then print the fields 1 (chrom), 2 (pos), i+1 (AA value)
                            #this is the format required by bcftools annotate to create a new field based on a tab delimited file
                                #--annotations: VCF file or tabix-indexed FILE with annotations: CHR\tPOS[\tVALUE]
                        #go to next row
                            #next is used to go to the next line once you have fulfilled the condition. It makes things faster if you have another if after, because if you already satisifed the condition, you do not need to do more stuff there.
                                #https://www.tecmint.com/use-next-command-with-awk-in-linux/
                                #https://www.biostars.org/p/304979/
            #save the result as a file
            #add a header to that file (tab separated names) but annotating that new line with "#" to avoid problems with tabix
                #https://unix.stackexchange.com/a/401673
            #compress using bgzip (from samtools) so the file is recognised by tabix, see below
                #https://github.com/samtools/bcftools/issues/668

    print_text("check whether the position in the new file is position sorted, as this is a requirement of tabix", header=4)
    print("extract the second column (POS) avoiding the header (NR=1)")
    list_pos = run_bash(" \
        gunzip \
            --stdout \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz | \
        awk \
            'BEGIN{FS=\"\t\"}{ \
                if(NR>1){print $2} \
            }'", return_value=True).split()
    print("convert to integer each position because we get it in the form of string from run_bash")
    list_pos = [int(pos) for pos in list_pos]
    print("sort in ascending order") 
    list_pos_sorted = sorted(list_pos)
        #sorted([1,2,3,5,4]) gives [1, 2, 3, 4, 5]
    print("check that both lists are the same and hence the original was already sorted")
    if list_pos==list_pos_sorted:
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! Coordinates are not sorted in ascending order")

    print_text("check also that the first column of chrom only have one unique chromosome", header=4)
    print("extract the unique chromosomes avoiding the header (NR=1)")
    uniq_chrom = run_bash(" \
        gunzip \
            --stdout \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz | \
        awk \
            'BEGIN{FS=\"\t\"}{ \
                if(NR>1 && !a[$1]++){print $1} \
            }'", return_value=True).strip()
        #print the unique chromosomes in the file avoiding the first row with the headers (i.e., NR=1)
            #a[$1] is the name of the array that holds $1 (chromosome) as keys.
            #uses the current value of $1 as key to the array "a", taking the value stored there. If this particular key was never referenced before, a[$2] evaluates to the empty string.
            #In "!a[$1]", the ! negates the value from before. If it was empty or zero (i.e., the key was never referenced before; false), we now have a true result. If it was non-zero (i.e., the key was referenced before; true in the original test), we have a false result. If the whole expression evaluated to true, meaning that a[$1] was not set to begin with, the whole line is printed as the default action. Therefore, it is only printing a value of $1 if was never references before, i.e., non-duplicate.
            #In other words:
                #a[$1]: look at the value of key $1, in associative array a. If it does not exist, automatically create it with an empty string.
                #!a[$1]++: negate the value of expression. If a[$1]++ returned 0 (a false value), the whole expression evaluates to true, and makes awk perform the default action print $1. Otherwise, if the whole expression evaluates to false, no further action is taken.
    print("check we have only one value of CHROM in the form of string and it is the correct one")
    if (type(uniq_chrom)==str) & (uniq_chrom=="chr"+selected_chromosome):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE MORE THAN 1 CHROMOSOME")

    print_text("check duplicated positions have the same AA because these are the different lines of the same multiallelic snp", header=4)
    print("calculate the number of unique AA_upcase values per POS value in the file that will be used to include AA_upcase into the VCF file")
    run_bash(" \
        gunzip \
            --stdout \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz | \
        awk \
            'BEGIN{FS=OFS=\"\t\"}{ \
                if(NR>1){ \
                    if(!found[$0]++){ \
                        val[$2]++ \
                    } \
                } \
            }END{ \
                print \"pos\tuniq_count\"; \
                for(i in val){ \
                    print i,val[i] \
                } \
            }' > check_unique_aa_chr" + selected_chromosome + ".tsv; \
            gzip \
                --force \
                check_unique_aa_chr" + selected_chromosome + ".tsv; \
            ls -l")
            #avoid the first row with the header, we only want to process actual positions and AA_upcase
            #found[$0]++
                #creates an array with the whole row ($0) as index and then sum 1 to the previous value of that index. We have two possibilities:
                    #the row (that specific combination of CHROM, POS and AA_upcase) has NOT been never included before.
                        #The value for that index is 0, so 0+1 gives 1
                        #Also, the expression gives a False, because found[$0] does NOT previously exists
                    #the row (that specific combination of CHROM, POS and AA_upcase) has been included before.
                        #If this combination was seen 1 time before, the value for that index is 1, so 1+1 gives 2
                        #Also, the expression gives True, because found[$0] does exists.
                    #you could also use specific columns to use them as index instead of the whole row (e.g., [$1, $2]). We use the whole row because we want to detect unique combinations of CHROM, POS and AA_upcase
            #!found[$0]++
                #"!found[$0]++" does the same than "found[$0]++", i,e., creates an array with the whole row as index and 1 as value if it has not been created before, or if it has been, just add 1 to the current value of that index.
                #the difference is that we negate, so if the row was not previously present as index in the array, you get now True, while you get False if it was previously included as an index.
                #using this, we are detecting two different cases
                    #CHROM-POS combination is present two times and have twoo different AA_upcase values. This makes that the same POS is associated with two different rows, thus it will be have a value of 2 in the "val" array instead of 1.
                    #POS-AA_upcase combination is present two times and have two different CHROM values. This makes that the same POS is associated with two different rows, thus it will be have a value of 2 in the "val" array instead of 1.
                    #we can detect two errors, as we should have only 1 chromosome and 1 AA_upcase value per POS.
                #therefore, you can do operations IF this is the first time the current row has appeared (see below).
            #if "!found[$0]++" gives True, this means that the row (the combination of CHROM, POS and AA_upcase) has NOT been previously included in "found" and the POS value of that row then has a new AA_upcase. If that is the case:
                #create an entry in the "val" array using "POS" ($2) as index and adding 1.
                #if this is the first time this POS appeared, its value will be 0+1=1, but if it is the second, then its val will be 1+1=2.
                #this will create an array with POS as index and then value as the number of times that POS appeared with a different combination of CHROM and AA_upcase
            #At the END
                #print a first row as header
                #run a loop across each index of "val", i.e., across each POS, 
                    #print
                        #the index, i.e., the position, and 
                        #the value, i.e., the number of times that POS appeared with a different AA_upcase value, i.e., the number of unique AA_upcase values per position
                        #index and value and separated by "\t" as this is the output field separator (OFS)
            #I have checked that
                #printing "found" at the END, gives you a list of the rows ($0) used as index and a value. I have seen that several rows with value=2 are rows whose position is repeated, meaning that the whole row is duplicated (i.e., multiallelic SNP).
                #printing val after doing "found[$0]++" instead of "!found[$0]++" gives a list of positions and some of them have 2 as value. These seem to be cases where the position is 3 times present. The first time is False for "found[$0]++", because the row is not present yet, and then the other two times is true, because the row is already present, so we add 2 to the value of val for the select POS index.
            #save as tsv file and then compress it
            #script based on:
                #https://stackoverflow.com/a/62467589/12772630
    print("load the tsv file into pandas")
    check_unique_aa = pd.read_csv(
        "check_unique_aa_chr" + selected_chromosome + ".tsv.gz", 
        sep="\t", 
        header=0, 
        low_memory=False)
    print(check_unique_aa)
    print("check that the position column has no duplicates, i.e., we have results per distinct value of pos AND these distinct position values are the same than those in list_pos (after selecting unique)")  
    if (sum(check_unique_aa["pos"].duplicated(keep=False))==0) & (sorted(check_unique_aa["pos"].unique().tolist()) == sorted(set(list_pos))):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WHEN CHECKING THE NUMBER OF UNIQUE CASES OF AA_UPCASE")
    print("Do now the actual important check here, which is the fact that the number of unique AA_upcase values per distinct position should be 1. Meaning that if two rows belong to the same multiallelic SNP, then both will have the same position and the same AA_upcase value")
    check_unique_aa_final_count = check_unique_aa["uniq_count"].unique()
    if (len(check_unique_aa_final_count)==1) & (check_unique_aa_final_count==1):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("FALSE! ERROR! WE HAVE A PROBLEM WHEN CHECKING THE NUMBER OF UNIQUE CASES OF AA_UPCASE")
    print("remove the file used to do the check")
    run_bash(" \
        rm check_unique_aa_chr" + selected_chromosome + ".tsv.gz; \
        ls -l")

    print_text("Then create an index for this tab-delimited file using tabix.", header=4)
    run_bash(" \
        tabix \
            --sequence 1 \
            --begin 2 \
            --end 2 \
            --force \
            --comment '#' \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz; \
        echo 'See first lines of the tab file with upper case ancestral alleles after creating the index'; \
        gunzip \
            --keep \
            --stdout \
            ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz | head -n 70")
            #create index of the tab delimited file with SNP positions and upper ancestral alleles
                #Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
                #After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)
                #tabix
                    #--sequence: column number for sequence names
                        #chromosome in our case (first column)
                    #--begin: column number for region start
                        #position of the SNP in our case (second column)
                        #I am using the coordinates from the VCF file, so they are 1-based, so I have not to use "--zero-based".
                    #--end: column number for region end (if no end, set INT to -b). The end column can be the same as the start column.
                        #again the SNP position in our case (second column)
                    #--force: overwrite existing index without asking
                    #--comment: skip comment lines starting with CHAR
                        #in our case we use "#" to comment the first line with the header
                    #http://www.htslib.org/doc/tabix.html

    print_text("take the indexed and tab-delimited file with the ancestral alleles in upper case and the position to create a new field with upper ancestral alleles, save the new VCF file and then take a look at it", header=4)
    run_bash(" \
        bcftools +split-vep \
            --annotation CSQ \
            --columns AA:String \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
        bcftools annotate \
            --remove INFO/CSQ | \
        bcftools view \
            --types snps | \
        bcftools annotate \
            --annotations ./results/00_vep_vcf_files/anc_alleles_uppercase_chr" + selected_chromosome + ".tsv.gz \
            --columns CHROM,POS,.AA_upcase \
            --header-line '##INFO=<ID=AA_upcase,Number=.,Type=String,Description=\"The AA field from INFO/CSQ after converting alleles to uppercase\">' \
            --output-type z\
            --output ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz; \
        bcftools query \
            --format '%TYPE %ID %CHROM %POS %REF %ALT AA:%AA AA_upcase:%AA_upcase\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        head -n 100")
            #From the CSQ field added by VEP, extract the tag "AA" as a string, which is the ancestral state.
            #remove the CSQ field
            #add a new INFO/Tag using the tab delimited file previously created
            #list the the columns from the tab file that will be used in the annotation
                #We do ".AA" because we want to include also missing values (i.e., 'AA=.') and overwrite existing tags with the same name as I am using a custom name for uppercase ancestral alleles
                    #.TAG 
                        #Add TAG even if the source value is missing. This can overwrite non-missing values with a missing value and can create empty VCF fields (TAG=.)
                    #TAG
                        #Add TAG if the source value is not missing (“.”). If TAG exists in the target file, it will be overwritten
                    #+TAG
                        #Add TAG if the source value is not missing and TAG is not present in the target file.
                    #.+TAG
                        #Add TAG even if the source value is missing but only if TAG does not exist in the target file; existing tags will not be overwritten.
                    #https://samtools.github.io/bcftools/howtos/annotate.html
            #add the header line for this new tag
            #save as a new compressed file
                #--output
                    #Write output to a file
                #--output-type
                    #u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    print("We can see how the new tag AA_upcase has ACTG in uppercase, while the rest of characters ('.', '-', 'N') remain the same")

    print_text("Check whether the new VCF file is exactly the same than the previous one but with the addition of AA_upcase. Only do it if NOT in debugging mode because we need the whole dataset to this check to work", header=4)
    if debugging==False:
        run_bash(" \
            bcftools +split-vep \
                --annotation CSQ \
                --columns AA:String \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
            bcftools annotate \
                --remove INFO/CSQ | \
            bcftools view \
                --types snps | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO %FORMAT\n' > chr" + selected_chromosome + "_file_check_1.txt; \
            bcftools annotate \
                --remove INFO/AA_upcase \
                ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
            bcftools query \
                -f '%TYPE %ID %CHROM %POS %REF %ALT %QUAL %FILTER %INFO %FORMAT\n' > chr" + selected_chromosome + "_file_check_2.txt; \
            echo '##See head first check file:'; head -n 1 chr" + selected_chromosome + "_file_check_1.txt; \
            echo '##See head second check file:'; head -n 1 chr" + selected_chromosome + "_file_check_2.txt; \
            check_status=$(cmp --silent chr" + selected_chromosome + "_file_check_1.txt chr" + selected_chromosome + "_file_check_2.txt; echo $?); \
            echo '##Do the check'; \
            if [[ $check_status -eq 0 ]]; then \
                echo 'TRUE'; \
            else \
                echo 'FALSE'; \
            fi; \
            rm chr" + selected_chromosome + "_file_check_1.txt; rm chr" + selected_chromosome + "_file_check_2.txt; \
            ls -l")
        #open the vep VCF file without AA_upcase using bcftools
            #extract AA as a new field and then remove CSQ. This is exactly what we did with VCF file before creating AA_upcase, so we need to recreate that.
        #open the new VCF with AA_upcase
            #remove the new AA_upcase field to have the same file than the previous VCF file
        #save both files and then check they are identical using cmp
    
    print_text("calculate the percentage of SNPs with low-confidence ancestral alleles, i.e., lower-case alleles", header=4)
    print("get the number of SNPs and SNPs with upper using awk")
    counts_case = run_bash(" \
        bcftools query \
            --format '%TYPE\t%AA\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        awk \
            'BEGIN{ \
                FS=\"\t\"} \
            { \
                if(index($2, \"N\")==0 && index($2, \"-\")==0 && index($2, \".\")==0){ \
                    if($1 == \"SNP\"){count_1++} \
                    if($1 == \"SNP\" && $2==toupper($2)){count_2++} \
                    if($1 == \"SNP\" && $2==tolower($2)){count_3++} \
                } \
            }END{ printf \"count_snps_acgt=%s,count_snps_acgt_upper_anc=%s,count_snps_acgt_lower_anc=%s\", count_1, count_2, count_3 }'", return_value=True).strip()
        #get the variant type and ancestral allele (with lower and upper case) from the VCF file for each variant, being separated by tab
            #very IMPORTANT to maintain the order (AA second), because we will use AWK to filter using columns indexes assuming that AA is the second
            #also use "\t" as separator so we can correctly load the result into awk with FS=\t
        #in awk
            #if the ancestral_allele NOT includes N, -, . 
                #We do not consider empty ("") because bcftools add "." to empty cases, so we should not have them. If this would be a problem, then the next check would fail, so no problem.
                    #see what index() do
                        #awk index(str1, str2) Function: This searches in the string str1 for the first occurrences of the string str2, and returns the position in characters where that occurrence begins in the string str1. String indices in awk starts from 1.
                        #For example:
                            #awk 'BEGIN{print index("Graphic", "ph"); print index("University", "abc")}'
                            #Gives 4 and 0, because "ph" is included in "Graphic" (position 4), while "abc" is not included in "University"
                        #therefore, if the result of index() is 0, means that the second string is NOT included in the first string
                        #https://www.geeksforgeeks.org/built-functions-awk/
                        #https://stackoverflow.com/a/25292338/12772630
                                #if the variant in the row is a SNP, then add one to count_1
                #if the variant in the row is a SNP
                    #add 1 to count_1
                    #these are  the SNPs with AA=ACGT or acgt, because we have already filtered out variants with AA equals to ".", "N" or "-"
                #if the variant in the row is a SNP and its ancestral allele is upper case
                    #add 1 to count_2
                #if the variant in the row is a SNP and its ancestral allele is lower case
                    #add 1 to count_3
            #at the END, print all counts
    print("extract the counts")
    count_snps_acgt = int(counts_case.split(",")[0].replace("count_snps_acgt=", ""))
    count_snps_acgt_upper_anc = int(counts_case.split(",")[1].replace("count_snps_acgt_upper_anc=", ""))
    count_snps_acgt_lower_anc = int(counts_case.split(",")[2].replace("count_snps_acgt_lower_anc=", ""))
    print("The number of SNPs with high-confidence ancestral alleles is " + str(count_snps_acgt_upper_anc))
    print("The number of SNPs with low-confidence ancestral alleles is " + str(count_snps_acgt_lower_anc))
    print("The total number of SNPs with lower or upper ancestral allele in this chromosome is " + str(count_snps_acgt))
    print("Calculate the percentage of SNPs with high-confidence ancestral alleles")
    percent_lower_anc = (count_snps_acgt_lower_anc/count_snps_acgt)*100
    if(count_snps_acgt == count_snps_acgt_lower_anc+count_snps_acgt_upper_anc):
        print("IMPORTANT RESULT: The percentage of SNPs with low-confidence ancestral allele is " + str(percent_lower_anc))
    else:
        raise ValueError("SERIOUS ERROR! We have not correctly calculated the number of snps with lower/upper case ancestral allele")

    print_text("calculate again the number of SNPs with low- and high-confidence ancestral alleles using bcftools this time", header=4)
    count_snps_acgt_lower_anc_2 = int(run_bash(" \
        bcftools query \
            --include 'TYPE=\"SNP\" && AA=\"a,c,g,t\"' \
            --format '%TYPE\t%AA\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}'", return_value=True).strip())
    count_snps_acgt_upper_anc_2 = int(run_bash(" \
        bcftools query \
            --include 'TYPE=\"SNP\" && AA=\"A,C,G,T\"' \
            --format '%TYPE\t%AA\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        awk \
            'BEGIN{FS=\"\t\"}END{print NR}'", return_value=True).strip())
        #select those variants that are SNPs and whose AA field is ACGT or acgt and obtain the type and AA
            #we have cases like A,A,A because every consequence has an ancestral allele.
            #AA=A includes also A,A,A. I checked that, with ACGT we are selecting any SNP that has ACGT in the AA field, even if you have several "A"s separated by comma.
            #This is ok because we know that every consequence has always the same ancestral allele for the same SNP (I have checked that). So if A is present, the rest of ancestral alleles for this SNP will be also A. That SNPs will not have "-", "N" or ".". Indeed, it will not have anything else than "A".
        #then calculate the number of rows after the filtering
    if(count_snps_acgt_lower_anc_2 == count_snps_acgt_lower_anc) & (count_snps_acgt_upper_anc_2 == count_snps_acgt_upper_anc):
        print("YES! GOOD TO GO!")
    else:
        raise ValueError("SERIOUS ERROR! We have not correctly calculated the number of snps with lower/upper case ancestral allele")

    print_text("Check the number of SNPs without ancestral allele respect to the total", header=4)
    #por aquii
    snps_no_ancestral = int(run_bash(" \
        bcftools query \
            --format '%TYPE\t%AA\n' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        awk \
            'BEGIN{FS=\"\t\"}{ \
                if($1==\"SNP\"){ \
                    if(index($2, \".\")!=0 || index($2, \"-\")!=0 || index($2, \"N\")!=0){ \
                        count++ \
                    } \
                } \
            }END{print count}'", return_value=True).strip())

    snps_total = int(run_bash(" \
        bcftools view \
            --drop-genotypes \
            --no-header \
            --include 'TYPE=\"SNP\"' \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz | \
        awk 'END{print NR}'", return_value=True).strip())
    
    snps_total == snps_no_ancestral + count_snps_acgt_upper_anc + count_snps_acgt_lower_anc








print_text("check that the new INFO/TAG with ancestral alleles in upper case is exactly the same than the original AA tag but in uppercase always", header=4)
print("calculate the number of SNPs for which AA is just AA_upcase but in lowercase")
count_aa = run_bash(" \
    bcftools view \
        ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            for(i=1;i<=NF;i++){ \
                if($i==\"AA\"){ \
                    if(toupper($(i+1))==$(i+3)){ \
                        count++; \
                        next \
                    } \
                } \
            } \
        }; \
        END{ \
            print count\
            }'", return_value=True).strip()
print("then check that this number is equal to the total number of SNPs")
check_aa = run_bash(" \
    n_variants=$( \
        bcftools view \
            ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.anc_up.vcf.gz \
            --drop-genotypes \
            --no-header | \
        awk \
            'END{print NR}'); \
    if [[ " + count_aa +  " -eq $n_variants ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi", return_value=True).strip()
    #load the VCF file without the header and the genotypes, having one row per SNP
    #count the number of rows at the end of the file with awk
    #if the total number of rows (i.e., SNPs) is equal to the number of SNPs for which AA is just AA_upcase in lowercase, then we are good, all SNPs are ok.
if (check_aa == "TRUE"):
    print("GOOD TO GO! The new tag with ancestral alleles is exactly the same than the original AA field but in upper case")
else:
    raise ValueError("SERIOUS ERROR! We have not correctly converted to upper case all ancestral alleles")

print_text("check that REF and ALT are always in upper case", header=4)
print("calculate the number of SNPs for which the REF and ALT alleles are in uppercase")
count_ref_alt_upper = run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
        --drop-genotypes \
        --no-header | \
    awk \
        'BEGIN{ \
            FS=\"\t|;|=\"; \
            OFS=\"\t\"}; \
        { \
            if(toupper($4)==$4 && toupper($5)==$5){ \
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
print("then check that this number is equal to the total number of SNPs")
check_ref_alt_upper = run_bash(" \
    n_variants=$( \
        bcftools view \
            ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
            --drop-genotypes \
            --no-header | \
        awk \
            'END{print NR}'); \
    if [[ " + count_ref_alt_upper +  " -eq $n_variants ]]; then \
        echo 'TRUE'; \
    else \
        echo 'FALSE'; \
    fi", return_value=True).strip()
    #load the VCF file without the header and the genotypes, having one row per SNP
    #count the number of rows at the end of the file with awk
    #if the total number of rows (i.e., SNPs) is equal to the number of SNPs for which AA is just AA_upcase in lowercase, then we are good, all SNPs are ok.
if (check_ref_alt_upper == "TRUE"):
    print("GOOD TO GO! The REF and ALT alleles are always in uppercase")
else:
    raise ValueError("SERIOUS ERROR! We do not have the all REF and ALT alleles in upper case so we cannot correctly select those SNPs whose REF is not AA, because our ancestral allele data is always upper case and bcftools filters are case sensitive")

print_text("Use the new AA_upcase tag to filter and check the behavior", header=4)
print("SNPs where REF is ancestral")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("SNPs where REF is NOT ancestral but ALT is")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF!=AA_upcase && ALT=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("SNPs where REF nor ALT is ancestral")
run_bash("\
    bcftools query \
        --include 'AA_upcase=\"A,C,G,T\" && REF!=AA_upcase && ALT!=AA_upcase' \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA_upcase\n' \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf")
print("We get as cases with AA not being the ALT those multiallelics where 1 of the ALT alleles are indeed the ancestral. For example, REF=A;ALT=T,G;AA=G is considered when filtering by ALT!=AA_upcase, when indeed G is the ancestral and it is one of the alternative alleles. Oddly enough, these SNPs are again filtered in when doing ALT=AA_upcase. In both cases makes sense, because we have ALT alleles that are the AA but other are not. For things like this, we should do ancestral filtering after dealing with multiallelic snps.")


#repeate this with the AA_upcase
run_bash("\
    bcftools view \
        --exclude 'INFO/AC=INFO/AN || INFO/AC=0' \
        ./results/00_vep_vcf_files/1kGP_high_coverage_Illumina.chr" + selected_chromosome + ".filtered.SNV_phased_panel.vep.vcf.gz | \
    bcftools norm \
        --multiallelic +snps | \
    bcftools view \
        --max-alleles 2 \
        --min-alleles 2 | \
    bcftools +split-vep \
        --annotation CSQ \
        --columns AA:String | \
    bcftools annotate \
        --remove INFO/CSQ | \
    bcftools view \
        --exclude 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT!=AA' | \
    bcftools view \
        --include 'AA=\"A,C,G,T,a,c,g,t\" && REF!=AA && ALT=AA' | \
    bcftools query \
        --format '%TYPE %ID %CHROM %POS %REF %ALT %AA\n' | \
    head -n 20")




print_text("See the dummy VCF file after all operations with the Ancestral Alleles as a new field with all bases in uppercase to match that of REF and ALT columns, so we will be able to apply filters in next steps", header=2)
run_bash(" \
    bcftools view \
        ./data/dummy_vcf_files/00_dummy_vcf_files_vep/dummy_example_vep_2_anc_up.vcf \
        --drop-genotypes")



    ###quick check if there is anything useful for this script from the previous version

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








    print("remove the file that was used as input for VEP")
    run_bash(" \
        rm " + input_vcf_file_vep + "; \
        ls -l")


    #restore sys.stdout using the previously saved reference to it
    #This is useful if you intend to use stdout for other things
    #only required if we are in production, as we changed stdout only in that case
    if debugging==False:
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




###########################################################################################
#### Number of SNPs with low- and high-confidence ancestral alleles across chromosomes ####
###########################################################################################
print_text("Number of SNPs with low- and high-confidence ancestral alleles across chromosomes", header=1)
print_text("empty lists to save the counts of high and low confidence ancestral alleles across chromosomes", header=4)
count_snps_acgt_list = []
count_snps_acgt_upper_anc_list = []
count_snps_acgt_lower_anc_list = []

print_text("run loop across chromosomes", header=4)
#chrom=1
for chrom in range(1,23):
    print("Doing chromosome " + str(chrom))

    print("from the output file of the selected chromosome, extract the row with the counts of high and low confidence ancestral alleles")
    row_results = run_bash(" \
        grep \
            'count_snps_acgt' \
            ./scripts/00_ancestral_calcs_outputs/chr" + str(chrom) + ".out", return_value=True).strip()

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
    count_snps_acgt_list.append(int(row_results_split[0].replace("count_snps_acgt=", "")))
    count_snps_acgt_upper_anc_list.append(int(row_results_split[1].replace("count_snps_acgt_upper_anc=", "")))
    count_snps_acgt_lower_anc_list.append(int(row_results_split[2].replace("count_snps_acgt_lower_anc=", "")))

print_text("process the results", header=4)
print("check we have all chromosomes")
#selected_list=list_results[0]
check_count = [len(selected_list) == 22 for selected_list in [count_snps_acgt_list, count_snps_acgt_upper_anc_list, count_snps_acgt_lower_anc_list]]
if sum(check_count) == len(check_count):
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs with low and high confidence ancestral alleles")

print("make the sum across chromosomes")
count_snps_acgt_total = sum(count_snps_acgt_list)
count_snps_acgt_upper_anc_total = sum(count_snps_acgt_upper_anc_list)
count_snps_acgt_lower_anc_total = sum(count_snps_acgt_lower_anc_list)

print("check the total is equal to the sum of high and low-confidence")
if count_snps_acgt_total == count_snps_acgt_upper_anc_total+count_snps_acgt_lower_anc_total:
    print("YES! GOOD TO GO!")
else:
    raise ValueError("FALSE! ERROR! We have a problem calculating the number of SNPs with low and high confidence ancestral alleles")

print("see the sums")
print("Total number of SNPs with ancestral allele: " + str(count_snps_acgt_total))
print("Total number of SNPs with high-confidence ancestral allele: " + str(count_snps_acgt_upper_anc_total))
print("Total number of SNPs with low-confidence ancestral allele: " + str(count_snps_acgt_lower_anc_total))

print("calculate the percentage of SNPs with low-confidence ancestral allele across all chromosomes")
print(count_snps_acgt_lower_anc_total/count_snps_acgt_total*100)








#por aquiii
    #run the script 
    #check the redirection of stdout
    #check the whole script in the meantime
    #while running
        #ask Jesus
            #about using low-confidence ancestral alleles
            #see email from Jesus about VEP installation, and answer once you have solved the upper case problem
            #https://mail.google.com/mail/u/0/?tab=rm&ogbl#drafts/QgrcJHsHpDRJdfjndBxlCjQHdCNBwJJqNSl 
        #ask David
            #about using low-confidence ancestral alleles
        #prepare next step


#for 01d_hap_map_calcs.py
    #You will do the polarization within each populatin. We need only biallelic snps (to easily exchange ref by alt if needed), and we need to do this within pop, porque an allele can be biallelic for one pop but not for other.
    #list of SNPs to for which we have to switch REF/ALT
        #WE HAVE TO EXCLUDE CASES WHERE REF NOR ALT ARE THE ANCESTRAL ALLELE
            #these can be multiallelic SNPs for which one of the ALTs is not present in the selected population and that very ALT is the ancestral. We need these SNPs OUT, if no ancestral, we cannot estimate selection.
            #you hav eto do it after removing truly multiallelci snps in the population
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
