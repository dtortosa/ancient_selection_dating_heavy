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

#script to check we have the correct files after calculating the ancestral allele

#go to the working directory
cd ./results/00_vep_vcf_files

#for each autosomal chromosome (sequence from 1 to 22)
#chrom=1
for chrom in $(seq 1 22) #https://stackoverflow.com/a/169517/12772630
do
	
	#list the files in the directory of the selected chromosome
	#then count the number of times we have files with specific names
	count_1=$(ls "./chr"$chrom | grep --count "1kGP_high_coverage_Illumina.chr"$chrom".filtered.SNV_phased_panel.vep.anc_up.vcf.gz")
	count_2=$(ls "./chr"$chrom | grep --count "1kGP_high_coverage_Illumina.chr"$chrom".filtered.SNV_phased_panel.vep.vcf.gz")
		#we should have two files matching this pattern
			#the VCF file
			#the html file
	count_3=$(ls "./chr"$chrom | grep --count "anc_alleles_uppercase_chr"$chrom".tsv.gz")
		#we should have two files matching this pattern
			#the file with the pos and AA allele
			#the tbi file
	
	#check we have the correct number of files with the corresponding names
	if [[ $count_1 -eq 1 && $count_2 -eq 2 && $count_3 -eq 2 ]]; then
		echo "TRUE";
	else
		echo "FALSE";
	fi
done
