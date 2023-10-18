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


####
echo -e "\n## check we have the correct files for each chromosome after calculating the ancestral allele ##"

#go to the working directory
cd ../
cd ./results/00_vep_vcf_files

#for each autosomal chromosome (sequence from 1 to 22)
#chrom=1
for chrom in $(seq 1 22) #https://stackoverflow.com/a/169517/12772630
do
	echo "chromosome "$chrom

	#list the files in the directory of the selected chromosome
	#then count the number of times we have files with specific names
	count_1=$(ls "./chr"$chrom | grep --count "1kGP_high_coverage_Illumina.chr"$chrom".filtered.SNV_phased_panel.vep.anc_up.vcf.gz")
	count_2=$(ls "./chr"$chrom | grep --count "1kGP_high_coverage_Illumina.chr"$chrom".filtered.SNV_phased_panel.vep.vcf.gz")
		#we should have one file matching this pattern
			#the html file
	count_3=$(ls "./chr"$chrom | grep --count "anc_alleles_uppercase_chr"$chrom".tsv.gz")
		#we should have two files matching this pattern
			#the file with the pos and AA allele
			#the tbi file
	
	#check we have the correct number of files with the corresponding names
	if [[ $count_1 -eq 1 && $count_2 -eq 1 && $count_3 -eq 2 ]]; then
		echo "TRUE";
	else
		echo "FALSE";
	fi
done


####
echo -e "\n## check the number of variants with '.', 'N' or '-' for ancestral allele ##"
#we could still lose more than these snps because some can have AA that is not REF nor ALT. See below.

#go to the working directory
cd ../..
cd ./scripts/00_ancestral_calcs_outputs/

#set variable with zero. Here we will be summing the number of SNPs without ancestral data across chromosomes
sum=0

#for each autosomal chromosome (sequence from 1 to 22)
#chrom=1
for chrom in $(seq 1 22)
do
	
	#calculate the number of snps without ancestral data
	number_to_sum=$( \
		grep "The number of SNPs without ancestral allele data is" "chr"$chrom".out" | \
		awk \
			'BEGIN{FS=" |\\."}{ \
				for(i=1;i<=NF;i++){ \
					if($i=="data" && $(i+1)=="is"){ \
						print $(i+2) \
					} \
				} \
			}')
		#first get, from the output file of the corresponding chromosome, the line with the percentage of SNPs without ancestral data
		#in awk
			#use space and "." as delimiters, so we have the "." in the middle of the sentence removed (we have to escape "." in order to be considered as "." and not as a special character)
			#for each field
				#if the field is "data" and the next one is "is", then print the next next field, which is the number of snps without ancestral data

	#save that number to the previous total
	sum=$(($sum + $number_to_sum))
		#https://unix.stackexchange.com/a/499036
done

#see the sum
echo "The total number of SNPs without ancestral allele is "$sum


####
echo -e "\n## check the number of SNPs for which the ancestral allele is not REF nor ALT ##"
#go to the working directory
cd ../..
cd ./scripts/00_ancestral_calcs_outputs/

#set variable with zero. Here we will be summing the number of SNPs for which AA is not REF nor ALT
sum_2=0

#for each autosomal chromosome (sequence from 1 to 22)
#chrom=1
for chrom in $(seq 1 22); do 
	number_to_sum_2=$( \
		awk '{ \
				if(/calculate the number of these problematic cases, to check this is not a problem/){ \
					getline; \
					print $0 \
				} \
			}' \
			"chr"$chrom".out")
		#if the text of the current line includes the string about this specific check, go to the next line (getline), which is the one with only the count of cases where AA is not REF nor ALT, and print it
		#https://stackoverflow.com/a/7451456/12772630
		#https://unix.stackexchange.com/a/127684
	
	#save that number to the previous total
	sum_2=$(($sum_2 + $number_to_sum_2))
		#https://unix.stackexchange.com/a/499036
done

#see the sum
echo "We have "$sum_2" SNPs for which the ancestral allele is NOT the REF nor the ALT allele"
