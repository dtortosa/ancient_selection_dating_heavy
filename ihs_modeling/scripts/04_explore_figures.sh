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

#script to copy aleplot images from different populations into a single folder for each gene set, ordering by superpopulations, i.e., all european first, then african....

#list of population codes
populations=("CEUD" "GBRD" "TSID" "IBSD" "FIND" "CDXD" "CHBD" "CHSD" "JPTD" "KHVD" "YRID" "LWKD" "GWDD" "MSLD" "ESND" "ASWD" "ACBD" "PURD" "MXLD" "CLMD" "PELD" "GIHD" "PJLD" "BEBD" "STUD" "ITUD")

#list of gene sets
gene_sets=("thermogenic" "bat")

#iterate over each gene set 
#gene_set="thermogenic"
for gene_set in "${gene_sets[@]}"; do

    #set name for a new folder for the gene set
    new_folder="./results/ihs_modeling_across_pops/00_all_aleplots/${gene_set}"

    #if the folder already exists, remove it
    if [[ -d "${new_folder}" ]]; then
        rm -rf "${new_folder}"
    fi

    #then create the new folder
    mkdir -p "${new_folder}"

    #initialize a counter to give a number to each image
    counter=0

    #iterate over each population
    for pop in "${populations[@]}"; do

        #construct the file name differnetly for each gene set
        if [[ "${gene_set}" == "thermogenic" ]]; then
            file_name="${pop}_${gene_set}_aleplot_${gene_set}_distance.png";
        else
            file_name="${pop}_${gene_set}_aleplot_${gene_set}_distance_percentile_1.png";
        fi

        #check if the file exists
        if [[ -f "./results/ihs_modeling_across_pops/${pop}/${gene_set}/aleplots/${file_name}" ]]; then
            
            #copy the file adding the counter to the file name
            cp "./results/ihs_modeling_across_pops/${pop}/${gene_set}/aleplots/${file_name}" "${new_folder}/${counter}_${file_name}"
            
            #increment the counter
            ((counter++))
        
            #show message of the copied file
            echo "File copied: ${file_name} to ${new_folder}/${counter}_${file_name}"

        else
            
            #print a warning if the file does not exist
            echo "ERROR! FALSE! File not found: ${file_name}"
        fi
    done
done

#to run the script:
#cd /home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating_heavy_analyses/dating_climate_adaptation/ihs_modeling/
#chmod +x ./scripts/04_explore_figures.sh
#./scripts/04_explore_figures.sh > ./scripts/04_explore_figures.out 2>&1
    #we use the container of previous steps because it has alibi installed. I have been unable to install alibi in the container of this step.
#grep -Ei 'error|false|fail' ./04_explore_figures.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
