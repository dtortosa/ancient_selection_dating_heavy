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

#script to copy process the results



#############################################################################################
### copy figures from different populations into a single folder per gene set and metric ####
#############################################################################################

#list of population codes
populations=("CEUD" "GBRD" "TSID" "IBSD" "FIND" "CDXD" "CHBD" "CHSD" "JPTD" "KHVD" "YRID" "LWKD" "GWDD" "MSLD" "ESND" "ASWD" "ACBD" "PURD" "MXLD" "CLMD" "PELD" "GIHD" "PJLD" "BEBD" "STUD" "ITUD")

#list of gene sets
gene_sets=("thermogenic" "bat" "smt" "all_thermogenic")

#list of metrics
metrics=("aleplots" "permutation")

#iterate over each gene set 
#gene_set="thermogenic"
for gene_set in "${gene_sets[@]}"; do

    #iterate over each metric
    #metric="aleplots"
    for metric in "${metrics[@]}"; do

        #skip if gene_set is "all_thermogenic" and metric is "aleplots"
        if [[ "${gene_set}" == "all_thermogenic" && "${metric}" == "aleplots" ]]; then
            echo "Skipping gene_set=${gene_set} and metric=${metric}"
            continue
        fi
        #continue: Skips the current iteration of the inner loop (for metric) and moves to the next iteration.

        #create a new folder for each metric
        new_folder="./results/ihs_modeling_across_pops/00_all_metrics/${gene_set}/${metric}"
        
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
            if [[ "${metric}" == "aleplots" ]]; then
                file_path="./results/ihs_modeling_across_pops/${pop}/${gene_set}/aleplots"
                if [[ "${gene_set}" == "thermogenic" ]]; then
                    file_name="${pop}_${gene_set}_aleplot_${gene_set}_distance.png";
                else
                    file_name="${pop}_${gene_set}_aleplot_${gene_set}_distance_percentile_1.png";
                fi
            elif [[ "${metric}" == "permutation" ]]; then
                file_path="./results/ihs_modeling_across_pops/${pop}/${gene_set}/permutation_importance"
                file_name="${pop}_0_${gene_set}_permutation_importance_all.png";
            fi

            #check if the file exists
            if [[ -f "${file_path}/${file_name}" ]]; then
                
                #copy the file adding the counter to the file name
                cp "${file_path}/${file_name}" "${new_folder}/${counter}_${file_name}"
                
                #increment the counter
                ((counter++))
            
                #show message of the copied file
                echo "File copied: ${file_path}/${file_name} to ${new_folder}/${counter}_${file_name}"

            else
                
                #print a warning if the file does not exist
                echo "ERROR! FALSE! File not found: ${file_name}"
            fi
        done
    done
done



#############################################################################################
### copy figures from different populations into a single folder per gene set and metric ####
#############################################################################################

#define the source folder
source_folder="./results/ihs_modeling_across_pops/"

#define the destination folder
destination_folder="./results/ihs_modeling_across_pops/01_all_permutation_importance/"

#remove the destination folder if it exists and then make a new one
if [[ -d "${destination_folder}" ]]; then
    rm -rf "${destination_folder}"
fi
mkdir -p "${destination_folder}"

#find and copy all files ending with *_permutation_importance_all.png
find "${source_folder}" -type f -name "*_permutation_importance_all.png" -exec cp {} "${destination_folder}" \;
    #find "${source_folder}" -type f -name "*_permutation_importance_all.png"
        #Searches for files (-type f) in the source_folder that match the pattern *_permutation_importance_all.png.
    #-exec cp {} "${destination_folder}" \;
        #Executes the cp command for each file found, copying it to the destination_folder.


#print a message indicating completion
echo "All *_permutation_importance_all.png files have been copied to ${destination_folder}."







#!/bin/bash

# Define the source folder
source_folder="./results/ihs_modeling_across_pops/"

# Define the output file
output_file="./results/ihs_modeling_across_pops/02_all_thermogenic_model_eval_combined.tsv"

# Get the header from the first file
header=$(find "${source_folder}" -type f -name "*_all_thermogenic_model_eval.tsv" | head -n 1 | xargs head -n 1)

# Write the header to the output file
echo "${header}" > "${output_file}"

# Concatenate all files, skipping their headers, without printing file names
find "${source_folder}" -type f -name "*_all_thermogenic_model_eval.tsv" -exec tail -n +2 {} \; >> "${output_file}"

# Print a message indicating completion
echo "All *_all_thermogenic_model_eval.tsv files have been concatenated into ${output_file} with a single header."





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
