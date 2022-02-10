#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



###########################################################################
########################### CALCULATE P-VALUES ############################
###########################################################################

#Calculation of the significance of the whole enrichment curve using randomized genomes (FDR approach in the pipeline of David)

#IMPORTANT: 
	#Note that this script was wrote and run in 2022 using the version of April 28th 2022 of David's pipeline.
		#The pipeline of David has been downloaded from 
			#https://github.com/DavidPierreEnard/Gene_Set_Enrichment_Pipeline
			#last commit done on Apr 28, 2021
				#7b755c0c23dd4d7c3f54c4b53e74366e4041ac8f 

	#Manual from David
		#/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/exdef_pipeline_manual.pdf



#################################################################
####################### PREVIOUS VERSIONS #######################
#################################################################



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################






#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



###########################################################################
########################### PREPARE RANKS #################################
###########################################################################

#create tables with the sweep rank files to test. The first column in a sweep file is the Ensembl gene IDs (if you are not using Ensembl human gene annotations, add ENSG to the start of your gene IDs). The following columns corresponding to gene ranks, for all the populations included, from the gene with the strongest sweep signal to the one with the lowest signal in each specific population. For example, the all_ihsfreqabs_ranks_”size” files have 27 columns with 26 corresponding to gene sweep ranks according to iHS in the 26 populations from 1000 Genomes phase 3.

#IMPORTANT: 
	#Note that this script was wrote and run in 2022 using the version of April 28th 2022 of David's pipeline.
		#The pipeline of David has been downloaded from 
			#https://github.com/DavidPierreEnard/Gene_Set_Enrichment_Pipeline
			#last commit done on Apr 28, 2021
				#7b755c0c23dd4d7c3f54c4b53e74366e4041ac8f 

	#Manual from David
		#/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/exdef_pipeline_manual.pdf



#################################################################
####################### PREVIOUS VERSIONS #######################
#################################################################



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(dplyr) #for full_join
require(plyr)
require(foreach) #for parallel
require(doParallel) #for parallel



#####################################
########### WRITE FUNCTION ##########
#####################################

#for debugging
#statistics=c("ihs"); pop_group="all"; window_sizes=c("50kb", "1000kb")
curve_significance = function(statistics, pop_group, window_sizes){

	pattern_to_search = paste(c(statistics, pop_group, window_sizes), collapse="|")

	paths_results = list.files("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/test_outputs", pattern=pattern_to_search, full=TRUE)


	paths_results_fake = paths_results[which(grepl("fake", paths_results))]
	

	#for debugging
	#path=paths_results_fake[2,]
	get_ids_fake = function(path){

		split_path = strsplit(as.character(path), split="_")[[1]]

		position_to_take = c(length(split_path), length(split_path)-1)
		id = paste(split_path[position_to_take], collapse="_")
		return(id)
	}

	unique_fake_ids = unlist(lapply(paths_results_fake, get_ids_fake))

	unique_fake_ids_paths = data.frame(ids=unique_fake_ids, path=paths_results_fake)

	length(unique_fake_ids_paths) == length(unique_fake_ids_paths)

	
	paths_results_real = paths_results[which(!grepl("fake", paths_results))]

	unique_real_ids_paths = data.frame(ids="real", path=paths_results_real)


	final_indexes = rbind.data.frame(unique_real_ids_paths, unique_fake_ids_paths)

	unique_indexes = unique(final_indexes$ids)


	#for debugging
	#unique_index=unique_indexes[2]
	calc_enrichment = function(index){


		#select the paths for all files of the corresponding real/fake genome
		paths_to_data = final_indexes[which(final_indexes$ids == unique_index),]$path
			#check that works with two paths
				#paths_to_data = final_indexes[2:3,]$path

		#WE HAVE TO READ THE FILES WITHPUT THE LAST ROW, BECAUSE IT GIVES PROBLEMS

		read.table(as.character(paths_to_data[1]), sep="\t", header=FALSE, nrow=length(count.fields(paths_to_data[1]))-1)

		lapply(paths_to_data, read.table, sep=


	}
}
