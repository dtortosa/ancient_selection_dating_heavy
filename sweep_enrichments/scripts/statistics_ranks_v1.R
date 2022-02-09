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


# IDs used in the pipeline
valid_file=read.table("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/valid_file.txt", sep="\t", header=FALSE)

#COPIAR IHS INTO DATA FOLDER?

#for debugging
#statistics=c("ihs"); pops=c("CDX", "CHB", "CHS", "JPT", "KHV"); window_sizes=c("50kb", "1000kb")
statistics_ranks = function(statistic, pops){

	paths_to_statistics=NULL
	if("ihs" %in% statistics){
		paths_to_statistics = append(paths_to_statistics, "/home/dftortosa/singularity/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows")

	}

	list_files_statistic_paths = list.files(paths_to_statistics, pattern="_mean_ihs_gene_windows_final_v1.txt.gz", full=TRUE)

	

	list_files_statistic_paths_selected = list_files_statistic_paths[which(grepl(paste(pops, collapse="|"), list_files_statistic_paths))]

	list_files_statistic = lapply(list_files_statistic_paths_selected, read.table, sep="\t", header=TRUE)

	#for debugging
	#window_size = window_sizes[1]
	ranks_cal_window = function(window_size){

		final_ranks = data.frame(gene_id=valid_file$V1)
		for(i in 1:length(list_files_statistic_paths_selected)){

			selected_path = list_files_statistic_paths_selected[[i]]

			pop_name_raw = strsplit(selected_path, split="/|_")

			pop_name = pop_name_raw[[1]][which(pop_name_raw[[1]] %in% paste(pops, "D", sep=""))]

			selected_pop = list_files_statistic[[i]]

			selected_pop_subset = selected_pop[, which(grepl(paste("gene_id|", window_size, sep=""), colnames(selected_pop)))]

			selected_pop_subset = selected_pop_subset[order(selected_pop_subset[,2], decreasing=TRUE),]
			selected_pop_subset[,3] = 1:nrow(selected_pop_subset)

			colnames(selected_pop_subset)[3] = pop_name

			selected_pop_subset = selected_pop_subset[which(!is.na(selected_pop_subset$mean_ihs_50kb)),]

			####QUITA NAS QUE SI NO TE CUENTAN!!

			final_ranks = full_join(final_ranks, selected_pop_subset[, c(1,3)], by="gene_id")
		}

		pops_names_file = paste(colnames(final_ranks)[which(!colnames(final_ranks) == "gene_id")], collapse="_")

		write.table(final_ranks, paste("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/", pops_names_file, "_", window_size, sep=""), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE) 

		#PONER NOMBRE POPS EN ORDERN Y TAMAÑO VENTANA
	}

	#RUN THE FUNCTION ACROSS WINDOWS SIZES W
}

