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

#Respect to V1:
	#Preparing the scripts for the container



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(dplyr) #for full_join and bind_row
require(plyr)
require(foreach) #for parallel
require(doParallel) #for parallel



##########################################
########## REMOVE PREVIOUS WORKSPACE #####
##########################################
remove(list=ls(all=TRUE))



##########################################
########## SET THE WORKING DIRECTORY #####
##########################################

#we do not need to specify path inside the container, starting with "/" and then opt or root, you get there
path_inside_container = "/"

#path of the starting folder in the image. This will let us to automatize the script, it will work both in my laptop and the HPC. We do not have to change the path. 
path_starting_folder = system("pwd", intern=TRUE) #intern: a logical (not 'NA') which indicates whether to capture the output

#set the path to data 
path_outside_results = paste(path_starting_folder, "/results", sep="")
path_outside_pipeline = paste(path_starting_folder, "/david_pipeline/exdef_folder", sep="")



########################################
########### WRITE FUNCTION #############
########################################

#for debugging
#pop_group="all"; statistics=c("ihs"); window_sizes=c("50kb", "1000kb"); pop_p_val=c("EUR", "EAS")
curve_significance = function(pop_group, statistics, window_sizes, pop_p_val){

	##Function for obtaining the significance of the enrichment of metabolic genes in rank thresholds considering specific summary statistics, window sizes and populations from the pipeline of David
		#pop_group: The population group considered to run the pipeline. For this group, the enrichment was calculated for different statistics and window sizes and subpopulations.
		#statistics: The summary statistics you want to consider to calculate the significance of the enrichment
		#window_sizes: The window sizes you want to consider to calculate the significance of the enrichment
		#pop_p_val: The subpopulations you want to consider to calculate the significance of the enrichment
			#If you have analyzed "all" populations, then you can select any possible subset of populations for calculating the p-value.

	#get the paths of all files in the output folder
	paths_results_raw = list.files(paste(path_outside_pipeline, "/test_outputs", sep=""), full=TRUE)

	#extract the name of each path, i.e., the statistic, window size and genome
	#we use paths_to_data, which is the input for getting all the outputs, so we have the same order
	files_names_raw = strsplit(as.character(paths_results_raw), split="/")
	files_names = sapply(files_names_raw, "[", length(files_names_raw[[1]]))
		#we can use the length of the first path split as reference because all paths should have the same length as split. We are splitting by "/", all the files are in the same folder with the same subfolders

	#select paths of those files belonging to the selected population AND window size AND statistic
	conditions_paths = grepl(pop_group, files_names, fixed=TRUE) & grepl(paste(window_sizes, collapse="|"), files_names, fixed=FALSE) & grepl(paste(statistics, collapse="|"), files_names, fixed=FALSE)
		#fixed=FALSE because we are using regular expressions.
		#fixed=TRUE for populations because we are only to use one population, either "all" or a subgroup of populations
		#we want only those having the selected window sizes ONLY for the selected statistics AND selected population

	#select only those belonging to the conditions
	paths_results = paths_results_raw[which(conditions_paths)]

	#extract the name of each final path, i.e., the statistic, window size and genome
	#we use paths_to_data, which is the input for getting all the outputs, so we have the same order
	final_files_names_raw = strsplit(as.character(paths_results), split="/")
	final_files_names = sapply(final_files_names_raw, "[", length(final_files_names_raw[[1]]))
		#we can use the length of the first path split as reference because all paths should have the same length as split. We are splitting by "/", all the files are in the same folder with the same subfolders

	#check we have selected the correct paths
	print("########################################################")
	print(paste(paste(pop_p_val, collapse="|"), ": CHECK WE HAVE SELECTED THE CORRECT PATHS", sep="")); print(length(which(!grepl(pop_group, final_files_names, fixed=TRUE) & !grepl(paste(window_sizes, collapse="|"), final_files_names, fixed=FALSE) & !grepl(paste(statistics, collapse="|"), final_files_names, fixed=FALSE))) == 0)
	print("########################################################")
		#we should have no paths without any of the selected populations, window sizes and statistics


	##get ids of fake genomes
	#extract the paths of the fake genomes
	paths_results_fake = paths_results[which(grepl("fake", final_files_names))]
	
	#function
	#for debugging
	#path=paths_results_fake[2]
	get_ids_fake = function(path){

		#split
		split_path = strsplit(as.character(path), split="_")[[1]]

		#positions of the indexes (two last positions)
		position_to_take = c(length(split_path)-1, length(split_path))
		
		#paste the id
		id = paste(split_path[position_to_take], collapse="_")
		
		#return
		return(id)
	}

	#apply the function across all paths to get the IDs of all fake genomes
	unique_fake_ids = unlist(lapply(paths_results_fake, get_ids_fake))

	#save the IDs with the corresponding paths
	unique_fake_ids_paths = data.frame(ids=unique_fake_ids, path=paths_results_fake)

	#check we have the paths for all fake genomes
	print("########################################################")
	print(paste(paste(pop_p_val, collapse="|"), ": CHECK WE HAVE THE PATHS FOR ALL FAKE GENOMES", sep="")); print(nrow(unique_fake_ids_paths) == length(paths_results_fake))
	print("########################################################")

	
	##get ids of real genome
	#extract the paths for the outputs of the real genome
	paths_results_real = paths_results[which(!grepl("fake", final_files_names))]

	#save the IDs with the corresponding paths
	unique_real_ids_paths = data.frame(ids="real", path=paths_results_real)

	#check we have the paths for the real genome
	print("########################################################")
	print(paste(paste(pop_p_val, collapse="|"), ": CHECK WE HAVE THE PATHS FOR THE REAL GENOME", sep="")); print(nrow(unique_real_ids_paths) == length(paths_results_real))
	print("########################################################")


	##bind all the paths and indexes
	#bind
	final_indexes = rbind.data.frame(unique_real_ids_paths, unique_fake_ids_paths)

	#get the unique indexes
	unique_indexes = unique(final_indexes$ids)
		#we can have several outputs for the same genome (real or fake), because for each window size and statistic a different file is created.


	##calculate the enrichment for each index, i.e., for the real genome and each of the fake genomes
	
	#write function
	#for debugging
	#unique_index=unique_indexes[2]
	calc_enrichment = function(unique_index){

		#select the paths for all files of the corresponding real/fake genome
		paths_to_data = final_indexes[which(final_indexes$ids == unique_index),]$path
			#we only have to match by index, but not by pop, statistic or window size as we have already matched by these conditions the paths

		#for debugging
		#paths_to_load=paths_to_data[2]
		read_fdr_outputs_func = function(paths_to_load){

			#convert the path to character to avoid problems with read.table and readLines
			paths_to_load = as.character(paths_to_load)


			##read the selected table avoiding the last row, if not, it gives and errors
			#open the connection to load use readLines
			con_file = file(paths_to_load)
			#read the table without the last line
			table_loaded = read.table(paths_to_load, sep=" ", header=FALSE, nrow=length(readLines(con_file))-1)
				#for that, we have to read the file with readLines (gives no error) and calculate the number of lines (the path has to be converted to "connection using "file"). Then, you can subtract 1 to get the total number of lines you want to read.
					#https://stat.ethz.ch/pipermail/r-help/2008-January/152807.html
					#https://stackoverflow.com/questions/25298840/con-not-a-connection-error-in-r-program
			#close the connection
			close(con_file)
				#https://stackoverflow.com/questions/6304073/warning-closing-unused-connection-n


			#change column names
			colnames(table_loaded) = c("threshold", "pop", "int_ctrl_ratio", "interest_genes", "control_genes", "unknown_A", "unknown_B", "unknown_C")

			#return the table
			return(table_loaded)
		}

		#load all the outputs for the selected index
		list_outputs = lapply(paths_to_data, read_fdr_outputs_func)

		#extract the name of each final path, i.e., the statistic, window size and genome
		#we use paths_to_data, which is the input for getting all the outputs, so we have the same order
		outputs_names_raw = strsplit(as.character(paths_to_data), split="/")
		outputs_names = sapply(outputs_names_raw, "[", length(outputs_names_raw[[1]]))
			#we can use the length of the first path split as reference because all paths should have the same length as split. We are splitting by "/", all the files are in the same folder with the same subfolders

		#change the column names
		names(list_outputs) = outputs_names

		#bind all the tables into a single data.frame
		all_outputs = bind_rows(list_outputs, .id="id")
			#When row-binding, columns are matched by name, and any missing columns will be filled with NA.
			#id. : Data frame identifier. When .id is supplied, a new column of identifiers is created to link each row to its original data frame. The labels are taken from the named arguments to bind_rows(). When a list of data frames is supplied, the labels are taken from the names of the list. If no names are found a numeric sequence is used instead.
			#https://dplyr.tidyverse.org/reference/bind.html

		#check we have the correct number of rows and all the required outputs
		print("########################################################")
		print(paste(paste(pop_p_val, collapse="|"), " - genome ", unique_index, ": CHECK WE HAVE THE CORRECT NUMBER OF ROWS AND ALL THE REQUIRED OUTPUTS", sep="")); print(nrow(all_outputs) == sum(sapply(list_outputs, nrow))); print(identical(unique(all_outputs$id), outputs_names))
		print("########################################################")

		
		##select the interest populations
		#remove ":" from the population names
		pop_no_dots = unlist(strsplit(as.character(all_outputs$pop), split=":"))
		
		#check we have removed correctly dots from pops names
		print("########################################################")
		print(paste(paste(pop_p_val, collapse="|"), " - genome ", unique_index, ": CHECK WE HAVE REMOVED CORRECTLY DOTS FROM POPS NAMES", sep="")); print(identical(as.character(all_outputs$pop), paste(pop_no_dots, ":", sep="")))
		print("########################################################")

		#save in the table the new pop variable
		all_outputs$pop = pop_no_dots

		#select the interest populations
		all_outputs = all_outputs[which(all_outputs$pop %in% pop_p_val),]

		#check we have selected the enrichemnt only for the p-value pops
		print("########################################################")
		print(paste(paste(pop_p_val, collapse="|"), " - genome ", unique_index, ": CHECK WE HAVE SELECTED THE ENRICHEMNT ONLY FOR THE P-VALUE POPS", sep="")); print(length(which(!all_outputs$pop %in% pop_p_val))==0)
		print("########################################################")
			#no output should belongs to a non-selected population

		#make an enrichment plot only for the real genome and if we have 200kb windows (intermediate size)
		outputs_plotting = which(grepl("200kb", all_outputs$id, fixed=TRUE))
		if(unique_index == "real" & length(outputs_plotting)>0){

			#select outputs for plotting
			all_outputs_plot = all_outputs[outputs_plotting,]

			#open the pdf
			pdf(paste(path_outside_results, "/figures/", paste(pop_p_val, collapse="_"), "_", paste(statistics, collapse="_"), "_200kb_fold_enrichment.pdf", sep=""))

			#plot enrichment of each threshold
			plot(x=all_outputs_plot$threshold, y=all_outputs_plot$int_ctrl_ratio, xlim=rev(range(all_outputs_plot$threshold)), xlab="Rank threshold", ylab="Fold enrichment", type="l", lwd=2, main=paste(paste(pop_p_val, collapse="|"), " - ", paste(statistics, collapse="|"), " - 200kb"))
				#we to use reverse in xlim for plotting first the bigger ranks
					#https://www.r-graph-gallery.com/77-turn-y-axis-upside-down.html

			#close the plot
			dev.off()
		}

		#sum the difference between interest genes and controls across the tops and windows
		int_ctrl_differ = sum(all_outputs$interest_genes - all_outputs$control_genes)
			#we are not interested in absolute value. If the controls are higher than interest genes we need the negative value to decrease the overall sum.

		#return
		return(cbind.data.frame(unique_index, int_ctrl_differ))
	}

	#apply the function
	results_enrichment_list = lapply(unique_indexes, calc_enrichment)

	#set names
	names(results_enrichment_list) = unique_indexes

	#bind all the tables into a single data.frame
	results_enrichment = bind_rows(results_enrichment_list, .id="id")
		#When row-binding, columns are matched by name, and any missing columns will be filled with NA.
		#id. : Data frame identifier. When .id is supplied, a new column of identifiers is created to link each row to its original data frame. The labels are taken from the named arguments to bind_rows(). When a list of data frames is supplied, the labels are taken from the names of the list. If no names are found a numeric sequence is used instead.
		#https://dplyr.tidyverse.org/reference/bind.html

	#check genome IDs
	if(identical(results_enrichment$id, as.character(results_enrichment$unique_index))){

		#remove the ID column
		results_enrichment$id = NULL
	} else { #if not

		#we have a problem, stop
		stop("WE HAVE A PROBLEM WITH THE GENOME IDs")
	}


	##check that we have the correct number of fake genomes
	#select only the fake genomes
	subset_fake_genomes = results_enrichment[which(results_enrichment$unique_index!="real"),]

	#open a connection to the file with the input parameters for the pipeline
	check_con=file("/opt/scripts/input_par_david_pipeline_v2.txt")

	#read the lines of the file
	lines_input_par_pipeline = readLines(check_con)

	#extract the row with FDR and avoid the lines with comments (i.e., "#")
	fdr_number_raw = lines_input_par_pipeline[which(grepl("FDR_number", lines_input_par_pipeline, fixed=TRUE) & !grepl("#", lines_input_par_pipeline, fixed=TRUE))]
		#fixed=TRUE because we are not using regular expressions
	fdr_number = strsplit(fdr_number_raw, split="FDR_number ")[[1]][2]

	#check the number indicated as FDR number is equal to the number of fake genomes we have
	print("########################################################")
	print(paste(paste(pop_p_val, collapse="|"), ": CHECK THE NUMBER INDICATED AS FDR NUMBER IS EQUAL TO THE NUMBER OF FAKE GENOMES WE HAVE", sep="")); print(as.character(nrow(subset_fake_genomes)) == fdr_number)
	print("########################################################")
		#remember that at the end of this script we have calculated the sum of the difference between interest genes and controls within each fake genome, so we should have one value per fake genome.

	#close the connection
	close(check_con)
		#https://stackoverflow.com/questions/6304073/warning-closing-unused-connection-n


	##calculate a p-value for the enrichment
	#extract the observed enrichment
	observed_enrichment = results_enrichment[which(results_enrichment$unique_index=="real"),]$int_ctrl_differ

	#calculate the number of fake genomes having a value of enrichment equal or higher than that of the real genome
	fake_equal_higher = subset_fake_genomes[which(subset_fake_genomes$int_ctrl_differ >= observed_enrichment),]
		#these are fake genomes where there is a larger enrichment of interest genes in rank thresholds compared to the real genome

	#check we have selected correct fake genomes to calculate p-values
	print("########################################################")
	print(paste(paste(pop_p_val, collapse="|"), ": CHECK WE HAVE SELECTED CORRECT FAKE GENOMES TO CALCULATE P-VALUES:", sep="")); print(length(which(fake_equal_higher$int_ctrl_differ<observed_enrichment))==0)
	print("########################################################")

	#calculate p.vale as the probability of a random genome having an enrichment of metabolic genes in rank thresholds equal or higher than the real genome
	p_value = prop.test(x=nrow(fake_equal_higher), n=nrow(subset_fake_genomes)) 
		#https://stats.stackexchange.com/questions/167164/test-of-equal-proportions-with-zero-successes

	#print the FINAL P-VALUE
	print("########################################################")
	print(paste("PVALUE of ", paste(pop_p_val, collapse="|"), " - ", paste(statistics, collapse="|"), " - ", paste(window_sizes, collapse="|"), " - ", nrow(subset_fake_genomes), " random_genomes", sep="")); print(p_value)
	print("########################################################")
}



#####################################
###### PARALLELIZE THE PROCESS ######
#####################################

#statistics
pops_to_get_results = list("EUR", "EAS", "SAS", "AMR", "AFR", c("EUR", "EAS"))
	#we can get the significance not only for a pre-defined group of pops but also for the combination of two groups like EUR and EAS

#set up cluster
clust <- makeCluster(length(pops_to_get_results), outfile="") #outfile let you to see the output in the terminal "https://blog.revolutionanalytics.com/2015/02/monitoring-progress-of-a-foreach-parallel-job.html"
registerDoParallel(clust)

#run the function for all populations
foreach(i=pops_to_get_results, .packages=c("plyr", "dplyr")) %dopar% {
	curve_significance(pop_group="all", statistics=c("ihs"), window_sizes=c("50kb", "100kb", "200kb", "500kb", "1000kb"), pop_p_val=i)
}

#stop the cluster 
stopCluster(clust)