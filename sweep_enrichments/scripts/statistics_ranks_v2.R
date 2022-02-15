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

#Respect to V1:
	#Preparing the scripts for the container



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(dplyr) #for full_join
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
path_outside_data = paste(path_starting_folder, "/data", sep="")
path_outside_pipeline = paste(path_starting_folder, "/david_pipeline/exdef_folder", sep="")



#######################################
############ FILES PREPARATION ########
#######################################

#use the valid_file as source for the gene IDs
valid_file=read.table(paste(path_outside_pipeline, "/valid_file.txt", sep=""), sep="\t", header=FALSE)
	#this file includes all the ensembl gene IDs used in the pipeline, which in our case come from our curated list of gene coordinates.

#make a new folder in data/ to save summary statistics data
#system("mkdir -p /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/summ_statistics")
	#mkdir -p because if you don't and build again, you will get an error that it exists. "p" flags is for "no error if existing, make parent directories as needed"


##add the summary statistics
#IMPORTANT: this script expects a value of the statistic per gene window (as many rows as genes) and window sizes separated in columns.

#iHS from our calculations and remove raw files and number iHS datapoints
#system("cp -avr /home/dftortosa/singularity/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/summ_statistics/; cd /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/summ_statistics/mean_ihs_gene_windows; rm *_raw_v1.txt.gz; rm *_n_ihs_gene_windows_final_v1.txt.gz")
	#a: preserve the specific attributes (e.g., mode)
	#v: verbose output
	#r: copy directories recursively
	#https://www.cyberciti.biz/faq/copy-folder-linux-command-line/



#####################################
########### WRITE FUNCTION ##########
#####################################

#for debugging
#statistics=c("ihs"); pop_group="eas"; window_sizes=c("50kb", "1000kb")
statistics_ranks = function(statistics, pop_group, window_sizes){

	##load the selected statistic
	#if ihs
	if("ihs" %in% statistics){

		#select the folder with this statistic
		paths_to_statistics = paste(path_outside_data, "/summ_statistics/mean_ihs_gene_windows", sep="")

		#load the files of this statistic across populations
		list_files_statistic_paths = list.files(paths_to_statistics, pattern="_mean_ihs_gene_windows_final_v1.txt.gz", full=TRUE)
	}


	##load the data for the selected population group
	#if all
	if(pop_group == "all"){
		pops = c("ACBD","ASWD","BEBD","CDXD","CEUD","CHBD","CHSD","CLMD","ESND","FIND","GBRD","GIHD","GWDD","IBSD","ITUD","JPTD","KHVD","LWKD","MSLD","MXLD","PELD","PJLD","PURD","STUD","TSID","YRID")
	}
	#if eur
	if(pop_group == "eur"){
		pops = c("CEUD", "FIND", "GBRD", "IBSD", "TSID")
	}
	#if eas
	if(pop_group == "eas"){
		pops = c("CDXD", "CHBD", "CHSD", "JPTD", "KHVD")
	}

	#extract the file names for each path
	list_files_statistic_names_raw = strsplit(list_files_statistic_paths, split="/")
	list_files_statistic_names = sapply(list_files_statistic_names_raw, "[", length(list_files_statistic_names_raw[[1]]))
		#we can use the length of the first path split as reference because all paths should have the same length as split. We are splitting by "/", all the files are in the same folder with the same subfolders

	#select only those paths belonging to the selected populations
	list_files_statistic_paths_selected = list_files_statistic_paths[which(grepl(paste(pops, collapse="|"), list_files_statistic_names, fixed=FALSE))]
		#we use fixed=FALSE so we can use a regular expression, which is matching for each pop name (separated by "|")

	#check that we have the paths for all selected pops
	print("########################################")
	print(paste(statistics, ": CHECK THAT WE HAVE THE PATHS FOR ALL SELECTED POPS", sep="")); print(length(list_files_statistic_paths_selected) == length(pops))
	print("########################################")

	#load the data for all the selected populations
	list_files_statistic = lapply(list_files_statistic_paths_selected, read.table, sep="\t", header=TRUE)


	##calculate the ranks per window size
	#write function
	#for debugging
	#window_size = window_sizes[2]
	ranks_cal_window = function(window_size){

		#open a data.frame to save the ranks for all populations as different columns. The first column will be gene_id (obtained from valid_file)
		final_ranks = data.frame(gene_id=valid_file$V1)

		#for each population
		for(i in 1:length(list_files_statistic_paths_selected)){

			#select the [i] path
			selected_path = list_files_statistic_paths_selected[[i]]

			##extract the population name
			#split the whole path
			pop_name_raw = strsplit(selected_path, split="/|_")
			#select the element with the population name
			pop_name = pop_name_raw[[1]][which(grepl(paste(pops, collapse="|"), pop_name_raw[[1]]))]
				#we use fixed=FALSE so we can use a regular expression, which is matching for each pop name (separated by "|")

			#extract the data for the selected population
			selected_pop = list_files_statistic[[i]]
				#remember that list_files_statistic comes from list_files_statistic_paths_selected, so we can use the same index

			#select the columns of gene id and the selected window size
			selected_pop_subset = selected_pop[, which(grepl(paste("gene_id|", window_size, sep=""), colnames(selected_pop), fixed=FALSE))]
				#we use fixed=FALSE so we can use a regular expression, which is matching for several conditions

			#change column name of the selected statistic and size
			colnames(selected_pop_subset)[which(grepl(window_size, colnames(selected_pop_subset), fixed=TRUE))] = "summary_statistic"


			##make the rank
			#calculate a new order of the rows based on the statistic value on decreasing order
			decreasing_order_statistic = order(selected_pop_subset$summary_statistic, decreasing=TRUE)
				#we select the column with the selected statistic and size
			#reorder
			selected_pop_subset_ordered = selected_pop_subset[decreasing_order_statistic, ]
			
			#make the rank
			selected_pop_subset_ordered$rank = 1:nrow(selected_pop_subset_ordered)
			#set the name of the column rank as the pop name
			colnames(selected_pop_subset_ordered)[which(colnames(selected_pop_subset_ordered) == "rank")] = pop_name


			##prepare the final file
			#remove NAs for the statistic
			selected_pop_subset_ordered = selected_pop_subset_ordered[which(!is.na(selected_pop_subset_ordered$summary_statistic)),]
				#select the column with the selected statistic and size

			#check the calculation of the rank with sort
			print("########################################")
			print(paste(statistics, " - ", window_size, " - ", pop_name, ": CHECK THE CALCULATION OF THE RANK WITH SORT", sep="")); print(identical(selected_pop_subset_ordered$summary_statistic, sort(selected_pop_subset$summary_statistic, decreasing=TRUE)))
			print("########################################")

			#merge the new columns used ensembl gene id
			final_ranks = full_join(final_ranks, selected_pop_subset_ordered[, c("gene_id", pop_name)], by="gene_id")
				#maintain all rows for both data.frame, adding NAs if there is mismatch
				#this can be done with full_join, which is similar to use merge with all=TRUE
					#https://stackoverflow.com/questions/21841146/is-there-an-r-dplyr-method-for-merge-with-all-true
					#https://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
					#https://www.datasciencemadesimple.com/join-in-r-merge-in-r/
		}

		#remove NAs, i.e., genes that do not have statistic data for all studied populations
		final_ranks = na.omit(final_ranks)

		#see the table
		summary(final_ranks)
		str(final_ranks)

		#check we have all the pops
		print("########################################")
		print(paste(statistics, " - ", window_size, ": CHECK WE HAVE ALL THE POPS", sep="")); print(ncol(final_ranks) == length(list_files_statistic_paths_selected) + 1)
			#we should a column for each pop plus the ensembl ID column
		print("########################################")

		#print the population names in the order of the table
		print("########################################")
		print(paste(statistics, " - ", window_size, ": PRINT THE POPULATION NAMES IN ORDER", sep="")); print(paste(colnames(final_ranks)[which(!colnames(final_ranks) == "gene_id")], collapse=","))
		print("########################################")

		#save the data
		write.table(final_ranks, paste(path_outside_pipeline, "/", pop_group, "_ranks_", statistics, "_", window_size, sep=""), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE) 
	}

	#apply the function across the window sizes
	lapply(window_sizes, ranks_cal_window)
}



#####################################
###### PARALLELIZE THE PROCESS ######
#####################################

#statistics
statistics_to_rank = c("ihs")

#set up cluster
clust <- makeCluster(length(statistics_to_rank), outfile="") #outfile let you to see the output in the terminal "https://blog.revolutionanalytics.com/2015/02/monitoring-progress-of-a-foreach-parallel-job.html"
registerDoParallel(clust)

#run the function for all populations
foreach(i=statistics_to_rank, .packages=c("plyr", "dplyr")) %dopar% {
    statistics_ranks(statistics=i, pop_group="all", window_sizes=c("50kb", "100kb", "200kb", "500kb", "1000kb"))
}

#run the function for East Asian populations
foreach(i=statistics_to_rank, .packages=c("plyr", "dplyr")) %dopar% {
    statistics_ranks(statistics=i, pop_group="eas", window_sizes=c("50kb", "100kb", "200kb", "500kb", "1000kb"))
} 

#run the function for East Asian populations
foreach(i=statistics_to_rank, .packages=c("plyr", "dplyr")) %dopar% {
    statistics_ranks(statistics=i, pop_group="eur", window_sizes=c("50kb", "100kb", "200kb", "500kb", "1000kb"))
}

#stop the cluster 
stopCluster(clust)

#you can use a second argument within foreach
if(FALSE){ #dummy example

	#write dummy function
	dummy_function = function(x, y){
		#for each element of y
		for(i in 1:length(y)){

			#sum it the value of x
			print(x+y[i])
		}
	}

	#dummy vector of Xs
	dumy_x = c(1,2)

	#run the function across X using in each case two Y values
	clust=makeCluster(length(dumy_x), outfile="")
	registerDoParallel(clust)
	foreach(i=dumy_x, .packages=c("plyr", "dplyr")) %dopar% {
	    dummy_function(x=i, y=c(1,2))
	} #the result is 2, 3, 3, 4, which makes sense. It has run X=1, and then sum each Y value (1+1 and 1+3), then repeat for X=2.
	stopCluster(clust)
}