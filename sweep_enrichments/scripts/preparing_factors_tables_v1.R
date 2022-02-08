#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



###################################################################################
########################### PREPARE FACTORS TABLE #################################
###################################################################################

#create a table file that contains all the confounding factors that the bootstrap test will control for. First column is Ensembl gene ID, next columns are each factor that will be matched. The example file is Factors_table.txt

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

require(data.table) #for fread
require(dplyr) #for full_join



###########################
######## SET WD ###########
###########################

#path
wd_path = "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments"

#set this path
setwd(wd_path)



##########################################################
########### CONFOUNDING FACTORS ACROSS WINDOWS ###########
##########################################################

#Factors_table.txt: create a table file that contains all the confounding factors that the bootstrap test will control for. First column is Ensembl gene ID, next columns are each factor that will be matched. The example file is Factors_table.txt

#for each confounding factor, we are going to consider two window sizes, i.e., average values in 50 and 500 kb windows, following the approach of David. Therefore, for each factor, we will have two columns, i.e., to variables for which the control genes have to be matched.

#get the path of all partitions
list_factors_across_windows_paths = list.files("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/predictors_gene_windows_v2/factors_across_windows", pattern=".txt", full.names=TRUE)

#create a function to read the factor
#selected_path=list_factors_across_windows_paths[1]
reading_factors_across_windows = function(selected_path){

	#read the table in the selected_path
	full_table = read.table(selected_path, sep="\t", header=TRUE)

	#select only columns corresponding to gene_id, along with 50 and 500kb
	columns_to_select = which(grepl("gene_id", colnames(full_table), fixed=TRUE) | grepl("50kb", colnames(full_table), fixed=TRUE) | grepl("500kb", colnames(full_table), fixed=TRUE))
		#If ‘TRUE’, ‘pattern’ is a string to be matched as is. This is our case, we are not using regular expressions, just select column names that include an string.

	#select these columns from the full table
	full_table_subset = full_table[, columns_to_select]

	#if the factors is related to transcription factor density (i.e., includes tbfs in any column name), change the column names, because all regulatory factors are named as tbfs
	if(TRUE %in% grepl("tfbs", colnames(full_table_subset))){

		#extract the name of the factor from the path
		raw_split = strsplit(selected_path, split="/|.txt")
			#we use "|" to apply two conditions for the split
		new_colname = raw_split[[1]][length(raw_split[[1]])]
			#select the last element, which is the factor name

		#change the names of the columns for 50 and 500kb
		colnames(full_table_subset)[which(colnames(full_table_subset) == "tfbs_50kb")] = paste(new_colname, "_50kb", sep="")
		colnames(full_table_subset)[which(colnames(full_table_subset) == "tfbs_500kb")] = paste(new_colname, "_500kb", sep="")
	}

	#return the result
	return(full_table_subset)
}

#load all of the tables into a list using the new function
list_factors_across_windows = lapply(list_factors_across_windows_paths, reading_factors_across_windows)

#see the list
str(list_factors_across_windows)

#merge all partitions maintaining all rows
factors_across_windows = Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by="gene_id"), list_factors_across_windows)
	#this can be done with full_join, which is similar to use merge with all=TRUE
		#https://stackoverflow.com/questions/21841146/is-there-an-r-dplyr-method-for-merge-with-all-true
		#https://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
		#https://www.datasciencemadesimple.com/join-in-r-merge-in-r/

#see the table
str(factors_across_windows)
head(factors_across_windows)
summary(factors_across_windows)

#check that we have the correct number of columns
print("################################################")
print("DO WE HAVE THE CORRECT NUMBER OF COLUMNS IN factors_across_windows?")
check_1 = (length(list_factors_across_windows_paths)*2)+1 == ncol(factors_across_windows)
	#the correct number of columns should be the number of paths multiplied by 2 and then summed 1 for the gene_id column
print(check_1)
print("################################################")



##########################################################
########### CONFOUNDING FACTORS CENTER WINDOWS ###########
##########################################################

#get the path of all partitions
list_factors_gene_center_paths = list.files("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/predictors_gene_windows_v2/factors_gene_center", pattern=".txt", full.names=TRUE)

#create a function to read the factor
#selected_path=list_factors_gene_center_paths[1]
reading_factors_gene_center = function(selected_path){

	#read the table in the selected_path
	full_table = read.table(selected_path, sep="\t", header=TRUE)

	#if the factors is related to transcription factor density (i.e., includes tbfs in any column name), change the column names, because all regulatory factors are named as tbfs
	if(TRUE %in% grepl("gene_expression", selected_path, fixed=TRUE)){

		#select only columns corresponding to gene_id, along average gene expression across all tissues, in testis and immune cells
		columns_to_select = which(grepl("gene_id", colnames(full_table), fixed=TRUE) | grepl("mean_expres_all_tissues", colnames(full_table), fixed=TRUE) | grepl("Testis", colnames(full_table), fixed=TRUE) | grepl("Cells...EBV.transformed.lymphocytes", colnames(full_table), fixed=TRUE))
			#If ‘TRUE’, ‘pattern’ is a string to be matched as is. This is our case, we are not using regular expressions, just select column names that include an string.

		#select these columns from the full table
		full_table_subset = full_table[, columns_to_select]

		#change the names of the columns for 50 and 500kb
		colnames(full_table_subset)[which(colnames(full_table_subset) == "Cells...EBV.transformed.lymphocytes")] = "immune_cells"
		colnames(full_table_subset)[which(colnames(full_table_subset) == "Testis")] = "testis"
		colnames(full_table_subset)[which(colnames(full_table_subset) == "mean_expres_all_tissues")] = "average_gene_expression"
	} else { #if not

		#save the table as it is
		full_table_subset = full_table
	}

	#return the result
	return(full_table_subset)
}

#load all of the tables into a list using the new function
list_factors_gene_center = lapply(list_factors_gene_center_paths, reading_factors_gene_center)

#see the list
str(list_factors_gene_center)

#merge all partitions maintaining all rows
factors_gene_center = Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by="gene_id"), list_factors_gene_center)
	#this can be done with full_join, which is similar to use merge with all=TRUE
		#https://stackoverflow.com/questions/21841146/is-there-an-r-dplyr-method-for-merge-with-all-true
		#https://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
		#https://www.datasciencemadesimple.com/join-in-r-merge-in-r/

#see the table
str(factors_gene_center)
head(factors_gene_center)
summary(factors_gene_center)

#check that we have the correct number of columns
print("################################################")
print("DO WE HAVE THE CORRECT NUMBER OF COLUMNS IN factors_gene_center?")
check_2 = length(list_factors_gene_center_paths)+1+2 == ncol(factors_gene_center)
	#the correct number of columns should be the number of paths plus 1 (gene_id column) plus 2 (gene expression in testis and immune cells)
print(check_2)
print("################################################")



###################################################
########### MERGE INTO ONE SINGLE TABLE ###########
###################################################

#merge maintaining all the rows
factors_table = full_join(factors_across_windows, factors_gene_center, by="gene_id")
	#this can be done with full_join, which is similar to use merge with all=TRUE
		#https://stackoverflow.com/questions/21841146/is-there-an-r-dplyr-method-for-merge-with-all-true
		#https://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
		#https://www.datasciencemadesimple.com/join-in-r-merge-in-r/

#see the table
str(factors_table)
head(factors_table)
summary(factors_table)



#####################################################
########### PREPARE THE TABLE TO BE SAVED ###########
#####################################################

#remove NAs
factors_table_final = na.omit(factors_table)
	#We not interested in any gene for which we do not have data across factors If a gene has no data for a factor, it cannot be matched for that factor. Maybe in some comparisons, that gene is not used, but it could be used in others. In order to use the same factors for matching always, we remove any gene with NA for any confounding factor.

#see the table
str(factors_table_final)
head(factors_table_final)

#save
write.table(factors_table_final, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/factors_table.txt", row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE)
	#separated with space, avoid column and row names, and remove quotes from the gene IDs to match the format used by David in his pipeline.