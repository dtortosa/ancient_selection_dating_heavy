#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



##########################################################################################
############ PREPARE VALID, GENE AND GENE DISTANCE FILES #################################
##########################################################################################

#IMPORTANT: 
	#Note that this script was wrote and run in 2022 using the version of April 28th 2022 of David's pipeline.
		#The pipeline of David has been downloaded from 
			#https://github.com/DavidPierreEnard/Gene_Set_Enrichment_Pipeline
			#last commit done on Apr 28, 2021
				#7b755c0c23dd4d7c3f54c4b53e74366e4041ac8f 

	#Manual from David
		#/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/exdef_pipeline_manual.pdf

#This script create three needed inputs for the pipeline of David:
	#valid_file: file with all the ensemble gene used
	#genes_set_file: file with the interest genes
	#distance_file.txt: file with the distance of each gene to the closest interest gene



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



########################################
######## CONFOUDNING FACTORS ###########
########################################

#Factors_table.txt: create a table file that contains all the confounding factors that the bootstrap test will control for. First column is Ensembl gene ID, next columns are each factor that will be matched. The example file is Factors_table.txt


#RUUUN SEPARATED WINDOW SIZE AND THE OTHER

#get the path of all partitions
list_factors_across_windows_paths = list.files("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/data/predictors_gene_windows_v2/factors_across_windows", pattern=".txt", full.names=TRUE)


#path=list_factors_across_windows_paths[1]
reading_function = function(path){

	full_table = read.table(path, sep="\t", header=TRUE)

	columns_to_select = which(grepl("gene_id", colnames(full_table)) | grepl("50kb", colnames(full_table)) | grepl("500kb", colnames(full_table)))

	full_table_subset = full_table[, columns_to_select]

	if(TRUE %in% c(grepl("tfbs", colnames(full_table_subset)))){

		raw_split = strsplit(path, split="/|.txt") #USA DOS CONDICIONES! LA RUTA Y TXT
		new_colname = raw_split[[1]][length(raw_split[[1]])]

		colnames(full_table_subset)[which(colnames(full_table_subset) == "tfbs_50kb")] = paste(new_colname, "_50kb", sep="")
		colnames(full_table_subset)[which(colnames(full_table_subset) == "tfbs_500kb")] = paste(new_colname, "_500kb", sep="")

	}

	return(full_table_subset)
}

#load all of them faster with fread
list_factors_across_windows = lapply(list_factors_across_windows_paths, reading_function)

#merge all partitions maintaining all rows
factors_across_windows = Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by="gene_id"), list_factors_across_windows)
	#this can be done with inner_join, which is similar to use merge with all=FALSE.
		#For ‘full_join()’, all ‘x’ rows, followed by unmatched ‘y’ rows.
		#https://stackoverflow.com/questions/21841146/is-there-an-r-dplyr-method-for-merge-with-all-true
		#https://stackoverflow.com/questions/28250948/how-to-dplyrinner-join-multi-tbls-or-data-frames-in-r
		#https://www.datasciencemadesimple.com/join-in-r-merge-in-r/
str(factors_across_windows)



###############################
######## VALID FILE ###########
###############################

#valid_file: Creation of file with the genes used for the test. The file contains the Ensembl IDs of the genes that are going to be used by the pipeline. One line in the file = one Ensembl gene ID. Use your own code to generate this file. In the provided folder the example file is valid_file.txt.

#load our curate list of ensembl genes
gene_coords = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", sep="\t", header=TRUE)
str(gene_coords) #this dataset include our curated list of genes with symbol and ensembl ID. It is well curated and has not repeated ensemble IDs.

#extract the unique ensemble IDs
unique_ensembl_id = unique(gene_coords$gene_id)

#save
valid_file = write.table(unique_ensembl_id, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/valid_file.txt", sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by space to follow the exact same format than David used.





genes_set_file = read.table("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/original_files/exdef_folder/genes_set_file.txt", sep=" ", header=FALSE)
colnames(genes_set_file)[which(colnames(genes_set_file) == "V1")] = "ensembl_gene"
colnames(genes_set_file)[which(colnames(genes_set_file) == "V2")] = "metabolic_gene"
genes_set_file$metabolic_gene = "no"
str(genes_set_file)
head(genes_set_file)


metabolic_gene_list = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating/data/metabolic_genes/metabolic_gene_list.txt.gz", sep="\t", header=TRUE)
str(metabolic_gene_list)
head(metabolic_gene_list)


metabolic_gene_list[which(!metabolic_gene_list$ensembl_gene %in% genes_set_file$ensembl_gene),]
	#three metabolic genes are not included, which are the genes rescued using gene symbol. makes sense because we are matching here by gene id...

genes_set_file[which(genes_set_file$ensembl_gene %in% metabolic_gene_list$ensembl_gene),]$metabolic_gene = "yes"


#check
length(which(genes_set_file[which(genes_set_file$metabolic_gene == "no"),]$ensembl_gene %in% metabolic_gene_list$ensembl_gene)) == 0


write.table(genes_set_file, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/genes_set_file.txt", row.names=FALSE, col.names=FALSE, sep=" ", quote = FALSE)
	#https://stackoverflow.com/questions/14846433/avoid-quotation-marks-in-column-and-row-names-when-using-write-table