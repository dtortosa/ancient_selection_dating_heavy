#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



####################################################################################################
########################### PREPARE VALID, GENE LIST AND DISTANCIE #################################
####################################################################################################

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

require(plyr) #for apply functions across lists and data.frames.



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



###############################
######## GENE COORD ###########
###############################

#we are going to use our curated list of gene coordinates for obtaining many of the files we need for the pipeline to run.

#load our curate list of ensembl genes (hg19)
gene_coords = read.table(paste(path_outside_data, "/gene_number_cds_coords.txt", sep=""), sep="\t", header=TRUE)
str(gene_coords) 
head(gene_coords) #this dataset include our curated list of genes with symbol and ensembl ID. It is well curated and has not repeated ensemble IDs.

#remove duplicated
gene_coords_no_duplicated = gene_coords[which(!duplicated(gene_coords$gene_id)),]
	#remember that in gene_coords we have several rows (i.e., exons) for the same gene, so the gene IDs are repeated.
str(gene_coords_no_duplicated) 
head(gene_coords_no_duplicated)

#check there are no duplicated ids in the gene coord file processed
print("##########################################################################")
print("CHECK THERE ARE NO DUPLICATED IDS IN THE GENE COORD FILE PROCESSED")
print(length(which(duplicated(gene_coords_no_duplicated$gene_id))) == 0)
print("##########################################################################")



###############################
######## VALID FILE ###########
###############################

#valid_file: Creation of file with the genes used for the test. The file contains the Ensembl IDs of the genes that are going to be used by the pipeline. One line in the file = one Ensembl gene ID. Use your own code to generate this file. In the provided folder the example file is valid_file.txt.

#we are going to use our list of genes with coordinates to get the list of gene ids used in the bootstrap. Some of these genes have no confounding factors data, but this is ok. The valid_file example of David also has IDs that are not present in factors_table and iHS ranks, and vice versa.

#extract the unique ensemble IDs
valid_file = gene_coords_no_duplicated$gene_id

#save
write.table(valid_file, paste(path_outside_pipeline, "/valid_file.txt", sep=""), sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by space to follow the exact same format than David used.



##############################
######## HGNC FILE ###########
##############################

#HGNC_file: file with HGNC names for Ensembl genes to be able to exclude HLA (extreme outliers that can bias results) and histone genes (notorious for abundant gene conversion). Replace this file with your own file if needed. This file has two columns: "Ensembl Gene ID" and "HGNC symbol".

#we are going to use our curated list of genes coordinates to get the gene id plus hgnc symbol. Our valid_files, i.e., the list of ensembl IDs considered in the pipeline, comes from here, so it makes sense to use it.

#from our curated list of genes select only the columns used by david
ensembl_hgnc_file = gene_coords_no_duplicated[, c("gene_id", "hgnc_symbol")]

#check the number of ensembl IDs without gene symbol
print("##########################################################################")
print("CHECK THE NUMBER OF ENSEMBL IDS WITHOUT GENE SYMBOL")
print(length(which(ensembl_hgnc_file$hgnc_symbol == "")))
print("##########################################################################")
	#From 19252 genes, 947 have no hgnc gene symbol. This is ok, because the example file of David has 1488 genes without gene symbol from 22080 genes: (1488/22080)*100=6.74, while (947/19252)*100=4.91. Therefore, we do not have a greater proportion of genes without gene symbol in our list.

#match names in the file of David
colnames(ensembl_hgnc_file)[which(colnames(ensembl_hgnc_file) == "gene_id")] = "Ensembl Gene ID"
colnames(ensembl_hgnc_file)[which(colnames(ensembl_hgnc_file) == "hgnc_symbol")] = "HGNC symbol"

#see the table
str(ensembl_hgnc_file)
head(ensembl_hgnc_file)

#save
write.table(ensembl_hgnc_file, paste(path_outside_pipeline, "/ensembl_hgnc_file.txt", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by tab to follow the exact same format than David used.



###############################
######## COORD FILE ###########
###############################

#create also a file with ensemble coordinates.

#from our curated list of genes select only the columns used by david
ensemble_gene_coords_v99 = gene_coords_no_duplicated[, c("gene_id", "chromosome_name", "gene_start", "gene_end")]
	#we are not using the strand
		#it is needed to include the strand? the variables about gene start and end I used referred only to the forward strand (see gene_coordinates_v10.r; line 1677)
		#CHECK

#match names in the file of David
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_id")] = "Ensembl Gene ID"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "chromosome_name")] = "Chromosome Name"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_start")] = "Gene Start (bp)"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_end")] = "Gene End (bp)"

#see the table
str(ensemble_gene_coords_v99)
head(ensemble_gene_coords_v99)

#save
write.table(ensemble_gene_coords_v99, paste(path_outside_pipeline, "/ensembl_gene_coords_v99.txt", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by tab to follow the exact same format than David used.



##################################
######## GENE SET FILE ###########
##################################

#create the file that defines the gene set of interest. Genes of interest (VIPs, disease genes, etc.) in the first column are flagged with “yes” in the second column. The other genes are flagged with “no” in the second column. Only the genes also present in the file of used genes (Step 1) will be used by the pipeline. In the folder the example file is genes_set_file.txt.

#load the list of metabolic genes
metabolic_gene_list = read.table(paste(path_outside_data, "/metabolic_gene_list.txt.gz", sep=""), sep="\t", header=TRUE)
str(metabolic_gene_list)
head(metabolic_gene_list)

#create a data.frame with the IDs used for valid file and then add a second column for indicating whether or not a gene belongs to our list
genes_set_file = data.frame(gene_id=valid_file, metabolic_gene="no")

#add "yes" as another level for metabolic_genes variable
levels(genes_set_file$metabolic_gene) = c("no", "yes")

#those genes with an ID included in the metabolic list are considered as metabolic genes
genes_set_file[which(genes_set_file$gene_id %in% metabolic_gene_list$ensembl_gene),]$metabolic_gene = "yes"

#check that genes with "no" for metabolic_gene are not present in the list of metabolic genes
print("##########################################################################")
print("CHECK THAT GENES WITH 'NO' FOR METABOLIC_GENE ARE NOT PRESENT IN THE LIST OF METABOLIC GENES")
length(which(genes_set_file[which(genes_set_file$metabolic_gene == "no"),]$gene_id %in% metabolic_gene_list$ensembl_gene)) == 0
print("##########################################################################")

#check that genes with "yes" for metabolic_gene are those present in the list of metabolic genes
print("##########################################################################")
print("CHECK THAT GENES WITH 'NO' FOR METABOLIC_GENE ARE NOT PRESENT IN THE LIST OF METABOLIC GENES")
length(which(!genes_set_file[which(genes_set_file$metabolic_gene == "yes"),]$gene_id %in% metabolic_gene_list$ensembl_gene)) == 0
print("##########################################################################")

#check we have all metabolic genes
print("##########################################################################")
print("CHECK WE HAVE ALL METABOLIC GENES")
length(which(genes_set_file$metabolic_gene == "yes")) == nrow(metabolic_gene_list)
print("##########################################################################")

#see the table
str(genes_set_file)
head(genes_set_file)
summary(genes_set_file)

#save
write.table(genes_set_file, paste(path_outside_pipeline, "/genes_set_file.txt", sep=""), row.names=FALSE, col.names=FALSE, sep=" ", quote = FALSE)
	#separated with space, avoid column and row names, and remove quotes from the gene IDs to match the format used by David in his pipeline.



##################################
######## GENE DIST FILE ##########
##################################

#Step 3: compute a file with the distance of every gene from the closest gene of interest (using gene genomic centers as reference points). If a gene is a gene of interest, then the distance is zero. This is done to be able to choose control genes far enough from the genes of interest. The example file is distance_file.txt


##write a function to do so
#for debugging
#selected_gene=genes_set_file[1,] #non-metabolic gene
#selected_gene=genes_set_file[94,] #metabolic gene
dist_close_met = function(selected_gene){

	#if the selected gene is NOT a metabolic gene
	if(selected_gene$metabolic_gene == "no"){

		#select the genomic center of the selected gene, we use our curated list of gene coordinates
		selected_gene_center = gene_coords_no_duplicated[which(gene_coords_no_duplicated$gene_id == selected_gene$gene_id),]$middle_point

		#select the genomic center of the metabolic genes, we use our curated list of gene coordinates
		metabolic_genes_center = gene_coords_no_duplicated[which(gene_coords_no_duplicated$gene_id %in%  genes_set_file[which(genes_set_file$metabolic_gene=="yes"),]$gene_id),]$middle_point

		#check we have all metabolic genes and just one position for the selected gene
		check_1 = length(metabolic_genes_center) == nrow(metabolic_gene_list) & length(selected_gene_center) == 1

		#calculate the distance en absolute value between the center of the selected gene and the rest of genes
		distances = abs(selected_gene_center - metabolic_genes_center)
			#we are not interested in the sense of the difference, but just the distance. We do not care if a metabolic gene is before or after a control gene, but the number of bases between them.

		#select the min distance between the selected gene and the metabolic genes
		min_dist_met_genes = min(distances)
	}

	#if the selected gene is a metabolic gene
	if(selected_gene$metabolic_gene == "yes"){

		#check_1
		check_1=NA

		#then the distance is zero
		min_dist_met_genes=0
	}

	#bind the results
	final_results = data.frame(selected_gene, min_dist_met_genes=min_dist_met_genes, check_1=check_1)

	#return the results
	return(final_results)
}

#apply the function
distance_file = ddply(.data=genes_set_file, .variables="gene_id", .fun=dist_close_met, .inform=TRUE, .parallel=FALSE, .paropts=NULL)
	#".inform=TRUE" generates and shows the errors. This increases the computation time, BUT is very useful to detect problems in your analyses.
	#".parallel" to paralelize with foreach. 
	#".paropts" is used to indicate additional arguments in for each, specially interesting for using the .export and .packages arguments to supply them so that all cluster nodes have the correct environment set up for computing. 
		#ADD PACKAGES USED INSIDE THE FUNCTION

#check we have calculated the distance for all genes
print("###################################################")
print("CHECK WE HAVE CALCULATED THE DISTANCE FOR ALL GENES")
nrow(distance_file) == nrow(genes_set_file)
print("###################################################")

#check we have all genes from valid file
print("###################################################")
print("CHECK WE HAVE ALL GENES FROM VALID FILE")
!FALSE %in% c(valid_file %in% distance_file$gene_id)
!FALSE %in% c(distance_file$gene_id %in% valid_file)
print("###################################################")

#check metabolic genes have zero distance and na for check_1
print("###################################################")
print("CHECK METABOLIC GENES HAVE ZERO DISTANCE AND NA FOR CHECK_1")
!FALSE %in% c(distance_file[which(distance_file$metabolic_gene=="yes"),]$min_dist_met_genes == 0)
!FALSE %in% c(is.na(distance_file[which(distance_file$metabolic_gene=="yes"),]$check_1))
print("###################################################")

#check metabolic genes have zero distance and na for check_1
print("###################################################")
print("CHECK NON-METABOLIC GENES HAVE DISTANCE AND NO NA FOR CHECK_1")
!FALSE %in% c(distance_file[which(distance_file$metabolic_gene=="no"),]$min_dist_met_genes != 0)
!c(FALSE, NA) %in% c(distance_file[which(distance_file$metabolic_gene=="no"),]$check_1)
print("###################################################")

#see the table
summary(distance_file)

#select only the gene_id and the distance columns
distance_file = distance_file[,which(colnames(distance_file) %in% c("gene_id", "min_dist_met_genes"))]

#save
write.table(distance_file, paste(path_outside_pipeline, "/distance_file.txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
	#separated with tabs, avoid column and row names, and remove quotes from the gene IDs to match the format used by David in his pipeline.