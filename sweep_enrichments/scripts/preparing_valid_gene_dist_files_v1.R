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



#################################################################
####################### REQUIRED PACKAGES #######################
#################################################################

require(plyr) #for apply functions across lists and data.frames.



###############################
######## VALID FILE ###########
###############################

#valid_file: Creation of file with the genes used for the test. The file contains the Ensembl IDs of the genes that are going to be used by the pipeline. One line in the file = one Ensembl gene ID. Use your own code to generate this file. In the provided folder the example file is valid_file.txt.

#we are going to use our list of genes with coordinates to get the list of gene ids used in the bootstrap. Some of these genes have no confounding factors data, but this is ok. The valid_file example of David also has IDs that are not present in Factors_table and iHS ranks, and vice versa.

#load our curate list of ensembl genes (hg19)
gene_coords = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", sep="\t", header=TRUE)
str(gene_coords) #this dataset include our curated list of genes with symbol and ensembl ID. It is well curated and has not repeated ensemble IDs.

#remove duplicated
gene_coords_no_duplicated = gene_coords[which(!duplicated(gene_coords$gene_id)),]
	#remember that in gene_coords we have several rows (i.e., exons) for the same gene, so the gene IDs are repeated.

#extract the unique ensemble IDs
valid_file = unique(gene_coords_no_duplicated$gene_id)

#save
write.table(valid_file, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/valid_file.txt", sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by space to follow the exact same format than David used.



###############################
######## COORD FILE ###########
###############################

#create also a file with ensemble coordinates

#from our curated list of genes select only the columns used by david
ensemble_gene_coords_v99 = gene_coords_no_duplicated[, c("chromosome_name", "gene_id", "gene_start", "gene_end")]
	#we are not using the strand
		#it is needed to include the strand? the variables about gene start and end I used referred only to the forward strand (see gene_coordinates_v10.r; line 1677)
		#CHECK

#reorder columns
ensemble_gene_coords_v99 = ensemble_gene_coords_v99[, c("gene_id", "chromosome_name", "gene_start", "gene_end")]

#match names in the file of David
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "chromosome_name")] = "Chromosome Name"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_id")] = "Ensembl Gene ID"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_start")] = "Gene Start (bp)"
colnames(ensemble_gene_coords_v99)[which(colnames(ensemble_gene_coords_v99) == "gene_end")] = "Gene End (bp)"

#save
write.table(ensemble_gene_coords_v99, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/ensemble_gene_coords_v99.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	#without quotes and separated by tab to follow the exact same format than David used.



##################################
######## GENE SET FILE ###########
##################################

#create the file that defines the gene set of interest. Genes of interest (VIPs, disease genes, etc.) in the first column are flagged with “yes” in the second column. The other genes are flagged with “no” in the second column. Only the genes also present in the file of used genes (Step 1) will be used by the pipeline. In the folder the example file is genes_set_file.txt.

#load the list of metabolic genes
metabolic_gene_list = read.table("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating/data/metabolic_genes/metabolic_gene_list.txt.gz", sep="\t", header=TRUE)
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
write.table(genes_set_file, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/genes_set_file.txt", row.names=FALSE, col.names=FALSE, sep=" ", quote = FALSE)
	#separated with space, avoid column and row names, and remove quotes from the gene IDs to match the format used by David in his pipeline.



##################################
######## GENE DIST FILE ##########
##################################

#Step 3: for every gene, compute a file with the distance of every gene from the closest gene of interest (using gene genomic centers as reference points). If a gene is a gene of interest, then the distance is zero. This is done to be able to choose control genes far enough from the genes of interest. The example file is distance_file.txt


#selected_gene=genes_set_file[94,]
dist_close_met = function(selected_gene){

	#if the selected gene is NOT a metabolic gene
	if(selected_gene$metabolic_gene == "no"){

		#select the genomic center of the selected gene, we use our curated list of gene coordinates
		selected_gene_center = gene_coords_no_duplicated[which(gene_coords_no_duplicated$gene_id == selected_gene$gene_id),]$middle_point

		#select the genomic center of the metabolic genes, we use our curated list of gene coordinates
		metabolic_genes_center = gene_coords_no_duplicated[which(gene_coords_no_duplicated$gene_id %in%  genes_set_file[which(genes_set_file$metabolic_gene=="yes"),]$gene_id),]$middle_point

		#check we have all metabolic genes
		check_1 = length(metabolic_genes_center) == nrow(metabolic_gene_list)

		#calculate the distance en absolute value between the center of the selected gene and the rest of genes
		distances = abs(selected_gene_center - metabolic_genes_center)

		#select the min distance between the selected gene and the metabolic genes
		min_dist_met_genes = distances[which(distances == min(distances))]
	}

	#if the selected gene is a metabolic gene
	if(selected_gene$metabolic_gene == "yes"){

		#check_1
		check_1=NA

		#then the distance is zero
		min_dist_met_genes=0
	}

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
write.table(distance_file, "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/distance_file.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
	#separated with tabs, avoid column and row names, and remove quotes from the gene IDs to match the format used by David in his pipeline.