#!/usr/bin/env Rscript

#This is done to have the possibility to run this script as an executable: 'chmod +x myscript.R' and then ' ./myscript.R'. If you run the script as 'R CMD BATCH myscript.R', i THINK this is not used, because it is annotated. 
	#https://www.jonzelner.net/statistics/make/docker/reproducibility/2016/05/31/script-is-a-program/

#In case you run this script as an executable, you can save the output without warnings "./myscript.R > myscript.Rout" or with errors "./myscript.R &> myscript.Rout"
	#https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file



##########################################################################################
####################### OBTAIN GENE LIST FROM MSIGDB #####################################
##########################################################################################

valid_file = read.table("/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/valid_file.txt", sep="\t", header=FALSE)
colnames(valid_file)[which(colnames(valid_file) == "V1")] = "ensembl_gene"
str(valid_file)
head(valid_file)

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