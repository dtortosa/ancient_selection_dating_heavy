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