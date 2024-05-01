#!/usr/bin/env Rscript
# coding: utf-8
#to run this script: chmod +x script.R; ./script.R
#!/bin/sh does not work with my terminal en msi of David.
#if you are using "$" to paste the path of the executable,
#you do not need to use "./" for running the executable.
#you can save the output and the errors
#./script.R > script.Rout #only output #nolint
#./script.R 2> error.Rout #only error
#./script.R > script.Rout 2> error.out #both in different files
#./script.R > script.Rout 2>&1 #both in the same file
#https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/



#########################################
######## OBTAIN GENE COORDINATES ########
#########################################

#From this script we will get the center of each
#coding gene along with several genomic factors
#This gene center will be used in other scripts to get the data
#about different confounding genomic



###################################################
###### CHANGES RESPECT TO PREVIOUS VERSIONS #######
###################################################

#This script comes from gene_coordinates_v10.r of the MDR paper



#################################################
###### MOVE FROM CONTAINER TO MAIN FOLDER #######
#################################################

#we move the WD to the parent folder
setwd("../")
getwd()



######################
###### LIBRARY #######
######################

#require package
library(biomaRt) #for connecting with ensemble
library(GenomicRanges) #for calculating overlap between exons of a gene
require(stringr) #for replacing "," in the total chromosome length




##################################################################################
###### SEARCH AND SELECTION OF THE MART AND DATABASE (SPECIES AND VERSION) #######
##################################################################################

#before run the final analyses we check if we have something in the cache.
#To save time and computing resources biomaRt will attempt to identify when you are re-running a query you have executed before. Each time a new query is run, the results will be saved to a cache on your computer. If a query is identified as having been run previously, rather than submitting the query to the server, the results will be loaded from the cache. You can get some information on the size and location of the cache using the function biomartCacheInfo():
biomartCacheInfo()

#The cache can be deleted using the command biomartCacheClear(). This will remove all cached files.
#biomartCacheClear()
#biomartCacheInfo()

#look for the list of marts currently included. At the time of writing this code (April 8th 2024), everything is version 111, which is the latest version released in Jan 2024.
listMarts()

#connect with the BioMart database we are interested. In this case, ensemble genes
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")

#see the list of datasets included in the database 'ensemble'
datasets <- listDatasets(ensembl)
head(datasets)

#search for the name of your dataset of interest, which should include 'hsapiens' in its name
searchDatasets(mart = ensembl, pattern = "hsapiens")

#We are interested in hsapiens_gene_ensembl
datasets[which(datasets$dataset == "hsapiens_gene_ensembl"), ]
#the version indicated is GRCh38.p14, which is the version we need, i.e., hg38.

#we want to view the available archived versions. listEnsemblArchives takes no arguments, and produces a table containing the names of the available archived versions, the date they were first available, and the URL where they can be accessed.
ensemble_archives <- listEnsemblArchives()
ensemble_archives

#we select the last version at the moment of writting this code, which is the Ensembl Release 111 (January 2024), this version has "*" as value in the column current_release, which means it is the current release, because the rest of releases have an empty value.
#https://jan2024.archive.ensembl.org
current_version <- ensemble_archives[ which(ensemble_archives$current_release == "*"), "version"]
if(current_version != "111") {
	stop("The current version is not 111. You need to change the host to 'https://jan2024.archive.ensembl.org' so we can get version 111 and get reproducible results")
}


#take a look for the marts
listMarts(host = "https://www.ensembl.org")
#REMEMBER: If you run this and the current version is no longer 111, then you have to change the mirror/host to https://jan2024.archive.ensembl.org
#we see ensemble genes, variation and regulation. In all cases release 111, which is the lastest in general (this is not the same than patch release, the last of which is 14 in GRCh38).

#extract the mart of ensemble genes for hg38, release Jan 2024
grch38_datasets <- listDatasets(useMart(host = "https://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL"))
#REMEMBER: If you run this and the current version is no longer 111, then you have to change the mirror/host to https://jan2024.archive.ensembl.org
head(grch38_datasets)

#We are interested in 'hsapiens_gene_ensembl' ()Human genes. 
grch38_datasets[which(grch38_datasets$dataset == "hsapiens_gene_ensembl"), ]
#The version indicated is GRCh38.p14, which is the last patch release of GRCh38 in the first quarter of 2024.
#this patch seems to be the latest in ncbi
#https://www.ncbi.nlm.nih.gov/assembly/?term=GRCh38

#We are going to use GCRh38.p14 (release of Jan 2024; see above).
#it is useful to use the specific URL for avoiding changes in the results if you compile the code again in several months and the current version (default) has changed. We select the url of the specific version we are interested, which is the release of feb 2014, last patch release of GRCh37

#Select the mart with ensemble genes of humans and the GRCh38.p14, Jan 2024 release (111)
grch38_human <- useMart(host = "https://www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#REMEMBER: If you run this and the current version is no longer 111, then you have to change the mirror/host to https://jan2024.archive.ensembl.org
str(grch38_human)

#check we have the correct assembly and patch in the human dataset
list_datasets_hg38 <- listDatasets(grch38_human)
if (list_datasets_hg38[which(list_datasets_hg38$dataset == "hsapiens_gene_ensembl"), "version"] != "GRCh38.p14") {
	stop("The version of the assembly and patch is not GRCh38.p14")
}

#See the attributes
grch38_human_attributes <- listAttributes(grch38_human)

#see for attributes related to exons
searchAttributes(mart = grch38_human, pattern = c("exon"))
#ensembl_exon_id is the only attribute within the list of feature_page type that is related to exons.
#exon_chrom_start and exon_chrom_end are attributes for sequences. You have to use getSequence with the id of a gene to get the beinign and end of the exons. For example: getSequence(id="ENSG00000198793",type="ensembl_gene_id", seqType="gene_exon",mart=ensembl). I think that adding the id of a gene in the argument value of getBM will alos work
#exon_chrom_start  and exon_chrom_end could alos work, but they are attributes for pages type structure.
#YOU CANNOT ask for atribbute of type 'pages' when you are filtering by a gene and a the same time get info of the exons.., so for most of the exon characteristics (except ensembl_exon_id), you must filter by gene. Previously you will have to get gen id in a different query.

#see for attributes related to coding density length
searchAttributes(mart = grch38_human, pattern = c("cds"))
searchAttributes(mart = grch38_human, pattern = c('genomic'))

#see for attributes related to the sense of the strand
searchAttributes(mart = grch38_human, pattern = c('strand')) 
#we have attraibutes for the stran in sequences and feature pages

#see for attributes related to entrez
searchAttributes(mart = grch38_human, pattern = c("entrez"))

#See what filter we can apply
grch38_human_filters <- listFilters(grch38_human)

#see for filters related to exons
searchFilters(mart = grch38_human, pattern = c("exon"))
#Only ensembl_exon_id Exon ID(s)

#see the values of one the filter we will apply, the chromosome name
listFilterOptions(mart = grch38_human, filter = "chromosome_name") 
#we are only interested in 'normal' data within the 1:22 chromosomes
#we will use only chromsome 1:22 removing all the highly variable regions like those related to the MHC. In these regions the chromosome name is not the actual chromosome where the sequence is, but a different notation for example if it is a highly varaible region with several versions included in the esemble

#extract the ids of ALL human genes. This will be used to split the query of all data into each gene (see below)
all_gene_ids_grch38_human <- getBM(attributes = c("ensembl_gene_id"), mart = grch38_human)




######################################################
###### CALCUALTE THE LENGTH OF EACH CHROMOSOME #######
######################################################

#This will be used to stop the windows that reach the start or end of the chromosomes.

##data from ncbi
#load chromosome length for GRCh38.p14 copied into a txt file from "https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38.p14"
chrom_length_ncbi_hg38 <- read.table("./data/gene_coordinates/chromosome_length/chromosome_length_hg38_from_ncbi.csv", sep = ",", header = TRUE) #chromosome lengths are calculated by summing the length of the placed scaffolds and estimated gaps.

#change col names
colnames(chrom_length_ncbi_hg38) <- c("chromosome", "length_bp", "GenBank.accession", "RefSeq.accession")

#see structure
str(chrom_length_ncbi_hg38)

#select only the autosomal chromosomes
chrom_length_ncbi_hg38 <- chrom_length_ncbi_hg38[which(chrom_length_ncbi_hg38$chromosome %in% 1:22),]

#check that the order of the chromosomes are ok
if (FALSE %in% c(chrom_length_ncbi_hg38$chromosome == 1:22)) {
	stop("The type of the chromosomes is not correct")
}

##data from uscs
#load chromosome length for GRCh38 as the chromInfo.txt file. This directory contains a dump of the UCSC genome annotation database for the Dec. 2013 (GRCh38/hg38) assembly of the human genome (hg38, GRCh38 Genome Reference Consortium Human Reference 38
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
chrom_length_ucsc_hg38 <- read.table("./data/gene_coordinates/chromosome_length/chromInfo.txt.gz", sep = "\t", header = FALSE)

#change colnames
colnames(chrom_length_ucsc_hg38) <- c("chromosome", "length_bp", "path_not_sure")

#select only autosomal
chrom_length_ucsc_hg38 <- chrom_length_ucsc_hg38[which(chrom_length_ucsc_hg38$chromosome %in% paste("chr", 1:22, sep = "")), ]

#see structure
str(chrom_length_ucsc_hg38)

#reorder the rows to have the chromosome number in increasing order
chrom_length_ucsc_hg38 <- chrom_length_ucsc_hg38[match(paste("chr", 1:22, sep = ""), chrom_length_ucsc_hg38$chromosome), ]

#check rows order
if (FALSE %in% c(chrom_length_ucsc_hg38$chromosome == paste("chr", 1:22, sep = ""))) {
	stop("The type of the chromosomes is not correct")
}


##data from uscs golden path patch 14
#data downloaded using this command
#cd ./data/gene_coordinates/chromosome_length/; wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p14/hg38.p14.chrom.sizes
#this is the last patch according to ncbi, released in feb 2022, while the chromosome file was relesaed in oct 2022. The information in the folder (hg38/bigZips) says we have there data for hg38 humans, including the different patches. 
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

#load the data
chrom_length_ucsc_hg38_14 <- read.table(
	"./data/gene_coordinates/chromosome_length/hg38.p14.chrom.sizes",
	header=FALSE,
	sep="\t"
)
colnames(chrom_length_ucsc_hg38_14) <- c("chromosome", "length_bp")
head(chrom_length_ucsc_hg38_14)
str(chrom_length_ucsc_hg38_14)

#select autosomals
chrom_length_ucsc_hg38_14 <- chrom_length_ucsc_hg38_14[which(chrom_length_ucsc_hg38_14$chromosome %in% paste("chr", 1:22, sep = "")), ]

#reorder by chromosome name
chrom_length_ucsc_hg38_14 <- chrom_length_ucsc_hg38_14[match(paste("chr", 1:22, sep = ""), chrom_length_ucsc_hg38_14$chromosome), ]


## compare the three sources
#calculate the difference of length between sources
ncbi_vs_uscs <- chrom_length_ncbi_hg38$length_bp - chrom_length_ucsc_hg38$length_bp
ncbi_vs_uscs_14 <- chrom_length_ncbi_hg38$length_bp - chrom_length_ucsc_hg38_14$length_bp

#check that three sources have the same chromosome length
if (FALSE %in% c(ncbi_vs_uscs == 0) || FALSE %in% c(ncbi_vs_uscs_14 == 0)) {
	stop("The chromosome length of the two sources is not the same")
} else {
	print("The chromosome length of the all sources is the same")
}

##save the chromosome length from ucsc
write.table(
	chrom_length_ucsc_hg38_14,
	"./data/gene_coordinates/chromosome_length/chrom_length_final_v1.txt",
	col.names = TRUE,
	row.names = FALSE,
	sep = "\t"
)
#I have checked that these lengths match those showed in ensembl hg38
#https://www.ensembl.org/Homo_sapiens/Location/Chromosome?r=1%3A1-1000
#https://www.ensembl.org/Homo_sapiens/Location/Chromosome?r=2%3A1-1000
#https://www.ensembl.org/Homo_sapiens/Location/Chromosome?r=3%3A1-1000
#and so on...



#######################################################################
###### CHECK BEGIN/END GENES REFERS TO REGIONS FLANKED BY EXONS #######
#######################################################################

#check that the beginning and end of genes refers to the beginning and end of the region with exons. You have to bear in mind that by definition, an intron is any nucleotide sequence within a gene that is removed by RNA splicing before translation. The word intro is derived from the them intragenic, i.e., a region inside a gene, thus I guess by definition, an intron must be between exons of a gene. Therefore, I don´t think that the end or start of the gene in ensemble would be an intronic region not flanked by exons of the corresponding gene. You can have non-coding sequences but inside a coding exon, like the 5´ UTR (see figure 11; the 3´ and 5´ below of the line i think that refer to the other strand). Said this, we are going to check it for two genes.

#load the chromosome name, gene id, transcript id, hgnc symbol (gene name), start and end position for genes in chromosome 1. We only use chromosome 1 because it is faster
chr_1_positions <- getBM(
    attributes = c("chromosome_name", "ensembl_gene_id", "hgnc_symbol","start_position", "end_position", "ensembl_transcript_id", "ensembl_exon_id"),
    mart = grch38_human,
    filter = "chromosome_name",
    values = "1"
)
colnames(chr_1_positions)
head(chr_1_positions)
tail(chr_1_positions)
nrow(chr_1_positions)

#SEMA4A
#manual checking with the ensembl browser
chr_1_positions[which(chr_1_positions$ensembl_gene_id == "ENSG00000196189"), ]
	#We check that the start (156147366) and end (156177752) correspond to the start and end in the GRCh38.p14 webpage ('https://jan2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000196189;r=1:156147366-156177752'), and that´s exactly right. The location indicated in the webpage is 'Chromosome 1: 156,147,366-156,177,752'

	#In figures 3,4 and 6 ('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego') you can see as the first exon of all transcripts matches the start position and the end position matches the last exon of all transcripts (figure 5).

#extract the structural info (structure characteristics) of the exons included in SEMA4A (ENSG00000196189) to make an additonal check
exons_sema4a <- getBM(
    attributes=c("exon_chrom_start", "exon_chrom_end", "is_constitutive", "ensembl_exon_id"),
    mart = grch38_human,
    filter="ensembl_gene_id",
    values="ENSG00000196189")
head(exons_sema4a)

#select the first exon, i.e., smaller pair base of start
start_first_exon_sema4a <- exons_sema4a[which(exons_sema4a$exon_chrom_start == min(exons_sema4a$exon_chrom_start)), ]$exon_chrom_start

#check that this exon begin just in the start of the gene
start_first_exon_sema4a == unique(chr_1_positions[which(chr_1_positions$ensembl_gene_id == "ENSG00000196189"), ]$start_position)

#select the last exon, i.e., bigger pair base of end
end_last_exon_sema4a <- exons_sema4a[which(exons_sema4a$exon_chrom_end == max(exons_sema4a$exon_chrom_end)), ]$exon_chrom_end

#check that this exon begin just in the start of the gene
end_last_exon_sema4a == unique(chr_1_positions[which(chr_1_positions$ensembl_gene_id == "ENSG00000196189"), ]$end_position)

#TMEM183A
#Done in the original gene coordinate script for the MDR paper

#We can conclude that the variables start and end positions correspond with the begining and end of the exons. However, an exon can included non-coding sequences:
#I have found that the start/stop position of each exon sum the transcript length, but this length is bigger than the coding sequence length variable. The transcript size includes non-coding sequences like the 5'UTR extreme, BUT NOT introns (Figure 9). I think we should use start and end positions of exons to calculate gene center, but for density of coding sequences, in that case i think the sum of coding sequences if better.
#David - Coding density. that is very true. You don want to include these sequences when calculating the coding density. You can have an exon, the half of which is translated into a protein, but the other half not (see figure 11; the 3´ and 5´ below of the line i think that refer to the other strand). We can use the attributes CDS length, cds start and cds end of a structure for that.	
#David - gene length: For gene length does not matter. Even include transcripts that are not translated to functional proteins (decay), because this only move some dozens of bases the gene position, while LD blocks in humans have a length of thousands of pair bases.




##############################################
###### LOAD GENE NAMES AND COORDINATES #######
##############################################

#load the chromosome name, gene id, transcript id, hgnc symbol (gene name), start and end position for all genes in ensembl
all_genes_grch38_human <- getBM(
    attributes = c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id", "strand", "ensembl_exon_id", "hgnc_symbol", "start_position", "end_position"),
    mart = grch38_human,
    filters = "ensembl_gene_id",
    values = all_gene_ids_grch38_human
)
#we apply the gene ids as a filter but setting as value the ids of ALL human genes, which were previously downloaded. I do the download in that way because if you download the whole data is an individual query with a wait limit of 5 minutes. This made the query be stopped. In contrast, if we apply the filter, you will download each gene as an independent query with 5 minutes each one BUT in one file. In addition, we get a bar of progress and an estimated time to the query be finished. See more information and other options for solving this problem here: https://github.com/grimbough/biomaRt/issues/20
colnames(all_genes_grch38_human)
head(all_genes_grch38_human)
tail(all_genes_grch38_human)
nrow(all_genes_grch38_human)
#I cannot obtain the start and end position of each gene (feature_page characteristics) and at the same time obtain info of each exon (structure or sequences characteristics) because we have to filter by gene for that. We need this information to calculate some confounding factors like the density of coding sequence density, but I think we will need to obtain the list of genes first, and the make a loop for gene calculating the confounding factors we can obtain.

#check differences between total rows and the number of uniques values of ensembl_gene_id
nrow(all_genes_grch38_human) #1791149
length(unique(all_genes_grch38_human$ensembl_gene_id)) #70711 unique gene id.
#this difference is because for each gene we can have several transcripts, and within each transcript we can have several exons. In addition, note that you can have several gene id for the same hgnc symbol (gene name), I think these are cases with high variability, so different versions of the gene are included, BUT I am not sure. In any case, i think after applying the filter for coding genes, these redundant gene id will be remove.

#remove all the rows that have a chromosome_name different than 1:22. In that way, we remove the exons of genes that are in  mithocondrial DNA or messy sequences like those of the MHC (highly variable regions for which the assemble include several versions, in that cases chromosome name is some like V_MHC...). In addition, we remove the sexual chromosomes.
#We don´t take data from the Y because is very small, it has very few genes, so the estimates of the slope of confounding factors only considering this chromosome would have very wide confidence interval (low sample size).
#We remove X because we do not have iHS for that chromosome right now, and it would entail a great amount of effort to calculate that. In addition, much of the analyses i could do for example with climate would be done within the autosomal only.
unique(all_genes_grch38_human$chromosome_name)
all_genes_grch38_human_filtered <- all_genes_grch38_human[which(all_genes_grch38_human$chromosome_name %in% c(1:22)), ]

#check that there is not row with other chromosome_name rather than 1:22
if(FALSE %in% c(unique(all_genes_grch38_human_filtered$chromosome_name) %in% c(1:22))) {
	stop("The type of the chromosomes is not correct")
}

#check that each gene name has only one gene id. I can not select those cases without gene name
#extract the gene names without empty cases ("")
gene_names_no_null <- unique(all_genes_grch38_human_filtered$hgnc_symbol)
gene_names_no_null <- gene_names_no_null[which(gene_names_no_null != "")]
!"" %in% gene_names_no_null
#for each gene name
test_duplicated_gene_id <- data.frame(selected_gene_name = NA, test_result = NA)
#i=1
for(i in 1:length(gene_names_no_null)){

	#select the [i] gene name
	selected_gene_name <- gene_names_no_null[i]

	#extract the unique cases of gene id for the [i] hgnc symbol
	unique_gene_id <- unique(all_genes_grch38_human_filtered[which(all_genes_grch38_human_filtered$hgnc_symbol == selected_gene_name), ]$ensembl_gene_id)

	#test that the number of unique cases is 1 (only 1 gene id for each gene name)
	test_result <- length(unique_gene_id) == 1

	#save
	test_duplicated_gene_id <- rbind.data.frame(test_duplicated_gene_id, cbind.data.frame(selected_gene_name, test_result))
}
#remove first row with NAs
test_duplicated_gene_id <- test_duplicated_gene_id[-1, ]
#take a look for genes with repeated gene id
summary(test_duplicated_gene_id) #we have several cases
#extract information of these genes
genes_with_duplicated_id <- test_duplicated_gene_id[which(test_duplicated_gene_id$test_result == FALSE), ]$selected_gene_name
all_genes_grch38_human_filtered[which(all_genes_grch38_human_filtered$hgnc_symbol %in% genes_with_duplicated_id), ]
#there are a lot of micro RNA with repeated gene names, but they should be removed after filtering by biotype (check)
#there is one case of non microRNAs: UGT2A1 has two gene ids, ENSG00000270386 (https://jan2024.archive.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000173610;r=4:69589194-69653249;t=ENST00000514019) and ENSG00000173610 (https://jan2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000173610;r=4:70454135-70518965). It is like ENSG00000270386 is a transcript of ENSG00000270386 that has been separated as an independent gene. They are in the same region, and has the same hgcn symbol.
#I do not know what to do with the latter case. Should I removed the gene id with only one transcript, if not we will have the start of two windows in that region. Moreover, I don´t know how detect more cases like this in genes that have no hgcn symbol.



###################################################################
###### LOAD CODING SEQUENCE LENGTH AND PREPARE CALCULATIONS #######
###################################################################

#we are going to load the positions of all exons of the genome
full_exon_data <- getBM(
    attributes = c("chromosome_name", "ensembl_gene_id", "gene_biotype", "ensembl_transcript_id", "transcript_biotype", "strand", "transcript_start", "transcript_end", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "genomic_coding_start", "genomic_coding_end", "cds_length", "cds_start", "cds_end"),
    mart = grch38_human,
    filters = "ensembl_gene_id",
    values = all_gene_ids_grch38_human)
#we apply the gene ids as a filter but setting as value the ids of ALL human genes, which were previously downloaded. I do the download in that way because if you download the whole data is an individual query with a wait limit of 5 minutes. This made the query be stopped. In contrast, if we apply the filter, you will download each gene as an independent query with 5 minutes each one BUT in one file. In addition, we get a bar of progress and an estimated time to the query be finished. See more information and other options for solving this problem here: https://github.com/grimbough/biomaRt/issues/20
head(full_exon_data)
#cds_start and cds_end work fine. For example with ENSG00000196189 (SEMA4A), I have checked the exons in the exon visualizer of the first transcript (ENST00000435124). The first exon is non-coding, so cds start and end is NA. This starts at 156147366 and ends at 156147430. Then an intron starts at 156147431 and ends at 156154549. At 156154550 begins the next exon, first coding. It has as cds start 1, because it is the first one. The first exon have a value of CDS length because it is a global value for the whole gene, but that exon does not contribute to coding length.
sema4a_full_first_transcript <- full_exon_data[which(full_exon_data$ensembl_gene_id == "ENSG00000196189" & full_exon_data$ensembl_transcript_id == "ENST00000435124"), ]
sema4a_full_first_transcript
#BUT the most interesting variables are genomic_coding_start and genomic_coding_end. For example, this first transcript of SEMA4A (ENST00000435124) has 139 coding basis in the first coding exon, as cds_start is 1 and cds_end is 139. The range of sequences included between genomic_coding_start and genomic_coding_end (156154717 - 156154579 + 1) is 139. The same goes for the last coding exon. The range according to genomic coding start and end is also 153 (156161498 - 156161346 + 1), also according to the data of cds_start and end (963-811+1).
#https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000196189;r=1:156147366-156177752;t=ENST00000435124


##check that exons from the whole data and from individual queries are the same
#from full_exon_data extract the data of SEMA4A
sema4a_exon_data_full_query <- full_exon_data[which(full_exon_data$ensembl_gene_id == "ENSG00000196189"), ]

#extract exon data of SEMA4A through and individual query
sema4a_exon_data_indv_query <- getBM(
	attributes = c("ensembl_transcript_id", "transcript_biotype", "strand", "transcript_start", "transcript_end", "exon_chrom_start", "exon_chrom_end", "ensembl_exon_id", "genomic_coding_start", "genomic_coding_end", "cds_length", "cds_start", "cds_end"),
	mart = grch38_human,
	filter = "ensembl_gene_id",
	values = "ENSG00000196189")

#order rows of the full query according to the exon start position
sema4a_exon_data_full_query_ordered <- sema4a_exon_data_full_query[order(sema4a_exon_data_full_query$exon_chrom_start), ]
sort(sema4a_exon_data_full_query_ordered$exon_chrom_start) == sema4a_exon_data_full_query_ordered$exon_chrom_start #check the order is correct

#order rows of the individual query according to the exon start position
sema4a_exon_data_indv_query_ordered <- sema4a_exon_data_indv_query[order(sema4a_exon_data_indv_query$exon_chrom_start), ]
sort(sema4a_exon_data_indv_query_ordered$exon_chrom_start) == sema4a_exon_data_indv_query_ordered$exon_chrom_start

#check that the transcript_id, start/end of the exons and the start/end and length of the coding sequence are the same in both queries
for (column in c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end", "cds_start", "cds_end", "cds_length", "genomic_coding_start", "genomic_coding_end", "strand")) {

	#check equality
	check_columns <- identical(
		sema4a_exon_data_full_query_ordered[, column], sema4a_exon_data_indv_query_ordered[, column]
	)

	#stop if FALSE
	if (check_columns == FALSE) {
		stop(paste("The column", column, "is not the same in both queries"))
	}
}


##check that the exons obtained from sema4a_exon_data_full_query are the same than those showed in the webpage 

#looking in the web
#https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000196189;r=1:156147366-156177752;t=ENST00000435124

#I have checked that the exons of the first transcript (ENST00000435124) have the same end and start in my data than in the webpage. The introns are not included between exon_chrom_start and exon_chrom_end (Figure 9 in '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego').
sema4a_full_first_transcript

#I have also checked the size of the transcript and the coding sequence. If we calculate the difference between the start and the end of each exon plus 1, i.e., calculate the length of each exon, and then sum all of them, we get 1057 for the length of all exons. This is exactly the length showed in the transcript table for this transcript. BUT the coding sequence length is slightly smaller, 963. I guess that transcript length include UTR extreme and other non-coding sequences. You can see these sequences in the server page.
exon_sizes_first_transcript <- (sema4a_full_first_transcript$exon_chrom_end - sema4a_full_first_transcript$exon_chrom_start) + 1 
#calculate the difference between and start of the exons in the first transcript. We sum 1 because we want to include the first and last base. If you calculate the difference between 3 and 5 (5-3), you get 2. 2 consider the middle number (4) and one of the extremes (3 or 5) but not bot of them, and we want the complete distance between 3 and 5 including extremes, so we have to sum 1, 2+1=3. This includes 3, 4, and 5.
sum(exon_sizes_first_transcript) 
#sum of the exon sizes of the first transcript is equal to the size of the transcript (1057). This size is bigger than coding sequence length so, both transcript and end/start exon positions included some non-coding sequences as the 5' UTR.
	
#calculate the length of the whole transcript including introns
whole_transcript_length <- unique(sema4a_full_first_transcript$transcript_end) - unique(sema4a_full_first_transcript$transcript_start) + 1

#calculate the length of introns
intron_lengths <- NULL
#i=1
for (i in 1:(nrow(sema4a_full_first_transcript) - 1)) { #we don´t want the last exon, because the transcript ends there, no more intros after that (only the UTR sequence)

	#calculate the difference between end of the [i] exon and the beginning of the [i+1] exon (i.e., the next one) without including the extremes. Because of this we subtract 1. 2:4 has 3 number including the extremes but only 1 without them (4-2-1=1). This will calculate the introns between each pair of exons. 
	calculated_intron_length <- sema4a_full_first_transcript[i + 1, ]$exon_chrom_start - sema4a_full_first_transcript[i, ]$exon_chrom_end - 1 #the start of the second exon has a higher position than the end of the previous one, because of this the former has to be first.

	#save
	intron_lengths <- append(intron_lengths, calculated_intron_length)
}

#The whole transcript including introns (but UTR before the first coding exon) should be equal to the sum of the exon sizes plus the sum of intron lengths
if( whole_transcript_length != sum(exon_sizes_first_transcript) + sum(intron_lengths)) {
	stop("The length of the whole transcript is not equal to the sum of the exon sizes plus the sum of intron lengths")
}



#################################################################
### CHECK CDS START-END MATCHES BETWEEN ENSEMBLE AND MY DATA ####
#################################################################

#cds start and end can be used to calculate when the coding regions start and end. Below you have support for this. HOWEVER, it much easier to use genomic coding start and genomic coding end (see below). These variables give the exact position in which the coding region start and end in each exon.
sema4a_full_first_transcript

#cds_start seems to indicate the beginning of the coding region in a exon, what it is the coding sequence length at that point. 1 for example would indicate that this is the beginning of the exon. If that exon have a total transcript length of 60, the coding sequences would end at 60. See several examples
#The transcript (ENST00000435124) of SEMA4A (ENST00000435124.5; https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000196189;r=1:156147366-156161498;t=ENST00000435124; FIGURE 14):
#The first exon (ENSE00001737906; 156147366-156147430) has no coding sequence, so cds_start is NA.
#The second exon (ENSE00003474176; 156154550-156154717) has the first coding sequence, so cds_start is 1, and the end is 139. The transcript length according to the webpage is 168, but according to the map, the start of the exon is not translated. If we count the splice region (removed in splicing and hence not translated; in orange) and the blue-greenish region (UTR) at the beginning we get 29 bases (I have seen this pasting the exon sequence into libre office, selecting the orange and blue bases at the beginning, then tools and word count to check the number of characters). 168 - 29 makes 139 which is exactly the value of cds_end.
#The third exon (ENSE00003784807; 156156414-156156574) has a transcript length of 161 and coding sequences are at the beginning and end of the exon, according the webpage. According to my data, cds_start and end are 140 and 300, which makes 300-140+1=161 (we sum 1 because we want to include both 300 and 140 as they are start and end of coding sequence in this exon), exactly the length of the transcript according to the webpage.

#the fourth exon (ENSE00003623497; 156158070-156158132) has a transcript length of 63 and coding sequences are at the beginning and end of the exon, according the webpage. It has a transcript length of 63 according to the webpage, and according to my data the cds start and end are 301 and 363, respectively. 363-301+1=63, which is exactly the length of transcripts

#the fifth exon (ENSE00003487912; 156158388-156158486) it has a transcript length of 99 and the beginning and end of the exon has coding basis according to the webpage. According to my data, the cds start and end is 364 and 462, so the length of the coding sequence is 462-364+1=99, same transcript length.

#For the last exon occurs the same, 811 and 963 for start and end of coding sequence, which makes 963-811+1=153, which is the transcript length according the webpage. Again, the first and last base is coding in this exon according to the webpage.

#The transcript (ENST00000368282.1):
sema4a_full_second_transcript <- full_exon_data[which(full_exon_data$ensembl_gene_id == "ENSG00000196189" & full_exon_data$ensembl_transcript_id == "ENST00000368282"), ]
#reorder by exon position
sema4a_full_second_transcript <- sema4a_full_second_transcript[match(sort(sema4a_full_second_transcript$exon_chrom_start), sema4a_full_second_transcript$exon_chrom_start), ]
sema4a_full_second_transcript
#The first exon (ENSE00001446766; 156154371-156154717) has as start and end cds 1 and 139, which makes 139-1+1=139 cosing bases. According the webpage, the last 139 bases are coding bases (removing the initial bases in green).

#The rest of exons are complete coding except the last one (ENSE00001446765; 156176405-156177752). The start and end coding sequence is 1694 and 2286, respectively. 2286 - 1694 + 1 = 593, which is the exact number of coding basis in the last exon according to the webpage (I have seen this pasting the exon sequence into libre office, selecting all bases except the orange tail, then tools and word count to check the number of characters).

#Again, with start and end cds we can calculate the position of the starting and end coding sequences

#CDS start = 1 indicates the first coding exon, without the UTR 5' sequence. The last cds end of the transcript has the same value of the cds length. Therefore we can use cds_start = 1 and cds_end=cds_length to select the data of start (first) and end (last) exon. We can count from the end of the first exon with coding sequence, as many bases as the cds end. For example, cds start and end of the first exon with coding sequence is 1 and 100, the end of that exon is at 110, we should rest 100 from 110 plus 1 to get the range of coding basis. For the end, we should look for the exon whose cds end is equal to cds length. In that exon, we calculate the range of coding baes (cds_end - cds_start) and sum to the exon start. We do not need to rest 1 to the range because the left extreme would be inclued when summing the range to the begining of the exon. For example, 500-400=100, and the first exon start at 1000, so 1000 + 100 makes 1100. The right extreme is included in 100 and the left in 1000.

#we are not going to calculate the number of coding bases using cds start/end. If you want a detailed example of how to calculate the number of coding bases using cds start/end, then you should go to the original script of gene coordinates for the MDR paper.



##############################################################################
#### CALC THE LENGTH OF CODING REGIONS USING GENOMIC CODING START AND END ####
##############################################################################

#Calculate coding density using genomic coding start/end instead of cds start/end and compare results

#cds start/end gives you the number of coding bases at the beginning and end of the coding sequence in each coding exon, but it does not give position of start and of these sequences. Genomic start and end gives you exactly the position of start and end of the coding sequence in each coding exon

#In addition, I think that cds start and end, only considers the order of the strand. If the strand is negative, the last coding exons in coordinates will be the first one according to cds start (cds_start=1). In contrast, genomic start considers the coordinates of the positive strand always.

#Therefore, we are going to use genomic start/end


## SEMA4A_1
#extract the exon data of the ENST00000435124 transcript of SEMA4A_1
test_genomic_start_end_sema4a_1 <- sema4a_exon_data_indv_query[which(sema4a_exon_data_indv_query$ensembl_transcript_id == "ENST00000435124"), ]
test_genomic_start_end_sema4a_1

#remove the exons without any coding sequence
test_genomic_start_end_sema4a_1 <- test_genomic_start_end_sema4a_1[which(!is.na(test_genomic_start_end_sema4a_1$genomic_coding_start) | !is.na(test_genomic_start_end_sema4a_1$genomic_coding_end)), ]

#reorder en basis on position
test_genomic_start_end_sema4a_1 <- test_genomic_start_end_sema4a_1[order(test_genomic_start_end_sema4a_1$cds_start), ]

#calculate the complete length of each coding region according to genomic coding start/end
coding_length_by_genom_start_end_sema4a_1 <- test_genomic_start_end_sema4a_1$genomic_coding_end - test_genomic_start_end_sema4a_1$genomic_coding_start + 1 #we sum 1 because we want the complete range

#First coding exon (ENSE00003474176) has a coding sequence of 139 bases, which is exactly the number of blue bases indicated in the webpage of ensemble. The rest has exactly the same coding length than that indicated in the webpage (https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000196189;r=1:156147366-156161498;t=ENST00000435124).
coding_length_by_genom_start_end_sema4a_1

#check that the the complete range of coding sequences according to genomic start/start is similar to that calculated with cds start/end
summary(coding_length_by_genom_start_end_sema4a_1 == (test_genomic_start_end_sema4a_1$cds_end - test_genomic_start_end_sema4a_1$cds_start + 1)) #we sum 1 because we want the complete range

#check that the genomic starts and ends are the same than the ranges I have calculated using cds start/end
#Check the original script for the MDR paper to see this check


## SEMA4A_2
#extract the exon data of the ENST00000414683 transcript of SEMA4A_2
test_genomic_start_end_sema4a_2 <- sema4a_exon_data_indv_query[which(sema4a_exon_data_indv_query$ensembl_transcript_id == 'ENST00000414683'), ]
test_genomic_start_end_sema4a_2

#remove the exons without any coding sequence
test_genomic_start_end_sema4a_2 <- test_genomic_start_end_sema4a_2[which(!is.na(test_genomic_start_end_sema4a_2$genomic_coding_start) | !is.na(test_genomic_start_end_sema4a_2$genomic_coding_end)), ]

#reorder en basis on position
test_genomic_start_end_sema4a_2 <- test_genomic_start_end_sema4a_2[order(test_genomic_start_end_sema4a_2$cds_start), ]

#calculate the complete length of each coding region according to genomic coding start/end
coding_length_by_genom_start_end_sema4a_2 = test_genomic_start_end_sema4a_2$genomic_coding_end - test_genomic_start_end_sema4a_2$genomic_coding_start + 1 #we sum 1 because we want the complete range
coding_length_by_genom_start_end_sema4a_2 

#First coding exon (ENSE00003589564) has a coding sequence of 3 bases, which is exactly the number of blue bases indicated in the webpage of ensemble. The rest has exactly the same coding length than that indicated in the webpage (https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000196189;r=1:156147366-156161498;t=ENST00000414683). 

#check that the the complete range of coding sequences accroding to genomic start/start is similar to that calculated with cds start/end
summary(coding_length_by_genom_start_end_sema4a_2 == (test_genomic_start_end_sema4a_2$cds_end - test_genomic_start_end_sema4a_2$cds_start + 1)) #we sum 1 because we want the complete range

#check that the genomic starts and ends are the same than the ranges I have calculated
#Check the original script for the MDR paper to see this check

#IN SUMMARY:
#The number of coding bases in each exon is the same calculated with cds start/end and genomic_coding_start/end.
#Importantly, genomic_coding_start/end give the coordinate position of the start and end of coding basis, which is much more useful to calculate the number of coding bases within a genomic window.
#Therefore, we will use genomic_coding_start/end for calculating coding density.



#############################################################
######  GENES AND TRANSCRIPTS THAT ARE NOT TRANSLATED #######
#############################################################

#extract the unique cases of gene and transcript type
unique_gene_types <- unique(full_exon_data$gene_biotype)
unique_transcript_types <- unique(full_exon_data$transcript_biotype)

#gene and transcript types are not exactly the same. All gene types are included in the transcript types, but some transcript types are not included in the gene types
unique_gene_types %in% unique_transcript_types
unique_transcript_types %in% unique_gene_types


## see what types of genes have only NA for cds_length ###
#for each gene types
test_gene_type_cds_length <- data.frame(selected_gene_type = NA, result_test = NA)
#i=1
for(i in 1:length(unique_gene_types)){

	#select the [i] gene type
	selected_gene_type <- unique_gene_types[i]

	#extract the exons belonging to gene of the [i] gene type
	cds_lengths <- full_exon_data[which(full_exon_data$gene_biotype == selected_gene_type), ]$cds_length

	#test if all lengths of coding sequences are equal to NA for this type
	result_test <- all(is.na(cds_lengths))

	#save the results
	test_gene_type_cds_length <- rbind.data.frame(test_gene_type_cds_length, cbind.data.frame(selected_gene_type, result_test))
}
#remove first row with NAs
test_gene_type_cds_length <- test_gene_type_cds_length[-1, ]

#take a look to case with ALL NA for coding sequence length
test_gene_type_cds_length[which(test_gene_type_cds_length$result_test == TRUE), ]


## see what types of transcripts have only NA for cds_length ###
#for each type of transcript
test_transcript_type_cds_length <- data.frame(selected_transcript_type = NA, result_test = NA)
#i=1
for(i in 1:length(unique_transcript_types)){

	#select the [i] transcript type
	selected_transcript_type <- unique_transcript_types[i]

	#extract the exons belonging to transcript of the [i] transcript type
	cds_lengths <- full_exon_data[which(full_exon_data$transcript_biotype == selected_transcript_type), ]$cds_length

	#test if all lengths of coding sequences are equal to NA for this type
	result_test <- all(is.na(cds_lengths))

	#save the results
	test_transcript_type_cds_length	<- rbind.data.frame(test_transcript_type_cds_length, cbind.data.frame(selected_transcript_type, result_test))
}

#remove first row with NAs
test_transcript_type_cds_length <- test_transcript_type_cds_length[-1, ]

#take a look to case with ALL NA for coding sequence length
test_transcript_type_cds_length[which(test_transcript_type_cds_length$result_test == TRUE), ]

#I have found that many gene and transcript types have only NA data for coding sequence length. These sequences are not translated. For example, ENST00000485575 of SEMA4A-008, which has no translation length (https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000196189;r=1:156147366-156177752;t=ENST00000485575)

#Removing these gene types would lead to consider genes that are not coding
#Non-coding genes will be removed.

#Removing these transcript types would not affect to coding density withing a coding gene if we use coding sequence length, as these transcript types have NA for that variable. The problem I see is that one of this transcript is at the end or beginning of a gene, then we would be considering non-coding regions for gene size. This would affect to gene length and the calculation of the center. I would REMOVE all the rows without cds_length data. Line 200.

#SOLUTION: David said that for gene center we can use the full gene length according gene start and end, but for coding density we need to avoid these transcripts.

#gene and transcript biotypes removing those cases with cds_length all NA
unique(full_exon_data[which(!is.na(full_exon_data$cds_length)), ]$gene_biotype)
unique(full_exon_data[which(!is.na(full_exon_data$cds_length)), ]$transcript_biotype)


##check for non-functional transcripts
#I have found some genes with coding sequences, but that are not functional, like a inmunglobulin gene (IGHV1OR15-9; ENSG00000188403; http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000188403;r=15:20169919-20170354;t=ENST00000338912).

#Do you think we should get ride of them? I am not very sure if there a filter or attribute to remove non functional transcripts. I have seen it in some gene types that include in their name the label variable (i guess too much variability in the sequence) like TR_V_gene (T cell receptor gamma variable) or IG_V_gene (immunoglobulin variable), but not in other genes of the variable category (ENST00000390633). In other types like T cell receptor gamma joining (TR_J_gene) or T cell receptor delta diversity (TR_D_gene) i have found functional proteins

#I don´t understand why these cases have a different gene type from protein coding...
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'IG_V_gene'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'TR_V_gene'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'TR_C_gene'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'TR_J_gene'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'TR_D_gene'), ], 100)

#transcripts with decay and pseudogenes
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'non_stop_decay'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'nonsense_mediated_decay'), ], 100)
head(full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$transcript_biotype == 'polymorphic_pseudogene'), ], 100)
#transcript with decay I guess means that these transcript are destroyed because have not stop codon or the stop codon is too soon. Am I right? We should consider this in our analyses?
#Same goes for pseudogenes. A gene that has homology to known protein-coding genes but contain a frame shift and/or stop codon(s) which disrupts the open reading frame. Thought to have arisen through duplication followed by loss of function.

#This is the list of all transcript biotypes:
unique(full_exon_data$transcript_biotype)
#"protein_coding"
#"nonsense_mediated_decay"
#"processed_transcript"
#"miRNA"
#"processed_pseudogene"
#"misc_RNA"
#"lincRNA"
#"snRNA"
#"snoRNA"
#"antisense"
#"rRNA"
#"transcribed_unprocessed_pseudogene
#"transcribed_processed_pseudogene"
#"unprocessed_pseudogene"
#"sense_intronic"
#"pseudogene"
#"retained_intron"
#"unitary_pseudogene"
#"sense_overlapping"
#"IG_V_gene"
#"IG_D_gene"
#"3prime_overlapping_ncrna"
#"IG_V_pseudogene"
#"non_stop_decay"
#"polymorphic_pseudogene"
#"Mt_tRNA"
#"Mt_rRNA"
#"TR_V_pseudogene"
#"IG_C_pseudogene"
#"TR_C_gene"
#"TR_J_gene"
#"TR_V_gene"
#"TR_J_pseudogene"
#"translated_processed_pseudogene"
#"IG_J_pseudogene"
#"IG_J_gene"
#"IG_C_gene"
#"TR_D_gene"
#def of gene types: https://uswest.ensembl.org/info/genome/genebuild/biotypes.html

#Should we remove all this types and only use protein coding? I know that between these categories there are non-functional proteins, but not sure if all of them.. should I revise each type to ensure that these types of sequences have no influence on phenotype?
#YES, we only want coding genes. For transcripts withing coding-functional genes, we will also remove this type of transcripts. We do not want them in coding density calculation. For gene length, we will use gene start and end, so it could include them, but it is not problematic (small changes of length and center...)

#select those rows without NA for coding sequence length and with gene and transcript type equal to protein coding
full_exon_data_filtered <- full_exon_data[which(!is.na(full_exon_data$cds_length) & full_exon_data$gene_biotype == "protein_coding" & full_exon_data$transcript_biotype == "protein_coding"), ]
head(full_exon_data_filtered)

#check the difference in size
length(unique(full_exon_data$ensembl_gene_id))
length(unique(full_exon_data_filtered$ensembl_gene_id))
length(unique(full_exon_data$ensembl_gene_id)) - length(unique(full_exon_data_filtered$ensembl_gene_id))

#check
if( 
	all(!is.na(full_exon_data_filtered$cds_length)) &
	unique(full_exon_data_filtered$gene_biotype) == "protein_coding" &
	unique(full_exon_data_filtered$transcript_biotype) == "protein_coding"
){
	print("All rows have coding sequence length and gene and transcript type equal to protein coding")
} else {
	stop("There are rows without coding sequence length or with gene or transcript type different to protein coding")
}

#Given that we remove all sequences without data for cds, if a transcript have at least one coding exon, it will have cds value and it will be included. This cds value will not consider the non-coding bases inside the coding exon like the 5´ UTR extreme, etc.. See figures 11 and 15 (the 3´ and 5´ below of the line i think that refer to the other strand).

#For gene length and center, we count these non-coding sequences inside the coding exons along with non-functional transcripts as we use gene start and end. According to David, this is not a problem. The difference in size and center position will be small in relation to the LD blocks in the human genome (thousands of bases).



##############################################################################
###### FILTER GENE LIST USING THE PREVIOUS FILTER APPLIES ON EXON DATA #######
##############################################################################

## select only genes that are included in the gene list and the exon data
#We are going to filter the gene list all_genes_grch38_human_filtered removing those id genes that are not included in our exon data, i.e., non coding genes or genes without protein coding transcripts.

#note that you can have exon data from genes in the sexual chromosomes, but the all_genes_grch38_human_filtered list is already filtered for that. Therefore, the loop for calculating windows and coding density will only use autosomal genes. Only coding density of autosomal genes will be calculated.

#REMEMBER: "all_genes_grch38_human_filtered" determine the genes that will be used in the analyses.

#extract from all_genes_grch38_human_filtered, i.e., dataframe with the gene names data, those rows with an gene id included in the filtered dataset with exon data.
all_genes_grch38_exon_filter <- all_genes_grch38_human_filtered[which(all_genes_grch38_human_filtered$ensembl_gene_id %in% unique(full_exon_data_filtered$ensembl_gene_id)), ]

#check that the filtering occurred well
summary(unique(all_genes_grch38_exon_filter$ensembl_gene_id) %in% unique(full_exon_data_filtered$ensembl_gene_id))
summary(unique(full_exon_data_filtered$ensembl_gene_id) %in% unique(all_genes_grch38_exon_filter$ensembl_gene_id))

#some genes of the exon data are not included in the list of genes
#extract the list of lost genes in exon data, i.e., genes that are includqed in the exon data but not in the gene list
lost_genes <- full_exon_data_filtered[which(!full_exon_data_filtered$ensembl_gene_id %in% unique(all_genes_grch38_exon_filter$ensembl_gene_id)), ]$ensembl_gene_id
#extract the exon data of these genes
exon_data_lost_genes <- full_exon_data_filtered[which(full_exon_data_filtered$ensembl_gene_id %in% lost_genes), ]
#take a look
head(exon_data_lost_genes, 100)
exon_data_lost_genes[1,] #the first one is from Y chromosome.
unique(exon_data_lost_genes$chromosome_name) #but there are also from Y chromosome and patches from areas with high variability.

#check no lost gene is included in the original autosomals
if (TRUE %in% c(1:22 %in% unique(exon_data_lost_genes$chromosome_name))) {
	stop("There are genes of the Y chromosome in the lost genes")
}

#These genes id are not included in all_genes_grch38_exon_filter because we removed those that do not have 1:22 chromosome data. There are some regions with great variability like the MHC for which different patches are included. Also fixing patches are created. similar, genes of the Y and X chromosome are not included. We have removed all of this in all_genes_grch38_exon_filter. Therefore, they will be not used in the analyses!

#In hg19, we also have a gene id of SLC25A26 (ENSG00000261657). This is id was not included the filtered gene list but the hgcs symbol does!!
exon_data_lost_genes[which(exon_data_lost_genes$ensembl_gene_id == "ENSG00000144741"), ]

#If you look for SLC25A26 in the filtered gene list, you will obtain a gene id that is already included exon data.
all_genes_grch38_exon_filter[which(all_genes_grch38_exon_filter$hgnc_symbol == "SLC25A26"), ]$ensembl_gene_id #gene id is ENSG00000144741

#According to the biomart online server, ENSG00000144741 is SLC25A26, a gene of the chromosome 3, while ENSG00000261657 is also SLC25A26, but it was removed in hg38 and merged with ENSG00000144741
#https://www.ensembl.org/Homo_sapiens/Gene/Idhistory?g=ENSG00000261657

#THEREFORE, if we use the gene list as a guide for searching genes, we will have these genes but NOT the patches

#remove these patches and the X-Y chromosome from exon data
full_exon_data_filtered_final <- full_exon_data_filtered[which(full_exon_data_filtered$ensembl_gene_id %in% unique(all_genes_grch38_exon_filter$ensembl_gene_id)), ]

#check that the X and Y chromosome are not included
if (TRUE %in% c(c("X", "Y") %in% unique(full_exon_data_filtered_final$chromosome_name))) {
	stop("WE HAVE NOT CORRECTLY FILTERED OUT X AND Y")
}

#check that only genes included in the gene list are present in the exon data
summary(full_exon_data_filtered_final$ensembl_gene_id %in% all_genes_grch38_exon_filter$ensembl_gene_id)
summary(all_genes_grch38_exon_filter$ensembl_gene_id %in% full_exon_data_filtered_final$ensembl_gene_id)


## Remove the genes and transcript from gene list that not produce coding protein ###
#select those transcript that belong to coding genes and produce proteins, i.e., gene and transcript biotype equal to 'protein_coding'
coding_transcripts <- full_exon_data_filtered_final[
	which(
		full_exon_data_filtered_final$gene_biotype == "protein_coding" &
		full_exon_data_filtered_final$transcript_biotype == "protein_coding"
	), ]$ensembl_transcript_id
length(coding_transcripts) == nrow(full_exon_data_filtered_final)
#all the rows of full_exon_data_filtered_final satisfy the condition of having gene and transcript biotype as 'protein_coding'

#select from the gene list those rows belonging to these coding transcripts
all_genes_grch38_exon_filter_final <- all_genes_grch38_exon_filter[which(all_genes_grch38_exon_filter$ensembl_transcript_id %in% unique(coding_transcripts)), ] #If a gene has an coding transcript, it will be included because that transcript will have as biotype 'protein_coding' and the gene biotype will be 'protein_coding'

#check 
summary(all_genes_grch38_exon_filter_final$ensembl_transcript_id %in% full_exon_data_filtered_final$ensembl_transcript_id)
summary(full_exon_data_filtered_final$ensembl_transcript_id %in% all_genes_grch38_exon_filter_final$ensembl_transcript_id)

#check for no name genes (no hgnc symbol) that remained after applying the filter
nrow(all_genes_grch38_exon_filter_final[which(all_genes_grch38_exon_filter_final$hgnc_symbol == ""), ])
#There are cases of uncharacterized proteins like AP000350.10 (ENSG00000251357)
#These are proteins that are not very well studied and we don´t know their function, but according to David we should include them. They are translated, and can be under selection...

#Summary:
#I have selected coding data data for coding density. This does NOT include non-coding genes, non-coding transcripts, non-functional transcripts (biotype different from coding protein), and the non-coding regions inside the coding exons.
#BUT for gene length I am using gene start and end, which can include in the gene length some areas of non-coding or non-functional transcripts. I think you told me that this only changes the gene length and center position of a few bases, which is nothing compare to the bigger size of LD blocks in human genome.
#That´s right? #YES, that it is perfect.
#Also we have removed genes whose chromosome names is not 1:22.



#######################################################
###### CHECKS FOR FILTERS IN THE FINAL DATASETS #######
#######################################################

##check that there is no row with other chromosome_name rather than 1:22 ###
#for all_genes_grch38_exon_filter_final
#check that no row has as a chromosome name different from 1:22
nrow(all_genes_grch38_exon_filter_final[which(!(all_genes_grch38_exon_filter_final$chromosome_name %in% c(1:22))), ]) == 0

#check that 1:22 are included as chromosome names
c(1:22) %in% unique(all_genes_grch38_exon_filter_final$chromosome_name)

#check that all chromosome names are included in 1:22
unique(all_genes_grch38_exon_filter_final$chromosome_name) %in% c(1:22)

#for full_exon_data_filtered_final
#chromosome names of the genes included in exon data
chromosome_names_exon_data <- all_genes_grch38_exon_filter_final[which(all_genes_grch38_exon_filter_final$ensembl_gene_id %in% unique(full_exon_data_filtered_final$ensembl_gene_id)), ]$chromosome_name

#check that 1:22 are included as chromosome names
c(1:22) %in% unique(chromosome_names_exon_data)

#check that all chromosome names are included in c(1:22)
unique(chromosome_names_exon_data) %in% c(1:22)


##check that all transcripts included has 'protein_coding' as biotype and cds_length is not NA ###
#check that the gene and transcript biotype in exon data is protein coding
unique(full_exon_data_filtered_final$gene_biotype) == "protein_coding"
unique(full_exon_data_filtered_final$transcript_biotype) == "protein_coding"

#extract the from the exon data those rows with gene id included in the gene list, and see if the gene biotype is protein coding
unique(full_exon_data_filtered_final[which(full_exon_data_filtered_final$ensembl_gene_id %in% unique(all_genes_grch38_exon_filter_final$ensembl_gene_id)), ]$gene_biotype) == "protein_coding"

#extract the from the exon data those rows with transcript id included in the gene list, and see if the transcript biotype is protein coding
unique(full_exon_data_filtered_final[which(full_exon_data_filtered_final$ensembl_transcript_id %in% unique(all_genes_grch38_exon_filter_final$ensembl_transcript_id)),]$transcript_biotype) == "protein_coding"


##check that each gene name has only one gene id. I can not select those cases without gene name
#extract the gene names without empty cases ('')
gene_names_no_null_v2 <- unique(all_genes_grch38_exon_filter_final$hgnc_symbol)
gene_names_no_null_v2 <- gene_names_no_null_v2[-which(gene_names_no_null_v2 == "")]
!"" %in% gene_names_no_null_v2

#for each gene name
test_duplicated_gene_id_v2 <- data.frame(selected_gene_name = NA, test_result = NA)
#i=1
for (i in 1:length(gene_names_no_null_v2)) {

	#select the [i] gene name
	selected_gene_name <- gene_names_no_null_v2[i]

	#extract the unique cases of gene id for the [i] hgnc symbol
	unique_gene_id <- unique(all_genes_grch38_exon_filter_final[which(all_genes_grch38_exon_filter_final$hgnc_symbol == selected_gene_name), ]$ensembl_gene_id)

	#test that the number of unique cases is 1 (only 1 gene id for each gene name)
	test_result <- length(unique_gene_id) == 1

	#save
	test_duplicated_gene_id_v2 <- rbind.data.frame(test_duplicated_gene_id_v2, cbind.data.frame(selected_gene_name, test_result))
}
#remove first row with NAs
test_duplicated_gene_id_v2 <- test_duplicated_gene_id_v2[-1, ]
#take a look for genes with repeated gene id
summary(test_duplicated_gene_id_v2) #we have several cases
#extract information of these genes
genes_with_duplicated_id_v2 <- test_duplicated_gene_id_v2[which(test_duplicated_gene_id_v2$test_result == FALSE), ]$selected_gene_name

#Problematic cases
problematic_cases <- all_genes_grch38_exon_filter_final[which(all_genes_grch38_exon_filter_final$hgnc_symbol %in% genes_with_duplicated_id_v2), ]

#HERC3 has two gene ids, ENSG00000287542 (https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000138641;r=4:88592434-88708541) and ENSG00000138641 (https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000138641;r=4:88592434-88708541). It is like ENSG00000287542 is a transcript of ENSG00000138641 that has been separated as an independent gene (it only has 1 transcript). They are in the same region, and has the same hgnc symbol.

#POLR2J3 has two gene ids: ENSG00000168255 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000168255;r=7:102537918-102572653) and ENSG00000285437 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000285437;r=7:102562133-102572583). Both IDs have a lot of transcripts, but they are in the same region. It is like they have split the transcripts in two different genes.

#ZNF724 has two gene ids, ENSG00000196081 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000196081;r=19:23221599-23250394) and ENSG00000283201 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000283201;r=19:23222755-23274221;t=ENST00000611392). It is like ENSG00000283201 is a transcript of ENSG00000138641 that has been separated as an independent gene (it only has 1 transcript). They are in the same region, and has the same hgnc symbol.

#PINX1 has two gene ids, ENSG00000258724 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000258724;r=8:10725399-10839847;t=ENST00000554914) and ENSG00000254093 (https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000254093;r=8:10764961-10839884). It is like they have splitted the transcripts of the same gene. They are relatively close.

#According to David in the MDR paper, cases like this are very strange, so he told me to remove each pair of cases like this.

#These are only 8 genes out 18K...
if (length(genes_with_duplicated_id_v2)>4) {
	stop("We have more gene symbols with several gene IDs than expected")
}

#remove these genes
all_genes_grch38_exon_filter_final <- all_genes_grch38_exon_filter_final[which(!all_genes_grch38_exon_filter_final$ensembl_gene_id %in% problematic_cases$ensembl_gene_id), ]
full_exon_data_filtered_final <- full_exon_data_filtered_final[which(!full_exon_data_filtered_final$ensembl_gene_id %in% problematic_cases$ensembl_gene_id), ]


##check AGAIN that each gene name has only one gene id. I can not select those cases without gene name
#extract the gene names without empty cases ('')
gene_names_no_null_v2 <- unique(all_genes_grch38_exon_filter_final$hgnc_symbol)
gene_names_no_null_v2 <- gene_names_no_null_v2[-which(gene_names_no_null_v2 == "")]
!"" %in% gene_names_no_null_v2
#for each gene name
test_duplicated_gene_id_v2 <- data.frame(selected_gene_name = NA, test_result = NA)
#i=1
for (i in 1:length(gene_names_no_null_v2)) {

	#select the [i] gene name
	selected_gene_name <- gene_names_no_null_v2[i]

	#extract the unique cases of gene id for the [i] hgnc symbol
	unique_gene_id <- unique(all_genes_grch38_exon_filter_final[which(all_genes_grch38_exon_filter_final$hgnc_symbol == selected_gene_name), ]$ensembl_gene_id)

	#test that the number of unique cases is 1 (only 1 gene id for each gene name)
	test_result <- length(unique_gene_id) == 1

	#save
	test_duplicated_gene_id_v2 <- rbind.data.frame(test_duplicated_gene_id_v2, cbind.data.frame(selected_gene_name, test_result))
}
#remove first row with NAs
test_duplicated_gene_id_v2 <- test_duplicated_gene_id_v2[-1, ]
#take a look for genes with repeated gene id
summary(test_duplicated_gene_id_v2) #No several cases
#extract information of these genes
genes_with_duplicated_id_v2 <- test_duplicated_gene_id_v2[which(test_duplicated_gene_id_v2$test_result == FALSE), ]$selected_gene_name
#check
if (length(genes_with_duplicated_id_v2) != 0) {
	stop("We have gene symbols with several gene IDs")
}


###############################################
###### SAVE THE GENE LIST AND EXON DATA #######
###############################################

#save the list of genes
write.table(
	x = all_genes_grch38_exon_filter_final,
	file = gzfile("./data/gene_coordinates/all_genes_grch38_exon_filter_final.txt.gz"),
	sep = "\t",
	row.names = FALSE)

#save the exon data
write.table(
	x = full_exon_data_filtered_final,
	file = gzfile("./data/gene_coordinates/full_exon_data_filtered_final.txt.gz"),
	sep = "\t",
	row.names = FALSE)



#############################################################################
#### FUNCTION TO EXTRACT GENE POSITIONS AND SOME FACTORS PER CHROMOSOME #####
#############################################################################

##load all the necessary data and do minimum preparations, so you can just run the script from here if you want to. 

#initial preps (wd and packages)
setwd("../")
getwd()
library(biomaRt) #for connecting with ensemble
library(GenomicRanges) #for calculating overlap between exons of a gene
require(stringr) #for replacing "," in the total chromosome length

#load the list of genes filtered
all_genes_grch38_exon_filter_final <- read.table(
	"./data/gene_coordinates/all_genes_grch38_exon_filter_final.txt.gz",
	sep = "\t",
	header = TRUE)
str(all_genes_grch38_exon_filter_final)
head(all_genes_grch38_exon_filter_final)

#load the exon data filtered
full_exon_data_filtered_final <- read.table(
	'./data/gene_coordinates/full_exon_data_filtered_final.txt.gz',
	sep = "\t",
	header = TRUE)
str(full_exon_data_filtered_final)
head(full_exon_data_filtered_final)

#load the chromosome length data to cut those windows that surpass the end of the chromosome (I will also cut those that surpass the start of the chromosome, i.e., negative bases)
chrom_length_ucsc_hg38 <- read.table(
	"./data/gene_coordinates/chromosome_length/chrom_length_final_v1.txt",
	sep = "\t",
	header = TRUE)
print(chrom_length_ucsc_hg38)
#you cannot have window where no chromosome exists. This can lead to problems, for example in the calculation of the recombination rate, where a window can be discarded if no data points are close to the extremes. Therefore, a window with one end outside of the chromosome could be removed even if recombination data exists within the boundaries of the chromosome.



##por aquiii




#load the gap data to subtract gap length from the window length when coding density is calculated. 
#IMPORTANT!!! REVISE YOUR ARE USING THE LAST VERSION OF GAPS!!
gap_results_extracted_final = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gap_length_final_v3.txt', sep='\t', header=TRUE) #the number of coding bases has to be divided by the real length of bases inside the window. If you have gaps of no data inside the window, you do not know if these regions are coding or not. Therefore, you have to remove them when calculating the proportion of coding bases respect to the total of bases. You have to count only those bases for which you have data. 
#check that all the gene ids for which we have gaps are included in the gene list and exon data.
#gene list
!FALSE %in% c(gap_results_extracted_final$gene_id %in% all_genes_grch37_exon_filter_final$ensembl_gene_id)
!FALSE %in% c(all_genes_grch37_exon_filter_final$ensembl_gene_id %in% gap_results_extracted_final$gene_id)
#exon data
!FALSE %in% c(gap_results_extracted_final$gene_id %in% full_exon_data_filtered_final$ensembl_gene_id)
!FALSE %in% c(full_exon_data_filtered_final$ensembl_gene_id %in% gap_results_extracted_final$gene_id)


#VERY IMPORTANT!!! EVERY TIME you make changes in the window calculation, you have to run this code, then run the script of gaps to calculate the gaps inside each window and then run again this script to calculate the coding density with the gaps calculated in the new windows. 

#write the function
#selected_chr=1
#selected_id_gene="ENSG00000058453" #gen with gaps close.
#selected_id_gene="ENSG00000185220" #gen very close to the end of the chromosome 1
#selected_id_gene="ENSG00000186092" #gen very close to the start of the chromosome 1
genomic_coords_coding_density = function(selected_chr){ #each chromosome will be analyzed independently. It does not make sense to create windows between chromosomes, linkage and other factors that affect close regions act within chromosomes
	#from all_genes_grch37_exon_filter_final, select the rows of the selected chromosome and then get the gene id, selecting the unique cases
	unique_gene_ids_chr = unique(all_genes_grch37_exon_filter_final[which(all_genes_grch37_exon_filter_final$chromosome_name == selected_chr),]$ensembl_gene_id)

	#for each id gene
	final_positions = data.frame(chromosome_name=NA, hgnc_symbol=NA, gene_id=NA, test_rows_initial_datasets_1=NA, test_rows_initial_datasets_2=NA, test_gene_id=NA, gene_biotype=NA, transcript_id=NA,transcript_biotype=NA, exon_id=NA, test_gene_transcripts_biotype=NA, gene_start=NA, gene_end=NA, gene_length=NA, middle_point=NA, test_gene_length_center=NA, windows_check_1=NA, windows_check_2=NA, windows_check_3=NA, windows_check_4=NA, windows_check_5=NA, lower_end_window_50kb=NA, upper_end_window_50kb=NA, lower_end_window_100kb=NA, upper_end_window_100kb=NA, lower_end_window_200kb=NA, upper_end_window_200kb=NA, lower_end_window_500kb=NA, upper_end_window_500kb=NA, lower_end_window_1000kb=NA, upper_end_window_1000kb=NA, n_genes_50kb=NA, n_genes_100kb=NA, n_genes_200kb=NA, n_genes_500kb=NA, n_genes_1000kb=NA, check_n_genes_0=NA, check_n_genes_1=NA, check_n_genes_2=NA, check_n_genes_3=NA, check_n_genes_4=NA, check_coding_density_1=NA, check_coding_density_2=NA, check_coding_density_4=NA, check_coding_density_5=NA, check_coding_density_6=NA, check_coding_density_7=NA, test_all_ranges_included=NA, test_iranges=NA, test_start_end=NA, test_ranges_overlap=NA, test_ranges_overlap_2=NA, coding_density_50kb=NA, coding_density_100kb=NA, coding_density_200kb=NA, coding_density_500kb=NA, coding_density_1000kb=NA, test_position_gene_exons=NA, test_chr_name=NA, test_na_cds=NA)
	for(i in 1:length(unique_gene_ids_chr)){

		#select the [i] gene
		selected_id_gene = unique_gene_ids_chr[i] #I am going to comment the results for the gene SCYL3 (ENSG00000000457; http://grch37.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000000457;r=1:169818772-169863093;t=ENST00000367771). This gene has all its transcripts in reverse sense
			#this is [i] = 1

		##########################################
		###### SELECT DATA OF THE [i] GENE #######
		##########################################

		#extract gene characteristics from the gene list
		gene_characteristics = all_genes_grch37_exon_filter_final[which(all_genes_grch37_exon_filter_final$ensembl_gene_id == selected_id_gene),]

		#from full_exon_data select those rows with exon data for the [i] gene
		exon_data = full_exon_data_filtered_final[which(full_exon_data_filtered_final$ensembl_gene_id == selected_id_gene),] #4 transcripts are included for SCYL3. In the webpage we have 5, but the last one is not protein coding, because of this we don´t have it.

		#reorder the rows of gene_characteristics to match exon data.
		gene_characteristics = gene_characteristics[match(paste(exon_data$ensembl_transcript_id, exon_data$ensembl_exon_id), paste(gene_characteristics$ensembl_transcript_id, gene_characteristics$ensembl_exon_id)),] #we match using both exon id and transcript id because the same exon id can be in different transcripts, thus we need two conditions for matching. The paste function create strings using the transcript id and the exon id, you can create these strings using start_end_coding_regions and exon_data, then use match to match the rows of exon data to those of start_end_coding_regions. See 'https://stackoverflow.com/questions/47404477/match-with-multiple-criteria-without-loop-in-r' for further details. 

		#remove those rows without coding sequence from exon data
		exon_data = exon_data[which(!is.na(exon_data$genomic_coding_start) & !is.na(exon_data$genomic_coding_end)),]

		#then we remove that rows from gene_characteristics, we select those rows of gene_characteristics that transcripts and exon id combinations included in exon_data. Bear in mind that the same exon ids can be in different transcripts
		gene_characteristics = gene_characteristics[which(paste(gene_characteristics$ensembl_transcript_id, gene_characteristics$ensembl_exon_id) %in% paste(exon_data$ensembl_transcript_id, exon_data$ensembl_exon_id)),]
		#check
		test_rows_initial_datasets_1 = gene_characteristics$ensembl_exon_id %in% exon_data$ensembl_exon_id & gene_characteristics$ensembl_transcript_id %in% exon_data$ensembl_transcript_id
		test_rows_initial_datasets_2 = exon_data$ensembl_exon_id %in% gene_characteristics$ensembl_exon_id & exon_data$ensembl_transcript_id %in% gene_characteristics$ensembl_transcript_id


		#######################################################
		###### EXTRACT CHROMOSOME, GENE SYMBOL, AND IDs #######
		#######################################################

		#extract the chromosome name
		#if the chromosome name matches with the chromosome is currently working the code
		if(all(selected_chr == gene_characteristics$chromosome_name)){

			#save the chromosome name
			chromosome_name = rep(unique(gene_characteristics$chromosome_name), nrow(exon_data))
		} else {

			#set an error
			chromosome_name = rep(NA,  nrow(exon_data))
		} #chromosome name 1 like in the webpage for SCYL3

		#extract the gene name
		hgnc_symbol = rep(unique(gene_characteristics$hgnc_symbol), nrow(exon_data)) #matches the data from the webpage for SCYL3

		#extract the id gene
		gene_id =  exon_data$ensembl_gene_id #matches the data from the webpage for SCYL3

		#check that the number of rows in each dataset (gene_characteristics and exon_data) for the [i] gene are the same
		test_gene_id = rep(nrow(gene_characteristics) == nrow(exon_data), nrow(exon_data))

		#extract the id transcript
		transcript_id =  exon_data$ensembl_transcript_id

		#extract the id exon
		exon_id =  exon_data$ensembl_exon_id #the same exon_id is repeated in different transcripts. I have checked in the webpage and it occurs the same. For example, ENSE00000789668 is included in three different transcripts of the same gene (SCYL3) and it has exactly the same sequence. You can see in the exon map exactly the same region for that exon in the three transcripts (is coding exon number 12 in the transcript ENST00000367772). Because of this, when we match rows between datasets, we cannot use only the exon id, but also the transcript id.



		###################################
		###### EXTRACT GENE BIOTYPE #######
		###################################

		#extract the gene biotype
		#if we only have 1 biotype of gene
		if(length(unique(exon_data$gene_biotype)) == 1 & length(unique(exon_data$transcript_biotype)) == 1){

			#save the biotype
			gene_biotype = exon_data$gene_biotype
			transcript_biotype = exon_data$transcript_biotype #both biotypes are 'protein_coding', so we have removed the last transcript for SCYL3, which is 'Processed transcript'

			#check that everything is ok with biotypes
			#if the gene
			if(unique(exon_data$gene_biotype) == 'protein_coding' & unique(exon_data$transcript_biotype) == 'protein_coding'){

				#everything seems to be ok
				test_gene_transcripts_biotype = rep(TRUE, nrow(exon_data))
			} else { #if not

				#we have a problem
				test_gene_transcripts_biotype = rep(FALSE, nrow(exon_data))
			}
		} else { #if not

			#we have a problem
			gene_biotype = rep(NA, nrow(exon_data))
			transcript_biotype = rep(NA, nrow(exon_data))			
			test_gene_transcripts_biotype = rep(NA, nrow(exon_data))
		}



		#####################################################################################################
		###### CALCULATE GENE LENGTH, MIDDLE POINT AND WINDOWS ACCORDING THE WHOLE REGION OF THE GENE #######
		#####################################################################################################

		#if the start and end position is unique for the whole gene, i.e., we only have one start and end for the whole gene
		if(length(unique(gene_characteristics$end_position)) == 1 & length(unique(gene_characteristics$start_position)) == 1){
		
			#save the gene start and end
			gene_start = gene_characteristics$start_position
			gene_end = gene_characteristics$end_position

			#calculate middle point and gen length
			#set the middle point calculating the mean of the start and stop positions
			middle_point = rep((unique(gene_characteristics$start_position) + unique(gene_characteristics$end_position))/2, nrow(exon_data)) #the length calculated with start and end of the webpage gives the same result for SCYL3 ((169863408+169818772)/2 = 169841090)

			#calculate the difference between start/end and adding 1 to calculate the length between them including both extremes
			gene_length = rep(unique(gene_characteristics$end_position) - unique(gene_characteristics$start_position) + 1, nrow(exon_data)) #the length calculated with start and end of the webpage gives the same result for SCYL3 (169863408-169818772 + 1 = 44637).
				#even if the strand is reversed, the start position occurs before the end position, i.e., the start position coordinate is smaller than the coordinate of the end position. The same goes for the start and end of the exons. The start coordinate is smaller than the end coordinate even if the strand is reverse. 
					#I have checked in SCYL3, which reverse strand: "https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000000457;r=1:169818772-169863408".

			#check that the middle point and the gene length are correctly calculated
			#extract the range or values beginning in the start position and ending in the end position of the gene
			gene_range = unique(gene_characteristics$start_position):unique(gene_characteristics$end_position)
			#if gene length is equal to the length of the gene range and the distance from the middle point to end and start of the gene are the same
			if(unique(gene_length) == length(gene_range) & unique(middle_point) - gene_range[1] == gene_range[length(gene_range)] - unique(middle_point)){
				
				#everything seems perfet
				test_gene_length_center = TRUE
			} else {
				#if not, we have a problem
				test_gene_length_center = FALSE
			}

			
			## calculate windows
			#set the length of the windows: 50, 100, 200, 500 and 1000 kb
			window_length = c(50000, 100000, 200000, 500000, 1000000)

			#create an empty data.frame to save results
			windows_coordinates = data.frame(selected_length=NA, lower_end=NA, upper_end=NA)

			#for each width
			for(j in 1:length(window_length)){

				#select the [j] width
				selected_length = window_length[j]

				#calculate the half of the total length, which will be used to calculate both ends around the gene center
				selected_length_half = (selected_length/2) - 0.5 #we have to subtract 0.5 because we have to count the point in the middle. You have to let 0.5 at each side to have 'space' for the middle point. For example, you have the middle of a gene in 20, and we want a window of 10, if you do 20+5 and 20-5 you will get 15 and 25, 25-15+1 gives 11. You need to remove 0.5 at each side of the center to have exactly the window size you want. 20+4.5=24.5 and 20-4.5=15.5; 24.5-15.5+1=10. In addition, length(15.5:24.5)=10

				#calculate the end and start position of the [j] window
				lower_end = unique(middle_point) - selected_length_half
				upper_end = unique(middle_point) + selected_length_half
				
				#save the results
				windows_coordinates = rbind.data.frame(windows_coordinates, cbind.data.frame(selected_length, lower_end, upper_end))
			}
			#remove the first row with NA
			windows_coordinates = windows_coordinates[-1,]

			#check that the width of each window is adequate
			windows_check_1 = rep(!FALSE %in% c(windows_coordinates$upper_end - windows_coordinates$lower_end + 1 == windows_coordinates$selected_length), nrow(exon_data)) #we calculate the complete length between the upper and lower end, because of this we sum 1, to include both extremes of the interval. We then check that the length of each interval correspond with the selected length used for calculating this interval.

			#check that the difference between each windows size for the lower and the upper extremes matches those differences in general of width size
			check_lower_end=NULL
			check_upper_end=NULL
			for(j in 1:(nrow(windows_coordinates)-1)){ #for each row of the dataset with the windows coordinate... we avoid the last row because there is no additional row to check differences
				
				#select the [j] and [j+1] rows
				selected_row_1 = windows_coordinates[j,]
				selected_row_2 = windows_coordinates[j+1,]

				#check that the lower end is correct
				check_lower_end = append(check_lower_end, selected_row_1$lower_end - selected_row_2$lower_end == (selected_row_2$selected_length - selected_row_1$selected_length)/2) #the difference in length between windows has to be divided by 2, because that difference is divided between both sides of the middle gene. For example, if you have 20 as middle gene, the lower end of a 10 size window would be 20-4.5=15.5, and for a 20 size window would be 20-9.5=10.5. The total difference between the windows is 10 (20-10=10), but the difference at each side of between the center is 5 (9.5-4.5=5)
					#we put selected_row_1$lower_end before the selected_row_2 because the lower end decrees with bigger window lengths, so it would give negative values if you put selected_row_2 first 

				#check that the upper end is correct
				check_upper_end = append(check_upper_end, selected_row_2$upper_end - selected_row_1$upper_end == (selected_row_2$selected_length - selected_row_1$selected_length)/2) #the difference in length between windows has to be divided by 2, because that difference is divided between both sides of the middle gene. For example, if you have 20 as middle gene, the upper end of a 10 window would be 20-4.5=15.5, and for a 20 window would be 20-9.5=10.5. The total difference between the windows is 10 (20-10=10), but the difference at each side of the center is 5 (9.5-4.5=5)
					#here we put the upper end of selected row 2 before, because the upper end is increasing with increased window size
			}

			#save the checks
			windows_check_2 = rep(!FALSE %in% c(check_upper_end), nrow(exon_data))
			windows_check_3 = rep(!FALSE %in% c(check_lower_end), nrow(exon_data))


			## remove decimals
			#after we have checked that the windows are properly calculated, we will floor both ends, in that way we maintain the size of the window we want. The window is not exactly exactly centered in the gene, there is a difference of 1 base if decimals exist, but it is only one base. This is not relevant for windows of 50000 bases or more. 
			windows_coordinates$lower_end = floor(windows_coordinates$lower_end)
			windows_coordinates$upper_end = floor(windows_coordinates$upper_end) #we have to floor both ends, if we floor one and apply ceiling to another, we will increase the size of the window. We want to remain the exact size we want, even if it is not exactly centered in the gene. For example, floor(9.5) would be 9 and ceiling(15.5) would be 16, so we are increasing the size of the window (from 15.5-9.5+1=7 to 16-9+1=8). If we do 15-9+1 we get 7.

			#RESPECT TO THE MIDDLE POSITION: we left the middle position with decimals to have the original middle position used to calculate the lower and the upper ends

			#final check that the window lengths matches what we expect
			windows_check_4 = !FALSE %in% c(windows_coordinates$upper_end - windows_coordinates$lower_end + 1 == windows_coordinates$selected_length) #we want the complete length of the window so we have to sum 1 to include both extremes, if not as the lower extreme is included in both lower and upper_end, we would left out one extreme.
		

			##check if the window surpass the limits of the chromosome
			#extract the length of the selected chromosome
			end_chromosome = chrom_length_ucsc_hg19[which(chrom_length_ucsc_hg19$chromosome == paste("chr", selected_chr, sep="")), which(colnames(chrom_length_ucsc_hg19) == "length_bp")]

			#if the any of the window ends surpass the chromosome limits: the lower end of any window is lower than 1 and hence negative OR the upper end of any window is higher than the length of the chromosome.
			if( TRUE %in% c(windows_coordinates$lower_end < 1) | TRUE %in% c(windows_coordinates$upper_end > end_chromosome) ){

				#if the lower end of any window is lower than 1 and hence negative
				if(TRUE %in% c(windows_coordinates$lower_end < 1)){

					#we set the window coordinates as NA
					windows_coordinates[which(windows_coordinates$lower_end < 1), c("lower_end", "upper_end")] <- NA #we have to remove them instead of trimming, because if you trim, you could have different windows with the same size, and we do not want that. We will have a little bit more of smaller windows, but this should not a problem because this is happening for around 400 genes.			
				}

				#if the upper end of any window is higher than the length of the chromosome.
				if(TRUE %in% c(windows_coordinates$upper_end > end_chromosome)){

					#we set the window coordinates as NA
					windows_coordinates[which(windows_coordinates$upper_end > end_chromosome), c("lower_end", "upper_end")] <- NA #we have to remove them instead of trimming, because if you trim, you could have different windows with the same size, and we do not want that. We will have a little bit more of smaller windows, but this should not a problem because this is happening for around 400 genes.			
				}
			}

			#check that no window surpass the chromosome limits
			windows_check_5 = !TRUE %in% c(windows_coordinates$lower_end < 1) & !TRUE %in% c(windows_coordinates$upper_end > end_chromosome) 
				#no TRUE should exist because no lower_end should be lower than 1 AND should exist because no upper_end should be higher than the length of the chromosome. 
				#you can have NAs because the window is not calculated, but NOT TRUE. TRUE would mean that the window surpass the chromosome limits.


			##save the windows
			#50kb
			lower_end_window_50kb = windows_coordinates[which(windows_coordinates$selected_length == 50000),]$lower_end
			upper_end_window_50kb = windows_coordinates[which(windows_coordinates$selected_length == 50000),]$upper_end

			#100kb
			lower_end_window_100kb = windows_coordinates[which(windows_coordinates$selected_length == 100000),]$lower_end
			upper_end_window_100kb = windows_coordinates[which(windows_coordinates$selected_length == 100000),]$upper_end

			#200kb
			lower_end_window_200kb = windows_coordinates[which(windows_coordinates$selected_length == 200000),]$lower_end
			upper_end_window_200kb = windows_coordinates[which(windows_coordinates$selected_length == 200000),]$upper_end

			#500kb
			lower_end_window_500kb = windows_coordinates[which(windows_coordinates$selected_length == 500000),]$lower_end
			upper_end_window_500kb = windows_coordinates[which(windows_coordinates$selected_length == 500000),]$upper_end

			#1000kb
			lower_end_window_1000kb = windows_coordinates[which(windows_coordinates$selected_length == 1000000),]$lower_end
			upper_end_window_1000kb = windows_coordinates[which(windows_coordinates$selected_length == 1000000),]$upper_end
		} else { #if we have different starts or ends, set NA

			#set NA
			gene_start = rep(NA, nrow(exon_data))
			gene_end = rep(NA, nrow(exon_data))
			middle_point = rep(NA, nrow(exon_data))
			gene_length = rep(NA, nrow(exon_data))
			test_gene_length_center = rep(NA, nrow(exon_data))
			windows_check_1 = rep(NA, nrow(exon_data))
			windows_check_2 = rep(NA, nrow(exon_data))
			windows_check_3 = rep(NA, nrow(exon_data))
			windows_check_4 = rep(NA, nrow(exon_data))
			windows_check_5 = rep(NA, nrow(exon_data))		
			lower_end_window_50kb = rep(NA, nrow(exon_data))
			upper_end_window_50kb = rep(NA, nrow(exon_data))
			lower_end_window_100kb = rep(NA, nrow(exon_data))
			upper_end_window_100kb = rep(NA, nrow(exon_data))
			lower_end_window_200kb = rep(NA, nrow(exon_data))
			upper_end_window_200kb = rep(NA, nrow(exon_data))
			lower_end_window_500kb = rep(NA, nrow(exon_data))
			upper_end_window_500kb = rep(NA, nrow(exon_data))
			lower_end_window_1000kb = rep(NA, nrow(exon_data))
			upper_end_window_1000kb = rep(NA, nrow(exon_data))
		}



		####################################
		###### CALCULATE GENE NUMBER #######
		####################################

		#extract the gene_characteristics_chromosome data for the chromosome (chromosome_name) to get all genes that are included in the window only for the selected chromosome. Note that the same position is in different chromosomes, so we could make an error here
		subset_gene_characteristics_chromosome = all_genes_grch37_exon_filter_final[which(all_genes_grch37_exon_filter_final$chromosome_name==unique(chromosome_name)),]
		
		#remove the gene ids not included in the corresponding chromosome
		subset_gene_characteristics_chromosome$ensembl_gene_id <- droplevels(subset_gene_characteristics_chromosome$ensembl_gene_id) #we remove the levels not used, i.e., genes not included in the corresponding chromosome, because in previous versions we used dplyr to check whether the sequence of each gene is included in the window. This was performed per gene applying the inrange function in each group of rows with the same gene id. Dplyr search for all levels of gene_id even if the gene id not included. We have to removed these unused levels.

		#open an empty data.frame to save the results 
		n_genes_windows = data.frame(selected_length=NA, lower_end=NA, upper_end=NA, n_genes=NA, check_n_genes_0=NA, check_n_genes_1=NA, check_n_genes_2=NA, check_n_genes_3=NA, check_n_genes_4=NA)

		#for each window
		for(j in 1:nrow(windows_coordinates)){

			#select the [j] window
			selected_window = windows_coordinates[j,]

			#we do not have NA in any of the window limits we perform the calculations
			if(!is.na(selected_window$lower_end) & !is.na(selected_window$upper_end)){

				#select those genes whose start OR end is equal/higher than the start of the window AND equal/smaller than the end of the window. 
				genes_inside_window = subset_gene_characteristics_chromosome[which(subset_gene_characteristics_chromosome$start_position >= selected_window$lower_end & subset_gene_characteristics_chromosome$start_position <= selected_window$upper_end | subset_gene_characteristics_chromosome$end_position >= selected_window$lower_end & subset_gene_characteristics_chromosome$end_position <= selected_window$upper_end | subset_gene_characteristics_chromosome$start_position < selected_window$lower_end & subset_gene_characteristics_chromosome$end_position > selected_window$upper_end),] #in that way we can include all the genes within the window, even if only a part of them (end or start) is included. The first two conditions select those genes whose start is within the window but the end could be included or not, whilst the next two conditions select genes that have their end inside the window but the start could be or not included. See figure 18 for further details.
					#we add an additional condition, which is to include genes with a start position lower than the lower end and a end position higher than the upper end. This is made for genes that are too big, for example C1orf112, when this gene is the center, the range of 50 kb will be smaller than the whole gene, so the start and the end are outside of the range, but because the windows is in the middle of the windows sides. I have checked that the first window when C1orf112 is the central gene does not include C1orf112 if the third condition is not included. C1orf112 starts before the end of the 50kb window and ends after the end of that window. the problem is solved adding the third condition.
					#NOTE: We should not have any problem with the reverse strand, because gen start/end (and exon start/end and genomic start/end) follow the forward values. For example, in the first transcript (ENST00000367771) of SCYL3 (ENSG00000000457), the last coding exon according to cds_start/end and the webpage, has the lowest values of exon_chrom_start and end, the same for genmomic_codin_start and end gene start/end. Therefore, start of gene, exon and coding region follow the coordinates of the positive strand.
					#FOR CODING DENSITY WOULD BE THE SAME, BUT USING EXON DATA!!
				
				#extract the genes outside the window to make a check
				genes_outside_window = subset_gene_characteristics_chromosome[which(!subset_gene_characteristics_chromosome$ensembl_gene_id %in% genes_inside_window$ensembl_gene_id),] #we don´t use the opposite selection with higher-lower etc... because some genes can be part outside and part inside of the gene window, they will be considered out and in. We directly select those genes not included inside the window. 
				check_n_genes_0 = nrow(genes_inside_window) + nrow(genes_outside_window) == nrow(subset_gene_characteristics_chromosome) #same number.

				#extract the number inside the window, whcih are unique cases of gene id
				n_genes = length(unique(genes_inside_window$ensembl_gene_id)) #we are considering any gene with end or start inside the window even the complete gene is not included in the window

				#checks n_genes
				#we only make the check if case with no hgnc symbol are present. For example, in the windows of SCYL3, we have 40 gene ids, but 37 hgnc symbols. This difference is caused because 4 genes has no hgnc symbol, hgnc_symbol variable is empty. Therefore, we have 36 hgnc real symbols, them a value of '' for the first gene without hgnc symbol (37) and the three more ids without symbol, summing up to 40.
				if(!'' %in% genes_inside_window$hgnc_symbol){
					check_n_genes_1 = length(unique(genes_inside_window$hgnc_symbol)) == n_genes #check that the number of genes matches the number of unique hgnc symbols
						#NOTE that we cannot make a check with the exon data dataset because in full_exon_data_filtered_final, non coding transcripts are already removed. So the first/last exon can be not included if it is non-coding and hence the corresponding gene would not be included, whilst it would be included considering the full gene length (gene start/end)
				} else {

					check_n_genes_1 = NA
				}

				#other check
				check_n_genes_2 = all(genes_inside_window$start_position - selected_window$lower_end >= 0 & selected_window$upper_end - genes_inside_window$start_position >= 0 | genes_inside_window$end_position - selected_window$lower_end >= 0 & selected_window$upper_end - genes_inside_window$end_position >= 0 | genes_inside_window$start_position - selected_window$lower_end < 0 & selected_window$upper_end - genes_inside_window$end_position < 0) #check that within the genes selected, the start position is equal or higher to the lower end of the window and AND equal or lower to the upper end OR the end position is equal or higher than the lower end AND equal or lower to the upper end. We check this making a rest between the start and end positions and the lower/upper ends: 
					#The gene start is inside the window if these two conditions does not occur:
						#a) if the difference between gen start position and lower end is negative, then the lower end is higher than the start position, i.e., the gene begins before the window does it; 
						#b) if the the difference between the upper end and the gene start position is negative, then the start is higher than the end, and hence the gene begins after the windows ends; 
					#The gene end is inside the window if these two conditions does not occur:
						#c) if the difference between the end of the gene and the lower end of the window is negative, then the lower end is higher than the end of the gene, and hence the gene ends after the window starts; 
						#d) if the difference between the upper end and the end of the gene is negative, then the end position of the gene is higher and then that gene ends after the window does it;
					#The gene starts before the start of the window and the gene ends after the window is finished. Therefore, even the start/end of the gene is not included in the window, the whole window is occupied by the gene.
						#f) if the difference between the start position of the gene and the  lower end of the window is negative, then the lower end is higher that the gene start and hence the gene starts before the window does it.
						#g) if the difference between the upper end of the window and the end position of the gene is negative, then the end of the gene is higher than the end of the window and hence the gene ends after the window does it.
					#ALL genes have to satisfy the following: 
						#That the start is between the upper and lower extremes, i.e., it is equal or higher than the lower limit and equal or lower than the upper limit
							#OR
						#That the end is between the upper and lower extremes, i.e., it is equal or higher than the lower limit and equal or lower than the upper limit. 
							#OR
						#The start occurs before the window and the end after the window. Therefore, the exon occupy the whole length of the window.

				#check that there is no gene outside the window that satisfy the conditions used in check_n_genes_2. Genes with part outside and inside has been removed from genes_outside_window
				check_n_genes_3_v1 = genes_outside_window$start_position - selected_window$lower_end >= 0 & selected_window$upper_end - genes_outside_window$start_position >= 0 | genes_outside_window$end_position - selected_window$lower_end >= 0 & selected_window$upper_end - genes_outside_window$end_position >= 0 | genes_outside_window$start_position - selected_window$lower_end < 0 & selected_window$upper_end - genes_outside_window$end_position < 0
				#check that all is FALSE, i.e., not TRUE is included
				check_n_genes_3 = ifelse(!TRUE %in% check_n_genes_3_v1 & all(!is.na(check_n_genes_3_v1)), TRUE, FALSE) #if all checks are FALSE and with have no NA, give me TRUE


				##Final check
				#we will calculate the number of genes inside the [j] window using genomic ranges

				#create the reference range using the window start and end
				window_range_check_gene_number = IRanges(start=selected_window$lower_end, end=selected_window$upper_end) #create a IRange object with the start and end of the [j] window

				#convert to a GRanges object. We need a GRanges object to select those exons in the the subject GRange file that included in the reference range. 
				window_range_check_gene_number_gr = GRanges(seqnames=unique(chromosome_name), ranges=window_range_check_gene_number, strand = '*') #seqnames are the names of the chromosome names. This is very important, because if we select as seq names the exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping using disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified in the exon data (I have checked the no impact of strand in this exact example for gap length).
					#the information was obtained from the script of the coding density (see below).

				#remove the duplicate rows for the same gene
				subset_gene_characteristics_chromosome_no_duplicated = subset_gene_characteristics_chromosome[which(!duplicated(subset_gene_characteristics_chromosome$ensembl_gene_id)),] #if you have several rows with the same gene, IRanges will consider them as different sequences.
				#check
				nrow(subset_gene_characteristics_chromosome_no_duplicated) == length(unique(subset_gene_characteristics_chromosome$ensembl_gene_id))

				#select the gene ranges
				gene_ranges = IRanges(start=subset_gene_characteristics_chromosome_no_duplicated$start_position, end=subset_gene_characteristics_chromosome_no_duplicated$end_position) #create a IRange object with the start and end of genes in the selected chromosome. 

				#convert the IRange file into a genomic range file
				gene_ranges_gr = GRanges(seqnames=unique(chromosome_name), ranges=gene_ranges, strand = '*') #sequnames are the names of the chromosome names. This is very important, because if we select as seq names de exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping using disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified in the exon data (I have checked the no impact of strand in this exact example for gap length).
						#this explanation is taken from gene_coordinates_v8.r. But for the 1000kb window of "ENSG00000058453" (CROCC) I have checked that putting the strand of the all genes ("subset_gene_characteristics_chromosome_no_duplicated") does not change anything. These are genes with positive and negative strand, but no change is produced.

				#we do not remove overlapping between the tested ranges because we want to count every gen that is included within the window. If you remove overlapping regions, for a given region you lost the different genes that are present there and get only 1 count. We do not want that. 

				#find what of the tested ranges and are overlapped with the reference range
				hits_gene_number = findOverlaps(window_range_check_gene_number_gr, gene_ranges_gr)

				#extract those tested ranges that hit and their number
				n_genes_2 = length(gene_ranges_gr[subjectHits(hits_gene_number)]@ranges)
					#subjectHits extracts what of the tested ranges are overlapped with the reference ranges.

				#check if the result is the same than the original gene number
				check_n_genes_4 = n_genes == n_genes_2
			} else { #if you have NA for the window limits

				#Set NA for everything, except selected window
				n_genes = NA
				check_n_genes_0 = NA
				check_n_genes_1 = NA
				check_n_genes_2 = NA
				check_n_genes_3 = NA
				check_n_genes_4 = NA
			}
			
			#save the results in a data.frame
			n_genes_windows = rbind.data.frame(n_genes_windows, cbind.data.frame(selected_window, n_genes, check_n_genes_0, check_n_genes_1, check_n_genes_2, check_n_genes_3, check_n_genes_4))
		}
		#remove first row with NAs
		n_genes_windows = n_genes_windows[-1,] #I have checked for SCYL3 (ENSG00000000457) the first two windows. The first window (50kb) includes two genes (SCYL3 and C1orf112), whilst the second one also includes KIFAP3. I have check the map around SCYL3 (http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000000457;r=1:169817880-169864300), and i have seen two genes closer to SCYL3 than C1orf112, which are METTL18 and SELL. However, I have checked that the end of these genes (169680839 and 169764107) occur before the start of the two first windows (169816090.5 and 169791090.5). In contrast, C1orf112 starts very far away from SCYL3, but its end (169823221) reached the start of these windows. In the next window, which is wider and hence starts before (169741090.5), METTL18 is included. The next windows (500kb) already include SELE, SELP and SEL. Finally, the widest window (100 kb) also includes BLZF1, CCDC181, SLC19A2, F5, METTL11B. In that way, we have included all the coding genes you can see in the window. If you scroll to both sides, you can see as no more coding genes are present between the boundaries of the last window (169341090.5 - 170341089.5). The rest are pseudogenes, RNa genes and processed transcripts.

		#save the number of genes per each window
		n_genes_50kb = rep(n_genes_windows[which(n_genes_windows$selected_length == 50000),]$n_genes, nrow(exon_data))
		n_genes_100kb = rep(n_genes_windows[which(n_genes_windows$selected_length == 100000),]$n_genes, nrow(exon_data))
		n_genes_200kb = rep(n_genes_windows[which(n_genes_windows$selected_length == 200000),]$n_genes, nrow(exon_data))
		n_genes_500kb = rep(n_genes_windows[which(n_genes_windows$selected_length == 500000),]$n_genes, nrow(exon_data))
		n_genes_1000kb = rep(n_genes_windows[which(n_genes_windows$selected_length == 1000000),]$n_genes, nrow(exon_data))


		##summarize all the gene_number_checks
		#create a vector with the check to be summarized
		gene_number_checks_to_summarize = c("check_n_genes_0","check_n_genes_1","check_n_genes_2","check_n_genes_3","check_n_genes_4")

		#for each check to be summarized
		for(j in 1:length(gene_number_checks_to_summarize)){

			#select the [j] check to summarize
			selected_gene_number_check = gene_number_checks_to_summarize[j]

			#extract the corresponding column
			selected_column_gene_number_check = n_genes_windows[,which(colnames(n_genes_windows) == selected_gene_number_check)]	

			#if no FALSE is present in the check
			if(!FALSE %in% selected_column_gene_number_check){

				#if there is no FALSE in is.na() and hence all cases are NA
				if(!FALSE %in% is.na(selected_column_gene_number_check)){

					#set the check as TRUE
					summary_gene_number_check = rep(NA, nrow(exon_data))
				} else { #if not, and hence not all cases are NA
					
					#if you have NAs and TRUEs
					if(TRUE %in% is.na(selected_column_gene_number_check) & TRUE %in% selected_column_gene_number_check){
						
						#set the check as TRUE/NA
						summary_gene_number_check = rep("TRUE/NA", nrow(exon_data))
					} else { #if not, we dot not have any NA

						#set the check as TRUE
						summary_gene_number_check = rep(TRUE, nrow(exon_data))
					}
				}
			} else { #if we have FALSE

				summary_gene_number_check = rep(FALSE, nrow(exon_data))
			}

			#assign the value of the summary check to the corresponding variable.
			assign(selected_gene_number_check, summary_gene_number_check) #see second response in "https://stackoverflow.com/questions/28909191/treat-string-as-object-name-in-a-loop-in-r" for further details about assign		
		}	



		#######################################
		###### CALCULATE CODING DENSITY #######
		#######################################

		#extract the exon data data for the chromosome (chromosome_name) to get all exons that are included in the window only for the selected chromosome. Note that the same position is in different chromosomes, so we could make an error here.
		subset_exon_data_chromosome = full_exon_data_filtered_final[which(full_exon_data_filtered_final$chromosome_name==unique(chromosome_name)),]
		
		#open an empty data.frame to save the results 
		coding_density_windows = data.frame(selected_length=NA, lower_end=NA, upper_end=NA, check_coding_density_1=NA, check_coding_density_2=NA, check_coding_density_4=NA, check_coding_density_5=NA, check_coding_density_6=NA, check_coding_density_7=NA, test_all_ranges_included=NA, test_iranges=NA, test_start_end=NA, test_ranges_overlap=NA, test_ranges_overlap_2=NA, coding_density=NA) #IMPORTANT: the check_coding_density_3 was in version 6, but in 7 was removed because it is very slow. I have maintained the order of the checks to remember that check 3 was that removed.

		#for each window
		for(j in 1:nrow(windows_coordinates)){

			#select the [j] window
			selected_window = windows_coordinates[j,]

			#we do not have NA in any of the window limits we perform the calculations. Few lines below there is a ifelse for separating cases with exons inside the window than those without exons. That ifelse is different from this one, because you have window but no exons, here we do not have window. 
			if(!is.na(selected_window$lower_end) & !is.na(selected_window$upper_end)){

				#select the gaps for the [i] gene and the [j] window
				selected_gap = gap_results_extracted_final[which(gap_results_extracted_final$gene_id == selected_id_gene), which( colnames(gap_results_extracted_final) == paste("length_gaps_", (selected_window$selected_length / 1000), "kb", sep="")) ,] #select the gaps (row) for the [i] gene. Select the gap (column) for the [k] window. We divided the length of the window by 1000 to have the length in kb and match the colnames in gap_results_extracted_final. 

				#select those exons whose start of the coding sequence OR its end is equal/higher than the start of the window AND equal/smaller than the end of the window. 
				exons_inside_window = subset_exon_data_chromosome[which(subset_exon_data_chromosome$genomic_coding_start >= selected_window$lower_end & subset_exon_data_chromosome$genomic_coding_start <= selected_window$upper_end | subset_exon_data_chromosome$genomic_coding_end >= selected_window$lower_end & subset_exon_data_chromosome$genomic_coding_end <= selected_window$upper_end | subset_exon_data_chromosome$genomic_coding_start < selected_window$lower_end & subset_exon_data_chromosome$genomic_coding_end > selected_window$upper_end),]
					#We select all exons that satisfy the following conditions:
						#Genomic start is equal/higher than the lower limit of the window AND equal/lower than the upper limit. Therefore, it is included between the window.
						#Genomic end is equal/higher than the lower limit of the window AND equal/lower than the upper limit. Therefore, it is included between the window.
						#Genomic start is lower than the lower limit AND higher than the upper limit. Therefore, it is included between the window. 
					#In that way we can include all the genes within the window, even if only a part of them (end or start) is included. The first two conditions we include those exons whose genomic start is within the window but their end could not be, whilst the next two conditions select genes that have their end inside the window but the start could be outside. See figure 19 for further details.
					#We add an additional condition, which is to include exons with a genomic start position lower than the lower end of the window and a genomic end position higher than the upper end of the window. This is made just in case a gene could have very big exon, bigger than the window centered in that gene. I do not know right know about this case, but just in case. In the case of big genes I have the case of C1orf112....
					#NOTE: We should not have any problem with the reverse strand, because gen start/end (and exon start/end and genomic start/end) follow the forward values. For example, in the first transcript (ENST00000367771) of SCYL3 (ENSG00000000457), the last coding exon according to cds_start/end and the webpage, has the lowest values of exon_chrom_start and end, the same for genmomic_codin_start and end gene start/end. Therefore, start of gene, exon and coding region follow the coordinates of the positive strand.
					
					#see figure 19 for further details

					#we use genomic start and end for getting ALL exons within the window. The same exon can be select several times, as the same exon id can be in different transcripts, bu this should not be a problem as we will calculate the coding sequence inside the window with not overlapped ranges (with genomic ranges), so two equal exons would share the same area, but that area will not be included two times. Even I have seen the same region shared by two different genes, like the beginning of SCYL3 and the end of C1orf112. Again, if the overlap is calculated, this should not be a problem.

				#if we have exons inside the window
				if(nrow(exons_inside_window) > 0){ #It is possible that no coding exon exits inside the window. For example when RAB27B (ENSG00000041353: chromosome 18, i=2) is the center, there is no other gene inside the 50kb window. Moreover, that gene has a very big intron in the middle, leading that the 50kb window has no coding exon inside: http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000041353;r=18:52381538-52566300. In the case of gene number, we do not need this, because if you center the window in the middle of a gene, you are for sure in the middle of a gene, so that counts for 1, even if the gene is very very small.

					#extract the genes outside the window to make a check
					exons_outside_window = subset_exon_data_chromosome[which(!subset_exon_data_chromosome$ensembl_exon_id %in% exons_inside_window$ensembl_exon_id),] #we don´t use the opposite selection with higher-lower etc... because some exons can be part outside and part inside of the window, they will be considered out and in. We directly select those exons not included inside the window. 
					nrow(exons_inside_window) + nrow(exons_outside_window) == nrow(subset_exon_data_chromosome) #same number.

					#other check_coding_density_1
					check_coding_density_1 = all(exons_inside_window$genomic_coding_start - selected_window$lower_end >= 0 & selected_window$upper_end - exons_inside_window$genomic_coding_start >= 0 | exons_inside_window$genomic_coding_end - selected_window$lower_end >= 0 & selected_window$upper_end - exons_inside_window$genomic_coding_end >= 0 | exons_inside_window$genomic_coding_start - selected_window$lower_end < 0 & selected_window$upper_end - exons_inside_window$genomic_coding_end < 0) #check that within the exons selected, the genomic start position is equal or higher to the lower end of the window and AND equal or lower to the upper end OR the genomic end position is equal or higher than the lower end AND equal or lower to the upper end. We check this making a rest between the genomic start and end positions and the lower/upper ends: 
						#The genomic start is inside the window if these two conditions does not occur:
							#a) if the difference between genomic start position and lower end is negative, then the lower end is higher than the start position, i.e., the coding region begins before the window does it; 
							#b) if the the difference between the upper end and the gene start position is negative, then the start is higher than the end of the window, and hence the gene begins after the windows ends; 
						#The end is outside the window if these two conditions does not occur:
							#c) if the difference between the end of the coding region of the exon and the lower end of the window is negative, then the lower end is higher than the end of the coding region, and hence this region ends after the window starts; 
							#d) if the difference between the upper end and the end of the genomic region is negative, then the end position of the genomic region is higher and then that region ends after the window does it;
						#The coding region starts before the start of the window and that region ends after the window is finished. Therefore, even the start/end of the gene is not included in the window, the whole window is occupied by the gene.
							#f) if the difference between the genomic start position of the exon and the  lower end of the window is negative, then the lower end is higher that the genomic start and hence that coding region starts before the window does it.
							#g) if the difference between the upper end of the window and the genomic end position of the exon is negative, then the end of the genomic region is higher than the end of the window and hence that region ends after the window does it.
						#ALL exons have to satisfy the following: 
							#That the genomic coding start is between the upper and lower extremes, i.e., it is equal or higher than the lower limit and equal or lower than the upper limit
								#OR
							#That the genomic coding end is between the upper and lower extremes, i.e., it is equal or higher than the lower limit and equal or lower than the upper limit. 
								#OR
							#The genomic starts occurs before the window and the genomic end after the window. Therefore, the exon occupy the whole length of the window.

					#check that there is no gene outside the window that satisfy the conditions used in check_coding_density_1. Genes with part outside and inside has been removed from exons_outside_window
					check_coding_density_2_v1 = exons_outside_window$genomic_coding_start - selected_window$lower_end >= 0 & selected_window$upper_end - exons_outside_window$genomic_coding_start >= 0 | exons_outside_window$genomic_coding_end - selected_window$lower_end >= 0 & selected_window$upper_end - exons_outside_window$genomic_coding_end >= 0 | exons_outside_window$genomic_coding_start - selected_window$lower_end < 0 & selected_window$upper_end - exons_outside_window$genomic_coding_end < 0
					#check that all is FALSE, i.e., not TRUE is included
					check_coding_density_2 = ifelse(!TRUE %in% check_coding_density_2_v1, TRUE, FALSE) #if all checks are FALSE give me TRUE. We remove the filter of the NAs respect to the code of gene number, because we can have exons without coding sequence (UTR tails). Even if some of these exons should be inside the window, they are not coding, so we are not interested in them for calculating coding density. In the case of gene number, all cases have a gene id, so NAs were not acceptable, and hence we have to check that they were not present.


					## from exons with part of the coding region outside and inside of the window, we remove the part outside

					#copy the exons inside the window to do some operations (remove parts outside of the window)
					exons_inside_window_copy <- exons_inside_window

					#In case you want to check this selection work, I have created code for generating exons that satisfy these conditions, i.e., that are in the frontier of the window
						#the genomic start is 100 bases before the start of the window
							#exons_inside_window_copy[3,]$genomic_coding_start <- selected_window$lower_end - 100
					
						#NOT RUN. the genomic start is 100 bases after the end of the window. We DO NOT use this because this scenario should not occur. If the start occurs after the window, the whole exon is outside the window and should be removed in the previous filter
							##exons_inside_window_copy[7,]$genomic_coding_start <- selected_window$upper_end + 100
					
						#the genomic end is 100 bases after the end of the window
							#exons_inside_window_copy[14,]$genomic_coding_end <- selected_window$upper_end + 100
					
						#NOT RUN. the genomic end is 100 bases before the start of the window. We DO NOT use this because this scenario should not occur. If the end occurs before the window, the whole exon is outside the window and should be removed in the previous filter
							##exons_inside_window_copy[9,]$genomic_coding_end <- selected_window$lower_end - 100

						#the genomic start is 100 bases before the start of the window AND the genomic end is 100 bases after the end of the window
							#exons_inside_window_copy[17,]$genomic_coding_start <- selected_window$lower_end - 100
							#exons_inside_window_copy[17,]$genomic_coding_end <- selected_window$upper_end + 100

					#select those exons whose genomic coding start occurs before the start of the window
					exons_inside_window_start_out = which(exons_inside_window_copy$genomic_coding_start < selected_window$lower_end & exons_inside_window_copy$genomic_coding_end >= selected_window$lower_end & exons_inside_window_copy$genomic_coding_end <= selected_window$upper_end) #we select those exons whose genomic start occurs before the start of the window but their genomic end is included in the window. In that way, we removed those with start before and end after. Those will be included in exons_inside_window_start_end_out.

					#select those exons whose genomic coding end occurs after the end of the window
					exons_inside_window_end_out = which(exons_inside_window_copy$genomic_coding_end > selected_window$upper_end & exons_inside_window_copy$genomic_coding_start >= selected_window$lower_end & exons_inside_window_copy$genomic_coding_start <= selected_window$upper_end) #we select those exons whose genomic end occurs after the end of the window but their genomic start is included in the window. In that way, we removed those with start before and end after. Those will be included in exons_inside_window_start_end_out.

					#select those exons whose genomic coding start occurs before the start of the window and genomic coding end occurs after the end of the window
					exons_inside_window_start_end_out = which(exons_inside_window_copy$genomic_coding_start < selected_window$lower_end & exons_inside_window_copy$genomic_coding_end > selected_window$upper_end)

					#if we have exons with the start before the beginning of the window
					#open check v4
					check_coding_density_4 = NULL
					if(length(exons_inside_window_start_out)>0){
				
						##for those exons with the start before the window but the end is inside the window, we set the start coding region at the start of the window, while the end of the coding region remains without change, as it is being included in the window.
						exons_inside_window_copy[exons_inside_window_start_out,]$genomic_coding_start <- selected_window$lower_end
					
						#check
						check_coding_density_4 = append(check_coding_density_4, exons_inside_window_copy[exons_inside_window_start_out,]$genomic_coding_start == selected_window$lower_end)
					}

					#if we have exons with the end after the beginning of the window
					if(length(exons_inside_window_end_out)>0){
					
						##for those exons with the end after the window but the start is inside the window, we set the end of the coding region at the end of the window, while the start of the coding region remains without change, as it is being included in the window.
						exons_inside_window_copy[exons_inside_window_end_out,]$genomic_coding_end <- selected_window$upper_end
					
						#check
						check_coding_density_4 = append(check_coding_density_4, exons_inside_window_copy[exons_inside_window_end_out,]$genomic_coding_end == selected_window$upper_end)
					}

					#if we have exons with the start before the beginning of the window the end after the beginning of the window
					if(length(exons_inside_window_start_end_out)>0){

						##for those exons with the start before the window and the end after the window, we set the start coding region at the start of the window and the end of the coding region at the end of the window. As the exon occupy the whole window, all the window is coding
						exons_inside_window_copy[exons_inside_window_start_end_out,]$genomic_coding_start <- selected_window$lower_end
						exons_inside_window_copy[exons_inside_window_start_end_out,]$genomic_coding_end <- selected_window$upper_end
					
						#check
						check_coding_density_4 = append(check_coding_density_4, exons_inside_window_copy[exons_inside_window_start_end_out,]$genomic_coding_start == selected_window$lower_end)
						check_coding_density_4 = append(check_coding_density_4, exons_inside_window_copy[exons_inside_window_start_end_out,]$genomic_coding_end == selected_window$upper_end)
					}

					#summarize checks
					if(length(check_coding_density_4) > 0){ #if we have checks for the check_4 (i.e., some exons are in part outside and inside of the window)
					
						#if all checks are TRUE, then TRUE for the v4 check
						check_coding_density_4 = all(check_coding_density_4)
					} else {

						#save the checks as zero
						check_coding_density_4 = NA
					}

					#Now that we have applied filters to removed those areas outside the windows in exons that have part in part out of the window, we can check that we do not have any exon with start and/or end outside of the windows.
					check_coding_density_5 = !TRUE %in% c(exons_inside_window_copy$genomic_coding_start < selected_window$lower_end | exons_inside_window_copy$genomic_coding_end < selected_window$lower_end | exons_inside_window_copy$genomic_coding_start > selected_window$upper_end | exons_inside_window_copy$genomic_coding_end > selected_window$upper_end)


					## calculate the not overlapping ranges of coding sequences considering all transcripts
					#create a IRange object with the start and end of each of the exon sequences of the [i] gene
					#require(IRanges)
					exons_ranges = IRanges(start=exons_inside_window_copy$genomic_coding_start, end=exons_inside_window_copy$genomic_coding_end) #we can see the start and end of each exon sequence along with the width of each one.
			
					#check that the row number is ok
					test_all_ranges_included = length(exons_ranges) == nrow(exons_inside_window_copy)

					#check that the start of each exon in the IRange object is the same than in the original data.frame. Also check that the sum of start more width less 1 (the start point is already included) is equal to the end point in the original data frame. Also check the exons end less the exons starts plus 1 (because we want the length of the exons including both extremes) is equal to the with calculated in iranges.
					test_iranges = all(exons_ranges@start == exons_inside_window_copy$genomic_coding_start & (exons_ranges@start + exons_ranges@width - 1) == exons_inside_window_copy$genomic_coding_end & exons_ranges@width == (exons_inside_window_copy$genomic_coding_end - exons_inside_window_copy$genomic_coding_start + 1))

					#create a GRanges object with the ranges of all exon sequences of the [i] gene consindering the all the overlapping ranges
					gr = GRanges(seqnames = exons_inside_window_copy$chromosome_name, ranges = exons_ranges, strand = '*') #seqnames are the names of the chromosome names. This is very important, because if we select as seq names the exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have check that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified.
					overlap_sequences = reduce(gr) #‘reduce’ returns an object of the same type as ‘x’ containing reduced ranges for each distinct (seqname, strand) pairing. So you can calculate the overlapped ranges for each chromosome or strand.\

					#check that the start of the first range matches with the start of the first exon within the window. We also check that the end of the last range matches the end of the last exon

					#extract first the exon with the first coding region inside the window and the exon with the last coding region within the window
					start_first_exon = unique(exons_inside_window_copy[which(exons_inside_window_copy$genomic_coding_start == min(exons_inside_window_copy$genomic_coding_start)),]$genomic_coding_start)
					end_last_exon = unique(exons_inside_window_copy[which(exons_inside_window_copy$genomic_coding_end == max(exons_inside_window_copy$genomic_coding_end)),]$genomic_coding_end) #we use unique because we can have the same exon in different transcripts, and hence the same exon can be the first or the last one
				
					#make the test
					test_start_end = min(overlap_sequences[,]@ranges@start) == start_first_exon & max(overlap_sequences[,]@ranges@start + overlap_sequences[,]@ranges@width - 1) == end_last_exon #the end of each range is calculated in overlaped_sequences using the start and the width of each range (start + width subtracting 1 because the start is already included in width, so you have it two times)

					#additional checks
					#check that no exons is not outside of the not overlapped ranges
					#if we have more than 1 range
					if(!length(overlap_sequences@ranges) == 1){
				
						#check no exon has been left outside of the ranges, i.e., is not included between three extremes of the created ranges
						test_ranges_overlap = NULL
						for(k in 1:length(overlap_sequences@ranges)){

							#if [k] is NOT the last row
							if(!k==length(overlap_sequences@ranges)){

								#select the [k] row of the range data and the next one
								selected_row_1 = overlap_sequences[k,]@ranges
								selected_row_2 = overlap_sequences[k+1,]@ranges

								#test that no exon is excluded form the ranges
								result_test_overlap = c(
									#any exon start position is higher than the end of the [k] created range and lower than the start of the [k+1] created range?
									exons_inside_window_copy$genomic_coding_start > (selected_row_1@start + selected_row_1@width - 1) & exons_inside_window_copy$genomic_coding_start < selected_row_2@start,
									#any exon end position is higher than the end of the [k] created range and lower than the start of the [k+1] created range?
									exons_inside_window_copy$genomic_coding_end > (selected_row_1@start + selected_row_1@width - 1) & exons_inside_window_copy$genomic_coding_end < selected_row_2@start)
							} else { #if not and hence this is the last row

								#select the [k] row of the range data and the next one
								selected_row = overlap_sequences[k,]@ranges

								#test that no exon is excluded form the ranges
								result_test_overlap = c(
									#the start positions of exons are higher than the end of the last range
										exons_inside_window_copy$genomic_coding_start > (selected_row@start + selected_row@width - 1),
										#the end positions of exons are higher than the end of the last range
										exons_inside_window_copy$genomic_coding_end > (selected_row@start + selected_row@width - 1))
							}

							#are false? i.e., the start and end positions are not included in the regions between created ranges?
							result_test_overlap_final = result_test_overlap == FALSE

							#save the results
							test_ranges_overlap = append(test_ranges_overlap, result_test_overlap_final)
						}

						#if all the test have been made (for each range we have 2 tests of TRUE/FALSE for each exon)
						if(length(test_ranges_overlap) == (length(overlap_sequences@ranges) * nrow(exons_inside_window_copy) * 2)){

							#calculate if all are true
							test_ranges_overlap = all(test_ranges_overlap)
						} else {

							#if not NA
							test_ranges_overlap = NA
						}	
					} else { #if not and hence we only have 1 row, we check that all start and end exons are included within that unique range

						#select the unique created range we have of the range data and the next one
						selected_ranges = overlap_sequences[,]@ranges

						#test that no exon is excluded form the ranges
						test_ranges_overlap = c(
								#the start position of the unique coding exon is lower than the start of the unique created range OR higher than the end of that unique created range?
								exons_inside_window_copy$genomic_coding_start < selected_ranges@start | exons_inside_window_copy$genomic_coding_start > (selected_ranges@start + selected_ranges@width - 1),
								#the end position of the unique coding exon is lower than the start of the unique created range OR higher than the end of that unique created range?	
								exons_inside_window_copy$genomic_coding_end < selected_ranges@start | exons_inside_window_copy$genomic_coding_end > (selected_ranges@start + selected_ranges@width - 1))

						#are false? i.e., the start and end positions are not included in the regions outside the created range?
						test_ranges_overlap = test_ranges_overlap == FALSE

						#if all the test have been made (for each range we have 2 tests of TRUE/FALSE for each exon)
						if(length(test_ranges_overlap) == (length(overlap_sequences@ranges) * nrow(exons_inside_window_copy) * 2)){

							#calculate if all are true
							test_ranges_overlap = all(test_ranges_overlap)
						} else {

							#if not NA
							test_ranges_overlap = NA
						}
					}

					#now check that any of the not overlapped range is not overlapped with the rest
					if(!length(overlap_sequences@ranges) == 1){ #if we have more than one range
						#require(data.table)
						#for each range we make another test
						test_ranges_overlap_2_v1 = NULL
						for(k in 1:length(overlap_sequences)){

							#select the [k] range to be tested
							range_to_test = overlap_sequences[k,]@ranges

							#extract the start and end of range [k] to be tested
							start_range_to_test = range_to_test@start
							end_range_to_test = (range_to_test@start + range_to_test@width - 1) #we have to subtract 1 because stars is included in width, so you have the start two time, and you have to remove it one time the get the last point of the range

							#copy the overlap_sequences and remove the range to be tested
							overlap_sequences_without_tested_range = overlap_sequences[-k,]@ranges #these will be the range where we will check if the [k] range is included

							#extract the starts and ends of the ranges that are not the [k] range
							starts_to_compare = overlap_sequences_without_tested_range@start
							ends_to_compare = (overlap_sequences_without_tested_range@start + overlap_sequences_without_tested_range@width - 1)


							##we are using exactly the same approach to extract the genes/exons with part of their sequence inside each window. In this case we check if the sequence of each range is included in the sequence of the other ranges. If the ranges are not overlapped, we would expect that this is not the case.

							#check if the start of the range to test ([k] range) is included within the start and the end of any the other ranges (ranges to compare)
							start_within_ranges_to_compare = start_range_to_test >= starts_to_compare & start_range_to_test <= ends_to_compare

							#check if the end of the range to test ([k] range) is included within the start and the end of any the other ranges (ranges to compare)
							end_within_ranges_to_compare = end_range_to_test >= starts_to_compare & end_range_to_test <= ends_to_compare

							#check if the start of the range to test ([k] range) is before the start of any the other ranges (ranges to compare) and the end the that range to test is after the end of any of the ranges to compare
							start_end_both_sides_ranges_to_compare = start_range_to_test < starts_to_compare & end_range_to_test > ends_to_compare #here we are considering the possibility that the range to be tested [k] is very big and completely includes other ranges within it

								#this approach of comparing the start or end of one range with a vector of ranges works. Se for example this, we have a three ranges starting at 4, 5, and 6 respectively, while they end at 8, 9, 10. You can check if the following formula works: X >= c(4,5,6) & X <= c(4,5,6). 
									#3 >= c(4,5,6) & 3 <= c(8,9,10)
										#FALSE FALSE FALSE

									#4 >= c(4,5,6) & 4 <= c(8,9,10)
										#TRUE FALSE FALSE

									#5 >= c(4,5,6) & 5 <= c(8,9,10)
										#TRUE  TRUE FALSE
									
									#6 >= c(4,5,6) & 6 <= c(8,9,10)
										#TRUE TRUE TRUE

									#7 >= c(4,5,6) & 7 <= c(8,9,10)
										#TRUE TRUE TRUE

									#8 >= c(4,5,6) & 8 <= c(8,9,10)
										#TRUE TRUE TRUE

									#9 >= c(4,5,6) & 9 <= c(8,9,10)
										#FALSE  TRUE  TRUE

									#10 >= c(4,5,6) & 10 <= c(8,9,10)
										#FALSE FALSE  TRUE

									#11 >= c(4,5,6) & 11 <= c(8,9,10)
										#FALSE FALSE FALSE

								#You can see as 3 and 11 are false for all the ranges. 4 and 10 are true only for the for the first and last range respectively. 5 and 9 are TRUE for first-second and second-last ranges respectively. Finally, 6, 7, and 8 are TRUE for all ranges, because all these number are between 4 and 10. For example, all TRUE disappear when we move from the limit of the first range (8), so 9 is not included between 4 and 8 because is higher than 4 but not smaller than 8.

							#check that we don´t have any TRUE in the three checks, and save
							test_ranges_overlap_2_v1 = append(test_ranges_overlap_2_v1, !TRUE %in% start_within_ranges_to_compare & !TRUE %in% end_within_ranges_to_compare & !TRUE %in% start_end_both_sides_ranges_to_compare)
						}

						#save the check result if all the ranges were tested
						if(length(test_ranges_overlap_2_v1) == length(overlap_sequences)){
							test_ranges_overlap_2 = all(test_ranges_overlap_2_v1)
						} else { #if not save NA
							test_ranges_overlap_2 = NA
						}	
					} else { #if not and hence we have only one range

						#save the check as NA
						test_ranges_overlap_2 = NA
					}

					#sum the width of all non-overlapped ranges
					coding_sequence_length = sum(overlap_sequences@ranges@width)


					##check that the number of coding bases are the same being calculated with Genomic ranges
					#create the reference range using the window start and end
					window_range = IRanges(start=selected_window$lower_end, end=selected_window$upper_end) #create a IRange object with the start and end of the [k] window

					#convert to a GRanges object. We need a GRanges object to select those exons in the the subject GRange file that included in the reference range. 
					window_range_gr = GRanges(seqnames=unique(chromosome_name), ranges=window_range, strand = '*') #sequnames are the names of the chromosome names. This is very important, because if we select as seq names de exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping using disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified in the exon data (I have checked the no impact of strand in this exact example for gap length).

					#remove the rows without coding region
					subset_exon_data_chromosome_check_coding_length = subset_exon_data_chromosome[which( !is.na(subset_exon_data_chromosome$genomic_coding_start) | !is.na(subset_exon_data_chromosome$genomic_coding_end) ),]

					#select the gap ranges
					exon_ranges = IRanges(start=subset_exon_data_chromosome_check_coding_length$genomic_coding_start, end=subset_exon_data_chromosome_check_coding_length$genomic_coding_end) #create a IRange object with the start and end of [k] window. 

					#convert the IRange file into a genomic range file
					exon_ranges_gr = GRanges(seqnames=unique(chromosome_name), ranges=exon_ranges, strand = '*') #sequnames are the names of the chromosome names. This is very important, because if we select as seq names de exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping using disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified in the exon data (I have checked the no impact of strand in this exact example for gap length).
						#this explanation is taken from gene_coordinates_v8.r

					#calculate the not overlapping ranges
					not_overlapped_exon_ranges_gr = disjoin(x=exon_ranges_gr, with.revmap=TRUE) #‘disjoin’ returns an object of the same type as ‘x’ containing disjoint ranges for each distinct (seqname, strand) pairing. Split all the segments to have no overlap. If for example you have 1-10 and 5-7, the first range would be 1-4, then 5-7, and finally 8-10. You have the second range included within the first one. If ‘with.revmap=TRUE’, a metadata column that maps the output ranges to the input ranges is added to the returned object. This is basically a map, tells you what initial ranges are included in each new non-overlapped range.

					#check that the disjoined version of the gap ranges (without the revision map) is similar to the initial ranges. If true, that means that we have no overlapping in the gap ranges
					check_coding_density_6 = identical(exon_ranges_gr, disjoin(exon_ranges_gr)) #in this case we should have overlapping because several transcripts of the same and different genes can be overlapped.

					#find what of the tested ranges and are overlapped with the reference range
					hits_2 = findOverlaps(window_range_gr, not_overlapped_exon_ranges_gr)

					#calculate the exact overlap, i.e., the gap regions that included within the window. 
					overlaps_2 <- pintersect(window_range_gr[queryHits(hits_2)], not_overlapped_exon_ranges_gr[subjectHits(hits_2)])
						#queryHits extracts what of the reference ranges are overlapped with the tested ranges
						#subjectHits extracts what of the tested ranges are overlapped with the reference ranges. For example, the first tested range is not overlapped, so the first integer in the column of subjectHits in hits_2 is 2, not 1. 
						#You associate each reference range with the overlapped tested ranges. That is the input needed for pintersect.
						#the result is for each of the tested ranges overlapped with the reference ranges, the bases not overlapped are removed.
							#obtained form here "https://support.bioconductor.org/p/72656/"
					
					#calculate the sum of all not-overlapping coding regions inside the window				
					coding_sequence_length_2 = sum(overlaps_2@ranges@width)

					#check that new calculation of coding region length is similar to the previous one
					check_coding_density_7 = coding_sequence_length == coding_sequence_length_2
		

					## calculate the coding density for the [i] gene considering all the transcripts
					#calculate the exact length of the window. The length should be equal to selected length, but in cases when the window surpass the chromosome limit, the window is smaller. Therefore, we need to calculate the window length using the window coordinates. 
					exact_window_length = (selected_window$upper_end - selected_window$lower_end)+1

					#remove the gap length from the length of the window
					exact_window_length_final = exact_window_length - selected_gap

					#coding density, length of coding sequences divided by the total length of the window
					coding_density = coding_sequence_length/exact_window_length_final
				} else { #if we do not have exons inside the window, set all NA, except coding density that will be zero
					check_coding_density_1 = NA
					check_coding_density_2 = NA
					check_coding_density_4 = NA
					check_coding_density_5 = NA
					check_coding_density_6 = NA
					check_coding_density_7 = NA
					test_all_ranges_included = NA
					test_iranges = NA
					test_start_end = NA
					test_ranges_overlap = NA
					test_ranges_overlap_2 = NA
					coding_density = 0
				}
			} else {#if you have NA for the window limits

				#Set NA for everything, except selected window. We cannot set coding density as zero!! zero is a value, NA is no calculation.
				check_coding_density_1 = NA
				check_coding_density_2 = NA
				check_coding_density_4 = NA
				check_coding_density_5 = NA
				check_coding_density_6 = NA
				check_coding_density_7 = NA
				test_all_ranges_included = NA
				test_iranges = NA
				test_start_end = NA
				test_ranges_overlap = NA
				test_ranges_overlap_2 = NA
				coding_density = NA
			}

			#save the results in a data.frame
			coding_density_windows = rbind.data.frame(coding_density_windows, cbind.data.frame(selected_window, check_coding_density_1, check_coding_density_2, check_coding_density_4, check_coding_density_5, check_coding_density_6, check_coding_density_7, test_all_ranges_included, test_iranges, test_start_end, test_ranges_overlap, test_ranges_overlap_2, coding_density))			
		}
		#remove first row with NAs
		coding_density_windows = coding_density_windows[-1,] #I have checked for SCYL3 (ENSG00000000457) the first two windows. The first window (50kb) includes two genes (SCYL3 and C1orf112), whilst the second one also includes KIFAP3. I have check the map around SCYL3 (http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000000457;r=1:169817880-169864300), and i have seen two genes closer to SCYL3 than C1orf112, which are METTL18 and SELL. However, I have checked that the end of these genes (169680839 and 169764107) occur before the start of the two first windows (169816090.5 and 169791090.5). In contrast, C1orf112 starts very far away from SCYL3, but its end (169823221) reached the start of these windows. In the next window, which is wider and hence starts before (169741090.5), METTL18 is included. The next windows (500kb) already include SELE, SELP and SEL. Finally, the widest window (100 kb) also includes BLZF1, CCDC181, SLC19A2, F5, METTL11B. In that way, we have included all the coding genes you can see in the window. If you scroll to both sides, you can see as no more coding genes are present between teh boundaries of the last window (169341090.5 - 170341089.5). The rest are pseudogenes, RNa genes and processed transcripts.

		#save the number of genes per each window
		coding_density_50kb = rep(coding_density_windows[which(coding_density_windows$selected_length == 50000),]$coding_density, nrow(exon_data))		
		coding_density_100kb = rep(coding_density_windows[which(coding_density_windows$selected_length == 100000),]$coding_density, nrow(exon_data))		
		coding_density_200kb = rep(coding_density_windows[which(coding_density_windows$selected_length == 200000),]$coding_density, nrow(exon_data))		
		coding_density_500kb = rep(coding_density_windows[which(coding_density_windows$selected_length == 500000),]$coding_density, nrow(exon_data))		
		coding_density_1000kb = rep(coding_density_windows[which(coding_density_windows$selected_length == 1000000),]$coding_density, nrow(exon_data))		


		##summarize all the coding_checks
		#create a vector with the check to be summarized
		coding_checks_to_summarize = c("check_coding_density_1","check_coding_density_2","check_coding_density_4","check_coding_density_5","check_coding_density_6","check_coding_density_7","test_all_ranges_included","test_iranges","test_start_end","test_ranges_overlap","test_ranges_overlap_2")

		#for each check to be summarized
		for(j in 1:length(coding_checks_to_summarize)){

			#select the [j] check to summarize
			selected_coding_check = coding_checks_to_summarize[j]

			#extract the corresponding column
			selected_column_coding_check = coding_density_windows[,which(colnames(coding_density_windows) == selected_coding_check)]	
			#if no FALSE is present in the check
			if(!FALSE %in% selected_column_coding_check){

				#if there is no FALSE in is.na() and hence all cases are NA
				if(!FALSE %in% is.na(selected_column_coding_check)){

					#set the check as TRUE
					summary_coding_check = rep(NA, nrow(exon_data))
				} else { #if not, and hence not all cases are NA
					
					#if you have NAs and TRUEs
					if(TRUE %in% is.na(selected_column_coding_check) & TRUE %in% selected_column_coding_check){
						
						#set the check as TRUE/NA
						summary_coding_check = rep("TRUE/NA", nrow(exon_data))
					} else { #if not, we dot not have any NA

						#set the check as TRUE
						summary_coding_check = rep(TRUE, nrow(exon_data))
					}
				}
			} else { #if we have FALSE

				summary_coding_check = rep(FALSE, nrow(exon_data))
			}

			#assign the value of the summary check to the corresponding variable.
			assign(selected_coding_check, summary_coding_check) #see second response in "https://stackoverflow.com/questions/28909191/treat-string-as-object-name-in-a-loop-in-r" for further details about assign		
		}


		#reorder data of the [i] gene from the gene list to get gene name (hgnc_symbol) and check that chromosome name and positions are correct 
		#reorder rows following exon_data again. We reorder rows of exon data following start_end_coding_regions
		gene_characteristics_ordered = gene_characteristics[match(paste(exon_data$ensembl_transcript_id, exon_data$ensembl_exon_id), paste(gene_characteristics$ensembl_transcript_id, gene_characteristics$ensembl_exon_id)),] #we match using both exon id and transcript id because the same exon id can be in different transcripts, thus we need two conditions for matching. The paste function create strings using the transcript id and the exon id, you can create these strings using start_end_coding_regions and exon_data, then use match to match the rows of exon data to those of start_end_coding_regions. See 'https://stackoverflow.com/questions/47404477/match-with-multiple-criteria-without-loop-in-r' for further details.

		#check that the first and the last exons (and the rest of them) are the same in gene list and exon data. I cannot compare start position because the position in gene list (all_genes_grch37_exon_filter_final) is calculate considering all exons, coding and non coding. we are only using coding exons. Also check that the chromosome name is ok. Also check that all rows have data for cds
		test_position_gene_exons = gene_characteristics_ordered$ensembl_exon_id == exon_data$ensembl_exon_id
		test_chr_name = unique(gene_characteristics_ordered$chromosome_name) == selected_chr
		test_na_cds = all(!is.na(exon_data$cds_length))

		#bind results in a row of dataframe
		binding_results = cbind.data.frame(chromosome_name, hgnc_symbol, gene_id, test_rows_initial_datasets_1, test_rows_initial_datasets_2, test_gene_id, gene_biotype, transcript_id, transcript_biotype, exon_id, test_gene_transcripts_biotype, gene_start, gene_end, gene_length, middle_point, test_gene_length_center, windows_check_1, windows_check_2, windows_check_3, windows_check_4, windows_check_5, lower_end_window_50kb, upper_end_window_50kb, lower_end_window_100kb, upper_end_window_100kb, lower_end_window_200kb, upper_end_window_200kb, lower_end_window_500kb, upper_end_window_500kb, lower_end_window_1000kb, upper_end_window_1000kb, n_genes_50kb, n_genes_100kb, n_genes_200kb, n_genes_500kb, n_genes_1000kb, check_n_genes_0, check_n_genes_1, check_n_genes_2, check_n_genes_3, check_n_genes_4, check_coding_density_1, check_coding_density_2, check_coding_density_4, check_coding_density_5, check_coding_density_6, check_coding_density_7, test_all_ranges_included, test_iranges, test_start_end, test_ranges_overlap, test_ranges_overlap_2, coding_density_50kb, coding_density_100kb, coding_density_200kb, coding_density_500kb, coding_density_1000kb, test_position_gene_exons, test_chr_name, test_na_cds)

		#check that we have the same number of rows than in exon data
		if(!nrow(binding_results) == nrow(exon_data)){

			#print the error
			print('ERROR! ROW NUMBER NOT MATCHES')

			#break the loop
			break			
		}

		#save the results
		final_positions = rbind.data.frame(final_positions, binding_results)
	}
	
	#remove NA row. We select those rows for which the sum of NAs is equal to the number of columns.
	final_positions = final_positions[-which( rowSums(is.na(final_positions)) == ncol(final_positions) ),]

	#return the final dataset
	return(final_positions)
}

##Parallelize the process
require(foreach)
require(doParallel) #for parallel

#write the list of chromosomes
list_chromosomes = c(1:22)

# set up cluster
clust <- makeCluster(7) 
registerDoParallel(clust)

#run the function
final_positions = foreach(i = list_chromosomes, .packages=c("GenomicRanges"), .combine="rbind") %dopar% { #we load genomic ranges for the overlapping ranges of coding sequences
    genomic_coords_coding_density(selected_chr = i)
} 

#stop the cluster 
stopCluster(clust)



################################
########## IMPORTANT ###########
################################

#the windows lengths were calculated using the middle points without removing the decimals, this having them also decimals. I make all the corresponding checks, and the I removed the decimals of upper and lower limits of the windows. These are the ends saved (without decimals) and used for calculating gene number and coding density. 



########################################
########## CHECK THE RESULTS ###########
########################################

#check head and tails
head(final_positions)
tail(final_positions)

#look for NAs or FALSE in the checks
summary(final_positions) #NO FALSES.Exceptions:
	#NAs in check_n_genes_1. 
	#FALSE check_coding_density_6.

#check check_n_genes_1.
summary(final_positions$check_n_genes_1) 
unique(final_positions$check_n_genes_1) #We have TRUE, TRUE/NA and NA. For this check we avoid doing it when some genes had not hgnc symbol, so you cannot calculate the number of genes using this variable
final_positions[which(is.na(final_positions$check_n_genes_1)),]$hgnc_symbol #there are cases with hgnc symbl. The gene at the center can have hgnc symbol, but any of the genes included in the window not, and then the test fail. In addition, there are NA for those windows with NA for coordinates (see below).

#check_coding_density_6
summary(final_positions$check_coding_density_6) #here we checked the overlap between coding exons inside a window. TRUE means that no overlap exist. We can have FALSE, because different transcript of the same gene can share an exon. In addition, different genes can overlap in their coding regions. Indeed all cases are NA (from NA windows, i.e., no coordinate) or FALSE (overlap within the gene).


##DIFFERENCES RESPECT TO V9

##NAs caused by removed windows because they surpass the chromosome limit. these NAs are new in v10
#IMPORTANT: You have thousands of NA ins the coordinates of big windows. This does no mean that there are thousand of genes without that window, because final_positions has a rows for EACH EXON! Therefore, we have removed thousand of exons when removing windows that surpass the chromosome limit, not genes. The real number of genes implicated can be obtained from the dataset after the duplicates have been removed.  

#In check_n_genes_1 and other checks we can have NAs because of the removed windows
#remove duplicates to do some checks. We need to have only one row per gene, all the rows of a gene should have the same values
final_positions_no_dupli = final_positions[-which(duplicated(final_positions$gene_id)),]
#check
length(unique(final_positions$gene_id)) == nrow(final_positions_no_dupli)

#check that the genes with NA for the coordinate of one end of a window, are the same than those with NA for the coordinate of the other extreme of that window 
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$upper_end_window_50kb)),]$gene_id == final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_50kb)),]$gene_id)
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$upper_end_window_100kb)),]$gene_id == final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_100kb)),]$gene_id)
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$upper_end_window_200kb)),]$gene_id == final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_200kb)),]$gene_id)
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$upper_end_window_500kb)),]$gene_id == final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_500kb)),]$gene_id)
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$upper_end_window_1000kb)),]$gene_id == final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_1000kb)),]$gene_id)

#extract the genes with at least one removed windows
genes_with_removed_windows = final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_50kb) | is.na(final_positions_no_dupli$lower_end_window_100kb) | is.na(final_positions_no_dupli$lower_end_window_200kb) | is.na(final_positions_no_dupli$lower_end_window_500kb) | is.na(final_positions_no_dupli$lower_end_window_1000kb)),] #which rows have NA for the coordinate of any of the windows? We use only the lower end. If one end has NA, the other will have NA (this was revised in the previous check).
#check that these genes are all included in those without coordinate for the 1000kb window
summary(final_positions_no_dupli[which(is.na(final_positions_no_dupli$lower_end_window_1000kb)),]$gene_id == genes_with_removed_windows$gene_id) #If 6 50kb windows are removed, then for the 100kb windows, you will have these 6 windows and new ones removed. If the 50 kb window surpass the chromosome limit, obviously the 100kb windows of these genes will also surpass the chromosome limits. Each size category includes all the previous genes. Finally, the 1000kb window size includes ALL the previous cases and some new ones. Because of this, the genes without 1000kb window are equal to the genes without at least one window.  
#see the number of cases. This is the total number of genes that suffered losses:
nrow(genes_with_removed_windows) #293 genes
#extract the total number of windows removed
length(which(is.na(final_positions_no_dupli$lower_end_window_50kb))) + length(which(is.na(final_positions_no_dupli$lower_end_window_100kb))) + length(which(is.na(final_positions_no_dupli$lower_end_window_200kb))) + length(which(is.na(final_positions_no_dupli$lower_end_window_500kb))) + length(which(is.na(final_positions_no_dupli$lower_end_window_1000kb))) #443 windows
	#There is 1 gene for which 50b windows have been removed. 6 in the case of 100kb windows, 22 for 200kb, 121 for 500kb and 293 for 1000kb. 
	 	#Therefore, 443 windows from 283 genes have been removed.  

#select those genes for which there is NA in the gene number of at least one window 
genes_with_na_gene_number = final_positions_no_dupli[which(is.na(final_positions_no_dupli$n_genes_50kb) | is.na(final_positions_no_dupli$n_genes_100kb) | is.na(final_positions_no_dupli$n_genes_200kb) | is.na(final_positions_no_dupli$n_genes_500kb) | is.na(final_positions_no_dupli$n_genes_1000kb)),]

#check if genes with at least one NA for gene coordinate are the same than those with at least one NA per gene number
summary(genes_with_removed_windows$gene_id == genes_with_na_gene_number$gene_id)

#check that all checks from gene number calculations that have NA are those with NA in windows coordinates
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_n_genes_0 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_n_genes_1 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id) #this check can have NA because of empty hgnc symbol. Therefore we can have different genes.
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_n_genes_2 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_n_genes_3 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_n_genes_4 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)

#calculate the number of genes for which there is NA in the coding density of at least one window 
genes_with_na_coding_density = final_positions_no_dupli[which(is.na(final_positions_no_dupli$coding_density_50kb) | is.na(final_positions_no_dupli$coding_density_100kb) | is.na(final_positions_no_dupli$coding_density_200kb) | is.na(final_positions_no_dupli$coding_density_500kb) | is.na(final_positions_no_dupli$coding_density_1000kb)),]

#check if genes with at least one NA for gene coordinate and at least one NA per gene number
summary(genes_with_removed_windows$gene_id == genes_with_na_coding_density$gene_id)

#check that all the cases with NA for coding density are those genes with NA for at least one window. WARNING: This can have have false without problem! See below.
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_1 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_2 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_4 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_5 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_6 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$check_coding_density_7 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$test_all_ranges_included %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$test_iranges %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$test_start_end %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$test_ranges_overlap %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
summary(final_positions_no_dupli[which(final_positions_no_dupli$test_ranges_overlap_2 %in% c(NA, "TRUE/NA")),]$gene_id == genes_with_removed_windows$gene_id)
	#coding_check_6 can give FALSE, because FALSE in this check is common. We are checking different coding regions within a gene are overlapped. This is very common, because you can have several transcripts in a gene. In general, if we find only one FALSE, everything is set as FALSE. This is like a red flag in general in the script, although in this cases is not very relevant.
	#In the case of check_coding_density_1, check_coding_density_2, check_coding_density_4, check_coding_density_5, check_coding_density_6, check_coding_density_7, test_all_ranges_included, test_iranges, test_start_end, test_ranges_overlap, test_ranges_overlap_2 you can obtain FALSE because you can have NAs if no coding region is included in the region.

#NOTE: The windows_check_5 does not have NAs. This was the first check after removing the windows that surpass the chromosome limit, but I simply check that no FALSE was present and hence no window surpass the limit. Therefore, the TRUEs includes windows with coordinate that do not surpass, and windows with NA in the coordinate. For the rest of check that were calculated for each window, we set NA if the windows have no coordinate. 


## DIFFERENCES RESPECT TO V8.
#There differences in the median start/end of the genes analyzed, the median gene length, median number of genes inside each windows size. This can be explained because we removed the whole chromosome X respect to version 8th. Therefore, the set of genes analyzed are different. The removal of the whole chromosome, and its genes affects the median number of genes you can find inside each window size. Similarly, the median length of the genes will different, because you have reduced your set of genes!! 
	
#There are also differences in the positions of the windows. We have removes windows that surpass the limits of the chromosome. Windows reaching the start of the chromosome (i.e., negative coordinates) and end of the chromosome. Indeed, you can see in the summary of the v8 that the min value for the lower end of all windows is negative!! That does not make any sense! This problem has been solved in the new version. The windows_check_5 check that no lower limit is < 1 and no upper limit is bigger than the corresponding chromosome length. We have all TRUE for this check.
#the minimum coordinate of all widows is equal or higher than 1
min(na.omit(final_positions$lower_end_window_50kb)) >= 1
min(na.omit(final_positions$lower_end_window_100kb)) >= 1
min(na.omit(final_positions$lower_end_window_200kb)) >= 1
min(na.omit(final_positions$lower_end_window_500kb)) >= 1
min(na.omit(final_positions$lower_end_window_1000kb)) >= 1 #we removed the NA for the windows that have been removed because of this problem of surpassing the chromosome limit.
#the maximum coordinate of all widows is equal or lower than the maximum chromosome length 
max(na.omit(final_positions$upper_end_window_50kb)) <= max(chrom_length_ucsc_hg19$length_bp)
max(na.omit(final_positions$upper_end_window_100kb)) <= max(chrom_length_ucsc_hg19$length_bp)
max(na.omit(final_positions$upper_end_window_200kb)) <= max(chrom_length_ucsc_hg19$length_bp)
max(na.omit(final_positions$upper_end_window_500kb)) <= max(chrom_length_ucsc_hg19$length_bp)
max(na.omit(final_positions$upper_end_window_1000kb)) <= max(chrom_length_ucsc_hg19$length_bp) #we removed the NA for the windows that have been removed because of this problem of surpassing the chromosome limit.

#in the case of coding density, we also made a change in the size of the window that affects. We removed the size corresponding with gaps inside the window. Therefore, the length of coding regions is divided by the length of the windows - gap length. 
 
 
#check that all exons included in final_positions are in full_exon_data_filtered_final. We have to calculate the intersection between transcript id and exon id, because exons from different transcripts can have the same exon id
summary(paste(final_positions$transcript_id, final_positions$exon_id) %in% paste(full_exon_data_filtered_final$ensembl_transcript_id, full_exon_data_filtered_final$ensembl_exon_id))
#the same for all_genes_grch37_exon_filter_final
summary(paste(final_positions$transcript_id, final_positions$exon_id) %in% paste(all_genes_grch37_exon_filter_final$ensembl_transcript_id, all_genes_grch37_exon_filter_final$ensembl_exon_id)) 
	#we did not do it in the inverse way (check if all exons of the input data are in the output) because the function removes those non-coding exons. 

#check that all genes per chromosome were included
#calculate the number of genes per chromosome
number_genes_per_chromosome = NA
for(i in 1:length(list_chromosomes)){ #for each chromosome

	#select the [i] chromosome
	selected_chromosome = list_chromosomes[i]

	#select those rows of the [i] chromosome in all_genes_grch37_exon_filter_final and extract the unique cases of gene id (the unique gene ids)
	chromosome_subset = unique(all_genes_grch37_exon_filter_final[which(all_genes_grch37_exon_filter_final$chromosome_name == selected_chromosome),]$ensembl_gene_id)

	#save the length of unique genes
	number_genes_per_chromosome = append(number_genes_per_chromosome, length(chromosome_subset))
}
#remove the first element of the vector, which is a NA
number_genes_per_chromosome = number_genes_per_chromosome[-1]

#calculate the total number of genes
total_number_genes = sum(number_genes_per_chromosome) 
#check this number is equal to the length of unique genes
length(unique(all_genes_grch37_exon_filter_final$ensembl_gene_id)) == total_number_genes

#check that you number of unique genes is equal to the real number of genes
length(unique(final_positions$gene_id)) == total_number_genes

#check that you number of unique exons is equal to the real number of exons
length(unique(final_positions$ensembl_exon_id)) == length(unique(all_genes_grch37_exon_filter_final[which(!is.na(all_genes_grch37_exon_filter_final$genomic_coding_start) & !is.na(all_genes_grch37_exon_filter_final$genomic_coding_end)),]$ensembl_exon_id))#we have to remove all exons without cds data of all_genes_grch37_exon_filter_final, because these were removed from final_positions 


##IMPORTANT: In gene coordinates, we do still have the three cases with window coordinated but completely overlapped with a gap indicated as coding density as zero, not NA. gene_coordinates has to be used only for window coordinates.
#save final_positions
write.table(final_positions, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt', sep='\t', row.names=FALSE)

#load the data to check that the writing worked ok
final_positions_2 = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt', sep='\t', header=TRUE)

#make the test
test=NULL
for(i in 1:ncol(final_positions)){ #for each column of the final positions

	#select that column in final_positions and final_positions_2
	selected_column_final_positions = final_positions[,i]
	selected_column_test = final_positions_2[,i]

	#if any of the selected columns has no NA
	if(all(!is.na(selected_column_final_positions)) & all(!is.na(selected_column_test))){

		#if the column is numeric we have to round to get exactly the same numbers
		if(is.numeric(selected_column_final_positions) & is.numeric(selected_column_test)){

			#round and check that all numbers are equal
			test = append(test, all(round(selected_column_final_positions) == round(selected_column_test)))
		} else { #if not, and hence we don´t have numeric data

			#make the test without doing anything else
			test = append(test, all(selected_column_final_positions == selected_column_test))
		}
	} else { #if not and hence we have NAs, we apply na.omit

		#if the column is numeric we have to round to get exactly the same numbers
		if(is.numeric(selected_column_final_positions) & is.numeric(selected_column_test)){

			#round and check that all numbers are equal
			test = append(test, all(round(na.omit(selected_column_final_positions)) == round(na.omit(selected_column_test))))
		} else { #if not, and hence we don´t have numeric data

			#make the test without rounding
			test = append(test, all(na.omit(selected_column_final_positions) == na.omit(selected_column_test)))
		}
	}
}
all(test)#all TRUE



####################################################
######## PREPARE THE FINAL FILES FOR DAVID #########
####################################################

##remove the duplicates
#remove those windows that are duplicated in gene_id. We are using the same code used in the other confounding factors to remove duplicated. It is tested. 
final_positions_no_duplicates = final_positions[-which(duplicated(final_positions$gene_id)),]
#check we have the correct number of rows
nrow(final_positions_no_duplicates) == length(unique(final_positions$gene_id))
#check we have the correct gene_ids
summary(final_positions_no_duplicates$gene_id == unique(final_positions$gene_id))


## gene length data
#select the gene_ids and the gene length
final_gene_length = final_positions_no_duplicates[,c("gene_id", "gene_length")]

#check we have the correct number of rows
nrow(final_gene_length) == length(unique(final_positions$gene_id))
#check we have the correct gene_ids
summary(final_gene_length$gene_id == unique(final_positions$gene_id))

#take a look
str(final_gene_length)
head(final_gene_length)
summary(final_gene_length)
	#the median and mean are 25811.00 and 65677.91.
	#No NA because in this case the windows surpassing the chromosome limit are not a problem. We are using the gene length according to ensemble, so the genes should not surpass the chromosome limit. In addition, if there are gaps inside the gene it is not a problem, we know where the gene starts and ends, that is the point here. 
	#According to wikipedia (https://en.wikipedia.org/wiki/Human_genome), the median size of a protein-coding gene is 26,288 bp (mean = 66,577 bp). In our case, we have a median value of 25,811.00 and a mean of 65,677.91. Values are very close. 

#check that maximum gene end in each chromosome is always lower than the chromosome length
check_gene_length = NULL
#for each chromosome
for(i in 1:length(unique(final_positions_no_duplicates$chromosome_name))){

	#select the [i] chromosome
	selected_chromosome = unique(final_positions_no_duplicates$chromosome_name)[i]
	
	#subset final_positions_no_duplicates and chrom_length_ucsc_hg19
	subset_final_gene_length = final_positions_no_duplicates[which(final_positions_no_duplicates$chromosome_name == selected_chromosome),]
	subset_chrom_length_ucsc_hg19 = chrom_length_ucsc_hg19[which(chrom_length_ucsc_hg19$chromosome == paste("chr", selected_chromosome, sep="")),]

	#save the check
	check_gene_length = append(check_gene_length, max(subset_final_gene_length$gene_end) < subset_chrom_length_ucsc_hg19$length_bp)
}
#take a look
summary(check_gene_length) #ALL TRUE, no gene end is bigger than the chromosome limit. 

#check that all gene start are higher than zero
min(final_positions_no_duplicates$gene_start) > 0 #Starts are higher than zero and ends are lower than the length of the corresponding chromosome.


## gene number data
# extract gene id and the number of gene inside each window
final_gene_number = final_positions_no_duplicates[,c("gene_id", "n_genes_50kb", "n_genes_100kb", "n_genes_200kb", "n_genes_500kb", "n_genes_1000kb")]

#check we have the correct number of rows
nrow(final_gene_number) == length(unique(final_positions$gene_id))
#check we have the correct gene_ids
summary(final_gene_number$gene_id == unique(final_positions$gene_id))

#take a look
str(final_gene_number)
head(final_gene_number)
summary(final_gene_number) 
	#NAs: There is 1 gene for which 50b windows have been removed. 6 in the case of 100kb windows, 22 for 200kb, 121 for 500kb and 293 for 1000kb. These are the cases removed because of the window surpass the chromosome limit.
	#The mean and median number of genes is increasing from 50kb windows to 1000kb windows.
	#The minimum is 1 and the maximum is 70, which sounds good.

#check that the rows with NA in gene number are the same than those with NA for any window coordinate in each window size
summary(which(is.na(final_gene_number$n_genes_50kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_50kb) | is.na(final_positions_no_duplicates$upper_end_window_50kb)))
summary(which(is.na(final_gene_number$n_genes_100kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_100kb) | is.na(final_positions_no_duplicates$upper_end_window_100kb)))
summary(which(is.na(final_gene_number$n_genes_200kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_200kb) | is.na(final_positions_no_duplicates$upper_end_window_200kb)))
summary(which(is.na(final_gene_number$n_genes_500kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_500kb) | is.na(final_positions_no_duplicates$upper_end_window_500kb)))
summary(which(is.na(final_gene_number$n_genes_1000kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_1000kb) | is.na(final_positions_no_duplicates$upper_end_window_1000kb)))

#load tbfs data to extract those rows with window coordinates but not tbfs data because the window is completely overlapped with a gap
tbfs_data = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/tbfs_density_raw_v3.txt', sep='\t', header=TRUE)

#from final_gene_number, select those rows with a gene_id that is included within the gene_id of the rows in final_tbfs that have window coordinates but the final window has a length of zero.
final_gene_number[which(final_gene_number$gene_id %in% tbfs_data[which(!is.na(tbfs_data$selected_lower_window) & !is.na(tbfs_data$selected_upper_window) & tbfs_data$window_length_no_gaps == 0),]$gene_id),]
	#The three cases with window coordinates but that the window is completely overlapped with a gap have a gene number of 1, which makes sense. We are counting the gene in the middle. 

	#I have checked that these 50kb windows are overlapped with only one gene
		#http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000148357;r=9:133046882-133309510
		#http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000168280;r=2:149733046-149783045
		#http://grch37.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000185974;r=13:114356089-114406088


## coding density data
# extract gene id and the coding density inside each window
final_coding_density = final_positions_no_duplicates[,c("gene_id", "coding_density_50kb", "coding_density_100kb", "coding_density_200kb", "coding_density_500kb", "coding_density_1000kb")]

#check we have the correct number of rows
nrow(final_coding_density) == length(unique(final_positions$gene_id))
#check we have the correct gene_ids
summary(final_coding_density$gene_id == unique(final_positions$gene_id))

#take a look
str(final_coding_density)
head(final_coding_density)
summary(final_coding_density)
	#NAs: There is 1 gene for which 50b windows have been removed. 6 in the case of 100kb windows, 22 for 200kb, 121 for 500kb and 293 for 1000kb. These are the cases removed because of the window surpass the chromosome limit.
	#The mean and median coding density is decreasing from 50kb windows to 1000kb windows. Bigger windows have more bases, therefore the coding bases have to be a lot of more to reach the similar density than smaller windows.
	#The minimum is zero and the maximum is 0.38, which sounds good.

#check that the rows with NA in gene number are the same than those with NA for any window coordinate in each window size
summary(which(is.na(final_coding_density$coding_density_50kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_50kb) | is.na(final_positions_no_duplicates$upper_end_window_50kb)))
summary(which(is.na(final_coding_density$coding_density_100kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_100kb) | is.na(final_positions_no_duplicates$upper_end_window_100kb)))
summary(which(is.na(final_coding_density$coding_density_200kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_200kb) | is.na(final_positions_no_duplicates$upper_end_window_200kb)))
summary(which(is.na(final_coding_density$coding_density_500kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_500kb) | is.na(final_positions_no_duplicates$upper_end_window_500kb)))
summary(which(is.na(final_coding_density$coding_density_1000kb)) == which(is.na(final_positions_no_duplicates$lower_end_window_1000kb) | is.na(final_positions_no_duplicates$upper_end_window_1000kb)))

#load tbfs data to extract those rows with window coordinates but not tbfs data because the window is completely overlapped with a gap
tbfs_data = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/tbfs_density_raw_v3.txt', sep='\t', header=TRUE)

#from coding density, select those rows with a gene_id that is included within the gene_id of the rows in final_tbfs that have window coordinates but the final window has a length of zero.
final_coding_density[which(final_coding_density$gene_id %in% tbfs_data[which(!is.na(tbfs_data$selected_lower_window) & !is.na(tbfs_data$selected_upper_window) & tbfs_data$window_length_no_gaps == 0),]$gene_id),]
	#coding_density_50kb is zero, but it should be NA. If the window is completely overlapped with a gap, we have no data about coding density inside the window, so we cannot calculate coding density, is NA, NOT zero. 

#set it as NA
final_coding_density[which(final_coding_density$gene_id %in% tbfs_data[which(!is.na(tbfs_data$selected_lower_window) & !is.na(tbfs_data$selected_upper_window) & tbfs_data$window_length_no_gaps == 0),]$gene_id),]$coding_density_50kb <- NA
	#note that the in the original gene_coordinate file, these three rows have zero in coding density but I set them as zero in the new versions.
#check
final_coding_density[which(final_coding_density$gene_id %in% tbfs_data[which(!is.na(tbfs_data$selected_lower_window) & !is.na(tbfs_data$selected_upper_window) & tbfs_data$window_length_no_gaps == 0),]$gene_id),]

#merge the data not duplicated and the gaps data to check that we are selecting the problematic genes
merged_data = merge(final_positions_no_duplicates, gap_results_extracted_final, by="gene_id")

#take those cases in which the window length of the 50kb window less the gap length is only zero
merged_data[which((merged_data$upper_end_window_50kb - merged_data$lower_end_window_50kb + 1) - merged_data$length_gaps_50kb == 0),]
	#These are the same cases and in gene_coordinates, these are coding_density_50kb equal to zero. The NA were only added to final_coding_density

#summary
summary(final_coding_density) #Now have three more NAs in coding_density_50kb. 


##save the files

#IMPORTANT: In gene coordinates, we do still have the three cases with window coordinated but completely overlapped with a gap indicated as coding density as zero, not NA. gene_coordinates has to be used only for window coordinates.

#save gene length
write.table(final_gene_length, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_length_v1.txt', sep='\t', row.names=FALSE)

#save gene length
write.table(final_gene_number, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_v1.txt', sep='\t', row.names=FALSE)

#save gene length
write.table(final_coding_density, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/coding_density_v1.txt', sep='\t', row.names=FALSE)



#####################################
######## SAVE THE WORKSPACE #########
#####################################
save.image('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/rdata/gene_coordinates_v10.RData')
library(biomaRt)
library(GenomicRanges)
library(data.table)
library(dplyr)
require(stringr)
require(D3GB)
