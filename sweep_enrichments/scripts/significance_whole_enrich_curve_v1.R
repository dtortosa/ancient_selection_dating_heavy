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