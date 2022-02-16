#Note that this script was wrote and run in 2022 using the version of April 28th 2022 of David's pipeline.
	#The pipeline of David has been downloaded from 
		#https://github.com/DavidPierreEnard/Gene_Set_Enrichment_Pipeline
		#last commit done on Apr 28, 2021
			#7b755c0c23dd4d7c3f54c4b53e74366e4041ac8f 

#Manual from David
	#/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/exdef_pipeline_manual.pdf

#remove the folder with the pipeline
#rm -rf /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder
	#r: removes directories and their content recursively
	#f: forces the removal of all files or directories
	#https://phoenixnap.com/kb/remove-directory-linux

#copy the folder with the original pipeline to the folder where we will run the pipeline
#cp -avr /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/original_files/exdef_folder /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/
	#a: preserve the specific attributes (e.g., mode)
	#v: verbose output
	#r: copy directories recursively
	#https://www.cyberciti.biz/faq/copy-folder-linux-command-line/

#run the scripts for preparing the files needed for the pipeline
R CMD BATCH /opt/scripts/preparing_valid_gene_dist_files_v2.R
R CMD BATCH /opt/scripts/preparing_factors_tables_v2.R
R CMD BATCH /opt/scripts/statistics_ranks_v2.R

#go to the folder with the pipeline
cd david_pipeline/exdef_folder
	#this step is MANDATORY in order to run the pipeline

#open the rights to use all the perl scripts
chmod +x *.pl

#run the pipeline using the input parameters file created by me in the script file
./exdef_pipeline.pl /opt/scripts/input_par_david_pipeline_v2.txt

#comeback to the original folder (two levels above the current folder)
cd ../..
	#https://askubuntu.com/questions/25813/commandline-shortcut-for-current-directory-similar-to-for-home-directory

#run the scripts for processing the output of the pipeline
R CMD BATCH /opt/scripts/significance_whole_enrich_curve_v2.R