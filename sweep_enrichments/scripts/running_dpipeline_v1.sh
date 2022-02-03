#Note that this script was wrote and run in 2022 using the version of April 28th 2022 of David's pipeline.
	#The pipeline of David has been downloaded from 
		#https://github.com/DavidPierreEnard/Gene_Set_Enrichment_Pipeline
		#last commit done on Apr 28, 2021
			#7b755c0c23dd4d7c3f54c4b53e74366e4041ac8f 

#Manual from David
	#/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder/exdef_pipeline_manual.pdf


#go to the folder with the pipeline
cd /home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/exdef_folder
	#the originals are located in "/home/dftortosa/singularity/dating_climate_adaptation/sweep_enrichments/david_pipeline/original_files"

#open the rights to use all the perl scripts
chmod +x ./*.pl
	#you have to add "./" first so the result of the wild card will be pasted to "./". In that way, you get all the perl script names along with "./", which is needed to use chmod and allow for running.

