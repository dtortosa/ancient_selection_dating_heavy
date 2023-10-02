- Ask jesus 
	- about VEP
		- finish writting email
		- show him output of installing VEP, the test is OK, but we have a few warnings like "The feature_type chromosome is being skipped". I have not found info in github about this. 
			- /home/dftortosa/singularity/dating_climate_adaptation/hg38_mig/output_installing_vep.txt 
			- when running VEP I get a warning, but it is not the same
				- "Smartmatch is experimental at /opt/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/File.pm line 472."
		- should I also use ancestral alleles that are lower case? i.e., they have low confidence?
			- check number of variants wiht low confidnce
    			- if only 2% is ok, half of the variants is not
	- about what accesibility masks is he using: L, Z and Q? also H?
		- These masks help to identify regions of the genome that are more or less accesible to next generation sequencing methods using short reads.
		- David say to avoid regions with a very low depth when it is not possible to do variant call.
			- Note that the average depth in past versions was 3X, but now the average is 10 times higher.
		- Info from README
			- masks
				- N - the base is an N in the reference genome GRCh37
				- L - depth of coverage is much lower than average
				- H - depth of coverage is much higher than average
				- Z - too many reads with zero mapping quality overlap this position
				- Q - the average mapping quality at the position is too low
				- P - the base passed all filters
				- 0 - an overlapping base was never observed in aligned reads 
			- Regions marked as L, H, Z or Q are less accessible to short reads. Although they can still be analyzed they are more prone to false positives.
				- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/
	- I guess the samples in the VCF files have the same order than in the pedigree right? So using that and the pedigree with population codes, we can subset the VCF file of each chromosome per population
		- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
		- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt


- calculate map file using genetic distance from decode hg38 original
	- check the original map is hg38
		- it says it is hg38 in the origina file, aau1043_DataS3.gz
	- check if the data is 1 or 0 based

	- DAVID: according to DAvid, you can check that with the assembly file, if the base the map is saying is in position X is indeed at position X, you are good.



- ask david about hgdp + 1000kgp
	- A new preprint has been posted about the harmonization of 1000KP high coverage and the HGDP. 
		- A harmonized public resource of deeply sequenced diverse human genomes
		- https://www.biorxiv.org/content/10.1101/2023.01.23.525248v2.full
	- we could gain populations at higher latitudes in eurasia, which could be useful for the climate adaptation paper. 
	- although I do not know if we would need new simulations for these pops...
	- do you think it is worth it?
		- For the future could be interesting, but FIRST we need to know if flex sweep works well with 1000GP and its number of indidividuals
		- Then, in the future, we can check whether Flexsweep works with smaller population sizes, as the HGDP has less samples in many populations.
		- Months later, David sent to you again this paper! He forgot that we talked about this, and he thought that it could be interesting for me! So in the future it would be worth it to mention it and check if flex-sweep works with this smaller sample size.