- ASK DAVID FIRST ABOUT hgdp + 1000kgp
	- A new preprint has been posted about the harmonization of 1000KP high coverage and the HGDP. 
		- A harmonized public resource of deeply sequenced diverse human genomes
		- https://www.biorxiv.org/content/10.1101/2023.01.23.525248v2.full
	- we could gain populations at higher latitudes in eurasia, which could be useful for the climate adaptation paper. 
	- although I do not know if we would need new simulations for these pops...
	- do you think it is worth it?

- calculate map file using genetic distance from decode hg38 original
	- check the original map is hg38

- Ask jesus 
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


