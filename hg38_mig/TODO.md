- calculate map file using genetic distance from decode hg38 original
	- check the original map is hg38

- Ask jesus about what accesibility masks is he using: L, Z and Q? also H?
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


