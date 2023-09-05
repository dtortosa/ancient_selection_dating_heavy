#!/usr/bin/env Rscript
# coding: utf-8
    #to run this script: chmod +x script.R; ./script.R
    #!/bin/sh does not work with my terminal en msi of David.
    #if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
    #you can save the output and the errors
        #./script.R > script.Rout #only output
        #./script.R 2> error.Rout #only error
        #./script.R > script.Rout 2> error.out #both in different files
        #./script.R > script.Rout 2>&1 #both in the same file
        #https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/



#########################################################
######## CHECK OVERLAPPING IN ACCESIBILITY MASKS ########
#########################################################

#I have downloaded masks about the accessbility of the genome to short reads sequencing (see 00_download_hg38_data.sh). These a BED files (0-based) with the intervals of the genome that are accessible according to the corresponding definition (pilot or strict). See 01_hap_map_calcs.py for further details about the masks.


#################
# load packages #
#################
require(GenomicRanges)



##############
# load masks #
##############

#pilot mask, less stringent
pilot_mask = read.table(
	"./data/masks/20160622.allChr.pilot_mask.bed",
	sep="\t",
	header=FALSE)
colnames(pilot_mask) = c("chrom", "chromStart", "chromEnd", "mask_type")
str(pilot_mask)

#strict mask, more stringent
strict_mask = read.table(
	"./data/masks/20160622.allChr.mask.bed.gz",
	sep="\t",
	header=FALSE)
colnames(strict_mask) = c("chrom", "chromStart", "chromEnd", "mask_type")
str(strict_mask)

#mask="pilot"
for(mask in c("pilot", "strict")){

	#
	print(paste("######## STARTING ", mask, " ########",sep=""))

	#select the mask file
	if(mask=="pilot"){
		selected_mask = pilot_mask
		print(paste("Do we have the correct mask? ", unique(selected_mask$mask_type) == "pilot", sep=""))
	}
	if(mask=="strict"){
		selected_mask = strict_mask
		print(paste("Do we have the correct mask? ", unique(selected_mask$mask_type) == "strict", sep=""))
	}

	#we should have 22 autosomals + 2 sexual chromosomes
	print(paste("Do we have 24 chromosomes? ", length(unique(selected_mask$chrom)) == 24, sep=""))

	#convert the segments with GC content in the [i] chromosome as a IRange object
	ranges_mask = IRanges(start=selected_mask$chromStart, end=selected_mask$chromEnd) 
		#create a IRange object with the start and end of each interval
	print("See the number of the intervals:")
	print(length(ranges_mask@width))
	print("See the median number of bases of the intervals:")
	print(summary(ranges_mask@width))
	print("Sum the total number of bases across all intervals:")
	print(sum(ranges_mask@width))
		#strict has more but smaller intervals leading to less total length

	#convert the IRange file into a genomic range file
	ranges_gr = GRanges(seqnames=selected_mask$chrom, ranges=ranges_mask, strand="*") 
		#sequnames are the names of the chromosome names. This is very important, because if we select as seq names de exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping with disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *.
	
	#calculate the not overlapping ranges with GC content
	ranges_gr_not_overlapped = disjoin(ranges_gr, with.revmap=TRUE) 
		#‘disjoin’ returns an object of the same type as ‘x’ containing disjoint ranges for each distinct (seqname, strand) pairing. Split all the segments to have no overlap. If for example you have 1-10 and 5-7, the first range would be 1-4, then 5-7, and finally 8-10. You have the second range included within the first one. If ‘with.revmap=TRUE’, a metadata column that maps the ouput ranges to the input ranges is added to the returned object. This is basically a map, tells you what initial ranges are included in each of the new non-overlapped range.
	
	#compare overlapped and not overlapped ranges
	check_overlapping = identical(
		ranges_gr@ranges, 
		ranges_gr_not_overlapped@ranges)
			#if there is no differences, then no overlap exists between the intervals

	#print
	print(paste(mask, " has NO intervals overlapped? ", check_overlapping, sep=""))
}
