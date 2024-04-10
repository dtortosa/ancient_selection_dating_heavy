########################################################
########## SCRIPTS FOR CALCULATING GC CONTENT ##########
########################################################

#In this script, we calculate the density of GC content inside windows of different sizes centered around each coding gene.



###################################################
###### CHANGES RESPECT TO PREVIOUS VERSIONS #######
###################################################

#Respect to v1, 
	#I have changed the way to calculate GC density. Now we sum the GC content of all segments inside a window and divide by the sum of the numbers of data points for all the segments inside the window.

	#I have added a check to review if the segments of GC content inside each chromosome are overlapped.

#Respect to v2, 
	#Rerun the script revising everything that could be wrong about the changes made about the coordinate windows in gene_coordinates_v10.r. We REMOVED those windows that surpass the chromosome limit (start or and). In the previous version we trimmed them. If you trimmed them, you can have 100kb and a 50 kb with similar size (less difference than it should be). We have to consider this. 



##########################################
########## REMOVE PREVIOUS WORKSPACE #####
##########################################
remove(list=ls(all=TRUE))



#################################
###### REQUIRE PACKAGES #########
#################################

require(plyr) #for apply functions across lists and data.frames. This is better than apply, because split rows of a data.frame without converting into matrix or array. In that way you can use "$" to call columns. In addition, you can save the output as a data frame or a list.



#############################################################
####### DESCRIPTION OF "GC PERCENT IN 5-BASE WINDOWS" #######
#############################################################

#From: "https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=801848801_WrOkvQQaHYn0IPPS2avBfoahfyXU&c=chr1&g=gc5Base"

#The GC percent track shows the percentage of G (guanine) and C (cytosine) bases in 5-base windows. High GC content is typically associated with gene-rich areas.



#######################################################
###### SCHEMA FOR "GC PERCENT IN 5-BASE WINDOWS" ######
#######################################################

#From: "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=gc5Base&hgta_table=gc5Base&hgta_doSchema=describe+table+schema".


##The file is in the WIG format. Info from here: "http://genome.ucsc.edu/goldenPath/help/wiggle.html"
	#The wiggle (WIG) format is an older format for display of dense, continuous data such as GC percent, probability scores, and transcriptome data. Wiggle data elements must be equally sized.

	#For speed and efficiency, wiggle data is compressed and stored internally in 128 unique bins. This compression means that there is a minor loss of precision when data is exported from a wiggle track (i.e., with output format "data points" or "bed format" within the Table Browser). The bedGraph format should be used if it is important to retain exact data when exporting.
		#we get the sum of percentages of GC content for 5 base windows, so you can calculate the mean but you do not have the GC-percentage of each window.

	#Wiggle format is line-oriented. For wiggle custom tracks, the first line must be a track definition line (i.e., track type=wiggle_0), which designates the track as a wiggle track and adds a number of options for controlling the default display.

	#Wiggle format is composed of declaration lines and data lines, and require a separate wiggle track definition line. There are two options for formatting wiggle data: variableStep and fixedStep. These formats were developed to allow the file to be written as compactly as possible. 

		#varaibleStep format: This format is used for data with irregular intervals between new data points, and is the more commonly used wiggle format. After the wiggle track definition line, variableStep begins with a declaration line and is followed by two columns containing chromosome positions and data values. The declaration line starts with the word variableStep and is followed by a specification for a chromosome. The optional span parameter (default: span=1) allows data composed of contiguous runs of bases with the same data value to be specified more succinctly. The span begins at each chromosome position specified and indicates the number of bases that data value should cover.

		#fixedStep format: This format is used for data with regular intervals between new data values and is the more compact wiggle format. After the wiggle track definition line, fixedStep begins with a declaration line and is followed by a single column of data values. The declaration line starts with the word fixedStep and includes specifications for chromosome, start coordinate, and step size. The span specification has the same meaning as in variableStep format

		#Note that for both variableStep and fixedStep formats, the same span must be used throughout the dataset. If no span is specified, the default span of 1 is used. As the name suggests, fixedStep wiggles require the same size step throughout the dataset. If not specified, a step size of 1 is used.

#load the data, This is Wig file.
raw_gc_data <- read.table(file="/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/gc_content/gc5Base.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
	#header=FALSE: because the file has not header, the first row already includes data.
	#sep="\t": The file is separated by tabs. I have checked if it works reading it with ";" or "," and it doesn´t. 
	#stringsAsFactors=FALSE: Not sure if including the 9th variable as a factor could give problems. For each row, they have a path. So we set stringsAsFactors as FALSE. 
	#Note: In case would have some problems reading the file ("EOF within quoted string" or "number of items read is not a multiple of the number of columns") you can add 'quote=""' ("https://www.biostars.org/p/170631/"; first response).
str(raw_gc_data)
head(raw_gc_data, 10)
nrow(raw_gc_data)
	#Download:
		#I did not find the link for this exact table, so from the path indicated in the David´s Cell paper for conserved elements ("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/"), I look for the GC track (5 window). I found the file called "gc5Base.txt.gz". I decompressed, and loaded it. 
			#It has exactly the same values than the 10 first rows showed in the scheme.
			#It has also exactly the same number of rows of the file described in the scheme of GC Percent - GC Percent in 5-Base Windows. 
	
		#There are other files that I have discarded, CHECK:
			#In the page for downloads ("http://hgdownload.soe.ucsc.edu/downloads.html#human"), I selected Genome sequence files and select annotations (2bit, GTF, GC-content, etc) in hg19, leading to "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/". There, I selected the file called "hg19.gc5Base.wig.gz". This is also a Wig file.
				#It has the same columns, and the first 10 rows are similar to those showed in the scheme, BUT it has more rows. This file was modified the last year, while the track explained in the scheme was last updated in 2009 (like the file included in the golden path for hg19).
					#More info of this file: wiggle database table for the GC Percent track this is an older standard alternative to the current bigWig format of the track, sometimes useful for analysis.

			#In the same golden Path ("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/"), we have the data to obtain this track in a bigWig format. This is a new format for this type of files, more compressed. I have not been able to open it, because require some operations.  
				#Info:
					#http://genome.ucsc.edu/goldenPath/help/wiggle.html
					#http://genomewiki.ucsc.edu/index.php/Using_hgWiggle_without_a_database

			#In the page for downloads ("http://hgdownload.soe.ucsc.edu/downloads.html#human"), I selected GC percent data in hg19, leading to "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/gc5Base/". There, I selected the file called "hg19.gc5Base.txt.gz".
				#This file is much heavier than the previous ones (1.5GB). The webpage says that " it is the raw data used to encode the gc5Base track on hg19." In addition, says that is not 0-based but 1-based, which does not match what the scheme says. 
				#It seems this data was used to create the final tracks.
				

##Scheme:
	
	#You can see as the database is hg19 (i.e., GRch37).

	#Definitions of each column:
		#bin: Indexing field to speed chromosome range queries.

		#chrom: The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
		unique(raw_gc_data$V2) #we have several names different from the autosomal and sexual chromosomes 
			#For what I have read, some cases can be sequences with not great reliability about the exact position in the chromosome (not finished) or different haplotypes from the Major Histocompatibility Complex. For example, COX is an haplotipe of MHC in chromosome 6 (https://www.biostars.org/p/7629/)

				#Starting with the April 2003 human assembly, these tables also include data for sequence that is not in a finished state, but whose location in the chromosome is known, in addition to the unordered sequence. Because this sequence is not quite finished, it could not be included in the main "finished" ordered and oriented section of the chromosome.

				#Also, in a very few cases in the April 2003 assembly, the random files contain data related to sequence for alternative haplotypes. This is present primarily in chr6, where we have included two alternative versions of the MHC region in chr6_random. There are a few clones in other chromosomes that also correspond to a different haplotype. Because the primary reference sequence can only display a single haplotype, these alternatives were included in random files. In subsequent assemblies, these regions have been moved into separate files (e.g. chr6_hla_hap1).

				#More info:
					#https://www.biostars.org/p/7629/
					#http://vega.archive.ensembl.org/info/data/MHC_Homo_sapiens.html
					#regions under review: https://www.ncbi.nlm.nih.gov/grc/human

			#We should remove all the rows that do not belong to the typical autosomal and X chromosome. Y has very few genes. The X will be used for another paper, but we will take the data. As we do not have gene ids for other chromosomes, no other chromosome will be considered when calculating density, because the windows will be calculated around genes of 1:22 and X.

		#chromStart: The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
		summary(raw_gc_data$V3) #min chromosome start is 0. In contrast with tbfs, no areas with a low score were removed. 

		#chromEnd: The ending position of the feature in the chromosome or scaffold.
		summary(raw_gc_data$V4) #the min (4260) is higher than the min of chromosome start, which makes sense.
			
		#Important: Like for tbfs, in this data: all start coordinates in our database are 0-based, not 1-based. Therefore, the chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), because the start is included instead of the ending ("0-based, not 1-based"), BUT note this will be represented as the position notation chr1:1-100.

			#Note that this is clearly stated in the scheme of "GC Percent - GC Percent in 5-Base Windows" as it is said "all start coordinates in our database are 0-based, not 1-based. See explanation here", referring to the same link that in the case of tbfs: "http://genome.ucsc.edu/FAQ/FAQtracks#tracks1"

			#I have checked the first row, start at 10000 (+1) and end at 15120. According to the browser, this region has 1024 data points of GC percentage, being the mean 58.5984. The mean calculate with the sum of squared of our table is 60360/1024=58.94. There is a small difference, I have checked other regions and difference is several decimals (0.3-0.4). 
				#"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A10001%2D15120&hgsid=801849049_39a6h5dF0SWu20mQNJOBVVVySYur"
				#detail "https://genome.ucsc.edu/cgi-bin/hgc?c=chr1&l=10000&r=15120&o=10000&t=15120&g=gc5Base&i=gc5Base&db=hg19"
				#For doing this query I had to select "full" in GC percent (Mapping and Sequencing).

			#THIS IS IMPORTANT: https://genome.ucsc.edu/FAQ/FAQtracks#tracks1
				#FAQ UCSC: I am confused about the start coordinates for items in the refGene table. It looks like you need to add "1" to the starting point in order to get the same start coordinate as is shown by the Genome Browser. Why is this the case?
					#Our internal database representations of coordinates always have a zero-based start (first base is zero) and a one-based end (the last number typed in the range is not really included). We add 1 to the start before displaying coordinates in the Genome Browser. Therefore, they appear as one-based start, one-based end in the graphical display. The refGene.txt file is a database file, and consequently is based on the internal representation.

					#We use this particular internal representation because it simplifies coordinate arithmetic, i.e. it eliminates the need to add or subtract 1 at every step. If you use a database dump file but would prefer to see the one-based start coordinates, you will always need to add 1 to each start coordinate.

					#If you submit data to the browser in position format (chr#:##-##), the browser assumes this information is 1-based. If you submit data in any other format (BED (chr# ## ##) or otherwise), the browser will assume it is 0-based. You can see this both in our liftOver utility and in our search bar, by entering the same numbers in position or BED format and observing the results. Similarly, any data returned by the browser in position format is 1-based, while data returned in BED format is 0-based.

					#For a detailed explanation, please see our blog entry for the UCSC Genome Browser coordinate counting systems: http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/

				#Following the previous example: 
					#The first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. 
					#We really want the 100 first bases to compare with ensemble, 1-100. If we sum 1 to the start and do nothing to the end, we would have chromStart=1 and chromEnd=100, which is exactly we want. All the positions are moved to the left, right? The end is the first base after the end of the segment, so we do not have to sum 1.

			#The regions are to be not overlapped, right? so a region covered by chromosome start and end in a given chromosome will have only 1 row and 1 score right?
				#I have checked the first 10 rows and seems not overlapped regions.
				#In addition, there is a check in the following lines of this scripts about this, just before loading gene coordinates.

		#name: Defines the name of the item. This label is displayed to the left of line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
		summary(raw_gc_data$V5)
		nrow(raw_gc_data) == length(unique(raw_gc_data$V5)) #No repeated names
		raw_gc_data[which(raw_gc_data$V5 == "chr1.0"),] #1 row. It is an identification of the data point.

		#span: each value spans this many bases
		summary(raw_gc_data$V6) #the window is ALWAYS 5 because we are working with GC percentage calculated in 5 base windows. 

		#count: number of values in this block
		summary(raw_gc_data$V7) #
		raw_gc_data[1,] #If you have 1024 values in the first row, and each value correspond with 5 bases, this block will have 1024*5=5120. The size og this block is 15120 - 10001 + 1 = 5120. MATCHES.
			#Note, we have 1024 data points of GC percentage, but this is summarized into one row by summing all of them (sumData column).

		#offset: offset in File to fetch data
		summary(raw_gc_data$V8) #it seems the total number of data point from the beginning to a given row. 1024+1024 for second 1024+1024+1024 for the third one.

		#file: path name to data file, one byte per value
		summary(raw_gc_data$V9)

		#lowerLimit: lowest data value in this block
		summary(raw_gc_data$V10)

		#dataRange: lowerLimit + dataRange = upperLimit
		summary(raw_gc_data$V11) #it is not higher than 100 (it is GC PERCENTAGE)

		#validCount: number of valid data values in this block
		summary(raw_gc_data$V12)

		#sumData: sum of the data points, for average and stddev calc
		summary(raw_gc_data$V13) #For each region of the genome, they take all the data points of GC percentage and sum them to get sumData. This sum divided by the total number of valid point will give us the average of GC percentage across all the 5 base - windows in a given region. 
			#IMPORTANT: We cannot use the sumData column as indicator of GC content. We have one value of sumData per row, higher GC content would lead to higher sum, but the problem is that not all rows refer to regions of the same size. If you regions have the same sum, but one has half of the size, the smaller has a higher GC content because it is the same content but in less space (see "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/figures_tentative_conf_text/figure_21.jpeg"). Because of this we have to use a metric that consider the N. Given they release the sumData in this table, we can calculate the mean.

		#sumSquares: sum of data points squared, for stddev calc
		summary(raw_gc_data$V14) 



#################################
###### PREPARE THE DATASET ######
#################################

##select only variables of interest
raw_gc_data_subset = raw_gc_data[,c(2,3,4,7,12,13)] #we are going to use only the chromosome, the start/end, the number and valid number of data points and the sum of percentages inside each region
head(raw_gc_data_subset)

#set the corresponding names of the columns
colnames(raw_gc_data_subset) <- c("chrom", "chromStart", "chromEnd", "count", "validCount", "sumData")
str(raw_gc_data_subset)


##remove chromosomes that are not the typical autosomal and sex chromosomes
#chromosome of interest
selected_chromosomes = paste("chr", 1:22, sep="") #Y has very few genes, thus it is not included. The X could be used for another paper, but we have no iHS data, so we will not use it right now. As we do not have gene ids for other chromosomes, no other chromosome will be considered when calculating density, because the windows will be calculated around genes of 1:22.

#select those rows whose chromosome name is included in selected_chromosomes
raw_gc_data_subset = raw_gc_data_subset[which(raw_gc_data_subset$chrom %in% selected_chromosomes),]

#check that everything works well
unique(raw_gc_data_subset$chrom)
all(unique(raw_gc_data_subset$chrom) %in% paste("chr", 1:22, sep=""))
all(paste("chr", 1:22, sep="") %in% unique(raw_gc_data_subset$chrom))
!c("chrY", "chrX") %in% unique(raw_gc_data_subset$chrom) #we have only the 22 autosomes, while we do not have the Y and X chromosome. 


##sum 1 to the start because the start begins at zero in bed files (see above)
raw_gc_data_subset$chromStart_new <- (raw_gc_data_subset$chromStart) + 1
#check 
summary(raw_gc_data_subset$chromStart_new-1 == raw_gc_data_subset$chromStart)

#remove the previous chromStart variable
raw_gc_data_subset$chromStart <- NULL

#reorder columns
raw_gc_data_subset = raw_gc_data_subset[,c("chrom", "chromStart_new", "chromEnd", "count", "validCount", "sumData")]

#see the new structure
str(raw_gc_data_subset)
head(raw_gc_data_subset)

#check that start is always smaller than the end
summary(raw_gc_data_subset$chromStart_new < raw_gc_data_subset$chromEnd)

#check if there are difference in the number of counts and valid counts
summary(raw_gc_data_subset$count == raw_gc_data_subset$validCount)
differences_count = raw_gc_data_subset[which(raw_gc_data_subset$count != raw_gc_data_subset$validCount),] #12 cases
summary(differences_count$count > differences_count$validCount) #in all cases count is bigger than validCount. ValidCount seems to indicate the total number of adequate counts. We will use validCount to calculate the mean.


##calculate the GC content multiplied by the total number of points inside each segment just in case want to use this variable
#we have to calculate the mean of data point. Remember, for each region there are several values of GC percentage, as many as 5 base - windows we have inside the region. The final total number of points is validCount, this is the N. All this data points of GC percentage are summarized into 1 value by summing all of them. We have their sum and the N
raw_gc_data_subset$gc_sample_size_product <- (raw_gc_data_subset$sumData * raw_gc_data_subset$validCount)
summary(raw_gc_data_subset$gc_sample_size_product) #IMPORTANT: David told me to multiply the sum of GC of each segment by the number of data points in that segment. Then sum all the products of all segments inside a window and divide it by the sum of the number of points of these segment. I agree with the second part, but not with the first, but I leave this multiplied just in case David tell me to do it. 


##check that the regions with GC content inside each chromosome are not overlapped. We will used Genomic Ranges for that
#load the required package
require(GenomicRanges)

#open an empty data.frame to save the results of the loop
check_overlapping_gc_content = data.frame(chromosome_name_gap_check=NA, check_overlapping_gc_content_1=NA, check_overlapping_gc_content_2=NA, check_overlapping_gc_content_3=NA)

#for each chromosome in raw_gc_data_subset (you can have the same coordinate in different chromosomes)
for(i in 1:length(unique(raw_gc_data_subset$chrom))){

	#extract the name of the [i] chromosome
	chromosome_name_gap_check = unique(raw_gc_data_subset$chrom)[i]

	#subset the gc content data for the [i] chromosome
	gc_content_subset_chromosome = raw_gc_data_subset[which(raw_gc_data_subset$chrom == chromosome_name_gap_check),]

	#check that you have only the chromosome [i]
	check_overlapping_gc_content_1 = unique(gc_content_subset_chromosome$chrom) == chromosome_name_gap_check

	#convert the segments with GC content in the [i] chromosome as a IRange object
	gc_ranges = IRanges(start=gc_content_subset_chromosome$chromStart_new, end=gc_content_subset_chromosome$chromEnd) #create a IRange object with the start and end of each segment with GC content.

	#convert the IRange file into a genomic range file
	gc_ranges_gr = GRanges(seqnames=gc_content_subset_chromosome$chrom, ranges=gc_ranges, strand = '*') #sequnames are the names of the chromosome names. This is very important, because if we select as seq names de exons ids, it consider each exon as independent and do not search for overlapping ranges within all exons; ranges are the start and end of each sequence, and are created with the function IRanges. I have checked that when an exon occupy the whole window, if you use exon id as sequence name, the complete length of the range is not calculated, I have to use the chromosome name. The same occurs when you try to get all the ranges without overlapping with disjoin, only are considered the ranges of the same seqname; strand indicate the sense of the sequence. We have the information of the strand in the strand variable, but as I have revised previously, genomic start/end have the same sense independently of the strand. Just in case, i will use *. In any case, I have seen no changes if this modified in the exon data (i have checked the no impact of strand in this exact example for gap length).
				#this explanation is taken from gene_coordinates_v8.r

	#check that you have all segments of the [i] chromosome inside the genomic range object
	check_overlapping_gc_content_2 = length(gc_ranges_gr@ranges) == nrow(gc_content_subset_chromosome)

	#calculate the not overlapping ranges with GC content
	gc_ranges_gr_not_overlapped = disjoin(gc_ranges_gr, with.revmap=TRUE) #‘disjoin’ returns an object of the same type as ‘x’ containing disjoint ranges for each distinct (seqname, strand) pairing. Split all the segments to have no overlap. If for example you have 1-10 and 5-7, the first range would be 1-4, then 5-7, and finally 8-10. You have the second range included within the first one. If ‘with.revmap=TRUE’, a metadata column that maps the ouput ranges to the input ranges is added to the returned object. This is basically a map, tells you what initial ranges are included in each of the new non-overlapped range.

	#compare overlapped and not overlapped ranges
	check_overlapping_gc_content_3 = identical(gc_ranges_gr@ranges, gc_ranges_gr_not_overlapped@ranges) #if there is no differences, then no overlap exists between segments with GC content.

	#save results
	check_overlapping_gc_content = rbind.data.frame(check_overlapping_gc_content, cbind.data.frame(chromosome_name_gap_check, check_overlapping_gc_content_1, check_overlapping_gc_content_2, check_overlapping_gc_content_3))
}

#remove the firs row with NAs
check_overlapping_gc_content = check_overlapping_gc_content[-which( rowSums(is.na(check_overlapping_gc_content)) == ncol(check_overlapping_gc_content) ),]

#check results
nrow(check_overlapping_gc_content) == 22 #we have all the autosomal chromosomes
summary(check_overlapping_gc_content) #no FALSE
summary(check_overlapping_gc_content$check_overlapping_gc_content_3) #Specially important the case of check 3. TRUE means that the no overlap exists between the segments with GC content in a given chromosome. 


##load coordinate data
gene_coordinates = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt', sep='\t', header=TRUE)
str(gene_coordinates)
head(gene_coordinates)
unique(gene_coordinates$chromosome) #Only autosomal. Y has very few genes. The X can be used for another paper, but we will take the data later as we do not have iHS data. As we do not have gene ids for other chromosomes, no other chromosome will be considered when calculating density, because the windows will be calculated around genes of 1:22.

#remove duplicates
gene_coordinates_no_duplicates = gene_coordinates[-which(duplicated(gene_coordinates$gene_id)),]
	#remember we have all exons of each gene, so we have several rows for the same gene. You have the same windows in all the rows of the same gene. 

#check that the remove of duplicates is correct
#coordinates_data_row_dupli_test = gene_coordinates[which(gene_coordinates$gene_id == "ENSG00000000457"),] #extract the rows of the first gene, with all the exons to check that the function will work in ddply. We will introduce a data.frame and "gene_id" as variable to split, so we will get a value of each window per gene.
#write a function to make the check
if(FALSE){ #we run it one time and worked. No more because it spends to much memory.
check_windows_duplicate = function(coordinates_data_row_dupli_test){

	#extract the gene id of the gene under study
	selected_gene_id = unique(coordinates_data_row_dupli_test$gene_id)

	#check that the number of cases for upper/lower limit in each size is exactly 1, i.e., for all the rows of a given gene, we have the same windows 
	lower_end_window_50kb_identical = length(unique(coordinates_data_row_dupli_test$lower_end_window_50kb)) == 1   
	upper_end_window_50kb_identical = length(unique(coordinates_data_row_dupli_test$upper_end_window_50kb)) == 1   
	lower_end_window_100kb_identical = length(unique(coordinates_data_row_dupli_test$lower_end_window_100kb)) == 1  
	upper_end_window_100kb_identical = length(unique(coordinates_data_row_dupli_test$upper_end_window_100kb)) == 1  
	lower_end_window_200kb_identical = length(unique(coordinates_data_row_dupli_test$lower_end_window_200kb)) == 1  
	upper_end_window_200kb_identical = length(unique(coordinates_data_row_dupli_test$upper_end_window_200kb)) == 1  
	lower_end_window_500kb_identical = length(unique(coordinates_data_row_dupli_test$lower_end_window_500kb)) == 1  
	upper_end_window_500kb_identical = length(unique(coordinates_data_row_dupli_test$upper_end_window_500kb)) == 1  
	lower_end_window_1000kb_identical = length(unique(coordinates_data_row_dupli_test$lower_end_window_1000kb)) == 1 
	upper_end_window_1000kb_identical = length(unique(coordinates_data_row_dupli_test$upper_end_window_1000kb)) == 1 

	#save the gene id and the checks
	results_check_coord_duplicate = cbind.data.frame(selected_gene_id, lower_end_window_50kb_identical, upper_end_window_50kb_identical, lower_end_window_100kb_identical, upper_end_window_100kb_identical, lower_end_window_200kb_identical, upper_end_window_200kb_identical, lower_end_window_500kb_identical, upper_end_window_500kb_identical, lower_end_window_1000kb_identical, upper_end_window_1000kb_identical)

	#return the results
	return(results_check_coord_duplicate)
}
#apply the function to gene_coordinates, splitting the data.frame for each gene id (all row of a gene together)
check_coords_dupli = ddply(.data=gene_coordinates, .variables="gene_id", .fun=check_windows_duplicate, .inform=TRUE, .parallel=FALSE, .paropts=NULL)
	#".inform=TRUE" generates and shows the errors. This increases the computation time, BUT is very useful to detect problems in your analyses.
	#".parallel" to paralelize with foreach. 
	#".paropts" is used to indicate additional arguments in for each, specially interesting for using the .export and .packages arguments to supply them so that all cluster nodes have the correct environment set up for computing. 
		#ADD PACKAGES USED INSIDE THE FUNCTION

#check we have the correct number of rows
nrow(check_coords_dupli) == length(unique(gene_coordinates$gene_id))
#check that all cases we have only TRUE, i.e., only 1 upper/lower limit for each window and gene. 
summary(check_coords_dupli)

#check the gene ids are correct
if(all(check_coords_dupli$gene_id == check_coords_dupli$selected_gene_id)){ #if all TRUE

	#remove the second column with id
	check_coords_dupli$selected_gene_id <- NULL
} else {

	#if not we have an error
	stop("ERROR!! We have a problem with the gene id in check_coords_dupli!!!!")
}
}



###############################################
###### CALCULATE GC CONTENT OF EACH GENE ######
###############################################

#write a function to calculate GC content for each gene and window 
#coordinates_data_row = gene_coordinates_no_duplicates[1,] #extract the first row to check the function will work with ddpply. We will do it for each gene id. The dataset used only have one row for each gene, so we do not need to extract all the rows of a given gene. The first row is for the first gene, the second row is for the second gene, and so on. Having one row for gene let us to obtain several outputs for each row (gen) and save them as different rows. So we will have 5 rows for each gene, one for each window. Each row will have its corresponding check
#coordinates_data_row = gene_coordinates_no_duplicates[which(gene_coordinates_no_duplicates$gene_id=="ENSG00000186092"),] #gen very close to the start of the chromosome 1. To check that the script goes well with windows without coordinates (removed because of the window surpass the chromosome limit)
gc_calc = function(coordinates_data_row){

	#extract the gene id of the row
	selected_gene_id = coordinates_data_row$gene_id

	#extract the chromosome name of the row. We add chr because this is the notation used in conserved elements data
	selected_chromosome = paste("chr", coordinates_data_row$chromosome_name, sep="")

	#subset conserved elements data for the selected chromosome
	gc_data = raw_gc_data_subset[which(raw_gc_data_subset$chrom == selected_chromosome),]

	#extract the opposite subset to make a check
	opposite_gc_data = raw_gc_data_subset[-which(raw_gc_data_subset$chrom == selected_chromosome),]

	#make a check. If the chromosome names of the working dataset are exactly the selected_chromosome of the row under study, and selected chromosome is not included in the opposite subset
	if(all(gc_data$chrom == selected_chromosome) & !selected_chromosome %in% opposite_gc_data$chrom){

		#everything is ok
		check_1 = TRUE
	} else {

		#problem
		check_1 = FALSE
	}

	#open an empty data frame to save the conserved elements density of each gene and window
	gc_content_final = data.frame(selected_gene_id=NA, selected_window=NA, selected_lower_window=NA, selected_upper_window=NA, gc_content=NA, check_1=NA, check_2=NA, check_3=NA, check_4=NA)

	#create a vector with the windows sizes
	window_size = c("50", "100", "200", "500", "1000")

	#for each window size
	for(k in 1:length(window_size)){

		#select the [k] window
		selected_window = window_size[k]

		#extract the values of the lower and upper limit of each window for the row under analysis
		selected_lower_window = eval(parse(text=paste("coordinates_data_row$lower_end_window_", selected_window, "kb", sep="")))
		selected_upper_window = eval(parse(text=paste("coordinates_data_row$upper_end_window_", selected_window, "kb", sep=""))) #we have to add as.numeric because if we extract number and name of these vector elements, we cannot compare each one with ALL the starts and ends in the next lines of code

		#if we do not have NA in any of the window limits we perform the calculations
		if(!is.na(selected_lower_window) & !is.na(selected_upper_window)){
		 
			#extract the segments included between the lower and upper windows
			starts_included = which(gc_data$chromStart_new >= selected_lower_window & gc_data$chromStart_new <= selected_upper_window)
			ends_included = which(gc_data$chromEnd >= selected_lower_window & gc_data$chromEnd <= selected_upper_window) #we only want the segments completely included within the window, thus we only want the starts and ends that are within the window. They have to larger or equal than the lower end of the window, and smaller or equal than the upper end of the window. We avoid regions that begin before and end within the window, regions that begin within the window and end after the window, and regions that begin before the window and end after the window.
				#see "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/figures_tentative_conf_text/figure_21.jpeg"
				#Note that the extremes of the windows are included inside the window (see "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/figures_tentative_conf_text/figure_18" and "gene_coordinates_v8.r").

			#now, we select those rows with both starts and ends are included inside the window. This is done obtaining the intersection between starts and ends included, i.e., those rows that have both start and end inside the window.  
			segments_included = gc_data[intersect(starts_included, ends_included),]

			#make another check. calculate the segments (rows) inside the window combining all the conditions
			segments_included_way_2 = gc_data[which(gc_data$chromStart_new >= selected_lower_window & gc_data$chromStart_new <= selected_upper_window & gc_data$chromEnd >= selected_lower_window & gc_data$chromEnd <= selected_upper_window),]
			#both dataset are equal?
			check_2 = identical(segments_included, segments_included_way_2)

			##additional checks
			#check that start and end of the segments included are exactly within the windows
			check_3 = all(segments_included$chromStart_new >= selected_lower_window & segments_included$chromStart_new <= selected_upper_window & segments_included$chromEnd >= selected_lower_window & segments_included$chromEnd <= selected_upper_window)
			
			#check that the rest of the conserved data dataset not included contains only segments not included within the window
			segments_not_included = gc_data[-intersect(starts_included, ends_included),]
			check_4 = all(segments_not_included$chromStart_new < selected_lower_window | segments_not_included$chromStart_new > selected_upper_window & segments_not_included$chromEnd < selected_lower_window | segments_not_included$chromEnd > selected_upper_window)
				#if the start is higher than the upper end of the window, then the segment is for sure outside of the window. The segment does not reach the window.
				#if the start is lower that the lower end of the window, then the segment could end inside the window (end inside or after the upper end), but we have additional conditions
					#if the end occurs before the lower end of the window, the segment does not reach the window
					#if the end occurs after the upper end of the window, the segment overlap the full window but have areas outside of it. We only select segments that completely inside the window.
				#if the end is lower than the lower end of the window, then the segment is for sure outside of the window. The segment does not reach the window.
				#if the end is higher than upper end of the window, then the segment could end inside the window (start inside or before the lower end), but we have additional conditions. 
					#if the start is higher than the upper window, then the segment is completely outside of the window.
					#if the start is lower that the lower window, then the segment overlap the full window but have areas outside of it. We only select segments that completely inside the window. 
				#Under the current definition (only segments fully included inside the window), the segments that follow any of these conditions should be included in our analyses.

			#calculate the gc content as the sum of GC of all segments inside the window divided by the sum of their length (number of data points)
			gc_content = sum(segments_included$sumData) / sum(segments_included$validCount)
			    #You told me I should multiply the GC content by the number of points in the corresponding segment. Then sum all the products for the segments inside the window and then divided by the sum of all lengths.
    				#I see the rationale of summing the number of point of all segments inside the window and use it as denominator. For the same GC content, windows with small number of data points, i.e., lower number of 5 base windows inside the segments, should be overweighted. It is more difficult for smaller region to have the same GC content than the bigger ones, but the density could be similar. We have to account for that. 
    				#but why we have to multiply GC by the number of points in the numerator. If you have a higher number of segments, then it is more likely to have more GC but not more density. Regions with more points will be overweighted. We do not want that.
    				#For now, I have sum all the GC content inside the window, and divided by sum of all lengths inside the window.
				#This is like a mean. Imagine you have 2 intervals, with 2 data points each one. Remember that each data point is the GC percentage in 5 bases. 
					#Scenario A: The sum of GC percentages is 100+100 for the first interval and 100+100 in the second interval. Therefore, we have 400, and this will be divided by 4, which is the number of points. 400/4=100. 100% of GC content, which is correct.
					#Scenario B: The sum of GC percentages is 100+100 for the first interval and 0+0 in the second interval. Therefore, we have 200, and this will be divided by 4, which is the number of points. 200/4=50. 50% of GC content, which is correct. Half of the regions included have GC content. 
		} else { #if not, and hence we have Na in any of the [j] window coordinates

			#set all results as NA
			gc_content = NA
			check_1 = NA
			check_2 = NA
			check_3 = NA
			check_4 = NA
		}

		#save conserved elements density data
		gc_content_final = rbind.data.frame(gc_content_final, cbind.data.frame(selected_gene_id, selected_window, selected_lower_window, selected_upper_window, gc_content, check_1, check_2, check_3, check_4))
	}

	#remove the first row with NAs (all entries are NA)
	gc_content_final = gc_content_final[-which(rowSums(is.na(gc_content_final)) == ncol(gc_content_final)),] #remove rows for which the sum of entries with NA is the same than the total number of columns, i.e., all entries have a NA.

	#return the final dataset
	return(gc_content_final)
}

#apply the function for each row of gene coordinates
final_gc = ddply(.data=gene_coordinates_no_duplicates, .variables="gene_id", .fun=gc_calc, .inform=TRUE, .parallel=FALSE, .paropts=NULL)
	#".inform=TRUE" generates and shows the errors. This increases the computation time, BUT is very useful to detect problems in your analyses.
	#".parallel" to paralelize with foreach. 
	#".paropts" is used to indicate additional arguments in for each, specially interesting for using the .export and .packages arguments to supply them so that all cluster nodes have the correct environment set up for computing. 
		#ADD PACKAGES USED INSIDE THE FUNCTION

	#NOTE: ddply reorder the output rows in basis on the variable used to split ("gene_id"), so the first gene id is the lowest number, which correspond with a gene of the chromosome X (ENSG00000000003). It should not be a problem, because we have the exact ID in each row, so the different datasets can be merged. When we will have to remove X genes, we have to select the gene_id of the autosomal chromosomes. 
		#I could use gene_id and chromosome_name as variables to split in ddply. In that way, I would get the rows ordered by chromosome name, and additional column with the chromosome name. I did not know this at the beginning, so I will let it in the way it is for this and conserved elements and tbfs, because it is not problematic. 

#check the gene ids are correct
if(all(final_gc$gene_id == final_gc$selected_gene_id)){ #if all TRUE

	#remove the second column with id
	final_gc$selected_gene_id <- NULL
} else {

	#if not we have an error
	stop("ERROR!! We have a problem with the gene id in final_gc!!!!")
}

#check we have all the genes
nrow(final_gc) == length(unique(gene_coordinates$gene_id))*5 #5 is the number of windows, for each gene we have 5 windows, so the number of rows should five times the number of genes

#take a look to the checks
summary(final_gc) #here are some windows without data

## NAs in gc content calculation
#take a look to these windows without gc content
genes_no_gc_content = final_gc[which(is.na(final_gc$gc_content)),]
nrow(genes_no_gc_content)

#check that the windows without median calculate are really windows without GC data inside
#extract the gene coordinates of the problematic genes
problematic_genes = gene_coordinates_no_duplicates[which(gene_coordinates_no_duplicates$gene_id %in% genes_no_gc_content$gene_id),]

#make a loop to make the check
test_na_genes = data.frame(gene_id=NA, selected_window=NA, gc_content=NA, chrom=NA, chromStart_new=NA, chromEnd=NA, count=NA, validCount=NA, sumData=NA, gc_sample_size_product=NA)
for(i in 1:nrow(problematic_genes)){

	#select the [i] gene
	selected_gene = problematic_genes[i,]

	#extract gene_id
	selected_gene_id = selected_gene$gene_id

	#extract the chromosome
	selected_chrom = selected_gene$chrom

	#extract row of the [i] gene from GC results
	gc_results = genes_no_gc_content[which(genes_no_gc_content$gene_id == selected_gene_id),] #we select from here the problematic window for the [i] gene

	#for each gc results check the number of GC data included in the corresponding window used
	for(k in 1:nrow(gc_results)){

		#selected gc result
		selected_gc_result = gc_results[k,]

		#extract the lower and upper end of the problematic window for the [i] gene
		lower_window = eval(parse(text=paste("selected_gene$lower_end_window_", selected_gc_result$selected_window, "kb", sep="")))
		upper_window = eval(parse(text=paste("selected_gene$upper_end_window_", selected_gc_result$selected_window, "kb", sep="")))

		#subset the GC data for the chromosome of the [i] gene
		gc_data = raw_gc_data_subset[which(raw_gc_data_subset$chrom == paste("chr", selected_chrom, sep="")),]

		#extract rows of GC data included in the window
		gc_window = gc_data[which(gc_data$chromStart_new >= lower_window & gc_data$chromStart_new <= upper_window & gc_data$chromEnd >= lower_window & gc_data$chromEnd <= upper_window),]

		#add gene_id
		if(nrow(gc_window) != 0){ #if we have GC data inside the window simply add the gene id and other data of this window
			gc_window = cbind.data.frame(selected_gc_result[,c("gene_id", "selected_window", "gc_content")], gc_window)
		} else { #if not and hence we do not have any row, set a NA row for the gc_window file, and then bind gene id and other data of this window
			gc_window[1,] <- NA
			gc_window = cbind.data.frame(selected_gc_result[,c("gene_id", "selected_window", "gc_content")], gc_window)
		}

		#save
		test_na_genes = rbind.data.frame(test_na_genes, gc_window)
	}
}
#remove the first row with NAs (all entries are NA)
test_na_genes = test_na_genes[-which(rowSums(is.na(test_na_genes)) == ncol(test_na_genes)),] #remove rows for which the sum of entries with NA is the same than the total number of columns, i.e., all entries have a NA.
nrow(test_na_genes) == nrow(genes_no_gc_content)

#take a look to the results
test_na_genes #EACH ROW comes from a segment with GC values, i.e., a region of the genome for which we have several values of GC per 5 base windows. With these values a mean and sd of GC percentage per region is calculated (mean_gc_percentage; sd_gc_percentage). For each gene and window that showed problems, we search all the regions that are inside the corresponding window. In all case, we have no data (NA), i.e., no region with GC data is included in the window. See inside the loop, when no segment is included inside the window I set NA for all the row.
	#in addition, we have cases with no window coordinate because they surpass the chromosome limit. 

#check that the median percentage per window is always between 0 and 100
min(na.omit(final_gc$gc_content)) >= 0 & max(na.omit(final_gc$gc_content)) <= 100


##take a look to those cases with window coordinate but the final window is full overlapped with a gap
#load tbfs data where we can see this
tbfs_data = read.table('/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/tbfs_density_raw_v3.txt', sep='\t', header=TRUE)

#from final_gc, select those rows with a gene_id that is included within the gene_id of the rows in tbfs_data that have window coordinates but the final window has a length of zero (i.e., window overlapped with a gap completely)
final_gc[which(final_gc$gene_id %in% tbfs_data[which(!is.na(tbfs_data$selected_lower_window) & !is.na(tbfs_data$selected_upper_window) & tbfs_data$window_length_no_gaps == 0),]$gene_id),] 
	#The three 50kb windows that are completely overlapped have no GC content, perfect!
		#We calculate the average percentage of GC by summing all the GC percentages within the window and divided by the number of valid points inside the window. If the window is completely overlap with a gap, we have no GC content, valid count os zero, and hence division by zero is NaN. If you have GC content data in the window, you will have a value, if not, you will have NA, even if we have windows coordinates. I have checked this with "ENSG00000185974".
	#In addition, for gene "ENSG00000185974" we have NA in the GC content 100kb window. 

#if you take a look on the 100kb window of ENSG00000185974, you have DNAse data, but not GC content.
tbfs_data[which(tbfs_data$gene_id == "ENSG00000185974"),] #the length of the window without gaps is only 5095, is very small, and indeed in the ucsc browser there is no GC content in the 100kb window (https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A114331088%2D114406088&hgsid=839993237_n8DYGiqfDekIq67wD1PS2AuSD39c), but we have GC content in the 200kb window (https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr13%3A114281088%2D114481088&hgsid=839993237_n8DYGiqfDekIq67wD1PS2AuSD39c).


#you can extract results of each gene
final_gc[which(final_gc$gene_id == "ENSG00000000457"),] #for each row, i.e., each gene because we removed duplicates, we extract calculate the gc content inside each window

#save the data
write.table(final_gc, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gc_content_raw_v3.txt', sep='\t', row.names=FALSE, col.names=TRUE)



####################################################################################################
###### OBTAIN A DATA.FRAME WITH ONE ROW PER GENE AND ONE COLUMN FOR GC CONTENT OF EACH WINDOW ######
####################################################################################################

#write a function to extract the gc content for each gene
#final_gc_row = final_gc[which(final_gc$gene_id == "ENSG00000000457"),] #extract the rows of the first gene, with all the window lengths to check that the function will work in ddply. We will introduce a dataframe and "gene_id" as variable to split, so we will get a value of each window per gene.
extract_gc = function(final_gc_row){

	#extract the selected gene id
	selected_gene_id = unique(final_gc_row$gene_id)

	#for the element of the list under study, extract the conserved element value for each window, from 50 to 1000kb
	gc_content_50kb = final_gc_row[which(final_gc_row$selected_window == 50),]$gc_content
	gc_content_100kb = final_gc_row[which(final_gc_row$selected_window == 100),]$gc_content
	gc_content_200kb = final_gc_row[which(final_gc_row$selected_window == 200),]$gc_content
	gc_content_500kb = final_gc_row[which(final_gc_row$selected_window == 500),]$gc_content
	gc_content_1000kb = final_gc_row[which(final_gc_row$selected_window == 1000),]$gc_content

	#return all conserved values values
	return(cbind.data.frame(selected_gene_id, gc_content_50kb, gc_content_100kb, gc_content_200kb, gc_content_500kb, gc_content_1000kb))
}

#we will use ddply because you can use a data.frame as an input and save it as an data.frame
gc_results_extracted_final = ddply(.data=final_gc, .variables="gene_id", .fun=extract_gc, .inform=TRUE, .parallel=FALSE, .paropts=NULL) 
	#".inform=TRUE" generates and shows the errors. This increases the computation time, BUT is very useful to detect problems in your analyses.
	#".parallel" to paralelize with foreach. 
	#".paropts" is used to indicate additional arguments in for each, specially interesting for using the .export and .packages arguments to supply them so that all cluster nodes have the correct environment set up for computing. 
		#ADD PACKAGES USED INSIDE THE FUNCTION

#check the gene ids are correct
if(all(gc_results_extracted_final$gene_id == gc_results_extracted_final$selected_gene_id)){ #if all TRUE

	#remove the second column with id
	gc_results_extracted_final$selected_gene_id <- NULL
} else {

	#if not we have an error
	stop("ERROR!! We have a problem with the gene id in gc_results_extracted_final!!!!")
}

#check we have all the genes
nrow(gc_results_extracted_final) == length(unique(gene_coordinates$gene_id))

#make checks for the extraction for each window with 50, 100, 200, 500, 10000
all(na.omit(final_gc[which(final_gc$selected_window == 50),]$gc_content) == na.omit(gc_results_extracted_final$gc_content_50kb))
all(na.omit(final_gc[which(final_gc$selected_window == 100),]$gc_content) == na.omit(gc_results_extracted_final$gc_content_100kb))
all(na.omit(final_gc[which(final_gc$selected_window == 200),]$gc_content) == na.omit(gc_results_extracted_final$gc_content_200kb))
all(na.omit(final_gc[which(final_gc$selected_window == 500),]$gc_content) == na.omit(gc_results_extracted_final$gc_content_500kb))
all(na.omit(final_gc[which(final_gc$selected_window == 1000),]$gc_content) == na.omit(gc_results_extracted_final$gc_content_1000kb))
	#we add NA because we have 10 rows without GC content for 50 and 100 kb windows

#save the data
write.table(gc_results_extracted_final, '/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gc_content_final_v3.txt', sep='\t', row.names=FALSE, col.names=TRUE)

#see summary
summary(gc_results_extracted_final)
	#NAs: There is 1 gene for which 50b windows have been removed. 6 in the case of 100kb windows, 22 for 200kb, 121 for 500kb and 293 for 1000kb. These are the cases removed because of the window surpass the chromosome limit.
		#All matches except in 50kb and 100kb. 
			#In 50kb, we have the gene without 50kb window coordinates (ENSG00000272636), the three with 50kb windows completely overlapped with a gap (ENSG00000148357, ENSG00000168280, ENSG00000185974) and then 5 more cases that are very close to gap, we have bases, thus the densities can be calculated (0 divided by a small number of bases, giving zero, which is a TRUE result), but no GC content data. If no GC value is found, we cannot have estimates of GC. GC is across the whole chromosome, while the other factors only have datapoint where a region is considered coding, conserved or DNAse binding site. In that cases, if we have bases but no one is considered any of these things, we can set zero. 
			#I know exactly the 100kb window that makes the difference (7 instead of 6). The 50kb of that gene is completely overlapped with a gap. The 100kb has areas outside the gap, but very few of them, therefore is normal if we do not have GC content there. 
	#The mean and median coding density is decreasing from 50kb windows to 1000kb windows. Bigger windows have GC data points, therefore the median GC content have to be a lot of more to reach the similar GC content than smaller windows.
	#The minimum is 30 and the maximum is 65, which sounds good for a percentage.

#save the workspace
save.image("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/rdata/gc_content_v3.RData")
require(plyr)
require(GenomicRanges)