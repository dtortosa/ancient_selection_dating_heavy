#!/usr/bin/env python3.9
# coding: utf-8
    #to run this script: chmod +x script.py; ./script.py
    #!/bin/sh does not work with my terminal en msi of David.
    #if you are using "$" to paste the path of the executable, you do not need to use "./" for running the executable.
    #you can save the output and the errors
        #./script.py > script.out #only output
        #./script.py 2> error.out #only error
        #./script.py > script.out 2> error.out #both in different files
        #./script.py > script.out 2>&1 #both in the same file
        #https://www.cyberciti.biz/faq/linux-redirect-error-output-to-file/




#####################################################
######## CALCULATE THE DISTANCE TO BAT GENES ########
#####################################################

#This script will calculate the distance to each coding gene to the closest gene related to the BAT connectome




##################
# import modules #
##################

import pandas as pd
import numpy as np




########################################
# define function to print text nicely #
########################################

#text="checking function to print nicely"
#header=1
def print_text(text, header=2):
    if header==1:
        print("\n#######################################\n#######################################")
        print(text)
        print("#######################################\n#######################################")
    elif header==2:
        print("\n###### " + text + " ######")
    elif header==3:
        print("\n## " + text + " ##")
    elif header==4:
        print("\n# " + text + " #")
print_text("checking function to print nicely: header 1", header=1)
print_text("checking function to print nicely: header 2", header=2)
print_text("checking function to print nicely: header 3", header=3)
print_text("checking function to print nicely: header 4", header=4)




########################################
# define function to run bash commands #
########################################

#create a wrapper for subprocess.run in order to define a set of arguments and avoid typing them each time. We will ensure that we are using bash always and not sh.
from subprocess import run, PIPE
#command="ls"
def run_bash(command, return_value=False):

    #run the command
    complete_process = run(
        command, 
        shell=True,
        executable="/bin/bash", 
        stdout=PIPE,
        stderr=PIPE, 
        text=True)
    #we have to use popen in order to ensure we use bash, os.system does not allow that
        #shell=True to execute the command through the shell. This means that the command is passed as a string, and shell-specific features, such as wildcard expansion and variable substitution, can be used.
            #THIS IS DANGEROUS IF UNTRUSTED DATA
        #executable="/bin/bash" to ensure the use of bash instead of sh
        #stdout=PIPE to capture the output into an python object. You can also capture the error doing stderr=PIPE. stdout and stderr are the standard output and error
            #you could also use capture_output=True to capture both stdout and stderr
        #text=True will return the stdout and stderr as string, otherwise as bytes
            #https://www.datacamp.com/tutorial/python-subprocess
            #https://docs.python.org/3/library/subprocess.html#subprocess.run

    #this generated a CompletedProcess instance where you can get
        #args: The command and arguments that were run.
        #returncode: The return code of the subprocess.
        #stdout: The standard output of the subprocess, as a bytes object.
        #stderr: The standard error of the subprocess, as a bytes object.

    #if the return code is 0, i.e., success 
        #https://askubuntu.com/a/892605
    if complete_process.returncode==0:

        #if stderr is empty
        if complete_process.stderr=="":

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
        else:

            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #print the standard error without stopping
            print(complete_process.stderr)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
    else:
        if(complete_process.returncode==1) & (complete_process.stderr=="warning [./data/all_human_gene-specific_connectomes_122015.zip]:  4294967296 extra bytes at beginning or within zipfile\n  (attempting to process anyway)\n"):
                #Error obtained from unzip. returncode=1 is only a warning according to unzip man.
                    #https://linux.die.net/man/1/unzip
            
            #print the standard output without "\n" and other characters
            print(complete_process.stdout)

            #print the standard error without stopping
            print(complete_process.stderr)

            #return also the value if required
            if return_value==True:
                return complete_process.stdout
        else:
            #print the standard error and stop
            raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print_text("check behaviour run_bash", header=1)
print_text("see working directory", header=2)
run_bash("pwd")
print_text("list files/folders there", header=2)
run_bash("ls")





#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
#we will use a bash script to run this python program two times instead of doing a function and parallelize that function. In this way, the same python script is run for each selective pressure and if we want to run again one selective pressure, we just need to modify the bash script, not the whole python program. This makes sense because we have only two selective pressured, and the parallelization will occur inside each selective pressure. Importantly, we can also have separated .out files, so we can look at the checks separately.
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--selective_pressure", type=str, default="smt", help="Name of the selective pressure used as input. Always string.")
parser.add_argument("--n_cores", type=int, default=5, help="Number of cores/threads requested. Integer always, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
pressure_name = args.selective_pressure
n_cores = args.n_cores




############################
# prepare folder structure #
############################
print_text("prepare folder structure", header=1)
run_bash(" \
    mkdir \
        --parents \
        ./data/" + pressure_name + "_distance; \
    ls -l ./data")




################
# prepare data #
################
print_text("prepare data", header=1)
print_text("load the required files", header=2)
print_text("load gene coordinates", header=3)
gene_coords = pd.read_table(\
	"/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", \
	sep="\t", \
	header=0, \
	low_memory=False)
print(gene_coords)


print_text("Select only one row per gene because in this file, each row is an exon, but all the rows of the same gene have the same middle gene position and gene windows as all the rows belongs to the same gene", header=4)
gene_coords_no_duplicated = gene_coords[ \
	~gene_coords["gene_id"].duplicated(keep="first")]
    #set as True all duplicates except the first occurrence.
    #then negate with "~", so all duplicates are False except the first occurrence, leading to select only that first occurrence
print(gene_coords_no_duplicated)

print("Check that the subset has as many rows as unique gene IDs in the original file. If True, remove the original gene coirdinate file")
check_gene_coords_subset = len(gene_coords["gene_id"].unique()) == gene_coords_no_duplicated.shape[0]
if(check_gene_coords_subset):
    del(gene_coords)
else:
    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE SUBSET OF THE GENE COORDINATE FILE")


print_text("select only the columns we need", header=4)
gene_coords_no_duplicated_subset = gene_coords_no_duplicated[[ 
    "chromosome_name", \
    "gene_id", \
    "hgnc_symbol", \
    "gene_start", \
    "gene_end", \
    "middle_point", \
    "lower_end_window_50kb", \
    "upper_end_window_50kb",  \
    "lower_end_window_100kb", \
    "upper_end_window_100kb",  \
    "lower_end_window_200kb", \
    "upper_end_window_200kb", \
    "lower_end_window_500kb", \
    "upper_end_window_500kb", \
    "lower_end_window_1000kb", \
    "upper_end_window_1000kb"]]
print(gene_coords_no_duplicated_subset)


print_text("reset the row index", header=4)
gene_coords_no_duplicated_subset = gene_coords_no_duplicated_subset \
    .reset_index(drop=True)
    #drop: Do not try to insert index into dataframe columns. This resets the index to the default integer index.
print(gene_coords_no_duplicated_subset)
    #https://stackoverflow.com/questions/20490274/how-to-reset-index-in-a-pandas-dataframe


print_text("make a depp copy of the dataset to add information about the gene of interest", header=4)
pressure_coords = gene_coords_no_duplicated_subset.copy(deep=True)
    #deep=True: 
        #a new object will be created with a copy of the calling object's data and indices. Modifications to the data or indices of the copy will not be reflected in the original object (see notes below).
print(pressure_coords)



print_text("load pressure data", header=3)
print_text("first extract the connectome of the interest gene", header=4)
if(pressure_name == "bat"):
    core_gene = "UCP1"
elif(pressure_name == "smt"):
    core_gene = "SLN"
run_bash(" \
    unzip \
        -p \
        ./data/all_human_gene-specific_connectomes_122015.zip \
        " + core_gene + ".txt > \
    ./data/" + pressure_name + "_distance/" + core_gene + ".txt")
        #extract the connectome of the gene from a Zip fila with ALL individual connectomes
            #downloaded from 
                #https://lab.rockefeller.edu/casanova/HGC
                #Human gene-specific connectomes - download specific genes
                #last file:
                    #all_human_gene-specific_connectomes_122015
        #unzip -p: 
            #extract files to pipe, no messages
            #we can get a specific file within the zip
                #https://unix.stackexchange.com/a/14125
        #We get a warning with both UCP1 and SNL
            #warning [./data/all_human_gene-specific_connectomes_122015.zip]:  4294967296 extra bytes at beginning or within zipfile (attempting to process anyway)
            #the returncode is equal to 1.
                #According to the man, this is not serious: 
                    #one or more warning errors were encountered, but processing completed successfully anyway. This includes zipfiles where one or more files was skipped due to unsupported compression method or encryption with an unknown password.
                        #https://linux.die.net/man/1/unzip
                #this is not error but warning, error would be returncode 2
                    #a generic error in the zipfile format was detected. Processing may have completed successfully anyway; some broken zipfiles created by other archivers have simple work-arounds.
            #We are unzipping just 1 file, and the unzipping was successful, also this is only a warning, so I think we are good here.
                #I am going to check below in the case of UCP1 is the file is identical to the original used for the BAT analyses. If so, the fact that we get a warning is not affecting.


print_text("load the connectome having " + core_gene + " as core gene", header=4)
pressure_conn = pd.read_csv( \
    "./data/" + pressure_name + "_distance/" + core_gene + ".txt", \
	sep="\t", \
	header=0, \
	low_memory=False)
print(pressure_conn)


print_text("check we have selected the correct connectome, i.e., the one with the correct core gene", header=4)
print(pressure_conn["Target"].iloc[0] == core_gene)


print_text("if UCP1 is the core gene, check that the connectome extracted from the zip is the same as the original used for BAT analyses", header=4)
if(pressure_name == "bat"):
    old_ucp1_conn = pd.read_csv( \
        "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/data/human_connectome/UCP1.txt", \
        sep="\t", \
        header=0, \
        low_memory=False)
    ucp1_2020 = pd.read_csv( \
        "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/data/human_connectome/UCP1_downloaded_2020.txt", \
        sep="\t", \
        header=0, \
        low_memory=False)
        #ucp1 connectome again but downloaded in 2020 (19/06/2020
    print(pressure_conn.equals(old_ucp1_conn))
    print(pressure_conn.equals(ucp1_2020))
        #both files are the same, suggesting that the warning obtained when zipping is not a problem
else:
    print("We are not working with BAT, but with " + pressure_name)


print_text("if UCP1 is the core gene, check that file with BAT relationships has the same genes than the connectome 1% obtained now", header=4)
if(pressure_name == "bat"):
    bat_relationship = pd.read_csv( \
        "/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/appendix_S1_ordered.csv", \
        sep=",", \
        header=0, \
        low_memory=False)
    from natsort import natsort_keygen
    print(pressure_conn \
        .loc[ \
            pressure_conn["Target_in_source_P-value(percentile)"] < 0.01, \
            "Target"] \
        .sort_values( \
            axis=0, 
            key=natsort_keygen(), \
            ignore_index=True) \
        .equals( \
            bat_relationship["Genes"]\
            .sort_values( \
                axis=0, 
                key=natsort_keygen(), \
                ignore_index=True)))
            #from the UCP1 connectome
                #select those rows for which the p-value percentile of the gene is below 1% and get the gene name
            #natural sort the gene names
                #"by" is not needed here because we only have 1 column
                #axis=0: rows
                #key: Apply the key function to the values before sorting.
                    #natsort_keygen()
                        #Generate a key to sort strings and numbers naturally.
                #ignore_index:
                    #If True, the resulting axis will be labeled 0, 1, …, n - 1
                    #so you avoid the previous index
            #check that the resulting series is equals to the Genes included in BAT relationships after sorting in the same way
                #https://stackoverflow.com/a/63890954/12772630
else:
    print("We are not working with BAT, but with " + pressure_name)




print_text("start with the specific subset of genes in the connectome according to P-value percentile", header=2)
#p_value_percentile = 1
for p_value_percentile in [0.5, 1, 5]:

    print_text("Percentile " + str(p_value_percentile) + "% :  select the interest genes genes", header=3)
    selected_connectome_genes = pressure_conn \
        .loc[ \
            pressure_conn["Target_in_source_P-value(percentile)"] < (p_value_percentile/100), \
            "Target"]
    print(selected_connectome_genes)


    print_text("Percentile " + str(p_value_percentile) + "% :  create a variable about pressure status using this set of genes (those belonging to the set of interest are 'yes')", header=3)
    pressure_coords[pressure_name + "_status"] = ["yes" if gene_symbol in selected_connectome_genes.to_list() else "no" for gene_symbol in pressure_coords["hgnc_symbol"]]
    print(pressure_coords[pressure_name + "_status"])


    print_text("Percentile " + str(p_value_percentile) + "% :  check we have set as 'yes' all genes with coords that are included in the interest set of genes", header=4)
    print(pressure_coords \
        .loc[ \
            pressure_coords["hgnc_symbol"].isin(selected_connectome_genes), \
            pressure_name+"_status"] \
        .unique() == "yes")



    print_text("Percentile " + str(p_value_percentile) + "% :  explore pressure data", header=3)
    print_text("Percentile " + str(p_value_percentile) + "% :  look for genes without gene symbol in gene coords", header=4)
    missing_genes = selected_connectome_genes.loc[~selected_connectome_genes.isin(pressure_coords["hgnc_symbol"])]
    if (missing_genes.shape[0] <= np.round(selected_connectome_genes.shape[0]*0.1)):
        print(f"We have lost {missing_genes.shape[0]} genes from the {core_gene} connectome")
        print(missing_genes)
    else:
        raise ValueError("ERROR: FALSE! MORE THAN 10% OF GENES OF THE CONNECTOME ARE NOT INCLUDED IN GENE COORDS FILE")


    print_text("Percentile " + str(p_value_percentile) + "% :  I obtained gene coords from biomart hg19, so I understand that I have all names (IDs and gene names) for those genes that passed my filters. If a interest gene is not included in my gene coord set, then it should not be used for hg19. Indeed, the two genes included in the 168 BAT connectome that are not in gene coordinates, have their gene id NOT included in gene coords", header=4)
    

    #check the missing cases only for the percentile 1%, which is the one used for now in the analyses
    if (p_value_percentile == 1):
        print_text("Percentile " + str(p_value_percentile) + "% :  check that we do not have any of the synonyms of these missing genes included in gene coords", header=4)
        if(pressure_name == "bat"):
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["CD132", "CIDX", "IL-2RG", "IMD4", "P64", "SCIDX", "SCIDX1"]), :].shape[0] == 0)
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["NRU", "P2P", "P2Y4", "UNR"]), :].shape[0] == 0)
                #the two missing genes have no synonmious in the dataset
                    #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000147168;r=X:70327254-70331958
                    #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000186912;r=X:69478016-69479654;t=ENST00000374519
                #THIS CHECK WAS DONE WITH HG19! LOOK AGAIN FOR HG38!
        elif (pressure_name == "smt"):
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["CHORDC3", "ITGB1BP", "MELUSIN"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000147166;r=X:70521584-70525221
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["FMRP", "FRAXA", "MGC87458", "POF", "POF1"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000102081;r=X:146993469-147032645
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["PKX1"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000183943;r=X:3522411-3631649
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["ECTD1", "ED1", "ED1-A1", "ED1-A2", "EDA-A1", "EDA-A2", "EDA1", "EDA2", "HED", "HED1", "ODT1", "STHAGX1", "XHED", "XLHED"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000158813;r=X:68835911-69259319
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["PMCA3", "PMCA3a", "SCAX1"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000067842;r=X:152783134-152848397
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["BRICD4", "CHM1L", "ChM1L", "TEM", "myodulin", "tendin"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000000005;r=X:99839799-99854882
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["CD132", "CIDX", "IL-2RG", "IMD4", "P64", "SCIDX", "SCIDX1"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000147168;r=X:70327254-70331958
            print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isin(["IRS-4", "PY160"]), :].shape[0] == 0)
                #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000133124;r=X:107975712-107979651;t=ENST00000372129
    else:
        print("ERROR: FALSE! This is not an actual error, just a warning. We have extracted genes close to UCP1 and Sarcopilin, not only according a p-value percentile of 1%, but also 0.5 and 5%. In these cases, the number of connectome genes that are NOT included in the gene_coords file is different than in the case of 1% as the set of genes is different. For example, for sarcopilin we lose 27 genes from the 5% connectome. I have NOT check if synonyms of these genes are in gene_coord. It should not be the case as no missing gene of the 1% had any synonym in gene_coords, but take this in mind if you use these sets.")

    print_text("Percentile " + str(p_value_percentile) + "% :  In the same vein, there are some genes in 'gene coords' that do not have gene name, only gene id. I understand that these genes have no valid gene name in hg19. Indeed, I have checked two of these genes and have NO description in ensembl. So they should not be included in our interest genes", header=4)
    print(pressure_coords.loc[pressure_coords["hgnc_symbol"].isna(),:]
    )
        #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000116883;r=1:36789335-36794822;t=ENST00000373137
        #https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000268538;r=22:30814212-30814469;t=ENST00000598426



    print_text("Percentile " + str(p_value_percentile) + "% :  define function to calculate the distance of each coding to gene to the closest BAT gene", header=3)
    #pressure_df = pressure_coords
    #pressure_dist_name=pressure_name
    #gene_id="ENSG00000003987"
    def distance_calc(pressure_df, pressure_dist_name, gene_id):
        
        #select the row of the selected gene id
        selected_gene_row = pressure_df[pressure_df["gene_id"] == gene_id]
        
        #extract its chromosome
        selected_chromosome = int(selected_gene_row["chromosome_name"])

        #select the genes related to the selective pressure
        selected_genes = pressure_df.loc[pressure_df[pressure_dist_name+"_status"] == "yes", :]

        #if the gene is a interest gene
        if all(selected_gene_row[pressure_dist_name+"_status"] == "yes"):
            #we need all() to avoid boolean problems with pandas
            
            #dist is zero
            distance=0
        else: 
            #else, status is "no" or nan so we have to calculate the
            #distance to the closest interest gene
            
            #check we have interest genes in the chromosome of the selected gene
            if(selected_chromosome in selected_genes["chromosome_name"].to_list()):
    
            	#from the interest genes, select those belonging to the same chromosome and then get their center
            	center_pressures = selected_genes \
                    .loc[
            	       selected_genes["chromosome_name"] == int(selected_chromosome), 
            	       "middle_point"]
            	    #int to avoid comparing two pandas series. We directly convert chromosome name
            	    #to integer
            	
            	#get the min distance between the center of the gene and those of the 
            	#interest genes in absolute value
            	distance=np.min(np.abs(float(selected_gene_row["middle_point"]) - center_pressures))
            else: #else, we do not have interest genes in the selected chromosome, so NA
            	distance = np.nan
        
        #return the distance and its gene_id as a tuple
        return tuple([gene_id, distance])



    print_text("Percentile " + str(p_value_percentile) + "% :  see function in action with two genes, one that is included in the group of interest genes and other not included", header=3)
    non_interest_gene = distance_calc( \
        pressure_df=pressure_coords, \
        pressure_dist_name=pressure_name, \
        gene_id=pressure_coords.loc[pressure_coords[pressure_name+"_status"]=="no", "gene_id"].iloc[0])
    interest_gene = distance_calc( \
        pressure_df=pressure_coords, \
        pressure_dist_name=pressure_name, \
        gene_id=pressure_coords.loc[pressure_coords[pressure_name+"_status"]=="yes", "gene_id"].iloc[0])
    print(non_interest_gene)
    print(interest_gene)



    print_text("Percentile " + str(p_value_percentile) + "% :  check that the non-interest gene has distance different from zero while the interest gene has a distance equal to zero", header=3)
    print(non_interest_gene[1] != 0)
    print(interest_gene[1] == 0)
        #the result of my function is tuple with two elements: gene id and distance to the closest interest gene. The second element [1] is the distance.



    print_text("Percentile " + str(p_value_percentile) + "% :  run the function", header=3)
    print_text("Percentile " + str(p_value_percentile) + "% :  use functools.partial to add a fixed parameter to the function. In this way, we can apply the function in parallel to each gene (third argument) while having the same value for the other two arguments, i.e., the dataset with the coordinates and the name of the selective pressure", header=4)
    from functools import partial
    distance_calc_fixed = partial(distance_calc, pressure_coords, pressure_name)
        #first you have the function then you have the arguments that will be fixed the first argument will always take the value of pressure_coords this the data.frame with the data for the selective pressure so it can be accessed by the function the second argument is for the name of the selective pressure, so we access columns by column name in pandas
        #gene id is not included because we are going to iterate across gene ids
            #https://stackoverflow.com/questions/25553919/passing-multiple-parameters-to-pool-map-function-in-python
    

    print_text("Percentile " + str(p_value_percentile) + "% :  open the pool with the selected number of cores", header=4)
    import multiprocessing as mp
    pool = mp.Pool(n_cores)
    print(pool)


    print_text("Percentile " + str(p_value_percentile) + "% :  run the function to calculate distance of each gene to the closest selective pressure gene. This takes the gene IDs as inputs and the function will output a list of tuples having each tuple the gene id and the distance to the closest interest gene", header=4)
    results_map = pool.map( \
        func=distance_calc_fixed, \
        iterable=pressure_coords["gene_id"].values)
            #Apply `func` to each element in `iterable`, collecting the results in a list that is returned
            #gene id is the only argument not fixed in "distance_calc_fixed", so we can iterate across it
                #https://stackoverflow.com/questions/64763867/parallel-processing-of-each-row-in-pandas-iteration
    pool.close()
    print("see first 10 gene ids")
    print(results_map[0:10])


    print_text("Percentile " + str(p_value_percentile) + "% :  convert the list with results to DF. The second column with the distance will be named used the pressure name and the p-value percentile used", header=4)
    results_df = pd.DataFrame( \
        results_map, 
        columns=[ \
            "gene_id", 
            pressure_name + "_distance_percentile_" + str(p_value_percentile)])
    print(results_df)
    

    print_text("Percentile " + str(p_value_percentile) + "% :  count how many NANs we have for the distance. NANs are caused because the gene has no interest genes within its chromosome, so no distance can be calculated", header=4)
    n_nans = sum(results_df[pressure_name + "_distance_percentile_" + str(p_value_percentile)].isna())
    print(n_nans)
    if(n_nans <= (results_df.shape[0]*0.03)):
        print(f"We have {n_nans} NANs for the distance to the closest {pressure_name.upper()} gene")
    elif (pressure_name=="smt") & (p_value_percentile==0.5) & (n_nans <= (results_df.shape[0]*0.1)):
        print("ERROR: FALSE! This is not an actual error, just a warning. We have many NANs for SMT and percentile 0.5%, around 1500 genes lost (only below 10%). Have this mind if you use this percentile.")
        print(f"We have {n_nans} NANs for the distance to the closest {pressure_name.upper()} gene")
    else:
        raise ValueError("ERROR: FALSE! WE HAVE TOO MUCH NANs FOR DISTANCE TO INTEREST GENES")


    print_text("Percentile " + str(p_value_percentile) + "% :  save the results", header=4)
    results_df.to_csv( \
        "./data/" + pressure_name + "_distance/" + pressure_name + "_data_percentile_" + str(p_value_percentile) + ".tsv", \
        sep='\t', \
        header=True, \
        index=False)


    print_text("Percentile " + str(p_value_percentile) + "% :  FINISHED", header=4)
