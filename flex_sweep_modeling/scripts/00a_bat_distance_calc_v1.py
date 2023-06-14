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



########################################################
######## COMPARE DIFFERENT MODELS IN PREDICTION ########
########################################################

#This script will perform model comparison for analysing flex-sweep probabilities 



##################
# import modules #
##################

import pandas as pd



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
        #print the standard error and stop
        raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM RUNNING COMMAND: " + complete_process.stderr)

#test it
print_text("check behaviour run_bash", header=1)
print_text("see working directory", header=2)
run_bash("pwd")
print_text("list files/folders there", header=2)
run_bash("ls")





####################
# prepare folder structure #
####################
print_text("prepare folder structure", header=1)
run_bash(" \
    mkdir \
        --parents \
        ./data/bat_distance; \
    ls -l ./data")





import pandas as pd

gene_coords = pd.read_table(\
	"/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/method_deep/data/search_diego/results/gene_number_cds_coords.txt", \
	sep="\t", \
	header=0, \
	low_memory=False)
gene_coords

gene_coords_no_duplicated = gene_coords[ \
	(-gene_coords["gene_id"].duplicated(keep="first"))]
gene_coords_no_duplicated

import numpy as np

print("Check that the subsetted has as many rows as unique gene IDs in the original file: ")
print(len(np.unique(gene_coords["gene_id"])) == gene_coords_no_duplicated.shape[0])


gene_coords_no_duplicated_subset = gene_coords_no_duplicated[[ 
	'chromosome_name', 
	'gene_id',
	'gene_start',
	'gene_end',
	'middle_point',
	'lower_end_window_50kb',
	'upper_end_window_50kb', 
	'lower_end_window_100kb',
	'upper_end_window_100kb', 
	'lower_end_window_200kb',
	'upper_end_window_200kb',
	'lower_end_window_500kb',
	'upper_end_window_500kb',
	'lower_end_window_1000kb',
	'upper_end_window_1000kb']]
gene_coords_no_duplicated_subset


gene_coords_no_duplicated_subset = gene_coords_no_duplicated_subset.reset_index(drop=True)
gene_coords_no_duplicated_subset
    #https://stackoverflow.com/questions/20490274/how-to-reset-index-in-a-pandas-dataframe





bat_coords = gene_coords_no_duplicated.copy()

#missing genes
bat_relationship.loc[~bat_relationship["Genes"].isin(bat_coords["hgnc_symbol"]),:]
	#faltan 2 genes

bat_coords.loc[bat_coords["hgnc_symbol"].isin(["CD132", "CIDX", "IL-2RG", "IMD4", "P64", "SCIDX", "SCIDX1"]), :]

bat_coords.loc[bat_coords["hgnc_symbol"].isin(["NRU", "P2P", "P2Y4", "UNR"]), :]
	#the two missing genes have no synonmious in the dataset
		#https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000186912;r=X:69478016-69479654;t=ENST00000374519
		#https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000147168;r=X:70327254-70331958


#load the connectome with UCP1 as core gene
ucp1_conn = pd.read_csv("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/data/human_connectome/UCP1.txt", \
	sep="\t", \
	header=0, \
	low_memory=False)



#load the information about BAT relationships
bat_relationship = pd.read_csv("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/other_projects/human_genome_connectome/bat_connectome/results/connectome_results/tables/appendix_S1_ordered.csv", \
	sep=",", \
	header=0, \
	low_memory=False)
	#CHECK THIS FILE

#we can filter by p-vlaue percentage
sum(ucp1_conn.loc[ucp1_conn.iloc[:, 3]< 0.25, "Target"].isin(bat_relationship["Genes"]))

#10% is 0.1, 1% (used in the connectome paper) is 0.01
selected_ucp_connectome_genes = ucp1_conn.loc[ucp1_conn.iloc[:, 3]< 0.01, "Target"]

bat_coords["bat_status"] = ["yes" if gene_id in selected_ucp_connectome_genes.to_list() else "no" for gene_id in bat_coords["hgnc_symbol"]]
	#CHECK THIS
#bat_coords["bat_status"] = ["yes" if gene_id in bat_relationship.loc[bat_relationship["BAT relationship"].isin([0, 1]), "Genes"].to_list() else "no" for gene_id in bat_coords["hgnc_symbol"]]
	#CHECK THIS


	#you want all genes in BAT relationships, known and unknown BAT genes

bat_coords.columns



#pressure_df = bat_coords
#pressure_dist_name="bat"
#gene_id="ENSG00000003987"
def distance_calc(pressure_df, pressure_dist_name, gene_id):
    
    #select the row of the selected gene id
    selected_gene_row = pressure_df[pressure_df["gene_id"] == gene_id]
        
    #if the gene is a interest gene
    if all(selected_gene_row[pressure_dist_name+"_status"] == "yes"):
        #we need all() to avoid boolean problems with pandas
        
        #dist is zero
        distance=0
    else: 
        #else, status is "no" or nan so we have to calculate the
        #distance to the closest interest gene
        
        #select the genes related to the selective pressure
        selected_genes = pressure_df[pressure_df[pressure_dist_name+"_status"] == "yes"]
        
        #check we have interest genes in the chromosome of the selected gene
        selected_chromosome = int(selected_gene_row["chromosome_name"])

        if(selected_chromosome in selected_genes["chromosome_name"].to_list()):

        	#from the interest genes, select those belonging to the same chromosome
        	#and the get their center
        	center_pressures = selected_genes.loc[
        	    selected_genes["chromosome_name"] == selected_chromosome, 
        	    "middle_point"]
        	    #int to avoid comparing two pandas series. We directly convert chromosome name
        	    #to integer
        	
        	#get the min distance between the center of the gene and those of the 
        	#interest genes in absolute value
        	distance=np.min(np.abs(float(selected_gene_row["middle_point"]) - center_pressures))
        else:
        	distance = np.nan
    
    #return the distance and its gene_id as a tuple
    return tuple([gene_id, distance])


import multiprocessing as mp
from functools import partial

#use partial to add a fixed parameter to the function
distance_calc_fixed = partial(distance_calc, 
                              bat_coords, 
                              "bat")
    #first you have the function
    #then you have the arguments that will be fixed
    #the first argument will always take the value of vip_distance_calc
    #this the data.frame with the data for the selective pressure
    #so it can be accessed by the function
    #the second argument is for the name of the selective pressure, so we access
    #columns by column name in pandas
    #gene id is not included because we are going to iterate across gene ids
    #https://stackoverflow.com/questions/25553919/passing-multiple-parameters-to-pool-map-function-in-python

#open the pool with 5 cores
pool = mp.Pool(5)

#run the function to calculate distance of each gene to the closest
#selective pressure gene. This takes the gene IDs as inputs
#and the function will calculate the distance for each gene id
#https://stackoverflow.com/questions/64763867/parallel-processing-of-each-row-in-pandas-iteration
bat_distance_first_map = pool.map(distance_calc_fixed, 
                                  bat_coords["gene_id"].values)

#close the pool
pool.close()

bat_distance_first_map_df = pd.DataFrame(bat_distance_first_map, 
                                         columns=["gene_id", 
                                                  "bat_distance"])
bat_distance_first_map_df


bat_distance_first_map_df.to_csv("./data/bat_distance/bat_distance.tsv",
	sep='\t', 
	header=True, 
	index=False)