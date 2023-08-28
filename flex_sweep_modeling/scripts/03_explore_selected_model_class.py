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



##################################################
######## EXPLORE THE SELECTED MODEL CLASS ########
##################################################

#This script will perform the modeling of flex-sweep probabilities using the selected model class.



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
        print("\n#######################################")
        print(text)
        print("#######################################")
    elif header==3:
        print("\n###### " + text + " ######")
    elif header==4:
        print("\n## " + text + " ##")
    elif header==5:
        print("\n# " + text + " #")
print_text("checking function to print nicely: header 1", header=1)
print_text("checking function to print nicely: header 2", header=2)
print_text("checking function to print nicely: header 3", header=3)
print_text("checking function to print nicely: header 4", header=4)
print_text("checking function to print nicely: header 5", header=5)



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





##########################
# ensure reproducibility #
##########################
print_text("set seeds for reproducibility", header=1)
#No need to set random_state for any model after this
#https://stackoverflow.com/a/52897216/12772630


print_text("Seed value", header=2)
seed_value=0
print(seed_value)


print_text("Set the `PYTHONHASHSEED` environment variable at a fixed value", header=2)
import os
os.environ["PYTHONHASHSEED"]=str(seed_value)
print(os.environ["PYTHONHASHSEED"])


print_text("Set the `python` built-in pseudo-random generator at a fixed value", header=2)
import random
random.seed(seed_value)


print_text("Set the `numpy` pseudo-random generator at a fixed value", header=2)
import numpy as np
np.random.seed(seed_value)


print_text("Set the `tensorflow` pseudo-random generator at a fixed value", header=2)
import tensorflow as tf
tf.random.set_seed(seed_value)


print_text("Configure a new global `tensorflow` session", header=2)
session_conf = tf.compat.v1.ConfigProto( \
    intra_op_parallelism_threads=1,  \
    inter_op_parallelism_threads=1)
sess = tf.compat.v1.Session( \
    graph=tf.compat.v1.get_default_graph(), \
    config=session_conf)
tf.compat.v1.keras.backend.set_session(sess)





############################
# prepare folder structure #
############################
print_text("prepare folder structure", header=1)
run_bash(" \
    mkdir \
        --parents \
        ./results/selected_model_class; \
    ls -l ./results")





####################
# data preparation #
####################
print_text("data preparation", header=1)
print_text("Set the window size selected for those variables that are calculated within windows centered around genes", header=2)
gene_window_size = "1000kb"
print(f"The selected window size is {gene_window_size}")




print_text("load and clean input data", header=2)
print_text("load flex-sweep probability and most predictors into pandas", header=3)
import pandas as pd
final_data_yoruba_raw = pd.read_csv( \
    "./data/flex_sweep_closest_window_center.txt.gz", \
    sep=",", \
    header=0, \
    low_memory=False, \
    compression="gzip")
print(final_data_yoruba_raw)



print_text("load BAT and SMT data, only 1% percentile for now", header=3)
p_value_percentile = 1
bat_distance = pd.read_csv( \
    "./data/bat_distance/bat_data_percentile_" + str(p_value_percentile) + ".tsv",
    sep='\t', 
    header=0, 
    low_memory=False)
smt_distance = pd.read_csv( \
    "./data/smt_distance/smt_data_percentile_" + str(p_value_percentile) + ".tsv",
    sep='\t', 
    header=0, 
    low_memory=False)
print(bat_distance)
print(smt_distance)


print_text("merge all data.frames", header=3)
print_text("make list with all DFs", header=4)
list_dataframes = [final_data_yoruba_raw, bat_distance, smt_distance]
print(list_dataframes)


print_text("merge all of them with reduce", header=4)
from functools import reduce
final_data_yoruba = reduce( \
    lambda x, y: \
        pd.merge( \
            left=x, \
            right=y, \
            how="left", \
            on="gene_id"), \
    list_dataframes)
        #reduce applies a function of two arguments cumulatively to the items of a sequence, from left to right, so as to reduce the sequence to a single value. For example, reduce(lambda x, y: x+y, [1, 2, 3, 4, 5]) calculates ((((1+2)+3)+4)+5).
        #therefore, in the first iteration, we merge the first DF of the list (left) with the second (right). The resulting DF (left) is then merged with the third DF (right), and so on....
        #the first DF is going to be the original dataset with flex-sweep probability and most of predictors, therefore using left we ensure we use only those rows with gene_id in the original DF.
            #in any case, we should not have NA for distance to BAT and SMT as these have been calculated using the original gene coordinate file, which was in turn used to calculate many of our predictors.
print(final_data_yoruba)


print_text("count the number of NANs in distance BAT and SMT genes", header=4)
n_nans_bat = sum(final_data_yoruba["bat_distance_percentile_1"].isna())
n_nans_smt = sum(final_data_yoruba["smt_distance_percentile_1"].isna())
if(n_nans_bat <= (final_data_yoruba.shape[0]*0.03)) & (n_nans_smt <= (final_data_yoruba.shape[0]*0.03)):
    print(f"We have {n_nans_bat} NANs for the distance to the closest BAT gene")
    print(f"We have {n_nans_smt} NANs for the distance to the closest SMT gene")
else:
    raise ValueError("ERROR: FALSE! WE HAVE MORE THAN 3% OF NANs FOR DISANCE TO INTEREST GENES")


print_text("remove NANs", header=4)
final_data_yoruba = final_data_yoruba.dropna()
print(final_data_yoruba)




'''
print_text("clean predicted class", header=3)
if FALSE:
    #if you want to do classification, you have to calculate the number of sweeps based on the probability and then convert "predicted_class" to 0-1 integer
    final_data_yoruba["predicted_class"] = ["neutral" if prob < 0.5 else "sweep" for prob in final_data_yoruba["prob(sweep)"]]
    decode_response = {"predicted_class": {"neutral": 0, "sweep": 1}}
    final_data_yoruba = final_data_yoruba.replace(decode_response)
       #https://pbpython.com/categorical-encoding.html
    final_data_yoruba["predicted_class"]=final_data_yoruba["predicted_class"].astype("int")
        #https://stackoverflow.com/questions/41925157/logisticregression-unknown-label-type-continuous-using-sklearn-in-python
    #discrete is not bad! but regression seems to work better and we avoid issues related to dicomitize a continuous variable and using an arbitrary threshold.
        #results RandomForestRegressor in CV - 5 fold
            #accuracy: array([0.80921584, 0.81601294, 0.82531338, 0.82248281, 0.82207845]),
            #precision: array([0.78202677, 0.80529301, 0.82278481, 0.77163904, 0.81785714]),
            #recall: array([0.53324641, 0.54755784, 0.57667934, 0.57084469, 0.57537688]),
            #f1: array([0.63410853, 0.65187452, 0.6780924 , 0.65622553, 0.67551622])}
        #results RandomForestRegressor in the test set
            #accuracy = 0.8479948253557568
            #precision = 0.8280346820809249
            #recall = 0.6201298701298701
            #f1 = 0.7091584158415841
        #explanation about metrics
            #accuracy does not work very well for imbalanced datasets like ours, we have 1/3 of sweeps compared to non-sweeps.
            #precision minimice false positives
            #recall minimize false negative
            #f1 is the harmomic average of precision and recall, so if any of the two parameters is low, the average will be much lower than it was an aritmetic or geometric average. 
                #this is useful if you want to minimize both false negatives and positives
                #https://towardsdatascience.com/essential-things-you-need-to-know-about-f1-score-dbd973bf1a3
                #https://stephenallwright.com/good-f1-score/
        #Interpretation
            #Overall we are OK: f1 value around 0.6 is OK but not good (pasable)! 
                #https://stephenallwright.com/good-f1-score/
            #we do not have a higher f1 because recall is lower than precision, i.e., we are not minimizing false negatives as much as false positives.
            #in other words, we have more sweeps that are classified as non-sweeps than neutral classified as sweep.
                #you can also see this with the continuous variable because there is more error in the prediction of strong sweeps candidates (log probability equal to zero).
        #setting the threshold at 0.75 make recall worse. Also not sure if using the prediciton improvement is a good a idea to select the threshold.
            #https://stats.stackexchange.com/a/94044
        #this is not too bad, but regression gets better results and we avoid the use of arbitrary thresholds for classification
'''



print_text("clean the data", header=3)
print_text("exclude some columns we are not interested in", header=4)
columns_to_exclude = [ \
    "gene_id", \
    "predicted_class", \
    "number_thermogenic_1000kb", \
    "number_vips_1000kb"]
final_data_yoruba_subset = final_data_yoruba[[column for column in final_data_yoruba.columns if column not in columns_to_exclude]]
print(final_data_yoruba_subset)
print(f"Columns excluded: {columns_to_exclude}")


print_text("put the response as first column", header=4)
response_name = "prob(sweep)"
first_column = final_data_yoruba_subset.pop(response_name)
final_data_yoruba_subset.insert(0, response_name, first_column)
    #https://www.geeksforgeeks.org/how-to-move-a-column-to-first-position-in-pandas-dataframe/
print(final_data_yoruba_subset)


print_text("make deep copy of the data to do further operations", header=4)
modeling_data = final_data_yoruba_subset.copy(deep=True)
    #deep=True: 
        #a new object will be created with a copy of the calling object's data and indices. Modifications to the data or indices of the copy will not be reflected in the original object (see notes below).
print(modeling_data)


print_text("Apply log transformation to the target variable using the original DF as source", header=4)
import numpy as np
modeling_data["prob(sweep)"] = final_data_yoruba_subset["prob(sweep)"].apply(lambda x: np.log(x))
    #It is should be ok to apply the log before splitting the dataset. There is a problem if you use a transformation that requires learn something from the rest of the data. For example, if you scale the whole dataset, you are using the mean and sd of the whole dataset, influencing data that will be used for test. In other words, there is room for a data leak. In this case, however, log(1.5) is always 0.4, independently of the rest of the data, so I think no data leak is possible. You could do a pipeline with log but it is a little bit more complicated (see [link](https://stats.stackexchange.com/questions/402470/how-can-i-use-scaling-and-log-transforming-together)), so we leave it for now.
        #Indeed I have found people in stack exchange saying this: However, yours (i.e. np.log1p) is a simple transformation that doesn't use any learnable parameters, and it won't matter if you do it before or after the split. It's like dividing a feature by 1000. 
            #https://stats.stackexchange.com/a/456056
    #From all these follow that if you use other transformations like preprocessing.PowerTransformer or QuantileTransformer ([link](https://yashowardhanshinde.medium.com/what-is-skewness-in-data-how-to-fix-skewed-data-in-python-a792e98c0fa6)), it is possible to have data leaks, so be careful.
    #In previous versions I was not using scaling or log transform for deep learning, because I assumed that the DNNs can deal with that, but maybe that was too much and in any case, we are going to use here also more simpler models that can be helped by scaling
    #Update: If I apply the log transformation within the pipeline I get a much lower R2 both in the training and test datasets! Not sure what is going on, but given this transformation does not summarize anything from the whole dataset, I can use it before splitting in training and evaluation. If this transformation was helping the training model to learn from the test set, the R2 in the test would be higher, but we have the opposite scenario.
    #in case you want to apply the log transformation within the pipeline, but this make the model MUCH WORSE compared to just apply the log to the original DF. Do not know why.
        #you can do it with func=np.log and inverse_func=np.exp in TransformedTargetRegressor


print_text("Explanations about scaling in pipeline", header=4)
#We could apply the preprocessing to the initial dataset and then split into training and test but this is problematic. As previously explained, the scaling is done using the mean and sd of the sample, so if you do it in the whole dataset, the part will be used for test will be also influenced but training data. If you included your test data in the scaling, that means that your new data is treated differently from the training set, which defeats the purpose of the training set. In practice, this is unlikely to have a large impact, as computing mean and standard deviation is relatively stable on well-behaved datasets. However, I recommend to adhere to best practices, and split off the test set before doing any processing ([more info](https://amueller.github.io/aml/01-ml-workflow/03-preprocessing.html)).

#In other words, when you fit the standard scaler on the whole dataset, information from the test set is used to normalize the training set. This is a common case of "data leakage", which means that information from the test set is used while training the model. This often results in overestimates of the model's performance because the model can access to information of the test set during training, so the test set is no longer new data. ([link](https://stackoverflow.com/questions/63037248/is-it-correct-to-use-a-single-standardscaler-before-splitting-data?noredirect=1&lq=1)). 

#We can easily avoid this problem by doing CV using a pipeline. This pipeline includes a transformer and the regressor, so when fitting it does the scaling using the training data and then fit the model on that training data, there is not access to test data. Test data is used only when predicting, once we have fit our model based on the scaled training data.
    #Pipeline in sklearn enables you to set up trainable blocks that contain both your models and transformations in order, so that when fitted, it only uses training data.

#standard.scaler() is similar to preprocessing.scale(). In both cases, scaling means standardizing by removing the mean and scaling to unit variance. The standard score of a sample x is calculated as: z = (x - u) / s, where u is the mean of the training samples or zero if with_mean=False, and s is the standard deviation of the training samples or one if with_std=False.

#Instead of using preprocessing.scale, we will use preprocessing.StandardScaler(), which is a transformer. You can open an instance of this transformer, call the method fit to learn the mean and sd from the data and then apply the trasfomation
from sklearn import preprocessing
dummy_sample = np.array([1,1,2,2])

scaler = preprocessing.StandardScaler()
scaler.fit(dummy_sample.reshape(-1,1))
print(f"See result with preprocessing.StandardScaler")
print(scaler.transform(dummy_sample.reshape(-1,1)))
print("We get the exactly the same if we use preprocessing.scale")
print(preprocessing.scale(dummy_sample))

#The great adventage of StandardScaler() is that we can call an instance and use it within a pipeline that will be feed into a gridsearch. Internally, gridsearch scale the training dataset using the training mean-sd, fit the data and then predict in evaluation dataset but after scaling that evaluation dataset with the mean sd of the training set. In the next iteration (next training-eval set) the processes is repeated ([link](https://stackoverflow.com/questions/51459406/how-to-apply-standardscaler-in-pipeline-in-scikit-learn-sklearn)).

#In previous versions I was not using scaling or log transform for deep learning, because I assumed that the DNNs can deal with that, but maybe that was too much and in any case, we are going to use here also more simpler models that can be helped by scaling


print_text("Explanations about the non-independence of genes", header=4)
#This is not a problem for the factors we are just controlling for, the problem comes when you want to test for causality for a given factor like BAT distance:
#Imagine you have a greenhouse experiment
    #you have a greenhouse experiment where you apply three levels of watering and you measure growth.
    #there is heterogeneity in conditions across the greenhouse, so some areas have more light and temperature than others. Therefore, two plants under the same water treatment can have different growth just because they are in a different part of the greenhouse and, consequently, they are exposed to different levels of light. You have to control for the light in order to have power to detect differences between the water treatments.
    #The equivalence in our case would be
        #plant = gene
        #water level = VIP/BAT/SMT
            #the problem is HERE, you have not applied a treatment to each gene
        #growth = positive selection
        #position in greenhouse = position in the chromosome
        #factors varying across the greenhouse (e.g., light) = factors varying across the chromosome (e.g., recombination or GC-content)
    #The position in the greenhouse impacts growth because abiotic factors co-varies with position. Similarly, the position in the chromosome impacts positive selection because this physical position covaries with genomic factors influencing the detection/occurrence of positive selection.
    #Therefore, if you control for the ultimate cases, for these factors, you indeed do not need to control for the chromosomal position, because you are considering the actual drivers of positive selection instead of using just a proxy.
    #This means that we could safely use genes that are close as long as we have enough information about their genomic characteristics. If two genes are close and share genomic features, the model will predict the same probability of selection. In contrast, if you genes are close but differ in a important factor like recombination (e.g., one gene is in a recombination hotspot), the model will predict different probability of selection. He has learned that two close regions with very different recombination rates will still have different probability of positive selection, thus it is able to correctly predict selection even for close genes, as it consider the actual drivers.
#There is a problem for the selective pressure that we want to test
    #In the greenhouse experiment, we randomly selected the plants under different water treatments, while we did not selected genes to be BAT or not.
    #for the genomic factors like recombination etc should be ok, because we have evidence of their importance (like for light) and we are just including them to control if there is an impact from them.
        #if two genes are close to the same VIP, I want to increase the probability of selection for both genes. I do not care if this is an independent point. We know that VIPs are enriched in positive selection (previous evidence) and I want to control for that effect in order to detect the influence of other selective pressures. If it is relevant in my dataset, the model will learn that close to VIPs genes tend to have more selection and will consider when predicting the probability of selection of a given gene, so we can control for this if, for example, a BAT gene is also a VIP, expecting for this gene more selection just because it is a VIP.
    #but in the case of BAT, we could have two close to the same BAT gene. If this BAT gene was under selection, the other two genes are going to belong to the same sweep, so they can have high positive selection but only because 1 sweep, 1 independent point.
#THE SOLUTION
    #Because of this, it is very important that when we want to check the impact of a selective pressure, we perform the enrichment test of David.
    #We select control genes based on the probability of selection (obtained from our model) but also at enough distance from the BAT genes.
    #In this way, we are not counting twice the same sweep
            #maybe some BAT genes are physically close to each other?





print_text("prepare CV scheme", header=1)
#It is VERY important that we avoid overfitting given the use we are going to make of the models.
    #In the future, we will use the final model to obtain a probability of selection considering genomic factors and then use it to select genes with the same expected probability of selection (according to these factors) than our genes of interest. The idea is that the interest genes should have the same probability of selection based on genomic features (predicted probability), but if they are target of a selective pressure, their observed probability of selection should be higher. If predicted and observed are exactly the same, there is no room for enrichment of selection in interest genes after controling for confounding factors, and this would be an methodological artifact due to a model that just fit the observed data without any generalization. This can be a problem with algorithms like RF or in deep learning. In other words, more overfitting, less power to detect the impact of selective pressures.

#Held out part of the dataset in each CV round (i.e., test set) 
    #This part will not be used for parameter optimization, but for the final validation after the final model has been optimized. In this way, we avoid potential overfiting in the evaluation sets ([see link](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation)). If you use train in a set of the data and then evaluate in other, you can see if the model trained is not overfitting to the training data and is flexible enough to predict the evaluation data. The problem is that in parameter optimization, we select the best parameters based on the evaluation metrics in the evaluation sets, thus we could get a model that fit too much the evaluation dataset, loosing generalization and thus making the evaluation metrics no longer metrics of generalization. To avoid this, we leave out a set of the data for final evaluation, i.e., the test set. This set will be NOT used in parameter optimization. Once we have selected the best parameters to train a model in the training dataset and predict well in the evaluation datasets of the CV, we use these parameter to create a final model, fit to the whole training data and then predict in the final evaluation dataset, which was not used for anything before. If the model works well, it means it is generalizable and it is not overfitting the data, so, in our case, we can say that it is explaining the variance in selection that it is really explained by the genomic factors and the rest would be variance that could be explained by selective pressures. If there is overfitting, the model fit too much the data, there is not non-explained variance.



print_text("split train/test set", header=2)
from sklearn.model_selection import train_test_split
train, test = train_test_split(
    modeling_data,
    test_size=0.20, #train size is automatically calculated using this value 
    #random_state=54, #we have already set the seed, so this is not required
    shuffle=True)
        #this function uses internally ShuffleSplit
        #so it can randomly split a dataset in training and evaluation
        #but instead of getting an interable (like in ShuffleSplit), 
        #you directly get the datasets splitted
            #https://stackoverflow.com/questions/66757902/differnce-between-train-test-split-and-stratifiedshufflesplit
print("Do we have the correct shapes")
print(train.shape)
print(test.shape)
print((train.shape[1] == modeling_data.shape[1]) & (test.shape[1] == modeling_data.shape[1]))
print(train.shape[0]+test.shape[0] == modeling_data.shape[0])
    #we are using 20% of test instead of 10% like in the model class comparison. 
        #In that case, we needed to repeated HPs optimization several times under pre-defined train/eval/test sets so we have exposed different model classes to the same data partitions randomly selected. As we increase the number of splits (repetitions), we increase the number of folds as we do 10 folds, each time 9 folds for training/eval and 1 for test and you repeat 10 times. The test set is only 10%. This is basically a nested CV where we need to repeat several times in order to get a good sense of the model performance of different classes.
        #Now, we only have 1 model class that we are going to tune, so we can tune HPs in training/eval set with CV and then check the performance in the held-out test. As we get a larger test set, we should get a more robust measurement about the generalization ability of the model.


print_text("set the CV scheme for the training set", header=2)
from sklearn.model_selection import KFold
cv_scheme = KFold( \
    n_splits=5,  \
    shuffle=True)
    #random_state is not required as we have already ensured reproducibility by setting the seeds of python, numpy and tensorflow
print(cv_scheme)
#[(train, test) for fold, (train, test) in enumerate(cv_scheme.split(train)) if fold==0]

#we will use k=5 for now for speed
    #it has been seen in many datasets that 10 folds works relatively well to estimate model performance ([link](https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/#:~:text=configure%20the%20procedure.-,Sensitivity%20Analysis%20for%20k,evaluate%20models%20is%20k%3D10.)). Although there is some debate about it. We will stick to 10 folds for now
    #This is very well explained here (https://stackoverflow.com/questions/45969390/difference-between-stratifiedkfold-and-stratifiedshufflesplit-in-sklearn).

#In the future you could compare different K values and select the best 
    #you could do a sensitivity analysis
        #use default parameters
        #then run several CVs with different number of folds
        #calculate the average R2 and min-max in the folds
        #select the number of folds that give the greatest R2
    #This will automatically select the best number of folds, thus also selecting the selecting the size of the training/evaluation sets. For example, 10 outer folds means that 90% will be used for training and 10% for evaluation.
    #link
        #https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/

#Respect to the repetition because of the stochastic nature of the machine learning models
    #If a large number of machine learning algorithms and algorithm configurations are being evaluated systematically on a predictive modeling task, it can be a good idea to fix the random seed of the evaluation procedure. Any value will do.
    #The idea is that each candidate solution (each algorithm or configuration) will be evaluated in an identical manner. This ensures an apples-to-apples comparison. It also allows for the use of paired statistical hypothesis tests later, if needed, to check if differences between algorithms are statistically significant.
    #This what we are gonna do, as we just want to select a model class. We can consider stochasiticity when working with the selected class, running different models with different seeds and get the average prediction. We will likely have several good models, not only one.
        #https://machinelearningmastery.com/different-results-each-time-in-machine-learning/




print_text("notes about XGBoost", header=1)
    #general notes about gradient boosting
        #The idea of boosting came out of the idea of whether a weak learner can be modified to become better. Hypothesis boosting was the idea of filtering observations, leaving those observations that the weak learner can handle and focusing on developing new weak learns to handle the remaining difficult observations.
        #A first implementation of this was AdaBoost
            #AdaBoost works by weighting the observations, putting more weight on difficult to classify instances and less on those already handled well. New weak learners are added sequentially that focus their training on the more difficult patterns.
            #This means that samples that are difficult to classify receive increasing larger weights until the algorithm identifies a model that correctly classifies these samples
            #Predictions are made by majority vote of the weak learners’ predictions, weighted by their individual accuracy.
        #This was further developed into Gradient Boosting Machines
            #the objective is to minimize the loss of the model by adding weak learners using a gradient descent like procedure.
                #The job of the algorithm is to find a set of internal model parameters (predictors) that perform well against some performance measure such as logarithmic loss or mean squared error. 
                #Optimization is a type of searching process and you can think of this search as learning. The optimization algorithm is called “gradient descent“, where “gradient” refers to the calculation of an error gradient or slope of error and “descent” refers to the moving down along that slope towards some minimum level of error.
            #This class of algorithms were described as a stage-wise additive model. This is because one new weak learner is added at a time and existing weak learners in the model are frozen and left unchanged.
                #Note that this stagewise strategy is different from stepwise approaches that readjust previously entered terms when new ones are added
            #The generalization allowed arbitrary differentiable loss functions to be used, expanding the technique beyond binary classification problems to support regression, multi-class classification and more.
        #How gradient boost works
            #1. A loss function to be optimized, i.e., you calculate the loss based on the predictions of your first tree.
                #many standard loss functions are supported and you can define your own.
                #For example, regression may use a squared error and classification may use logarithmic loss.
                #A benefit of the gradient boosting framework is that a new boosting algorithm does not have to be derived for each loss function that may want to be used, instead, it is a generic enough framework that any differentiable loss function can be used.
            #2. A weak learner to make predictions.
                #Decision trees are used as the weak learner in gradient boosting.
                #We use regression trees. They output real values for splits, thus their output can be added together. This allows to add the output of subsequent models and “correct” the residuals in the predictions.
                #Trees are constructed in a greedy manner, choosing the best split points based on purity scores like Gini or to minimize the loss.
                #Initially, such as in the case of AdaBoost, very short decision trees were used that only had a single split, called a decision stump. Larger trees can be used generally with 4-to-8 levels.
                #It is common to constrain the weak learners in specific ways, such as a maximum number of layers, nodes, splits or leaf nodes. This is to ensure that the learners remain weak, but can still be constructed in a greedy manner. 
                    #THIS SEEMS IMPORTANT, as it seems to reduce overfitting when adding many learners.
            #3. An additive model to add weak learners to minimize the loss function.
                #Trees are added one at a time, and existing trees in the model are not changed.
                #A gradient descent procedure is used to minimize the loss when adding trees.
                #Traditionally, gradient descent is used to minimize a set of parameters, such as the coefficients in a regression equation or weights in a neural network. After calculating error or loss, the weights are updated to minimize that error.
                #Instead of parameters, we have weak learner sub-models or more specifically decision trees. After calculating the loss, to perform the gradient descent procedure, we must add a tree to the model that reduces the loss (i.e. follow the gradient). We do this by parameterizing the tree, then modify the parameters of the tree and move in the right direction by reducing the residual loss.
                    #Basically we train the new tree in order to fit the loss of the current set of trees.
                #The output for the new tree is then added to the output of the existing sequence of trees in an effort to correct or improve the final output of the model.
                    #you combine the existing ensemble of tree with the new tree.
                #The new ensemble is used to predict and the new loss is then fitted by a new tree.
                #A fixed number of trees are added and then we stop. We can also stop training stops once loss reaches an acceptable level or no longer improves on an external validation dataset.
        #summary
            #steps
                #you train a decision tree
                #calculate the loss, 
                #train a new tree using the LOSS, the error.
                #add the new trained tree to the ensemble with the previous tree, i.e., you are BOOSTING the previous tree with a new one. Each time you add a new tree is a boosting round.
                #make predictions again
                #calculate again the loss
                #use this loss to train a new tree
                #and so on...
                    #https://www.youtube.com/watch?v=yw-E__nDkKU&ab_channel=EmmaDing
            #I think this approach makes process to focus on the parts of the data that are more difficult to learn as you try to fit the residuals of the current model when adding a new tree.
        #Improvements to Basic Gradient Boosting
            #Gradient boosting is a greedy algorithm and can overfit a training dataset quickly.
            #It can benefit from regularization methods that penalize various parts of the algorithm and generally improve the performance of the algorithm by reducing overfitting.
                #1.Tree Constraints
                    #It is important that the weak learners have skill but remain weak. There are a number of ways that the trees can be constrained.
                    #A good general heuristic is that the more constrained tree creation is, the more trees you will need in the model, and the reverse, where less constrained individual trees, the fewer trees that will be required.
                    #Below are some constraints that can be imposed on the construction of decision trees:
                        #Number of trees, generally adding more trees to the model can be very slow to overfit. The advice is to keep adding trees until no further improvement is observed.
                            #This means there is a low-risk of overfitting by adding trees, so just add them until there is no improvement.
                        #Tree depth, deeper trees are more complex trees and shorter trees are preferred. Generally, better results are seen with 4-8 levels.
                        #Number of nodes or number of leaves, like depth, this can constrain the size of the tree, but is not constrained to a symmetrical structure if other constraints are used.
                        #Number of observations per split imposes a minimum constraint on the amount of training data at a training node before a split can be considered
                        #Minimum improvement to loss is a constraint on the improvement of any split added to a tree.
                #2.Shrinkage - Weighted Updates
                    #The predictions of each tree are added together sequentially.
                    #The contribution of each tree to this sum (previous ensemble) can be weighted to slow down the learning by the algorithm. This weighting is called a shrinkage or a learning rate.
                        #Similar to a learning rate in stochastic optimization, shrinkage reduces the influence of each individual tree and leaves space for future trees to improve the model.
                        #each individual tree has less impact on the global output because you limit how much change an individual tree can make in a single step.
                    #In other words
                        #The learning rate controls the amount of contribution that each model has on the ensemble prediction.
                    #The effect is that learning is slowed down, in turn require more trees to be added to the model, in turn taking longer to train, providing a configuration trade-off between the number of trees and learning rate.
                    #It is common to have small values in the range of 0.1 to 0.3, as well as values less than 0.1.
                #3.Random sampling
                    #A big insight into bagging ensembles and random forest was allowing trees to be greedily created from subsamples of the training dataset, so trees and errors are less correlated.
                    #This same benefit can be used to reduce the correlation between the trees in the sequence in gradient boosting models.
                    #This variation of boosting is called stochastic gradient boosting.
                    #at each iteration a subsample of the training data is drawn at random (without replacement) from the full training dataset. The randomly selected subsample is then used, instead of the full sample, to fit the base learner
                        #I understand the base learner is the weak learner that is being used for boosting.
                    #A few variants of stochastic boosting that can be used:
                        #Subsample rows before creating each tree.
                        #Subsample columns before creating each tree
                        #Subsample columns before considering each split.
                    #Generally, aggressive sub-sampling such as selecting only 50% of the data has shown to be beneficial.
                    #According to user feedback, using column sub-sampling prevents over-fitting even more so than the traditional row sub-sampling
                        #XGBoost: A Scalable Tree Boosting System, 2016
                #4.Penalized Learning
                    #Additional constraints can be imposed on the parameterized trees in addition to their structure.
                    #Classical decision trees like CART are not used as weak learners, instead a modified form called a regression tree is used that has numeric values in the leaf nodes (also called terminal nodes). The values in the leaves of the trees can be called weights in some literature.
                        #In classical decision trees you have a condition like predictor X1 is higher than 30. 
                        #I understand that in regression trees we have just a number of weight in the node instead of a condition.
                    #As such, the leaf weight values of the trees can be regularized using popular regularization functions, such as:
                        #L1 regularization of weights.
                        #L2 regularization of weights.
                    #The additional regularization term helps to smooth the final learnt weights to avoid over-fitting. Intuitively, the regularized objective will tend to select a model employing simple and predictive functions.
                        #These are the regularization methods of Ridge and Lasso.
        #links
            #https://machinelearningmastery.com/gentle-introduction-gradient-boosting-algorithm-machine-learning/
            #https://youtu.be/yw-E__nDkKU
    #general notes about Extreme Gradient Boosting
        #Gradient boosting refers to a class of ensemble machine learning algorithms that can be used for classification or regression predictive modeling problems.
        #Ensembles are constructed from decision tree models. Trees are added one at a time to the ensemble and fit to correct the prediction errors made by prior models. This is a type of ensemble machine learning model referred to as boosting.
        #Models are fit using any arbitrary differentiable loss function and gradient descent optimization algorithm. This gives the technique its name, “gradient boosting,” as the loss gradient is minimized as the model is fit, much like a neural network.
        #Extreme Gradient Boosting, or XGBoost for short, is an efficient open-source implementation of the gradient boosting algorithm. As such, XGBoost is an algorithm, an open-source project, and a Python library.
        #It is designed to be both computationally efficient (e.g. fast to execute) and highly effective, perhaps more effective than other open-source implementations.
        #The two main reasons to use XGBoost are execution speed and model performance.
        #XGBoost dominates structured or tabular datasets on classification and regression predictive modeling problems. The evidence is that it is the go-to algorithm for competition winners on the Kaggle competitive data science platform.
        #link
            #https://machinelearningmastery.com/xgboost-for-regression/
    #About the HPs
        #To fully understand their impact see general notes and, in particular, regularization methods used in order to avoid overfitting.
            #booster: 
                #The booster to use, you can use trees as base learner or linear models. We use trees as linear does not seem to improver over a simpler linear model
                #https://stats.stackexchange.com/questions/230388/how-does-linear-base-learner-works-in-boosting-and-how-does-it-works-in-the-xgb
            #n_estimators: 
                #An important hyperparameter for the XGBoost ensemble algorithm is the number of decision trees used in the ensemble.
                #Recall that decision trees are added to the model sequentially in an effort to correct and improve upon the predictions made by prior trees. As such, more trees is often better.
                #The number of trees can be set via the “n_estimators” argument and defaults to 100.
                #The number of trees in the ensemble, often increased until no further improvements are seen.
                #We will explore up to around 1000.
            #max_depth: 
                #Varying the depth of each tree added to the ensemble is another important hyperparameter for gradient boosting.
                #The tree depth controls how specialized each tree is to the training dataset: how general or overfit it might be. Trees are preferred that are not too shallow and general (like AdaBoost) and not too deep and specialized (like bootstrap aggregation [of RF?]).
                #Gradient boosting generally performs well with trees that have a modest depth, finding a balance between skill and generality.
                #Tree depth is controlled via the “max_depth” argument and defaults to 6.
                #The maximum depth of each tree, often values are between 1 and 10. We explore this range along with None (i.e., no depth limitation).
            #min_child_weight
                #If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. The larger min_child_weight is, the more conservative the algorithm will be.
                    #https://xgboost.readthedocs.io/en/latest/parameter.html#
                #For a regression task with squared loss (the default) min_child_weight is just the number of instances in a child.
                    #https://stackoverflow.com/questions/69786993/tuning-xgboost-hyperparameters-with-randomizedsearchcv/69830350#69830350
                    #https://xgboost.readthedocs.io/en/latest/parameter.html#
                #If, for example, you have 500000 observations, it will probably not make (much of) a difference wether 1, 2, 3 or 4 observations end up in a leaf. So you have to use a number of instances that can be relevant given your sample size.
                #Used to control over-fitting. Higher values prevent a model from learning relations that might be highly specific to the particular sample selected for a tree.
                #Too high values can lead to under-fitting; hence, it should be tuned using CV.
                    #https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
                #we will explore 1, 10 and 100 as we have thousands of samples. So I only see impact after setting a minimum of at least tens of samples.
            #eta:
                #Learning rate controls the amount of contribution that each model has on the ensemble prediction.
                    #Makes the model more robust by shrinking the weights on each step
                #Smaller rates may require more decision trees in the ensemble (because each individual tree has less impact).
                #The learning rate can be controlled via the “eta” argument and defaults to 0.3.
                #It is common to have small values in the range of 0.1 to 0.3, as well as values less than 0.1 (0.01-0.2). We have also added 0.001 to cover more space.
            #subsample
                #The number of samples used to fit each tree can be varied. This means that each tree is fit on a randomly selected subset of the training dataset.
                #Using fewer samples introduces more variance for each tree, although it can improve the overall performance of the model.
                #The number of samples used to fit each tree is specified by the “subsample” argument and can be set to a fraction of the training dataset size. By default, it is set to 1.0 to use the entire training dataset.
                #Typical values: 0.5-1. We will explore from 0.1 to 1 just in case, like in RF.
            #“colsample_bytree” and “colsample_bylevel”
                #The number of features used to fit each decision tree can be varied.
                #Like changing the number of samples, changing the number of features introduces additional variance into the model, which may improve performance, although it might require an increase in the number of trees.
                #The number of features used by each tree is taken as a random sample and is specified by the “colsample_bytree” argument and defaults to all features in the training dataset, e.g. 100 percent or a value of 1.0. You can also sample columns for each split, and this is controlled by the “colsample_bylevel” argument, but we will not look at this hyperparameter here.
                #We will explore from 0.1 to 1 like in RF.
            #links
                #https://machinelearningmastery.com/xgboost-for-regression/ 
                #https://stackoverflow.com/a/69830350/12772630
                #https://datascience.stackexchange.com/a/108242
        #Additional HPs considered
            #there are more HPs that can be useful to combat overfitting and improve performance like lamba or alpha (for regularization) or gamma. This can be relevant in our case given that we have a distribution with most of the data around the extreme values (close to 0 and close to 1). Therefore, I think we can make good use of a flexible model.
                #https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
            #more manual tuning if you want
                #You can also narrow the hyperparametric space by using a sequential approach: All default, just tune max_depth and min_child_weight. See what range of values give the best and select it. Repeat now with subsample and colsample leaving all default except max_depth and min_child_weight. Select the best range for subsample and colsample. Then tune eta and n_estimators.
                    #https://datascience.stackexchange.com/a/108242
                    #https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
            #New HPS
                #gamma:
                    #A node is split only when the resulting split gives a positive reduction in the loss function. Gamma specifies the minimum loss reduction required to make a split.
                    #Larger values make the algorithm conservative. The values can vary depending on the loss function and should be tuned.
                    #Range: [0, infinite]; default 0
                    #Following analysis vidha, I have explored values between 0 and 1, getting better results with 0.4, so I going to restrict the search to that range.
                #alpha
                    #L1 regularization term on weight (analogous to Lasso regression)
                    #It can be used in case of very high dimensionality so that the algorithm runs faster when implemented
                    #Increasing this value will make model more conservative.
                    #default 0
                    #I have used the range of values from Analytics Vidhya as template, which is very similar to the range of values I used for L1 in elastic net. I have also added 1e-6 because I got the best results with 1e-5.
                    #I have checked you get the same results with alpha and reg_alpha
                #lambda:
                    #L2 regularization term on weights (analogous to Ridge regression)
                    #This is used to handle the regularization part of XGBoost. Though many data scientists don’t use it often, it should be explored to reduce overfitting.
                    #Increasing this value will make model more conservative
                    #default 1
                    #i have used the same range than in elastic net for L2 but with less data points
                #colsample_by*
                    #The parameters specify the percentage of random columns to sample from total columns available. colsample_by parameters work cumulatively, as each tree has different levels which end in nodes. The sampling can be done at each tree, level, and/or node.If you set the sampling to 0.5, you will use half off your columns. For example, the combination {colsample_bytree:0.5, colsample_bylevel: 0.5, colsample_bynode:0.5} with 64 initial features will randomly halve features used by the tree to 32, subsequently halve to 16 at the tree levels, and finally halve to 8 features at the node.
                        #https://stackoverflow.com/questions/51022822/subsample-colsample-bytree-colsample-bylevel-in-xgbclassifier-python-3-x
                    #Summary
                        #colsample_bytree is the subsample ratio of columns when constructing each tree. Subsampling occurs once for every tree constructed.
                        #colsample_bylevel is the subsample ratio of columns for each level. Subsampling occurs once for every new depth level reached in a tree. Columns are subsampled from the set of columns chosen for the current tree.
                        #colsample_bynode is the subsample ratio of columns for each node (split). Subsampling occurs once every time a new split is evaluated. Columns are subsampled from the set of columns chosen for the current level.
                    #All colsample_by* parameters have a range of (0, 1], the default value of 1, and specify the fraction of columns to be subsampled.
                    #this seems to be specially useful when you have many features, as you can use maany different subsets of features across trees, depth levels and nodes. 
                    #I have used the same range between 0 and 1 for all HPs.
                #Maximum delta step
                    #this is useful for imbalanced class. I used it anyways because we have a continuous variables with a lot of values in the extremes of the distribution.
                    #If the value is set to 0, it means there is no constraint. If it is set to a positive value, it can help making the update step more conservative. Usually this parameter is not needed, but it might help in logistic regression when class is extremely imbalanced. Set it to value of 1-10 might help control the update.
                        #https://stats.stackexchange.com/questions/233248/max-delta-step-in-xgboost
                    #It gives me problems, I am not able to find best perofrmance around optimum. Given that and the fact is recommended for logistic regression, I am NOT going to use it.

#you can have slightly different results with XGBoost even setting the seeds
    #Changing subsample and colsample_bytree  to '1' and increasing early_stopping_rounds to '1000' (or whatever n_estimators is set to) should do the trick - let me know if this solves your problem or not. – 
        #https://stackoverflow.com/questions/61764057/how-to-get-reproducible-results-from-xgboostregressor-random-state-has-no-effec
    #I have checked that increasing the number of rounds in the final model makes things more stable.





print_text("Start with the tuning", header=1)
print_text("Fix the learning rate and number of estimators for tuning tree-based parameters", header=2)
#for this step I am going to use, not only the sklearn wrapper for XGBoost (XGBRegressor), but also the original xgb package in order to use its cv function.
#with the cv function of xgb, we can run XGBoost with a pre-define number of boosts, but stopping when the there is no improvement. In this way, we can easily tune the number of estimators.
#we will update the number of estimator after tunning other parameters.
import xgboost as xgb


print_text("define a function to run CV + early top within xgb package", header=3)
#alg=xgb.XGBRegressor(objective="reg:squarederror", eval_metric="rmse", nthread=10)
#train_data=train
#predictors=[x for x in train.columns if x not in ["prob(sweep)"]]
#useTrainCV=True
#early_stopping_rounds=50
#cv_schema=cv_scheme
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
def modelfit(alg, train_data, predictors, useTrainCV=True, early_stopping_rounds=50, cv_schema=cv_scheme):
    
    #check we do not have NANs in the training data
    non_na_data_check = train_data.dropna().equals(train_data)
    if non_na_data_check == False:
        raise ValueError("ERROR: FALSE! WE HAVE NANs IN THE TRAINING DATA, THIS SCRIPT WILL GIVE PROBLEMS WITH NAN BECAUSE OF THE CONVERSION OF XGBOOST DMATRIX TO NUMPY")

    #if you want CV
    if useTrainCV:

        #extract the parameters of the XGBoost instance we have opened
        xgb_param = alg.get_xgb_params()

        #convert data to the structure required by XGBoost
        xgtrain = xgb.DMatrix(train_data[predictors].values, label=train_data["prob(sweep)"].values)
            #use predictors and response (label) of training data in the form of numpy array as input

        #define function to preprocess training and evaluation data separately within the CV
            #this function should take (dtrain, dtest, param) and returns transformed versions of those.
            #if no change is required, just return without change
        #example_index = [(train_index, test_index) for split, (train_index, test_index) in enumerate(cv_scheme.split(X=train_data)) if split==0][0]
        #dtrain=xgtrain.slice(example_index[0])
        #dtest=xgtrain.slice(example_index[1])
            #get the index of the first split and use them to slice the DMatrix into training and test using slice
            #Slice the DMatrix and return a new DMatrix that only contains rindex.
        def fpreproc(dtrain, dtest, param):

            #convert the XGBoost DMatrix to array
            dtrain_X_array = dtrain.get_data().toarray()
            dtest_X_array = dtest.get_data().toarray()
                #get_data():
                    #Get the predictors from DMatrix as a CSR matrix.
                #toarray():
                    #Return a dense ndarray representation of this sparse array.
            dtrain_y_array = dtrain.get_label()
            dtest_y_array = dtest.get_label()
                #get_label()
                    #Get the label of the DMatrix.
                #These approaches seem to give problems when having NANs, as these are considered as zero instead of missing. Therefore, it is very important that we do not have NANs in our data
                    #https://stackoverflow.com/questions/37309096/convert-python-xgboost-dmatrix-to-numpy-ndarray-or-pandas-dataframe

            #check we have the correct number of predictor columns
            if (dtrain_X_array.shape[1] != modeling_data.shape[1]-1) | (dtest_X_array.shape[1] != modeling_data.shape[1]-1):
                raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM WITH THE NUMBER OF PREDICTOR COLUMNS IN THE CONVERSION FROM DMATRIX TO NUMPY")

            #use a second approach to do the conversion of predictors
            from dmatrix2np import dmatrix_to_numpy
                #https://github.com/aporia-ai/dmatrix2np
                #https://datascience.stackexchange.com/questions/20340/python-xgboost-dmatrix-get-feature-values-or-convert-to-np-array
            dtrain_X_array_check = dmatrix_to_numpy(dtrain)
            dtest_X_array_check = dmatrix_to_numpy(dtest)
            if (np.array_equal(dtrain_X_array, dtrain_X_array_check)==False) | (np.array_equal(dtest_X_array, dtest_X_array_check)==False):
                    raise ValueError("ERROR: FALSE! WE HAVE A PROBLEM CONVERTING THE XGBOOST DMATRIX TO NUMPY")

            #open two instances of scaler to scale X and y
            scaler_X = preprocessing.StandardScaler()
            scaler_y = preprocessing.StandardScaler()

            #let the scaler learn the mean and sd of the training data
            scaler_X.fit(dtrain_X_array)
            scaler_y.fit(dtrain_y_array.reshape(-1,1))

            #transform both the training and the eval sets using the training mean and SD
            train_X_array_scaled = scaler_X.transform(dtrain_X_array)
            eval_X_array_scaled = scaler_X.transform(dtest_X_array)
            train_y_array_scaled = scaler_y.transform(dtrain_y_array.reshape(-1,1))
            eval_y_array_scaled = scaler_y.transform(dtest_y_array.reshape(-1,1))
                #pipeline works exactly like this. It calculates the mean and sd of the training data and then use it to transform both training and evaluation data

            #return to DMatrix after preprocesing
            dtrain=xgb.DMatrix(train_X_array_scaled, label=train_y_array_scaled)
            dtest=xgb.DMatrix(eval_X_array_scaled, label=eval_y_array_scaled)

            #return the scaled dtrain and dtest, along with the parameters without changes
            return (dtrain, dtest, param)
                #https://xgboost.readthedocs.io/en/stable/python/examples/cross_validation.html

        #perform CV within XGBoost
        cvresult = xgb.cv( \
            params=xgb_param, \
            dtrain=xgtrain, \
            num_boost_round=alg.get_params()["n_estimators"], \
            folds=cv_schema, \
            metrics="rmse", \
            maximize=False,
            early_stopping_rounds=early_stopping_rounds, \
            fpreproc=fpreproc, \
            seed=0, \
            verbose_eval=100)
                #params : dict
                    #Booster params.
                #dtrain : DMatrix
                    #Data to be trained.
                #num_boost_round : int
                    #Number of boosting iterations. We use the number of estimators previously selected in the XGBoost instance
                #folds : a KFold or StratifiedKFold instance or list of fold indices
                #metrics : string or list of strings
                    #Evaluation metrics to be watched in CV.
                #maximize : bool
                    #Whether to maximize evaluation metric.
                    #False in our case, because we are using root mean squared error, so larger values means more error. We want to minimize.
                #early_stopping_rounds:
                    #Activates early stopping. Cross-Validation metric (average of validation metric computed over CV folds) needs to improve at least once in every **early_stopping_rounds** round(s) to continue training. The last entry in the evaluation history will represent the best iteration. If there's more than one metric in the **eval_metric** parameter given in **params**, the last metric will be used for early stopping.
                #fpreproc : function
                    #Preprocessing function that takes (dtrain, dtest, param) and returns transformed versions of those.
                #verbose_eval : bool, int, or None, default None 
                    #Whether to display the progress. If None, progress will be displayed when np.ndarray is returned. If True, progress will be displayed at boosting stage. If an integer is given, progress will be displayed at every given `verbose_eval` boosting stage.

        #get the number of boosters before stopping
        alg.set_params(n_estimators=cvresult.shape[0])
        print("\nNumber of trees")
        print(cvresult.shape[0])
    
    ##use the model with the selected number of boosters to predict the whole training data
    #Fit the algorithm on the data
    alg.fit(preprocessing.scale(train_data[predictors]), preprocessing.scale(train_data["prob(sweep)"]))
        #you have to use the same preprocesing used during CV
        #no problem to use the whole training data for this because we already did the CV

    #Predict training set:
    dtrain_predictions = alg.predict(preprocessing.scale(train_data[predictors]))

    #Print model report:
    print("\nModel Report")
    print("R2 (train): %.10g" % r2_score(preprocessing.scale(train_data["prob(sweep)"].values), dtrain_predictions))

    #return the number of estimators
    return cvresult.shape[0]

print_text("run the function", header=3)
predictors = [x for x in train.columns if x not in ["prob(sweep)"]]
n_estimators_1 = modelfit(
    xgb.XGBRegressor(
        learning_rate=0.1,
        n_estimators=5000,
        max_depth=13,
        min_child_weight=1,
        gamma=0,
        subsample=0.7,
        colsample_bytree=0.7,
        objective="reg:squarederror",
        nthread=10, 
        eval_metric="rmse", 
        seed=0),
            #there is no R2 in XGboost
            #seed is not required when running the whole script (we set the seed at the start), but it is when running interactively just parts
            #select eval metric to avoid warning
                #https://stackoverflow.com/questions/66097701/how-can-i-fix-this-warning-in-xgboost
        train, 
        predictors)



print_text("Tune max_depth and min_child_weight", header=2)
print_text("wide search", header=3)
gsearch_depth_child_1 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_1,
                max_depth=13,
                min_child_weight=1,
                gamma=0,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                seed=0))]), 
    param_grid={ \
        "regressor__max_depth": [i for i in range(1,24,4)] + [None], \
        "regressor__min_child_weight": [1,5,10,15,25,30,40,50,100]}, \
    scoring="neg_root_mean_squared_error", #there is no R2 in xgboost. We have here different values respect to XGBoost cv, but it seems to be caused by the 'neg'. I have used mae in both, getting again different results... \
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
        #Exhaustive search over specified parameter values for an estimator.
            #The parameters of the estimator used to apply these methods are optimized by cross-validated grid-search over a parameter grid
        #estimator
            #This is assumed to implement the scikit-learn estimator interface. Either estimator needs to provide a ``score`` function, or ``scoring`` must be passed.
        #param_grid
            #Dictionary with parameters names (`str`) as keys and lists of parameter settings to try as values, or a list of such dictionaries, in which case the grids spanned by each dictionary in the list are explored. This enables searching over any sequence of parameter settings.
        #scoring
            #Strategy to evaluate the performance of the cross-validated model on the test set.
            #It can be one or several (included in a list or dict...)
        #n_jobs
            #Number of jobs to run in parallel
            #we are using only 1 because we are parallelizing per fold*model class combination.
        #Refit=True
            #Refit an estimator using the best found parameters on the WHOLE dataset used for training
        #cv:
            #Determines the cross-validation splitting strategy.
        #pre_dispatch:
            #Controls the number of jobs that get dispatched during parallel execution. Reducing this number can be useful to avoid an explosion of memory consumption when more jobs get dispatched than CPUs can process.
            #If `n_jobs` was set to a value higher than one, the data is copied for each point in the grid (and not `n_jobs` times). This is done for efficiency reasons if individual jobs take very little time, but may raise errors if the dataset is large and not enough memory is available.  A workaround in this case is to set `pre_dispatch`. Then, the memory is copied only `pre_dispatch` many times. A reasonable value for `pre_dispatch` is `2 * n_jobs`.
            #if n_jobs>1, then pre_dispatch should be "1*n_jobs", if not, the dataset will be copied many times increasing a lot memory usage
gsearch_depth_child_1.fit(train[predictors],train["prob(sweep)"])
gsearch_depth_child_1.cv_results_, gsearch_depth_child_1.best_params_, gsearch_depth_child_1.best_score_


print_text("narrow search", header=3)
gsearch_depth_child_2 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
        ('scale', preprocessing.StandardScaler()), \
        ('regressor', xgb.XGBRegressor(
            learning_rate=0.1,
            n_estimators=n_estimators_1,
            max_depth=13,
            min_child_weight=1,
            gamma=0,
            subsample=0.7,
            colsample_bytree=0.7,
            objective="reg:squarederror",
            nthread=1, 
            eval_metric="rmse", 
            seed=0))]),
    param_grid={ \
        "regressor__max_depth": [11,12,13,14,15], \
        "regressor__min_child_weight": [23,24,25,26,27]}, \
    scoring="neg_root_mean_squared_error", \
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4, \
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_depth_child_2.fit(train[predictors],train["prob(sweep)"])
gsearch_depth_child_2.cv_results_, gsearch_depth_child_2.best_params_, gsearch_depth_child_2.best_score_


print_text("Tune gamma", header=2)
print_text("wide search", header=3)
gsearch_gamma_1 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_1,
                max_depth=14,
                min_child_weight=23,
                gamma=0,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__gamma": np.arange(0, 1, 0.2)}, 
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_gamma_1.fit(train[predictors],train["prob(sweep)"])
gsearch_gamma_1.cv_results_, gsearch_gamma_1.best_params_, gsearch_gamma_1.best_score_


print_text("narrow search", header=3)
gsearch_gamma_2 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_1,
                max_depth=14,
                min_child_weight=23,
                gamma=0,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__gamma": [0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]}, 
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_gamma_2.fit(train[predictors],train["prob(sweep)"])
gsearch_gamma_2.cv_results_, gsearch_gamma_2.best_params_, gsearch_gamma_2.best_score_


print_text("narrower search", header=3)
gsearch_gamma_3 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_1,
                max_depth=14,
                min_child_weight=23,
                gamma=0,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__gamma": [0.0, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175]}, 
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_gamma_3.fit(train[predictors],train["prob(sweep)"])
gsearch_gamma_3.cv_results_, gsearch_gamma_3.best_params_, gsearch_gamma_3.best_score_



print_text("re-calibrate the number of boosters", header=2)
n_estimators_2 = modelfit(
    xgb.XGBRegressor(
        learning_rate=0.1,
        n_estimators=5000,
        max_depth=14,
        min_child_weight=23,
        gamma=0.01,
        subsample=0.7,
        colsample_bytree=0.7,
        objective="reg:squarederror",
        nthread=10, 
        eval_metric="rmse", 
        seed=0),
        train, 
        predictors)



print_text("Tune subsample and colsample_bytree", header=2)
print_text("wide search", header=3)
gsearch_sample_1 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_2,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__subsample":  [i for i in np.arange(0.1,0.8,0.3)]+[1], \
        "regressor__colsample_bytree": [i for i in np.arange(0.1,0.8,0.3)]+[1], \
        "regressor__colsample_bylevel": [i for i in np.arange(0.1,0.8,0.3)]+[1], \
        "regressor__colsample_bynode": [i for i in np.arange(0.1,0.8,0.3)]+[1]},
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_sample_1.fit(train[predictors],train["prob(sweep)"])
gsearch_sample_1.cv_results_, gsearch_sample_1.best_params_, gsearch_sample_1.best_score_


print_text("narrow search", header=3)
gsearch_sample_2 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_2,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=0.7,
                colsample_bytree=0.7,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__subsample":  [0.85, 0.9, 0.95, 1.0], \
        "regressor__colsample_bytree": [0.85, 0.9, 0.95, 1.0], \
        "regressor__colsample_bylevel": [0.6, 0.65, 0.7, 0.75], \
        "regressor__colsample_bynode": [0.6, 0.65, 0.7, 0.75]},
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_sample_2.fit(train[predictors],train["prob(sweep)"])
gsearch_sample_2.cv_results_, gsearch_sample_2.best_params_, gsearch_sample_2.best_score_



print_text("Tune regularization parameters", header=2)
print_text("wide search", header=3)
gsearch_regu_1 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_2,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=1.0,
                colsample_bytree=0.9,
                colsample_bylevel=0.65,
                colsample_bynode=0.65,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__alpha": [0, 1e-8, 1e-6, 1e-5, 1e-2, 1, 10, 100], \
        "regressor__lambda": np.arange(0, 1.01, 0.2)},
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_regu_1.fit(train[predictors],train["prob(sweep)"])
gsearch_regu_1.cv_results_, gsearch_regu_1.best_params_, gsearch_regu_1.best_score_


print_text("narrow search", header=3)
gsearch_regu_2 = GridSearchCV( \
    estimator=Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.1,
                n_estimators=n_estimators_2,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=1.0,
                colsample_bytree=0.9,
                colsample_bylevel=0.65,
                colsample_bynode=0.65,
                objective="reg:squarederror",
                nthread=1, 
                eval_metric="rmse", 
                seed=0))]), 
    param_grid={
        "regressor__alpha": [2.5e-7, 5e-7, 7.5e-7, 1e-6, 2.5e-6, 5e-6, 7.5e-6], \
        "regressor__lambda": [0.15, 0.175, 0.2, 0.225, 0.25]},
    scoring="neg_root_mean_squared_error",
    n_jobs=10, \
    cv=cv_scheme, \
    verbose=4,
    refit=True, \
    pre_dispatch="1*n_jobs")
gsearch_regu_2.fit(train[predictors],train["prob(sweep)"])
gsearch_regu_2.cv_results_, gsearch_regu_2.best_params_, gsearch_regu_2.best_score_
    #the best parameters for lambda and alpha are sligly worse than the model without these parameters
    #not sure why because I have already included the default in the wide search. Maybe random differences due to sampling and low number of booster interations...


print_text("re-calibrate the number of boosters to see the impact", header=2)
n_estimators_3 = modelfit(
    xgb.XGBRegressor(
        learning_rate=0.1,
        n_estimators=5000,
        max_depth=14,
        min_child_weight=23,
        gamma=0.01,
        subsample=1.0,
        colsample_bytree=0.9,
        colsample_bylevel=0.65,
        colsample_bynode=0.65,
        reg_lambda=0.2,
        reg_alpha=1e-6,
        objective="reg:squarederror",
        nthread=10, 
        eval_metric="rmse", 
        seed=0),
        train, 
        predictors)



print_text("reduce the learning rate and re-calibrate the number of boosters", header=2)
n_estimators_4 = modelfit(
    xgb.XGBRegressor(
        learning_rate=0.01,
        n_estimators=5000,
        max_depth=14,
        min_child_weight=23,
        gamma=0.01,
        subsample=1.0,
        colsample_bytree=0.9,
        colsample_bylevel=0.65,
        colsample_bynode=0.65,
        reg_lambda=0.2,
        reg_alpha=1e-6,
        objective="reg:squarederror",
        nthread=10, 
        eval_metric="rmse", 
        seed=0),
        train, 
        predictors)
        #After 5000 improvements are very slow, from 0.608 to 0.606 after 14300 iterations. Probably not worth it.



print_text("check performance after tunning XGBoost parameters", header=2)
print_text("define the model", header=3)
final_model = Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.01,
                n_estimators=n_estimators_4,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=1.0,
                colsample_bytree=0.9,
                colsample_bylevel=0.65,
                colsample_bynode=0.65,
                reg_lambda=0.2,
                reg_alpha=1e-6,
                objective="reg:squarederror",
                nthread=10, 
                eval_metric="rmse", 
                seed=0))])
final_model.fit(train[predictors], train["prob(sweep)"])


print_text("predict on the test set", header=3)
y_pred = final_model.predict(test[predictors])


print_text("calculate evaluation score in the test set", header=3)
score = r2_score(test["prob(sweep)"], y_pred)
print(score)
#We get R2=0.66, instead of 0.69 like in the model class comparison. The difference is caused by the size of the test set, as in this case we are using 20% while in the first case we used the 10% (9 folds for training/evaluation and 1 for test). If we decrease the test size to 10%, then we got R2=0.705!! so we gained predictive power after doing all this manual tuning.





print_text("combine different instance of XGBoost with same best HP combination but different seeds", header=1)
#As Janizek et al 2023 did, we want to consider the variability associated with multiple good fitting models. Following their example, we will run several XBoost models with the best HP combination (already defined) but different seed, so the sampling of samples and features across iterations will be different between trees. 
    #After selecting XGBoost as the best-performing model class for the prediction of anti-AML drug synergy, we then wanted to account for the full diversity of possible good XGBoost models fit to the highly cor- related AML gene expression data. We therefore trained 100 models and explained the ensemble model. Each individual model had both row and column subsampling turned on for each additional tree fit, and the difference between the models in the ensemble was the random seed given to generate the subsampling
        #Note that you can have slightly different results with XGBoost even setting the seeds
            #Changing subsample and colsample_bytree  to '1' and increasing early_stopping_rounds to '1000' (or whatever n_estimators is set to) is used to increase stability        
                #https://stackoverflow.com/questions/61764057/how-to-get-reproducible-results-from-xgboostregressor-random-state-has-no-effec
        #so it makes sense to run several models with different seeds.
    #After identifying GBMs as the best-performing model class for our dataset, we ensembled individual models until the ensemble model attributions were stable, leading to a final ensemble of 100 XGBoost models (see Supplementary Fig. 2 and Methods). We then analysed the resultant ensemble model attributions to look for genes with ‘global’ importance for drug combination synergy, that is, genes whose expres- sion is related to synergy across many different drug pairs in our data- set18.
    #These steps are described in the script "ensemble_explanations.py" of this manuscript
#We are using just 5 XGBoost instance because of 
    #computation time limitations. We have to use around 5K booster iterations per XGBoost instance, so this makes things slower.
    #using a low number of boosters, I have checked what happens by adding 10 instead of 5 models (i.e., instances), and there is no improvement.
    #maybe check in the future
#This approach of combining models with different seeds for improving performance is not strange. See this comment form Analytics Vidhya
    #As I mentioned in the end, techniques like feature engineering and blending have a much greater impact than parameter tuning. For instance, I generally do some parameter tuning and then run 10 different models on same parameters but different seeds. Averaging their results generally gives a good boost to the performance of the model. Hope this helps. Please share your thoughts.
        #https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
#People says that this is usually useful when you have very different models, but as I said, other people and papers also use this to combine different instances of the same model buth with different seed. I have seen 1% of improvement in R2 using this approach.
    #Stacking is appropriate when multiple different machine learning models have skill on a dataset, but have skill in different ways. Another way to say this is that the predictions made by the models or the errors in predictions made by the models are uncorrelated or have a low correlation.
    #https://machinelearningmastery.com/stacking-ensemble-machine-learning-with-python/
    #https://www.thekerneltrip.com/python/stacking-an-introduction/
#Voting regressor vs stacking
#The fundamental difference between voting and stacking is how the final aggregation is done. In voting, user-specified weights are used to combine the classifiers whereas stacking performs this aggregation by using a blender/meta classifier. In other words, stacking uses a meta-estimator to model what weight you have to give to the prediction of each individual model in order to obtain a global prediction similar to the observed. You are modeling the response as a function of the individual predictions of each base estimator.
    #https://towardsdatascience.com/ensemble-methods-comparing-scikit-learns-voting-classifier-to-the-stacking-classifier-f5ab1ed1a29d
#I have found better results with stacking and I do not have to define the combinations of importance between estimators, so I am going to use it.
#Super-learner
#this is an special type of stacking, but it seems to be used when you have hundreds of models.
    #https://machinelearningmastery.com/super-learner-ensemble-in-python/



print_text("run GS for tuning the meta-estimator", header=2)
print_text("define the CV schema for stacking", header=3)
cv_stacking = KFold( \
    n_splits=5,  \
    shuffle=True)
#The final stacked model will be evaluated using the previous defined CV schema with 5 folds, but we need an additional CV schema in order to train the meta-estimator with the predictions of the base estimators. See explanations above for further details
    #https://machinelearningmastery.com/stacking-ensemble-machine-learning-with-python/

#it seems that the cv of StackingRegressor finds parameters of the meta-estimator, i.e., the importance given to each base estimator, so the HPs of the base estimators are not optimized, nor the HPs of the meta-estimator (in case it has, for example alpha for Ridge)
    #https://machinelearningmastery.com/stacking-ensemble-machine-learning-with-python/  

#info from scikit learn
#During training, the estimators are fitted on the whole training data X_train. They will be used when calling predict or predict_proba. To generalize and avoid over-fitting, the final_estimator is trained on out-samples using sklearn.model_selection.cross_val_predict internally.
    #https://scikit-learn.org/stable/modules/ensemble.html#stacking

#the prediction of a base estimator is done knowing the mean/sd of the whole training data, yes, but it would also know about the training data even if we do not scale because StackingRegressor fit the base estimators using the whole training data, while use CV fitting the meta-estimator on the predictions of the base estimators.
    #note that we are using another CV schema on top of the stackingregressor schema, so in each repetition, one fold is not used for any of the base estimators, so the prediction in that fold would give info about generalization ability. We are predicting exactly in the same fold we have predicted when doing HP tunning for XGBoost




print_text("run the stacked model and evaluate", header=2)
#LinearRegression as meta-estimator
#The meta-model is often simple, providing a smooth interpretation of the predictions made by the base models. As such, linear models are often used as the meta-model, such as linear regression for regression tasks (predicting a numeric value) and logistic regression for classification tasks (predicting a class label). Although this is common, it is not required.
#The level-1 model or meta-model is provided via the “final_estimator” argument. By default, this is set to LinearRegression for regression and LogisticRegression for classification, and these are sensible defaults that you probably do not want to change.
    #https://machinelearningmastery.com/stacking-ensemble-machine-learning-ith-python/
#Often simple ensembling models (e.g. simple average of predicted probabilities, weighted average of predicted probabilities or logistic regression of logits - possibly with regularization towards a simple average) perform a lot better than trying to fit fancier models on the second level
    #https://stats.stackexchange.com/questions/561584/why-is-my-stacking-meta-learning-not-outperforming-the-best-base-model
from sklearn.ensemble import StackingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
score_stacking = cross_val_score( \
    estimator=Pipeline( \
        steps=[ \
                ('scale', preprocessing.StandardScaler()), \
                ("regressor", StackingRegressor( \
                    estimators=[( \
                        "est_"+str(seed), \
                        xgb.XGBRegressor( \
                            learning_rate=0.01,
                            n_estimators=n_estimators_4,
                            max_depth=14,
                            min_child_weight=23,
                            gamma=0.01,
                            subsample=1.0,
                            colsample_bytree=0.9,
                            colsample_bylevel=0.65,
                            colsample_bynode=0.65,
                            reg_lambda=0.2,
                            reg_alpha=1e-6,
                            objective="reg:squarederror",
                            nthread=1, 
                            eval_metric="rmse",  
                            seed=seed)) \
                        for seed in range(0,5,1)], \
                    final_estimator=LinearRegression(), \
                    cv=cv_stacking, \
                    n_jobs=1))]), \
    X=train[predictors],
    y=train["prob(sweep)"], 
    scoring="neg_root_mean_squared_error", \
    cv=cv_scheme, \
    n_jobs=10, \
    verbose=4, \
    error_score='raise', \
    pre_dispatch="1*n_jobs")
        #StackingRegressor, passthrough=False as default
            #Sometimes, better performance can be achieved if the dataset prepared for the meta-model also includes inputs to the level-0 models, e.g. the input training data. This can be achieved by setting the “passthrough” argument to True and is not enabled by default.
print(f"Mean neg_root_mean_squared_error across repetitions {np.mean(score_stacking)}")
    #the lowest value obtained during HP tuning in XGBoost was -1.8744
    #we get here under the same CV scheme an average of -1.837 across 5 folds. Results per split are:
        #[CV] END ... score: (test=-1.800) total time=70.3min
        #[CV] END ... score: (test=-1.807) total time=70.3min
        #[CV] END ... score: (test=-1.860) total time=70.5min
        #[CV] END ... score: (test=-1.896) total time=70.6min
        #[CV] END ... score: (test=-1.824) total time=70.7min
        #Mean neg_root_mean_squared_error across repetitions -1.8372924401780832
            #it took a few hours (not 10), because each split was run in parallel (10 jobs)




print_text("check performance after stacking in test set", header=2)
print_text("define the model using the best estimator from the last GS", header=3)
final_model_stacking = Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ("regressor", StackingRegressor( \
                estimators=[( \
                    "est_"+str(seed), \
                    xgb.XGBRegressor( \
                        learning_rate=0.01,
                        n_estimators=n_estimators_4,
                        max_depth=14,
                        min_child_weight=23,
                        gamma=0.01,
                        subsample=1.0,
                        colsample_bytree=0.9,
                        colsample_bylevel=0.65,
                        colsample_bynode=0.65,
                        reg_lambda=0.2,
                        reg_alpha=1e-6,
                        objective="reg:squarederror",
                        nthread=1, 
                        eval_metric="rmse",  
                        seed=seed)) \
                    for seed in range(0,5,1)], \
                final_estimator=LinearRegression(), \
                cv=cv_stacking, \
                n_jobs=10))])
final_model_stacking.fit(train[predictors], train["prob(sweep)"])
import joblib
joblib.dump(final_model_stacking, "./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_training_dataset.pkl.gz", compress=("gzip", True))
    #https://stackoverflow.com/questions/34143829/sklearn-how-to-save-a-model-created-from-a-pipeline-and-gridsearchcv-using-jobli
#final_model_stacking=joblib.load("./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_training_dataset.pkl.gz")


print_text("predict on the training and test set", header=3)
y_pred_stacking_train = final_model_stacking.predict(train[predictors])
y_pred_stacking_test = final_model_stacking.predict(test[predictors])

print_text("calculate evaluation scores", header=3)
from sklearn.metrics import mean_squared_error
final_model_stacking_rmse_train = mean_squared_error(train["prob(sweep)"], y_pred_stacking_train, squared=False)
final_model_stacking_rmse_test = mean_squared_error(test["prob(sweep)"], y_pred_stacking_test, squared=False)
    #sklearn.metrics has a mean_squared_error function with a squared kwarg (defaults to True). Setting squared to False will return the RMSE.
        #y_true
            #Ground truth (correct) target values.
        #y_pred
            #Estimated target values.
        #sample_weight=None
            #Sample weights.
        #multioutput=’uniform_average’
            #Defines aggregating of multiple output values.
                #‘uniform_average’ :
                    #Errors of all outputs are averaged with uniform weight.
        #squared=True
            #If True returns MSE value, if False returns RMSE value.
        #https://stackoverflow.com/a/18623635/12772630
final_model_stacking_score_train = r2_score(train["prob(sweep)"], y_pred_stacking_train)
final_model_stacking_score_test = r2_score(test["prob(sweep)"], y_pred_stacking_test)
print(f"RMSE of final model on train set: {final_model_stacking_rmse_train}")
print(f"RMSE of final model on test set: {final_model_stacking_rmse_test}")
print(f"R2 of final model on train set: {final_model_stacking_score_train}")
print(f"R2 of final model on test set: {final_model_stacking_score_test}")
#RMSE of final model on train set: 0.2842961032821582
#RMSE of final model on test set: 1.7297735272800812
#R2 of final model on train set: 0.9913054963448299
#R2 of final model on test set: 0.6760527299906672


print_text("plot predicted vs observes in training and test sets", header=3)
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.scatter(train["prob(sweep)"], y_pred_stacking_train, s=0.5)
ax1.set_ylabel("Predicted log(Flex-sweep probability)")
ax1.yaxis.set_label_coords(-0.11, -0.1)
ax1.yaxis.label.set_size(13.5)
ax1.annotate("Root-mean-square error (RMSE; training): " + str(np.round(final_model_stacking_rmse_train, 4)), xy=(0.025, 0.85), xycoords='axes fraction')
ax1.annotate("$R^{2}$ (training): " + str(np.round(final_model_stacking_score_train, 4)), xy=(0.025, 0.75), xycoords='axes fraction')
ax2.scatter(test["prob(sweep)"], y_pred_stacking_test, s=0.5)
ax2.set_xlabel("Observed log(Flex-sweep probability)")
ax2.xaxis.label.set_size(13.5)
ax2.annotate("RMSE (test): " + str(np.round(final_model_stacking_rmse_test, 4)), xy=(0.025, 0.85), xycoords='axes fraction')
ax2.annotate("$R^{2}$ (test): " + str(np.round(final_model_stacking_score_test, 4)), xy=(0.025, 0.75), xycoords='axes fraction')
plt.savefig( \
    fname="./results/selected_model_class/yoruba_hg19_stack_5_xgboost_scatter.png", dpi=300)
plt.close()
    #Increased error at top sweep candidated
        #There is an increase in the residuals (error) as we get closer to high Flex-Sweep values, specifically, the residuals tend to be more negative, i.e., predicted values tend to be smaller than observed.
        #My intuition is the following: 
            #At low-medium flex-sweep values, genomic features like recombination rate or conserved elements density explain well the distribution of sweep signals. Lower recombination makes easier to detect sweep like higher conservation increase the probability of having regions with phenotypic relevance. 
            #However, as we get closer to top flex-sweep candidates, these features are not enough. We have here genes whose functions are likely target of different selective pressures, so they have more sweep signals compared to other genes with similar genomic features. We have not included all selective pressures influencing the human genome, thus it is not surprising that we have unexplained variability in strong sweep candidates.
    #predicted probability above 1
        #There are some genes with a predicted value above 0, meaning flex-sweep probability above 1. This is caused because the distribution of flexsweep is veery skewed, with a very high and narrow peak at 0.99-1.
        #Therefore, it is not surprising that the model smooth that peak making it a bit wider and then surpassing zero a bit. 
        #Note, anyways, that most of the genes at that level (top sweep candidates) are predicted to have less (not more) flex-sweep probabilities, supporting the explanation of the previous point.



print_text("prepare the TRUE final model", header=1)
print_text("write the final model from scratch and fit to the whole data", header=2)
final_model_stacking_full_data = Pipeline( \
        steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ("regressor", StackingRegressor( \
                estimators=[( \
                    "est_"+str(seed), \
                    xgb.XGBRegressor( \
                        learning_rate=0.01,
                        n_estimators=n_estimators_4,
                        max_depth=14,
                        min_child_weight=23,
                        gamma=0.01,
                        subsample=1.0,
                        colsample_bytree=0.9,
                        colsample_bylevel=0.65,
                        colsample_bynode=0.65,
                        reg_lambda=0.2,
                        reg_alpha=1e-6,
                        objective="reg:squarederror",
                        nthread=1, 
                        eval_metric="rmse",  
                        seed=seed)) \
                    for seed in range(0,5,1)], \
                final_estimator=LinearRegression(), \
                cv=cv_stacking, \
                n_jobs=10))])
final_model_stacking_full_data.fit(modeling_data[predictors], modeling_data["prob(sweep)"])
print_text("save the model", header=2)
import joblib
joblib.dump(final_model_stacking_full_data, "./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_full_dataset.pkl.gz", compress=("gzip", True))
    #https://stackoverflow.com/questions/34143829/sklearn-how-to-save-a-model-created-from-a-pipeline-and-gridsearchcv-using-jobli
#final_model_stacking_full_data=joblib.load("./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_full_dataset.pkl.gz")
