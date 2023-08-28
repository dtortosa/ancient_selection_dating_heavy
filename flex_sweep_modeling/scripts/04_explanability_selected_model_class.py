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
        ./results/selected_model_class/ale_plots/saved_explained; \
    mkdir \
        --parents \
        ./results/selected_model_class/permutation_importance; \
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





print_text("feature explanations in training/test/whole?", header=2)
#Permutation importance
    #https://christophm.github.io/interpretable-ml-book/feature-importance.html#feature-importance-data
    #the case for using test data
        #This is a simple case: Model error estimates based on training data are garbage -> feature importance relies on model error estimates -> feature importance based on training data is garbage. Really, it is one of the first things you learn in machine learning: If you measure the model error (or performance) on the same data on which the model was trained, the measurement is usually too optimistic, which means that the model seems to work much better than it does in reality. And since the permutation feature importance relies on measurements of the model error, we should use unseen test data. The feature importance based on training data makes us mistakenly believe that features are important for the predictions, when in reality the model was just overfitting and the features were not important at all.
    #the case for using training data
        #The arguments for using training data are somewhat more difficult to formulate, but are IMHO just as compelling as the arguments for using test data. We take another look at our garbage SVM. Based on the training data, the most important feature was X42. Let us look at a partial dependence plot of feature X42. The partial dependence plot shows how the model output changes based on changes of the feature and does not rely on the generalization error. It does not matter whether the PDP is computed with training or test data.
        #The plot clearly shows that the SVM has learned to rely on feature X42 for its predictions, but according to the feature importance based on the test data (1), it is not important. Based on the training data, the importance is 1.19, reflecting that the model has learned to use this feature. Feature importance based on the training data tells us which features are important for the model in the sense that it depends on them for making predictions.
        #As part of the case for using training data, I would like to introduce an argument against test data. In practice, you want to use all your data to train your model to get the best possible model in the end. This means no unused test data is left to compute the feature importance. You have the same problem when you want to estimate the generalization error of your model. If you would use (nested) cross-validation for the feature importance estimation, you would have the problem that the feature importance is not calculated on the final model with all the data, but on models with subsets of the data that might behave differently.
        #However, in the end I recommend to use test data for permutation feature importance. Because if you are interested in how much the model’s predictions are influenced by a feature, you should use other importance measures such as SHAP importance.
#ALE plots and PDPs
    #Christopher says that "The partial dependence plot shows how the model output changes based on changes of the feature and does not rely on the generalization error. It does not matter whether the PDP is computed with training or test data.". I think this also applies for ALE plots as they do not use the performence, i.e., an error metric from the model obtained by comparing observed and predicted, but they show changes in the prediction as the feature changes in certain data points. Christopher talks in length about the problem of calculating feature importance by perturbation because it requires to calculate changes in predictive power. See previous section
    #In addition, I saw the follwing comment "Partial dependence plots can be performed over either the training or validation set, and examples of both cases can be found. For instance, fastbook uses the validation set, whereas Interpretable Machine Learning: A Guide for Making Black Box Models Explainable, in chapter 8.1 4, uses the training set. Frankly, in my experience, it ultimately does not matter which strategy you choose, and there are only a couple of important considerations. First, if the training set is large, PDP may take excessively long, in which case you can use a small chunk of it or resort to the validation set. Second, if the validation set is too small, PDP’s results would understandably be not very reliable, so the training set might be the wiser option."
    #Therefore, we can use data previously used for training in ALE plots, but in permutation test, which is based in changes in predictive power, we should use the test set.





print_text("permutation importance", header=2)
#The concept is really straightforward: We measure the importance of a feature by calculating the increase in the model’s prediction error after permuting the feature. A feature is “important” if shuffling its values increases the model error, because in this case the model relied on the feature for the prediction. A feature is “unimportant” if shuffling its values leaves the model error unchanged, because in this case the model ignored the feature for the prediction. 
    #https://christophm.github.io/interpretable-ml-book/feature-importance.html#feature-importance-data
#It should be used in the test set (or CV see previous section)
#Important limitations for us
    #If features are correlated, the permutation feature importance can be biased by unrealistic data instances. The problem is the same as with partial dependence plots: The permutation of features produces unlikely data instances when two or more features are correlated. In other words, for the permutation feature importance of a correlated feature, we consider how much the model performance decreases when we exchange the feature with values we would never observe in reality. 
    #Another tricky thing: Adding a correlated feature can decrease the importance of the associated feature by splitting the importance between both features.
    #The permutation is random, so the results can be different across runs.
    #WE SOLVE THIS PROBLEMS WITH SHAP IN THE FUTURE



print_text("split train/test set", header=2)
#permutation importance should be done on test set! so we use the same test set we used for test after HP tuning.
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



print_text("load the final model fit to the training data", header=2)
import joblib
final_model_training=joblib.load("./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_training_dataset.pkl.gz")
print(final_model_training)



print_text("run permutation importance", header=2)
print_text("define function for using rmse as loss function", header=3)
#This is the same metric we have used for CV in the training set. R2 orders the features in the same way but strangely, we get a lot of variability (large error bar) for recombination rate. Therefore, we are using rmse, which does not have this problem and produce similar results for the rest. 
from sklearn.metrics import mean_squared_error
def loss_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return mean_squared_error(y_true=y_true, y_pred=y_pred, squared=False)
        #alibi requires to define the signature of the function. The signatures is number and type of input arguments the function takes and the type of the result the function returns.
            #inputs
                #y_true is the array of ground-truth values, 
                #y_pred is the output of the predictor used in the initialization of the explainer
            #outputs
                #a float
            #https://stackoverflow.com/a/72789090/12772630
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


print_text("get feature names in a nice format", header=3)
dict_nice_feature_names={
    "chip_density_1000kb": "Regulatory density (ChIP-seq)",
    "chip_immune_1000kb": "Immune regulatory density (ChIP-seq)",
    "chip_testis_1000kb": "Testis regulatory density (ChIP-seq)",
    "coding_density_1000kb": "Coding density",
    "cons_elements_1000kb": "Density of conserved elements",
    "gc_content_1000kb": "GC-content",
    "expression_all_tissues": "Gene expression across multiple tissues",
    "testis_expression": "Gene expression in testis",
    "lympho_expression": "Gene expression in lymphocytes",
    "gene_length": "Gene length",
    "gene_number_1000kb": "Gene number",
    "n_dis_1000kb": "Number of Mendelian disease variants",
    "pip_v3": "Number of protein-protein interactions",
    "recombination_final_1000kb": "Recombination rate",
    "tbfs_density_1000kb": "Regulatory density (DNAseI)",
    "thermogenic_distance": "Distance to thermogenic genes",
    "vip_distance": "Distance to VIPs",
    "bat_distance_percentile_1": "Distance to Brown Adipose Tissue (BAT) genes",
    "smt_distance_percentile_1": "Distance Skeletal Muscle Tissue (SMT) genes"}
print("check that the keys in the dict are the features in the same order than in the data")
predictors=[i for i in modeling_data.columns if i not in ["prob(sweep)"]]
predictors_nice=[dict_nice_feature_names[i] for i in modeling_data.columns if i not in ["prob(sweep)"]]
print(predictors_nice)


print_text("Initialize the permutation feature importance", header=3)
from alibi.explainers import PermutationImportance
#Implementation of the permutation feature importance for tabular datasets. The method measure the importance of a feature as the relative increase/decrease in the loss/score function when the feature values are permuted. Supports black-box models.
    #https://docs.seldon.io/projects/alibi/en/stable/examples/permutation_importance_classification_leave.html
    #https://docs.seldon.io/projects/alibi/en/stable/methods/PermutationImportance.html
    #https://docs.seldon.io/projects/alibi/en/stable/api/alibi.explainers.html#alibi.explainers.PermutationImportance
explainer_perm = PermutationImportance( \
    predictor=final_model_training.predict, \
    loss_fns={"rmse": loss_rmse}, \
    score_fns=None, \
    feature_names=predictors_nice, \
    verbose=True)
    #Initialize the permutation feature importance.
        #https://docs.seldon.io/projects/alibi/en/stable/api/alibi.explainers.html#alibi.explainers.PermutationImportance
        #predictor:
            #A PREDICTION FUNCTION which receives as input a `numpy` array of size `N x F`, and outputs a | `numpy` array of size `N` (i.e. `(N, )`) or `N x T`, where `N ` is the number of input instances, `F` is the number of features, and `T` is the number of targets.
        #loss_fns=None
            #loss function, or a dictionary of loss functions
            #options: 
                #[‘mean_absolute_error’, ‘mean_squared_error’, ‘mean_squared_log_error’, ‘mean_absolute_percentage_error’, ‘log_loss’],
        #score_fns=None
            #scorer function
            #options:
                #[‘mean_absolute_error’, ‘mean_squared_error’, ‘mean_squared_log_error’, ‘mean_absolute_percentage_error’, ‘log_loss’]
        #Note about the metrics
            #Note that based on the metric used, the importance of the features, and implicitly their ordering can differ.
            #I have checked that in our case RMSE and R2 give the same order
        #feature_names
            #A list of feature names used for displaying results.
        #verbose=False
            #Whether to print the progress of the explainer.
    #Note:
        #Remember that the PermutationImportance explainer measures the importance of a feature f as the degradation of the model when the feature values of f are permuted. The degradation of the model can thus be quantified as either the increase in the loss function or the decrease in the score function. Although one can transform a loss function into a score function an vice-versa (i.e., simply negate the value and optionally add an offset), the equivalent representation might not be always be natural to interpret (e.g., transforming mean squared error loss into the equivalent score given by the negative mean squared error). Thus, the alibi API allows the user to provide the suitable metric either as a loss or a score function.


print_text("compute permutation importance", header=3)
perm_exp = explainer_perm.explain( \
    X=test[predictors].to_numpy(),  \
    y=test["prob(sweep)"].to_numpy(),  \
    method="estimate",  \
    kind="ratio", \
    n_repeats=50)
        #Computes the permutation feature importance for each feature with respect to the given loss or score functions and the dataset (X, y).
            #X (ndarray)
                #A N x F input feature dataset used to calculate the permutation feature importance. This is typically the test dataset.
            #y (ndarray) 
                #Ground-truth labels array of size N (i.e. (N, )) corresponding the input feature X.
            #features 
                #An optional list of features or tuples of features for which to compute the permutation feature importance. If not provided, the permutation feature importance will be computed for every single features in the dataset. Some example of features would be: [0, 2]
            #method="estimate"
                #[‘estimate’, ‘exact’])
                #The method to be used to compute the feature importance. 
                #If set to 'exact', a “switch” operation is performed across all observed pairs, by excluding pairings that are actually observed in the original dataset. This operation is quadratic in the number of samples (N x (N - 1) samples) and thus can be computationally intensive. 
                #If set to 'estimate', the dataset will be divided in half. The values of the first half containing the ground-truth labels the rest of the features (i.e. features that are left intact) is matched with the values of the second half of the permuted features, and the other way around. This method is computationally lighter and provides estimate error bars given by the standard deviation. Note that for some specific loss and score functions, the estimate does not converge to the exact metric value.
                    #In other words, "the dataset is divided in half and the first and half values for Y and X2 (the predictor not permuted) are matched with the second half values of X1, and the other way around. Besides the light computation, this approach can provide confidence intervals by computing the estimates over multiple data splits."
                    #so you just make one split of the data and match Y from 1 half to X1 values of the other side. This should be faster and you can ran several iterations with different random sets, then obtaining an average and SD
                        #https://docs.seldon.io/projects/alibi/en/stable/methods/PermutationImportance.html 
            #kind=ratio
                #[‘ratio’, ‘difference’]
                #Whether to report the importance as the loss/score ratio or the loss/score difference.
                    #i.e., difference between the base value and the value after permutation, or calculate the ratio.
                #if ratio
                    #importance of 1 means that a feature not relevant for the model (i.e., we are using the ratio between the base and permuted score and a ratio close to 1 means that the original score and the permuted score are approximately the same).
            #n_repeats=50
                #Number of times to permute the feature values. Considered only when method='estimate'.
            #sample_weight
                #Optional weight for each sample instance.


print_text("save permutations with pickle", header=3)
import pickle
pickle.dump(perm_exp, open("./results/selected_model_class/permutation_importance/yoruba_hg19_stack_5_xgboost_permutation_importance_test_set.sav", 'wb'))
#import pickle; perm_exp = pickle.load(open("./results/selected_model_class/permutation_importance/yoruba_hg19_stack_5_xgboost_permutation_importance_test_set.sav", 'rb'))


print_text("plot results with Recombination", header=3)
import matplotlib.pyplot as plt
from alibi.explainers import plot_permutation_importance
ax=plot_permutation_importance( \
    exp=perm_exp, \
    features="all", \
    metric_names="all", \
    n_cols=1, \
    sort=True, \
    top_k=None, \
    ax=None, \
    bar_kw=None, \
    fig_kw={'figwidth': 14, 'figheight': 6})[0][0]
    #exp
        #An Explanation object produced by a call to the alibi.explainers.permutation_importance.PermutationImportance.explain() method.
    #features="all"
        #A list of feature entries provided in feature_names argument to the alibi.explainers.permutation_importance.PermutationImportance.explain() method, or 'all' to plot all the explained features. You have to use numbers.
    #metric_names="all"
        #A list of metric entries in the exp.data[‘metrics’] to plot the permutation feature importance for, or 'all' to plot the permutation feature importance for all metrics (i.e., loss and score functions). The ordering is given by the concatenation of the loss metrics followed by the score metrics.
    #n_cols
        #Number of columns to organize the resulting plot into.
    #sort=True
        #Boolean flag whether to sort the values in descending order.
    #top_k=None
        #Number of top k values to be displayed if the sort=True. If not provided, then all values will be displayed.
    #ax
        #A matplotlib axes object or a numpy array of matplotlib axes to plot on
        #you can add this figure to a axis already defined
    #bar_kw
        #Keyword arguments passed to the matplotlib.pyplot.barh function.
    #fig_kw –
        #Keyword arguments passed to the matplotlib.figure.set function.
ax.update({"xlim": (1, 1.8), "xlabel": "Permutation importance (original RMSE / permuted RMSE)"})
    #update the axis to change the xlim and focus on values above 1
    #https://www.geeksforgeeks.org/matplotlib-axes-axes-update-in-python/
ax.set_title("", fontsize=18)
ax.xaxis.label.set_size(13.5)
ax.tick_params(axis='y', labelsize=11)
plt.savefig( \
    fname="./results/selected_model_class/permutation_importance/permutation_importance_all.png", dpi=300)
plt.close()
    #RESULTS:
        #As expected, recombination rate is the most important feature by far. the we have the density of conserved elements.
        #BAT and SMT are the next more important features according to permutation and then we have mulitple features that are expected to be in the top like GC-content or regulatory density (DNAseI).
        #VIPs are relatively low, this can be explained by the correlation between features. If two features are correlated, their importance is split, maybe this makes VIPs to be reduced because we have VIPs in BAT and thermogenic genes.
            #404 out 923 thermogenic genes are VIPs
                #modeling_data.loc[(modeling_data["thermogenic_distance"]==0) & (modeling_data["vip_distance"]==0), ["vip_distance", "bat_distance_percentile_1"]]
            #53 out 149 BAT genes are VIPs
                #modeling_data.loc[(modeling_data["bat_distance_percentile_1"]==0) & (modeling_data["vip_distance"]==0), ["vip_distance", "bat_distance_percentile_1"]]
            #69 out 143 SMT genes are VIPs
                #modeling_data.loc[(modeling_data["smt_distance_percentile_1"]==0) & (modeling_data["vip_distance"]==0), ["vip_distance", "bat_distance_percentile_1"]]




print_text("ALE plots", header=1)
#Accumulated local effects describe how features influence the prediction of a machine learning model on average. ALE plots are a faster and unbiased alternative to partial dependence plots (PDPs). Given that we have correlated features (see above), we have to consider this bias. 
    #ALE explanations obtained from this book: https://christophm.github.io/interpretable-ml-book/ale.html
#If features of a machine learning model are correlated, the partial dependence plot cannot be trusted. The computation of a partial dependence plot for a feature that is strongly correlated with other features involves averaging predictions of artificial data instances (samples) that are unlikely in reality. This can greatly bias the estimated feature effect. Imagine calculating partial dependence plots for a machine learning model that predicts the value of a house depending on the number of rooms and the size of the living area. We are interested in the effect of the living area on the predicted value. As a reminder, the recipe for partial dependence plots is: 1) Select feature. 2) Define grid of values of the feature to be tested. 3) Per grid value: a) Replace feature with the selected grid value in each samples and b) average predictions across samples; repeate with another values of the grid... 4) Draw curve across grid values. In this way, you can see what happens with the prediction across a range of values of the feature of interest.
#For the calculation of the first grid value of the PDP – say 30 m2 – we replace the living area for all samples by 30 m2, even for houses with 10 rooms. Sounds to me like a very unusual house, 10 rooms but only 30 m2?. The partial dependence plot includes these unrealistic houses in the feature effect estimation and pretends that everything is fine.
#In figure 8.5 of the tutorial you can see this with two continuous variables: When analyzing the value of x1=0.75, all samples get a value of 0.75 for this feature, including those with a x2=0.2. If you see the relationship between the 2 features, there are no samples with x1=0.75 and x=0.2, this is a combination that does not exist in the data because when x1 is high, x2 is also high (both features are positively correlated). 
#In our case, for example, a very low coding density with very high number of genes is unlikely, if this is not present in the dataset, specially if this is the whole dataset, it does not make ansy sense to predict selection considering this combination to check the impact of high number of genes in the predictions.
#What can we do to get a feature effect estimate that respects the correlation of the features? We could average over the conditional distribution of the feature, meaning at a grid value of x1, we average the predictions of instances with a similar x1 value. The solution for calculating feature effects using the conditional distribution is called Marginal Plots, or M-Plots (confusing name, since they are based on the conditional, not the marginal distribution). Wait, did I not promise you to talk about ALE plots? M-Plots are not the solution we are looking for. Why do M-Plots not solve our problem? If we average the predictions of all houses of about 30 m2, we estimate the combined effect of living area and of number of rooms, because of their correlation. If number of rooms and area are correlated, houses with around 30 m2 will tend to have a similar number of rooms. Suppose that the living area has no effect on the predicted value of a house, only the number of rooms has. The M-Plot would still show that the size of the living area increases the predicted value, since the number of rooms increases with the living area.
#Figure 8.6 shows for two correlated features how M-Plots work. Now, we change to x1=0.75 for only those with x2 around 0.8, which are the ones actually having values close x1=0.75, so we avoid unlikely combinations of the features,like x1=0.75 and x2=0.2 but, again, the effect on the prediction is because x1=0.75 or because x2 is equal to 0.8?
#Just to clarify: We take samples with x1 around 0.75. For each one, we change x1 to 0.75 and make a prediction. Then average across all predicitons for these samples. Repeate with samples having x1 around 0.77, 0.79 and so on. Imagine that the average prediction is higher in samples around x1=0.79 compared to samples with x1 around 0.75. This means that the model has learned to use x1 as a positive predictor? Not necessarily, becuase the samples with x1 around 0.79 also tend to have higher values of x2 compared to the samples around x1=0.75, as we said, x1 and x2 are positively correlated. Therefore, what is causing the increase in the prediction, x1 or x2?
#M-Plots avoid averaging predictions of unlikely data instances, but they mix the effect of a feature with the effects of all correlated features. ALE plots solve this problem by calculating – also based on the conditional distribution of the features – differences in predictions instead of averages. For the effect of living area at 30 m2, the ALE method uses all houses with about 30 m2, gets the model predictions pretending these houses were 31 m2 minus the prediction pretending they were 29 m2. This gives us the pure effect of the living area and is not mixing the effect with the effects of correlated features. The use of differences blocks the effect of other features. The rest of features are fixed!
#Just to clarify: In each interval, we have houses with similar values for the feature of interest, i.e., area. We will only use these houses, their area and the corresponding number of rooms, which will be also similar. Therefore, we are selecting real combinations of the two features. Now, we change the area for each of the datapoints inside this intervals and see what happens with the prediction. It only make two changes, the lowest and the highest area of all points in the interval. We change the value of all sample to the lowest. Then, repeate with the highest value obtaining a new set of predictions. Calculate the difference between the lowest and the highest prediction for each sample, and average across all samples within the interval.
#These samples will also have similar values of the correlated second feature, i.e., number of rooms. BUT we ae not calculating the average prediction in these samples, instead, we are checking what happens if, within samples with similar area and hence similar number of rooms, we increase and decrease the area. So the impact on the prediction is only caused by the change in area, not by the number of rooms. In other words, we are actually modifying the value of the feature of interest, area, within a set of samples that have similar number of rooms and we are just looking at the average price in these data points. Therefore, the prediction change is caused only by area, as we have similar rooms for the analyzed datapoints. In addition, we are only looking at combinations of area and rooms that are observed in the data. In summary, we can see how THE MODEL HAS LEARNED to use area in order to predict the price, independently of the number of rooms. This is exacly what we want. Figure 8.7 shows this.
#But is it not interesting to see that our model behaves oddly at x1 > 0.7 and x2 < 0.3? Well, yes and no. Since these are data instances that might be physically impossible or at least extremely unlikely, it is usually irrelevant to look into these instances. But if you suspect that your test distribution might be slightly different and some instances are actually in that range, then it would be interesting to include this area in the calculation of feature effects. But it has to be a conscious decision to include areas where we have not observed data yet and it should not be a side-effect of the method of choice like PDP. If you suspect that the model will later be used with differently distributed data, I recommend to use ALE plots and simulate the distribution of data you are expecting.
#first-order effects
    #Therefore, the ALE value can be interpreted as the main effect of the feature at a certain value compared to the AVERAGE PREDICTION OF THE DATA. This interpretation is made possible because ALE plots are centered at zero so each point of the ALE curve represents the difference with mean prediction.[4] For instance, an ALE value of 5 for a feature value 12.2 would mean that if the feature of interest equals 12.2, the prediction it would yield is higher by 5 than the average prediction.
        #https://towardsdatascience.com/explainable-ai-xai-methods-part-3-accumulated-local-effects-ale-cf6ba3387fde
    #The number of instances is the same in ALL INTERVALS, thanks to the use of quantiles as grid, but in some areas there will be many short intervals and the ALE curve will consist of many more estimates. But for long intervals, which can make up a big part of the entire curve, there are comparatively fewer instances. This happened in the cervical cancer prediction ALE plot for high age for example.
    #ALE plots can become a bit shaky (many small ups and downs) with a high number of intervals. In this case, reducing the number of intervals makes the estimates more stable, but also smoothes out and hides some of the true complexity of the prediction model. There is no perfect solution for setting the number of intervals. If the number is too small, the ALE plots might not be very accurate. If the number is too high, the curve can become shaky.
    #Second order effects
        #The second-order effect is the additional interaction effect of the features after we have accounted for the main effects of the features.
        #Second-order effect plots can be a bit annoying to interpret, as you always have to keep the main effects in mind. It is tempting to read the heat maps as the total effect of the two features, but it is only the additional effect of the interaction. The pure second-order effect is interesting for discovering and exploring interactions, but for interpreting what the effect looks like, I think it makes more sense to integrate the main effects into the plot
#Important note about interpretation
    #The main interpretation of the ALE plot is qualitative—fixing the feature value and looking at the ALE plot as a function at that point, the tangent at that point (or THE SLOPE OF LINEAR INTERPOLATION BETWEEN THE CLOSEST BIN ENDPOINTS) shows how sensitive the target prediction is with respect to small changes of the feature value. Since we have a linear regression model, the tangent/slope is the same across the whole feature range so the feature sensitivity is identical at any point in the feature range.Using this, we can see the average effect for datapoints close to the feature value of interest (inside the bin).
    #The ALE value at a point is the relative feature effect with respect to the mean feature effect. The interpretation of the ALE plot is that, given a feature value, the ALE value corresponding to that feature value is the difference to the mean effect of that feature. Put differently, the ALE value is the relative feature effect on the prediction at that feature value.
    #Comparing the ALE plots of multiple models on the same axis should be done with care. In general, we can only make qualitative comparisons of the plots between different intervals of the feature values as we have done here (see example below).
        #https://docs.seldon.io/projects/alibi/en/stable/examples/ale_regression_california.html
    #an example of alibi
        #We can also say a few more quantitative things about this plot (https://docs.seldon.io/projects/alibi/en/stable/examples/ale_regression_california.html). The ALE value for the point MedInc(median income)=6 ($60,000) is ~1 which has the interpretation that for areas with this median income the model predicts an up-lift of ~$100,000 with respect to the average effect of MedInc (in y axis, 1 is 100,000). This is because the ALE plots are centered such that the average effect of the feature across the whole range of it is zero.
        #On the other hand, for neighbourhoods with MedInc=4 ($40,000), the ALE value is ~0 which indicates that the effect of the feature at this point is the same as the average effect of the feature. For even lower values of MedInc, below $40,000, the feature effect becomes less than the average effect, i.e. a smaller median income in the area brings the predicted house value down with respect to the average feature effect.
        #An additional feature of the ALE plot is that it shows feature deciles on the x-axis (https://docs.seldon.io/projects/alibi/en/stable/_images/examples_ale_regression_california_44_0.png). This helps understand in which regions there is low data density so the ALE plot is interpolating. For example, for the AveOccup feature (average number of household members), there appears to be an outlier in the data at over ~1,200 which causes the plot to linearly interpolate over a large range where there is no data. Note that this can also be seen by the lack of markers on the plot within that large range.
    #But there is a PROBLEM WITH CORRELATED FEATURES:
        #ALE plots do have their limitations as well. It is trickier to interpret than PDPs or ICE plots. The interpretation also is usually limited within the “window” or “interval” defined. If the features are highly correlated (which is often the case), interpreting an effect across intervals is impossible. Keep in mind that the effects are calculated within each interval and that each interval contains different sets of data points. These calculated effects are accumulated (hence the name “Accumulated” Local Effects) only for the sake smoothing the line.
            #info from Christopher
    #I do not fully understand this point. 
        #If within a given interval, we have points with the similar values of the correlated feature, and we change the values of the interest feature, then I am seeing the impact on the predictions of the first feature, not the second one (in contrast with M-plots; see above). 
        #If within a given interval I can quantify the impact on the predictions of my feature of interest, I could go interval by interval checking if higher values of my feature consistently do the same to the prediction.
        #I will do as alibi says, just look at the interpolation line between the interval of interest and the closest interval, to see how sensitive the target prediction is with respect to small changes of the feature value in that interval.
#Summary. To summarize how each type of plot (PDP, M, ALE) calculates the effect of a feature at a certain grid value v:
    #Partial Dependence Plots: “Let me show you what the model predicts on average when each data instance has the value v for that feature (i.e., all samples have v for the feature). I ignore whether the value v makes sense for all data instances.”
    #M-Plots: “Let me show you what the model predicts on average for data instances that have values close to v for that feature. The effect could be due to that feature, but also due to correlated features.”
    #ALE plots: “Let me show you how the model predictions change in a small ”window” of the feature around v for data instances in that window.”
#training vs test
    #At this point, I do not see any reason to not use training+test set for ALE plots. Remember that is mandatory to use the test set if the approach relies on error measurements, which is NOT the case for ALE plots.
    #Alibi docs say that X is "An N x F tabular dataset used to calculate the ALE curves. This is typically the training dataset or a representative sample"
    #See this thread about using training or set for shap values.
        #https://datascience.stackexchange.com/a/97522
        #They say that if you just want model explanations (I guess how the model has learnt to predict the response), you can just refit to the whole dataset, being this the final model used in "production".
        #In constrast, if you are interested in data explanations, you have to check how good proxy are model explanations for data explanations. In that case, you need "globally how importance is attributed to each feature, and how variable those importances are". 
            #We do that with permutation importance across different splits of the test set. We get how important is a feature for the model in order to generalize and how this changes as the data changes.
            #Remember we have done this in test set, not in the whole dataset, and making 50 random splits of the test set for each feature.
    #I think we can just use the whole dataset for ALE plots in order to get the maximum amount of information to learn how the model uses the features to predict flex sweep.
print_text("Load the final model fit to the whole dataset (training+set)", header=2)
import joblib
final_model_full_data=joblib.load("./results/selected_model_class/yoruba_hg19_stack_5_xgboost_fit_to_full_dataset.pkl.gz")
print(final_model_full_data)


print_text("define dict with nice axes names for ALE", header=2)
dict_nice_feature_names_ale={
    "chip_density_1000kb": "Regulatory density (ChIP-seq)",
    "chip_immune_1000kb": "Immune regulatory density (ChIP-seq)",
    "chip_testis_1000kb": "Testis regulatory density (ChIP-seq)",
    "coding_density_1000kb": "Coding density",
    "cons_elements_1000kb": "Density of conserved elements",
    "gc_content_1000kb": "GC-content (%)",
    "expression_all_tissues": "Gene expression across multiple tissues (log TPM)",
    "testis_expression": "Gene expression in testis (log TPM)",
    "lympho_expression": "Gene expression in lymphocytes (log TPM)",
    "gene_length": "Gene length (pair base)",
    "gene_number_1000kb": "Gene number",
    "n_dis_1000kb": "Number of Mendelian disease variants",
    "pip_v3": "Number of protein-protein interactions (log interactions)",
    "recombination_final_1000kb": "Recombination rate (cM/Mb)",
    "tbfs_density_1000kb": "Regulatory density (DNAseI)",
    "thermogenic_distance": "Distance to thermogenic genes (pb)",
    "vip_distance": "Distance to VIPs (pb)",
    "bat_distance_percentile_1": "Distance to BAT genes (pb)",
    "smt_distance_percentile_1": "Distance SMT genes (pb)"}
print("check that the keys in the dict are the features in the same order than in the data")
predictors_ale=[i for i in modeling_data.columns if i not in ["prob(sweep)"]]
predictors_nice_ale=[dict_nice_feature_names_ale[i] for i in modeling_data.columns if i not in ["prob(sweep)"]]
print(predictors_nice_ale)



print_text("initialize ALE plots using alibi", header=2)
#alibi is much more maintained than aleplot and it is much easier to install in container
#we have already used this package for permutation importance
    #https://docs.seldon.io/projects/alibi/en/stable/index.html
from alibi.explainers import ALE
#we have a warning when importing alibi: "NumbaDeprecationWarning: The 'nopython' keyword argument was not supplied to the 'numba.jit' decorator."
    #It seems that this is ok and it is related to shap. This is going to be solved in new versions of shap and we are not using shap anyways.
        #https://github.com/slundberg/shap/issues/2909
lr_ale = ALE( \
    predictor=final_model_full_data.predict, \
    feature_names=predictors_ale, \
    target_names=None, \
    check_feature_resolution=True, \
    low_resolution_threshold=10, \
    extrapolate_constant=True, \
    extrapolate_constant_perc=10.0, \
    extrapolate_constant_min=0.1)
    #Accumulated Local Effects for tabular datasets. Current implementation supports first order feature effects of numerical features.
        #https://docs.seldon.io/projects/alibi/en/stable/api/alibi.explainers.ale.html#alibi.explainers.ale.ALE
        #predictor
            #A callable that takes in an N x F array as input and outputs an N x T array (N - number of data points, F - number of features, T - number of outputs/targets (e.g. 1 for single output regression, >=2 for classification)).
        #feature_names=None
            #A list of feature names used for displaying results.
        #target_names=None
            #A list of target/output names used for displaying results.
        #check_feature_resolution=True
            #If True, the number of unique values is calculated for each feature and if it is less than low_resolution_threshold then the feature values are used for grid-points instead of quantiles. This may increase the runtime of the algorithm for large datasets. Only used for features without custom grid-points specified in alibi.explainers.ale.ALE.explain().
        #low_resolution_threshold=10 
            #If a feature has at most this many unique values, these are used as the grid points instead of quantiles. This is to avoid situations when the quantile algorithm returns quantiles between discrete values which can result in jumps in the ALE plot obscuring the true effect. Only used if check_feature_resolution is True and for features without custom grid-points specified in alibi.explainers.ale.ALE.explain().
            #this and the previous argument are used to determine whether to use the quantile algorithm or just use the unique values of the feature as grid points.
                #I do not fully understand the quantile algorithm, but given we have over 15K points, I think we will always have enough data points for each feature.
                #we will use the default.
        #extrapolate_constant=True
            #If a feature is constant, only one quantile exists where all the data points lie. In this case the ALE value at that point is zero, however this may be misleading if the feature does have an effect on the model. If this parameter is set to True, the ALE values are calculated on an interval surrounding the constant value. The interval length is controlled by the extrapolate_constant_perc and extrapolate_constant_min arguments.
                #I guess you use values of the feature around its constant value to see how the prediction change, so we have points outside zero if the feature has effect.
        #extrapolate_constant_perc
            #Percentage by which to extrapolate a constant feature value to create an interval for ALE calculation. If q is the constant feature value, creates an interval [q - q/extrapolate_constant_perc, q + q/extrapolate_constant_perc] for which ALE is calculated. Only relevant if extrapolate_constant is set to True.
        #extrapolate_constant_min
            #Controls the minimum extrapolation length for constant features. An interval constructed for constant features is guaranteed to be 2 x extrapolate_constant_min wide centered on the feature value. This allows for capturing model behaviour around constant features which have small value so that extrapolate_constant_perc is not so helpful. Only relevant if extrapolate_constant is set to True.



print_text("calculate ALE curves", header=2)
lr_exp = lr_ale.explain( \
    X=modeling_data[predictors].to_numpy(), \
    features=None, \
    min_bin_points=4, \
    grid_points=None) 
    #Calculate the ALE curves for each feature with respect to the dataset X.
        #X
            #An N x F tabular dataset used to calculate the ALE curves. This is typically the training dataset or a representative sample.
        #features
            #Features for which to calculate ALE.
        #min_bin_points
            #Minimum number of points each discretized interval should contain to ensure more precise ALE estimation. Only relevant for adaptive grid points (i.e., features without an entry in the grid_points dictionary).
            #This argument determines the number of bins the range of each feature is subdivided into so that the ALE estimate for each bin is made with at least min_bin_points. Smaller values can result in less accurate local estimates while larger values can also result in less accurate estimates by averaging across large parts of the feature range.
            #In general, the ALE for a non-linear model doesn’t have to be monotonic, although small departures from monotonicity which may be due to artifacts from the grid-size used to calculate the ALE. It may be useful to experiment with different resolutions of the grid size
            #as you increase this number, you get less intervals because more points are required per interval.
            #According to Christopher Molnar: "ALE plots can become a bit shaky (many small ups and downs) with a high number of intervals. In this case, reducing the number of intervals makes the estimates more stable, but also smoothes out and hides some of the true complexity of the prediction model. There is no perfect solution for setting the number of intervals. If the number is too small, the ALE plots might not be very accurate. If the number is too high, the curve can become shaky."
            #so we can check different number of minimum number of points per interval and hence different number of intervals and see if the pattern is clear.
        #grid_points
            #Custom grid points. Must be a dict where the keys are features indices and the values are monotonically increasing numpy arrays defining the grid points for each feature. See the Notes section for the default behavior when potential edge-cases arise when using grid-points. If no grid points are specified (i.e. the feature is missing from the grid_points dictionary), deciles discretization is used instead.



print_text("save the explanations into pickle", header=2)
import pickle
pickle.dump(lr_exp, open("./results/selected_model_class/ale_plots/saved_explained/yoruba_hg19_stack_5_xgboost_fit_to_full_dataset_ale_alibi_explanations.sav", 'wb'))
#import pickle; lr_exp = pickle.load(open("./results/selected_model_class/ale_plots/saved_explained/yoruba_hg19_stack_5_xgboost_fit_to_full_dataset_ale_alibi_explanations.sav", 'rb'))
    #https://machinelearningmastery.com/save-load-machine-learning-models-python-scikit-learn/



print_text("make ale plots", header=2)
import matplotlib.pyplot as plt
#Plotting ale_values against feature_values recovers the ALE curves. For convenience we include a plotting function plot_ale which automatically produces ALE plots using matplotlib:
#Note that the ALE is estimated for each interval edge and linearly interpolated in between, for real applications it is important to have a sufficiently fine grid but also one that has enough points into each interval for accurate estimates. 
#Important
    #The x-axis also shows feature deciles of the feature to help judge in which parts of the feature space the ALE plot is interpolating more and the estimate might be less trustworthy.
    #Therefore, first bar is percentile 10%, second bar is percentile 20%, and so on...
#feature=[feature for (pos,feature) in enumerate(dict_nice_feature_names_ale.keys()) if pos==4][0]
from alibi.explainers import plot_ale
for feature in dict_nice_feature_names_ale.keys():
    ax=plot_ale( \
        exp=lr_exp, \
        features=np.where(np.array(predictors_ale) == feature)[0].tolist(), \
        targets="all", \
        n_cols=1, \
        sharey="all", \
        constant=False, \
        ax=None,
        fig_kw={'figwidth':6, 'figheight':3}, \
        line_kw={"linewidth":0.2, "markersize":2, "markeredgewidth":0.25,"fillstyle":"none"})[0][0]
        #you can change arguments of plot using fig and line_kw
            #look posibilities in "plt.rcParams.keys()"
    ax.update({"xlabel": dict_nice_feature_names_ale[feature], "ylabel": "Change in prediction"})
        #https://www.geeksforgeeks.org/matplotlib-axes-axes-update-in-python/
        #The ALE on the y-axes of the plot above is in the units of the prediction variable.
    ax.xaxis.label.set_size(12)
    ax.yaxis.label.set_size(13.5)
    ax.get_legend().remove()
        #https://stackoverflow.com/questions/59352887/how-to-remove-legend-from-an-image-plot
    plt.savefig( \
        fname="./results/selected_model_class/ale_plots/aleplot_" + feature + ".png", dpi=300)
    plt.close()
        #exp
            #An Explanation object produced by a call to the alibi.explainers.ale.ALE.explain() method.
            #features
                #A list of features for which to plot the ALE curves or 'all' for all features. Can be a mix of integers denoting feature index or strings denoting entries in exp.feature_names. Defaults to 'all'.
        #targets
            #A list of targets for which to plot the ALE curves or 'all' for all targets. Can be a mix of integers denoting target index or strings denoting entries in exp.target_names. Defaults to 'all'.
        #n_cols
            #Number of columns to organize the resulting plot into.
        #sharey
            #A parameter specifying whether the y-axis of the ALE curves should be on the same scale for several features. Possible values are: 'all' | 'row' | None.
        #constant
            #A parameter specifying whether the constant zeroth order effects should be added to the ALE first order effects.
        #ax
            #A matplotlib axes object or a numpy array of matplotlib axes to plot on.
        #line_kw
            #Keyword arguments passed to the plt.plot function.
        #fig_kw
            #Keyword arguments passed to the fig.set function.
#Results
#the ALe plot of VIP makes sense, the prediction gets below the average as we move far wawy from VIPs, then increase but when the data starts to get very disperse
#The ALE plot of BAT shows that the predictions tend to be above average close to BAT genes. There is dispersion, but for most of the points is the case. Note that most of the data is in this part, supporting this pattern. Once we get away from BAT genes, the predictions tend to be in the average or below average. Then, as data is more disperse, we get more noisy pattern around the average. The lowest value is observed just before the percentile 90%. After that percentile, data is MUCH more disperse.
    #Yes, the pattern is noisy, but check the pattern of testis regulatory density. There is a decrease of selection but very noisy even the permutation importance is also not low. We know from our previous model that this variable is strongly and negatively associated with iHS, likely caused by recombination fine scale patterns. We are detecting here the same effect but with more noise. 
    #therefore, the pattern of BAT can be still legit.
#For thermogenic and SMT genes there is no clear pattern, predictions tend to be around the average. In the case of thermogenic, there is some pattern at the beginning but we have much less accumulation of points compared to the case of BAT genes.




print_text("plot flex-sweep distribution", header=1)
import matplotlib.pyplot as plt
plt.hist(modeling_data["prob(sweep)"], bins=50, color="blue", alpha=0.4, label="Observed log(Flex-sweep probability)")
plt.legend(loc="upper left")
plt.savefig( \
    fname="./results/selected_model_class/yoruba_hg19_observed_log_flex_sweep.png", dpi=300)
plt.close()



'''
#shap for the future

#####tree shap

#you cannot you stacking with tree shap!!
    #you could check what they did in the deep learning paper


        #The short answer to your question is yes, if you are taking the mean of the 10 XGBoost model outputs (margin outputs), then you can average the 10 SHAP values matrices and get the right answer. Just make sure you are averaging the margin outputs and not the probabilities if it’s for classification
            #https://github.com/slundberg/shap/issues/112

#only feature importance?
    #we need feature importance not based on error! so we can show BAT genes are important


#tree shap also avoids problems with correlated features, but it is only for tree-based methods, this is not the case for the other SHAP methods

#we can use treeshap that does not include unlikely instances!!!
#Feature dependencies can skew the approximations made by KernelSHAP. The algorithm estimates SHAP values by randomly sampling feature values. The issue is that, when features are correlated, the sampled values can be unlikely. This means that when using SHAP values we can put too much weight on unlikely observations.
#TreeSHAP does not have the problem. However, the algorithm has a different issue due to feature dependencies. That is a feature that has no impact on the prediction can get a SHAP value that is non-zero. This can happen when the feature is correlated with another feature that does impact the prediction. In this case, we can make the incorrect conclusion that a feature has contributed to a prediction.
    #that is ok, if a feature is correlated with recombination, I expect that features to show a signal, but in the case of GC, the positive effect would remain?
    #https://christophm.github.io/interpretable-ml-book/shap.html
    #https://towardsdatascience.com/kernelshap-vs-treeshap-e00f3b3a27db


#feature importance with shap is based on error? if that is the case, then we should think about training/test
    #https://datascience.stackexchange.com/questions/62913/shap-explanations-in-case-of-repeated-train-test-split
    #https://datascience.stackexchange.com/questions/61395/shap-value-analysis-gives-different-feature-importance-on-train-and-test-set


#check if shaps gives effect after controling for the rest of factors.
    #Note that GC content is also negatively associated with selection accoridng to ALE plots!!

import shap

#https://mail.google.com/mail/u/0/?tab=rm&ogbl#inbox/FMfcgzGtwCtShmTtxhBvnFxjKGtNtvbh
explainer = shap.TreeExplainer(final_model.named_steps['regressor'], feature_names=train[predictors].columns)
    #errors numpy
        #module 'numpy' has no attribute 'int'.
        #module 'numpy' has no attribute 'bool'.
        #https://stackoverflow.com/questions/74946845/attributeerror-module-numpy-has-no-attribute-int
    #solution
        #np.int=int
        #np.bool=bool
        #https://github.com/WongKinYiu/yolov7/issues/1280

shaps = explainer(final_model.named_steps['scale'].fit_transform(train[predictors]))
    #https://stackoverflow.com/questions/55867862/how-to-use-shap-with-a-linear-svc-model-from-sklearn-using-pipeline
#shap.plots.waterfall(shaps[0])

    #we are using scaling outside the pipeline, so the axes are NOT GOING TO BE RAW



shap.summary_plot(shaps, train[predictors])
import matplotlib.pyplot as plt
plt.savefig( \
    fname="./eso4.png")
plt.close()

shap.summary_plot(shaps, train[predictors], plot_type="bar")
    #VIPs are very low

#By default a SHAP bar plot will take the mean absolute value of each feature over all the instances (rows) of the dataset.
shap.plots.bar(shaps)



#plot per feature and instance

shap.plots.scatter(shaps[:,"recombination_final_1000kb"], color=shaps)
shap.plots.scatter(shaps[:,"cons_elements_1000kb"], color=shaps[:,"recombination_final_1000kb"])
shap.plots.scatter(shaps[:,"gc_content_1000kb"], color=shaps[:, "recombination_final_1000kb"])
    #negative association, this makes sense. If a feature is correlated with a causal feature, it can get non-zero shapely values. I understand then that the impact of a feature is not cleaned from the impact of other features. 
    #I think this should not be the case for ale plot, but CHECK
    #maybe combine feature importante of shap with ale plots?
shap.plots.scatter(shaps[:,"chip_testis_1000kb"], color=shaps[:, "recombination_final_1000kb"])
shap.plots.scatter(shaps[:,"chip_immune_1000kb"], color=shaps[:, "recombination_final_1000kb"])
shap.plots.scatter(shaps[:,"vip_distance"], color=shaps[:, "recombination_final_1000kb"])
shap.plots.scatter(shaps[:,"bat_distance_percentile_1"], color=shaps[:, "recombination_final_1000kb"])
shap.plots.scatter(shaps[:,"smt_distance_percentile_1"], color=shaps[:, "recombination_final_1000kb"])
    #As we go away from BAT and SMT genes, the shapely values tend to zero or negative for genes with high recombination
    #In constras, close to BAT and SMT genes, we find positive shapely values even for high-recombination genes, suggesting that proximity to these genes has a positive impact on the probability of selection independently of recombination rate.
    #There is, of course, noise with many genes being close to BAT/SMT genes but still having non-positive shapely values. So we need the bootstrap test.
'''
    