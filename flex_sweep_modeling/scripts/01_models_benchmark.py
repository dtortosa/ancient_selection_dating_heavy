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
os.environ['PYTHONHASHSEED']=str(seed_value)
print(os.environ['PYTHONHASHSEED'])


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
        ./results/model_comparison; \
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
final_data_yoruba_raw = pd.read_csv( \
    "data/flex_sweep_closest_window_center.txt.gz", \
    sep=",", \
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
if(n_nans_bat <= (final_data_yoruba.shape[0]*0.03)) | (n_nans_smt <= (final_data_yoruba.shape[0]*0.03)):
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





###por aquiii





print_text("define function to run across folds and model classes within a nested CV schema", header=1)


#we are going to use nested Cross-Validation, so we can tune hyperparameters and select the best model class without having too optismistic results about model performance. If the same dataset is used to tune hyperparameters and then evaluate the tuned models to compare between model classes, there is a risk of overfitting.
    #https://machinelearningmastery.com/nested-cross-validation-for-machine-learning-with-python/




##check this explanations that are from previous versions


#Held out part of the dataset. This part will not be used for parameter optimization, but for the final validation after the final model has been optimized. In this way, we avoid potential overfittin in the evaluation sets ([see link](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation)). If you use train in a set of the data and then evaluate in other, you can see if the model trained is not overfitting to the training data and flexible enough to predict the evaluation data. The problem is that in parameter optimization, we select the best parameters based on the evaluation metrics in the evaluation sets, thus we could get a model that fit too much the evaluation dataset, loosing generalization and thus making the evaluation metrics no longer metrics of generalization. To avoid this, we leave out a set of the data for final evaluation. This set will be NOT used in parameter optimization. Once we have selected the best parameters to train a model in the training dataset and predict well in the evaluation dataset, we use these parameter to create a final model, fit to the whole training data and then predict in the final evaluation dataset, which was not used for anything before. If the model works well, it means it is generalizable and it is not overfitting the data, so, in our case, we can say that it is explaining the variance in selection that it is really explained by the genomic factors and the rest would be variance that could be explained by selective pressures. If there is overfitting, the model fit too much the data, there is not non-explained variance.

#In the future, we may want to use the final model to obtain a probability of selection considering genomic factors and then use it to select genes with the same expected probability of selection (according to these factors) than our genes of interest. The idea is that the interest genes should have the same probability of selection based on genomic features (predicted probability), but if they are target of a selective pressure, their observed probability of selection should be higher. If predicted and observed are exactly the same, there is no room for enrichment of selection in interest genes after controling for confounding factors, and this would be an methodological artifact due to a model that just fit the observed data without any generalization. This can be a problem with algorithms like RF or in deep learning. In other words, more overfitting, less power to detect the impact of selective pressures.

#It is usually recommended a 80-20 ratio for training-test when the dataset is large, but for small datasets is better 70-30 ([link](https://www.researchgate.net/post/Is-there-an-ideal-ratio-between-a-training-set-and-validation-set-Which-trade-off-would-you-suggest)). 

#Note that the test dataset gives information about the generalization ability of the model, and this is important to us as a measure of overfitting. We need to leave enough amount of data to check this but without leaving the training dataset so small that we underfit.

#Because of this, we are going to use 70-30 even though we do not have a specially small dataset. We want to ensure our validation set covers well the variability found in the genome, it should be a difficult problem to the model so we can have more confidence we are not overfitting.

#Update: We are having problems to get good R2 with the probability closest at the center of the window. I am going to increase the sample size of the training set so maybe we can get better models, we have a relatively large sample size, so we can use 80-20.



print_text("prepare models, HPs and CV scheme", header=2)
print_text("define a dict with model instances and grids of HPs", header=3)
#define a dict
    #keys:
        #the name of the model
    #Values:
        #the code for creating an instance of the corresponding model with no arguments.
        #We will use "eval()" to run this line of code and create a new instance of the corresponding model
        #Note that the reproducibility is ensured as we have set the seeds of python, numpy and tensorflow before.
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
import xgboost as xgb
from tensorflow import keras
dict_models = { \
    "elastic_net": {"instance": "ElasticNet()"}, \
    "random_forest": {"instance": "RandomForestRegressor()"}, \
    "xgboost": {"instance": "xgb.XGBRegressor()"}, \
    "neural_nets": {"instance": "keras.Sequential()"}}
dict_models


print_text("add HPs for gridsearch", header=3)
#the regressor will be part of a pipeline, so we need to specify the name of the step in the pipeline (i.e., regressor) before the name of the hyperparameter
print_text("elastic net", header=4)
#general notes about elastic nets
    #Linear regression is the standard algorithm for regression that assumes a linear relationship between inputs and the target variable. An extension to linear regression involves adding penalties to the loss function during training that encourage simpler models that have smaller coefficient values. These extensions are referred to as regularized linear regression or penalized linear regression.
    #Elastic net is a popular type of regularized linear regression that combines two popular penalties, specifically the L1 and L2 penalty functions.
    #Linear regression refers to a model that assumes a linear relationship between input variables and the target variable. With a single input variable, this relationship is a line, and with higher dimensions, this relationship can be thought of as a hyperplane that connects the input variables to the target variable. The coefficients of the model are found via an optimization process that seeks to minimize the sum squared error between the predictions (yhat) and the expected target values (y).
    #A problem with linear regression is that estimated coefficients of the model can become large, making the model sensitive to inputs and possibly unstable. This is particularly true for problems with few observations (samples) or more samples (n) than input predictors (p) or variables (so-called p >> n problems).
    #One approach to addressing the stability of regression models is to change the loss function to include additional costs for a model that has large coefficients. Linear regression models that use these modified loss functions during training are referred to collectively as penalized linear regression.
    #One popular penalty is to penalize a model based on the sum of the squared coefficient values. This is called an L2 penalty. An L2 penalty minimizes the size of all coefficients, ALTHOUGH IT PREVENTS ANY COEFFICIENTS FROM BEING REMOVED FROM THE MODEL.
    #Another popular penalty is to penalize a model based on the sum of the absolute coefficient values. This is called the L1 penalty. An L1 penalty minimizes the size of all coefficients and allows some coefficients to be minimized to the value zero, WHICH REMOVES THE PREDICTOR FROM THE MODEL.
    #Elastic net is a penalized linear regression model that includes both the L1 and L2 penalties during training. The benefit is that elastic net allows a balance of both penalties, which can result in better performance than a model with either one or the other penalty on some problems.
        #https://machinelearningmastery.com/elastic-net-regression-in-python/
    #MY OPINION
        #with elastic net, we can check a wide range in the degree of penalties from maximum penalty to no penalty (usual linear model). We include the penalty of Lasso and Ridge, so we can covering the different options in out GridSearch.
            #The main difference between Lasso and Ridge is the penalty term they use. Ridge uses L2 penalty term which limits the size of the coefficient vector. Lasso uses L1 penalty which imposes sparsity among the coefficients and thus, makes the fitted model more interpretable. Elasticnet is introduced as a compromise between these two techniques, and has a penalty which is a mix of L1 and L2 norms.
            #https://stats.stackexchange.com/a/93195
dict_models["elastic_net"]["HPs"] = { \
    "regressor__l1_ratio": np.arange(0, 1, 0.01),
    "regressor__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]}
    #HPs
        #l1_ratio
            #alpha in elastic net theory (different from alpha in scikit learn) determines how much weight is given to each of the L1 and L2 penalties (see above). Alpha is a value between 0 and 1 and is used to weight the contribution of the L1 penalty and one minus the alpha value is used to weight the L2 penalty.
            #For example, an alpha of 0.5 would provide a 50 percent contribution of each penalty to the loss function. An alpha value of 0 gives all weight to the L2 penalty and a value of 1 gives all weight to the L1 penalty.
            #the alpha hyperparameter can be set via the “l1_ratio” argument that controls the contribution of the L1 and L2 penalties 
        #alpha
            #Another hyperparameter is provided called “lambda” (alpha in scikit learn) that controls the weighting of the sum of both penalties to the loss function. A default value of 1.0 is used to use the fully weighted penalty; a value of 0 excludes the penalty. Very small values of lambada, such as 1e-3 or smaller, are common.
            #In scikit learn, “alpha” argument controls the contribution of the sum of both penalties to the loss function
            #``alpha = 0`` is equivalent to an ordinary least square, solved by the :class:`LinearRegression` object
                #Therefore, typical linear model
        #Confusingly, as you can see the notation is different in scikit learn, as indicated in the help:
            #The parameter l1_ratio corresponds to alpha in the glmnet R package while alpha corresponds to the lambda parameter in glmnet.
    #grid
        #Janizek 2023
            #the ‘alpha’ parameter was tuned over values ranging from 0.1 to 100, whereas the ‘l1_ratio’ parameter was tuned from 0.25 to 0.75.
        #ML mastery
            #One approach would be to gird search l1_ratio values between 0 and 1 with a 0.1 or 0.01 separation and alpha values from perhaps 1e-5 to 100 on a log-10 scale and discover what works best for a dataset.
        #The grid of ML Mastery is wider but it takes less than 1.5 hours to run across all folds, so we should use that.
    #we get warning about not converging. ML Mastery says that is ok.



#check other parameters of elastic net



print_text("random forest", header=4)
#general notes about random forest
dict_models["random_forest"]["HPs"] = { \
    "regressor__l1_ratio": np.arange(0, 1, 0.01),
    "regressor__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]}





#X = modeling_data.iloc[:, 1:]

#think about the number of K, usually K in outer cv is larger

from sklearn.model_selection import KFold
cv_outer = KFold( \
    n_splits=10,  \
    shuffle=True)
    #random_state is not required as we have already ensured reproducibility by setting the seeds of python, numpy and tensorflow
print(cv_outer)




print_text("define a dict with model instances and grids of HPs", header=2)


indexes_cv_outer = []
#train_index, test_index = [(train_index, test_index) for fold, (train_index, test_index) in enumerate(cv_outer.split(modeling_data)) if fold==0]
    #select the first fold for debugging
for fold, (train_index, test_index) in enumerate(cv_outer.split(modeling_data)):

    for model_class in dict_models.keys():

        indexes_cv_outer.append((fold, train_index, test_index, model_class))
print(indexes_cv_outer)


len(indexes_cv_outer) == cv_outer.n_splits * len(dict_models.keys())



print_text("Define function and run it only for elastic net", header=2)
from sklearn.pipeline import Pipeline
from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score
from sklearn.base import clone
#fold, train_index, test_index, model_class = indexes_cv_outer[0]
def model_evaluation(fold, train_index, test_index, model_class):

    print_text(f"Starting with fold: {fold}, model: {model_class}", header=3)
    print_text("split the data", header=4)
    X_train, y_train = modeling_data.iloc[train_index, 1:], modeling_data.iloc[train_index, 0]
    X_test, y_test = modeling_data.iloc[test_index, 1:], modeling_data.iloc[test_index, 0]


    print_text("check correct shape datasets", header=4)
    print((X_train.shape[0] == len(train_index)) & (y_train.shape[0] == len(train_index)))
    print((X_test.shape[0] == len(test_index)) & (y_test.shape[0] == len(test_index)))


    print_text("set the inner CV", header=4)
    cv_inner = KFold( \
        n_splits=3,  \
        shuffle=True)
        #random_state is not required as we have already ensured reproducibility by setting the seeds of python, numpy and tensorflow
    print(cv_inner)


    print_text("open instance of the model", header=4)
    regressor = eval(dict_models[model_class]["instance"])
        #we have a dict the code to create a new instance of the model between '""', so we can evaluate it
    print(regressor)


    print_text("define pipeline with scaling and regressor", header=4)
    pipeline = Pipeline( \
        steps=[ \
        ('scale', preprocessing.StandardScaler()), \
        ('regressor', regressor)])
        #Pipeline:
            #Sequentially apply a list of transforms and a final estimator. Intermediate steps of the pipeline must be 'transforms', that is, they must implement `fit` and `transform` methods. The final estimator only needs to implement `fit`.
        #We obtain the model instance from the dict of models selecting the corresponding model class
    print(pipeline)


    print_text("extract dict with HPs for the GridSearch", header=4)
    space = dict_models[model_class]["HPs"]
    print(space)
    

    print_text("define the GridSearch", header=4)
    search = GridSearchCV( \
        estimator=pipeline, \
        param_grid=space, \
        scoring="r2", \
        n_jobs=1 , \
        cv=cv_inner, \
        refit=True, \
        pre_dispatch="2*n_jobs") #if n_jobs>1, then pre_dispatch should be "1*n_jobs", if not, the dataset will be copied many times increasing a lot memory usage
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
                #Refit an estimator using the best found parameters on the WHOLE dataset
            #cv:
                #Determines the cross-validation splitting strategy.
            #pre_dispatch:
                #Controls the number of jobs that get dispatched during parallel execution. Reducing this number can be useful to avoid an explosion of memory consumption when more jobs get dispatched than CPUs can process.
                #If `n_jobs` was set to a value higher than one, the data is copied for each point in the grid (and not `n_jobs` times). This is done for efficiency reasons if individual jobs take very little time, but may raise errors if the dataset is large and not enough memory is available.  A workaround in this case is to set `pre_dispatch`. Then, the memory is copied only `pre_dispatch` many times. A reasonable value for `pre_dispatch` is `2 * n_jobs`.


    print_text("run the GridSearch", header=4)
    search_results = search.fit(X_train, y_train)
    print(search_results)


    print_text("get the parameters of the best model", header=4)
    best_params = search_results.best_params_
        #Parameter setting that gave the best results on the hold out data
    print(best_params)


    print_text("extract the best regressor already fitted", header=4)
    best_model = search_results.best_estimator_
        #Estimator that was chosen by the search, i.e. estimator which gave highest score (or smallest loss if specified) on the left out data. Not available if ``refit=False``.
    print(best_model)


    print_text("predict on the test set", header=4)
    y_pred = best_model.predict(X_test)


    print_text("calculate evaluation score in the test set", header=4)
    score = r2_score(y_test, y_pred)


    print_text("return results as a tuple", header=4)
    tuple_results = (fold, model_class, best_params, score)
    print(tuple_results)
    return tuple_results

#run it across folds for elastic net only
[model_evaluation(*i) for i in indexes_cv_outer if i[3]=="elastic_net"]
    #"*" is used to unpack a tuple and use its elements as arguments
    #in our case, each tuple has as elements de fold, indexes and name of the model class
        #https://stackoverflow.com/a/1993732/12772630



print_text("parallelize the function", header=2)
print_text("open pool", header=3)
import multiprocessing as mp
pool = mp.Pool(5)
print(pool)


print_text("run the function using starmap, which is useful to apply function across iterable whose elements are in turn also iterables storing arguments (tuples in our case)", header=3)
results = pool.starmap(model_evaluation, indexes_cv_outer)
    #map: Apply `func` to each element in `iterable`, collecting the results in a list that is returned
    #starmap: Like `map()` method but the elements of the `iterable` are expected to be iterables as well and will be unpacked as arguments. Hence `func` and (a, b) becomes func(a, b).
        #this is our case, as we have a list of tuples. Therefore, each element of the iterable (list) is in turn another iterable (tuple).
        #elements of the tuple are unpacked as argument which is exactly what we want.
        #https://stackoverflow.com/a/47506842/12772630
        #https://discuss.python.org/t/differences-between-pool-map-pool-apply-and-pool-apply-async/6575
print(results)





pd.DataFrame(results)

    #you could do a sensitivity analysis for K
        #https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/


#gridsearch, random serach or optuna?
    #probably grdisearch with narrow list is the best option
#r2 or mse?
#apply scaling to prob(sweep) before like yifei?


from sklearn.model_selection import train_test_split
modeling_data_train, modeling_data_test = train_test_split(
    modeling_data, \
    test_size=0.20, \
    random_state=54, \
    shuffle=True)
    #train size is automatically calculated using "test_size" 
    #this function uses internally ShuffleSplit so it can randomly split a dataset in training and evaluation but instead of getting an interable (like in ShuffleSplit), you directly get the datasets splitted
        #https://stackoverflow.com/questions/66757902/differnce-between-train-test-split-and-stratifiedshufflesplit
print("see shapes of training and test sets")
print(modeling_data_train.shape)
print(modeling_data_test.shape)


print_text("Create create random splits for training and test sets", header=4)
#We are going to make 5 splits. Therefore, we will obtain the average score of 5 CV splits, to catch models that are overfitting. If we exposed the model to model different datasets, there are more possibilities to detect if the model has problems to generalize. 1/5 of 14000 is 1400 which I think it is ok.
n_folds = 5
sample_size_eval = modeling_data_train.shape[0]/n_folds
    #the training set is divided by the number of folds, and one of the folds is used for evaluation
if(sample_size_eval > 2000):
    print(f"We will use 5 folds. The size of the evaluation set in cross-validation is {sample_size_eval}")
else:
    raise ValueError("ERROR: FALSE! We have less than 2000 genes for evaluation in cross-validation")

#It is very important for us to avoid overfitting becuase we want to quantify the true impact of genomic factors on selection probability. If the model overfits, it will assume that genomic factors like recombination or GC-content explain much selection variability across the genome than they actually explain, reducing in this way the remaining variability that could be explained by other factors. 

#This is specially relevant in the context of botstrapping, where we compare genes related to a selective pressures with control genes that have a similar probability of selection according to confounding factors. The goal is to test whether interest genes have more selection than predicted, but if selection is already fully explained by confounding factors, there is no room for interest genes to be enriched or depleted in selection. They have the exact selection the confounding factors mark.

#In addition, it has been seen in many datasets that 10 folds works relatively well to estimate model performance ([link](https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/#:~:text=configure%20the%20procedure.-,Sensitivity%20Analysis%20for%20k,evaluate%20models%20is%20k%3D10.)). Although there is some debate about it. We will stick to 10 folds for now, if we have problems with the predictive power in the test set, we can re-evaluate.

#Update: I have problems to reach a R2 above 0.5 with deep nets, so I am using 5 fold to have more genes in each evaluation set.

#We will use KFold instead of ShuffleSplit because we want to have training and validation datasets that are not overlapped so we can use all the data. KFold with shuffle=True only shuflle one time at the beginning and then make the folds, while ShuffleSplit shuffles every time making possible to select as evaluation two times the same sample. In the case of 5 splits with KFold, one time 1/5 is the validation and the rest is for training, then another 1/5 and so on... until the 5 partitions have been used for validation, i.e., the whole dataset have been used for validation. 

#It is important to us to use the whole training dataset for CV because we need to obtain models that are generalizable enough so we do not overestimate the influence of genomic confounding factors on selection due to overfitting (see above).

#This is very well explained here(https://stackoverflow.com/questions/45969390/difference-between-stratifiedkfold-and-stratifiedshufflesplit-in-sklearn).

from sklearn.model_selection import KFold
shuffle_split = KFold( \
    n_splits=5,  \
    shuffle=True,  \
    random_state=61)
print(shuffle_split)


print_text("set also the number of jobs as only 2 due to memory usage errors (see below)", header=4)
number_jobs = 2
    #Using more jobs we get
    #"The exit codes of the workers are {SIGABRT(-6)}"
    #A lot of people is having the same problem even RAM seems to be ok. Some solve it by decreasing n_jobs and others by increase RAM usage.    
    #"Turns out allocating all CPUs can be unstable, specially when there are other independent programs running that can suddenly have an uncontrolled spike in memory usage."
    #it seems a problem of memory usage with scikit
        #https://github.com/scikit-learn-contrib/skope-rules/issues/18
    #I have detected that the decreasing the number of cores reduces the usage of memory even having the same number_jobs and optuna processes. Using for example 20 cores, I can run 10 optuna processes with number_jobs=2 using just 200GB per core. If I increase the number of cores to 40, things seems to seed up, but I get the workers error. So we are keeping things conservative with just 20 cores.


print_text("Convert to numpy array the training and test sets along with the whole set", header=4)
modeling_data_array = modeling_data.values
modeling_data_train_array = modeling_data_train.values
modeling_data_test_array = modeling_data_test.values
print(modeling_data_array)
print(modeling_data_train_array)
print(modeling_data_test_array)


print_text("Separate the response and predictors", header=4)
y = modeling_data_array[:, 0]
X = modeling_data_array[:, 1:]
y_train = modeling_data_train_array[:, 0]
X_train = modeling_data_train_array[:, 1:]
y_test = modeling_data_test_array[:, 0]
X_test = modeling_data_test_array[:, 1:]
print(f"see shape of the X and y training {(X_train, y_train)}")
print(f"see shape of the X and y test {(X_test, y_test)}")
print(f"Do we have the correct number of predictors in the X sets? {(X_train.shape[1] + 1 == modeling_data_train.shape[1], X_test.shape[1] + 1 == modeling_data_test.shape[1], X.shape[1] + 1 == modeling_data.shape[1])}")
print(f"Do we have the correct number of observations in the y sets? {(len(y_train)==modeling_data_train.shape[0], len(y_test)==modeling_data_test.shape[0], len(y)==modeling_data.shape[0])}")


print_text("Explanations about scaling", header=4)
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




print_text("comparison of multiple models", header=2)



##por aqui
#create function that run 1 optuna process per model and core??



#Create a pipe that includes the tranformer to do the scaling so it is fit (get mean and SD) to the training sets and not evaluation/validation sets. When using grid search ([link](https://stackoverflow.com/questions/51459406/how-to-apply-standardscaler-in-pipeline-in-scikit-learn-sklearn)) and cross_val_score ([link](https://stackoverflow.com/questions/44446501/how-to-standardize-data-with-sklearns-cross-val-score)), it uses the mean and SD of the training set of a given iteration and apply the transformation to the test set, then in the next iteration does the same but with the new training and test set.

#We need to first create a pipeline ([link](https://stackoverflow.com/questions/33091376/what-is-exactly-sklearn-pipeline-pipeline)) including scaling of predictors and the regressor. Then add it to the transformer of the response using TransformedTargetRegressor. We cannot just use StandardScaling on the response ([link](https://stackoverflow.com/questions/67824676/how-to-scale-both-x-and-y-data-in-sklearn-pipeline)).


from sklearn.pipeline import Pipeline
from sklearn import preprocessing

from sklearn.ensemble import RandomForestRegressor


#import xgboost as xgb
#model = xgb.XGBRegressor(random_state=23534, n_estimators=100, booster="gbtree")
#model_name="XGBoost"
    #booster: gbtree, gblinear or dart.

#from sklearn.ensemble import ExtraTreesRegressor
#model = ExtraTreesRegressor(random_state=23534, n_estimators=10)
#model_name = "extra_tree_regressor"

#from sklearn.ensemble import VotingRegressor, GradientBoostingRegressor
#from sklearn.linear_model import Ridge
#reg1 = GradientBoostingRegressor(random_state=1)
#reg2 = ExtraTreesRegressor(random_state=1)
#reg3 = Ridge()
#model = VotingRegressor(estimators=[('gb', reg1), ('et', reg2), ('rdg', reg3)])
#model_name = "voting_regressor"

model = RandomForestRegressor(random_state=23534, n_estimators=10)
model_name = "random_forest"
    #RF is 0.52 for VC and 0.86 for whole dataset

estim = Pipeline([ \
    ('scale', preprocessing.StandardScaler()), \
    ('regressor', model)])

#4 fives of models used by the paper: elastic net, RF, XGboost and DNNs
    #From my previous experience, RF and DNNs work well
    #XGboost is a new approach with very good results in regression
        #https://machinelearningmastery.com/extreme-gradient-boosting-ensemble-in-python/
    #elastic net is a linear approach that combines the strengths of ridge and lasso when controlling the impact of multiple features. We do not have maaany features (i.e., more features than samples), but still have a considerable amount of features.
        #https://machinelearningmastery.com/elastic-net-regression-in-python/
#5-fold CV
    #4 folds for training and 1 for test
        #this step is important because we ensure that the same data is not used for test twice. 
        #it is like we were repeating the same experiment 5 times, using 5 completely different test set
    #within training
        #3 folds for learning and 1 for validation
        #use different combinations of hyperparameters
        #select the one with the best performance in validation
        #use that combination to train on the 4 folds and then predict in the test set
    #we obtain one value type of model and test set
#Within each test set, compare each model with the other 3
#we will have 15 comparisons (compare each model with the other 3 x 5 tests = 15)
#select the model class being the best in the larger number of comparisons
#once you have the model class
    #try to ensemble different models of the same class because there are differences in the feature importance between highly-predicitve models (see paper)
    #they used the combination of shapely values across models, and they know that shapely fails with multicolinearty! check what the did to control for that


#Gene correlations
    #I am not sure this is a problem
    #Imagine you have a house pricing dataset and want to predict these prices using features of the houses
        #use 80% of houses for training and 20% for evaluation.
        #you could have several houses belonging to the same owner, one in training and other in evaluation. This owner tend to do good renovations and improve the condition. Therefore, the houses of this owner tend to be above the average of each area.
        #the model learns to predict price using features of the house and this includes some of the houses of this owner.
        #when evaluating, the model predict prices for some other houses of this owner. These houses are related to some houses in the training set as they share the owner, but the model does not know that!
        #The model just take the features of the house, and given the houses of this owner are in good condition, it predicts a high price.
        #the ultimate cause of their higher price is the owner, the proximate causes are the house features that owner has changed. 
        #The would go for genes that are physically close. They will likely have a more similar probability of selection than more distant genes. The ultimate cause is the proximity, but the proximate causes would be the recombination rate, GC content, etc... 
        #Therefore, I do not think there is any leak of data from training to evalution. The model learns patterns between genomic features and selection, then extrapolating to other genes.



#yifei applied preprocessing.StandardScaler() after the log to y, maybe we can applied woth the transformed in the pipeline before fitting, but I have seen does not improve THINK


'''
#in case you want to apply the log transformation within the pipeline, but this make the model MUCH WORSE compared to just apply the log to the original DF. Do not know why.
if FALSE:
    from sklearn.compose import TransformedTargetRegressor
    estim = TransformedTargetRegressor( \
        func=np.log, \
        inverse_func=np.exp, \
        regressor=Pipeline([ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', model)]), \
        check_inverse=True)
            #func: function applied to "y" before doing anything, i.e., even before fitting. 
                #This can be np.log or np.sqrt for example.
            #inverse_func: function applied to the prediction of the regressor, i.e., we are reverting the transformation done before fitting, so we can get back the predictions of response without transformation, being comparable with the raw values of the response.
                #If you do log, the inverse would be exp. For example np.log(1)=0.0, while np.exp(0.0)=1
            #regressor: the regressor used for fitting, it can be a pipeline, so you can add scaling (applied to Y and Xs) along with modeling
            #transformer: i understand this transformer is only aplied to "y", so it cannot be used if func/inverse_func are used.
            #check_inverse: Whether to check that `transform` followed by `inverse_transform` or `func` followed by `inverse_func` leads to the original targets 
    estim
'''




    #we can trandofm only y, the R2 is lowet in CV, but fit in the hist is good but not the scatter
        #Sometimes it can be beneficial to transform a highly exponential or multi-modal distribution to have a uniform distribution. This is especially useful for data with a large and sparse range of values, e.g. outliers that are common rather than rare.
            #https://machinelearningmastery.com/quantile-transforms-for-machine-learning/




#QuantileTransformer().fit(modeling_data_array[:, 0].reshape(-1, 1)).transform(modeling_data_array[:, 0].reshape(-1, 1))

#modeling_data_array[:, 0] = QuantileTransformer(n_quantiles=100, output_distribution="uniform").fit(modeling_data_array[:, 0].reshape(-1, 1)).transform(modeling_data_array[:, 0].reshape(-1, 1)).reshape(15458)
    #DATA LEAKAGE!!! if use this, you have to do the transformation of Y in the pipeline


#See the keys, you can see how you need to write two times regressor to reach alpha parameter of Ridge ("regressor__regressor__alpha"). Ridge is a regressor of the pipeline, but the pipeline is in turn a regressor of TransformedTargetRegressor.




#CHANGE X_train por modeling_data_train


###POR AQUIIII
#WE MODEL WITH LOG, BUT THEN ALEPLOT GIVE DIFFERNECES OF 0.2? THIS MAKES SENSE?
#THINK ABOUT THIS!!!



#Fit the model using the input data
estim.fit(X_train, y_train)

#See R2 in the whole dataset and in CV-subsets
from sklearn.metrics import r2_score
from sklearn.model_selection import cross_val_score

print("R2 in the whole dataset:")
r2_whole_dataset = r2_score(y, estim.predict(X))
print(r2_whole_dataset)
print("R2 across CV folds:")
r2_cv = np.mean(cross_val_score(estimator=estim, 
    X=X_train, \
    y=y_train, \
    cv=shuffle_split, \
    scoring="r2", \
    n_jobs=shuffle_split.n_splits))
print(r2_cv)

    #repeate with mse? like in uncovering expression signatures of synergistic drug response using an ensemble of explainable ai models?

print("R2 in the test dataset:")
r2_test_dataset = r2_score(y_test, estim.predict(X_test))
print(r2_test_dataset)


#Plot the observed sweep probability and the prediction in the whole dataset
import matplotlib.pyplot as plt

plt.hist(y, bins=50, color="green", alpha=0.4, label="Observed sweep probability")
plt.hist(estim.predict(X), bins=50, color="blue", alpha=0.4, label="Prediction")
plt.annotate("R2 whole dataset: " + str(np.round(r2_whole_dataset, 4)), 
             xy=(0.05, 0.7),
             xycoords='axes fraction')
plt.annotate("R2 CV: " + str(np.round(r2_cv, 4)), 
             xy=(0.05, 0.6),
             xycoords='axes fraction')
plt.annotate("R2 test set: " + str(np.round(r2_test_dataset, 4)), 
             xy=(0.05, 0.5),
             xycoords='axes fraction')
plt.legend(loc='upper left')
plt.savefig( \
    fname="./results/model_comparison/" + model_name + "_hist_pred_observed.png")
plt.close()


#When interpreting the R-Squared it is almost always a good idea to plot the data. That is, create a plot of the observed data and the predicted values of the data. This can reveal situations where R-Squared is highly misleading. For example, if the observed and predicted values do not appear as a cloud formed around a straight line, then the R-Squared, and the model itself, will be misleading. Similarly, outliers can make the R-Squared statistic be exaggerated or be much smaller than is appropriate to describe the overall pattern in the data.
    #https://www.displayr.com/8-tips-for-interpreting-r-squared/#:~:text=Don't%20use%20R%2DSquared%20to%20compare%20models&text=There%20are%20two%20different%20reasons,the%20variables%20are%20being%20transformed.
plt.scatter(y, estim.predict(X), s=0.5)
plt.xlabel("Observed log Flex-sweep probability")
plt.ylabel("Predicted log Flex-sweep probability")
plt.annotate( \
    "R2 whole dataset: " + str(np.round(r2_whole_dataset, 4)), \
    xy=(0.05, 0.9), \
    xycoords='axes fraction')
plt.annotate("R2 CV: " + str(np.round(r2_cv, 4)), \
    xy=(0.05, 0.8), \
    xycoords='axes fraction')
plt.annotate("R2 test set: " + str(np.round(r2_test_dataset, 4)), \
    xy=(0.05, 0.7), \
    xycoords='axes fraction')
plt.savefig( \
    fname="./results/model_comparison/" + model_name + "_scatter_pred_observed.png")
plt.close()

    #R2=0.5 in evaluation and test sets, which is usually considered OK
        #https://stephenallwright.com/good-r-squared-value/
    #the problem is that we have more error for strong sweep candidates (log probability close to zero). In classification occurs the same, as we have more false negatives than false positives (lower recall than precision). 
    #the MDR for iHS has around 0.7 for the whole dataset (we are above here) but fit very well the distribution of iHS, while here we have a problem with the right tail.
    #R2 or MSE/MAE for optimization?
        #https://machinelearningmastery.com/regression-metrics-for-machine-learning/

    #maybe this is a result itself, the genomic features considered cannot fully explain all sweeps predicted by flex-sweep. This makes sense because we do not have included all selective pressures affecting the human genome, so functions that have been targeted by these pressures would have an excess of sweep probability based on their genomic features.
        #we have strong sweep candidates that have less probability, but their probability is not zero. The model is predictiing some probability of sweep, but there is something else doing these genes enriched in sweep probability, maybe something related with their function.


    #mira vip distance and thermo

    #myabe using classification so we can optimize specifically recall?



    #log improve a bit prediction and the histogram is much better, we lose the big peak around zero probability having only one at 1
        #maybe you should apply log for DNN







from alepython import ale_plot


print_text("get the predictors used for training as a DF, which is required for ale_plot. We can just use the training DF that was converted to numpy array", header=4)
modeling_data_train_predictors = modeling_data_train.iloc[:, 1:]
    #the first column is the response, so we need to avoid it
print(modeling_data_train_predictors)


print_text("run ale_plot", header=4)
ale_plot(model=estim, 
    train_set=modeling_data_train_predictors,
        #pandas DF with values of predictors for training data
    features=["bat_distance_percentile_" + str(p_value_percentile)], 
    bins=10,
        #Number of bins used to split feature's space
        #I understand each bin has the same number of datapoints
        #so large intervals means very sparse data
        #Default 10
    monte_carlo=True,
        #Run or not a Monte Carlo, using only a sample of the 
        #dataset. Monte Carlo can help to detect data ranges
        #where data is very sparse and we cannot be confident
    monte_carlo_rep=50,
        #number of montecarlo replicas, default: 50
    monte_carlo_ratio=0.1,
        #Proportion of randomly selected samples from dataset
        #for each Monte-Carlo replica
        #default: 0.1
    rugplot_lim=1000)
        #A rug plot displays marks along an axis to visualise 
        #the distribution of the data.
        #The default is 1000, meaning that if you have more than 1000
        #samples, you do not need to see the bar where each sample is
        #set to None to plot to make always the rug plot
plt.close()


#bat 1% without no other selective pressure shows a clear decrease in the probability of selection as we move away from the BAT connectome genes. Note that there is a peak, an increase in the probability of selection at 717kb of distance from BAT genes. This is still relatively close. Remember that we consider the impact of genomic factors on the probability of selection of a gene that are up to 500kb from the center of the gene (500+500=1000kb windows) and then there is a clear decrease.
#whean adding vip, thermogenic and smt distance, the pattern is much more CLEAR with a decrease of flex-sweep probability. With specific combinations of these predictors, bAT pattern is not super clear, but it is with all of them. This is not a problem to me becuase alone this factor already shows pattern. and we adding other factors, i.e., controlling for then we see even clearer impact.
#note that the differences respect to the average prediction are similar to those of VIP distance! so I see potential in BAT set!!! We have this set validated, knowjing that is cohesive...!



#vip distance gives the expected result very well, in contrast with vip number. we are going to use ditance? think
    #we know the expected result because we already know VIPs are enrichd in positive selection so we should select the approach with power to detect this.
#thermogenic distance doe snot work, but this is Yoruba, we should check in non-african pops exposed to cold conditions
#check BAT? very good for climahealth. It is not comprehensive, we are losng genes important for thermo and BAT, but we know the genes included are important for thermo fiven what we found...

#there is some pattern of decrease with distance to the top 1% of the UCP1 connectome but it is noisy. The decrease occurs, but the pattern is not very clear. In contrast, the pattern is very clear for higher tops (0.05%) and lower tops from 2% to 8%, then from 9 to 50% the pattern tends to disappear in general.

#Maybe we can say that we know that genes biologically close to UCP1 tend to show cohesion and be enriched in associations with BAT, which is expected given the key role of UCP1 in BAT thermogenesis (cite our paper), but instead of just using 1% of genes, we can show ALE plots for several percentages from 0.5% to 10%, so we can show that in general there is a deficit of positive selection as we get away from genes biologically close to UCP1.

#check 15, 17, 22....
#explica readme in the data folder saying calc of dataset was done in the other folder
#create a dataset with bat distance at different p-value thresholds
#repasa codigo de BAT y sigue modeling, trying to decide best transformation and then model comparison
    #hink if run DNN optuna with log...



'''

##########
# optuna #
##########
print_text("optuna", header=1)


# We have to define a function that runs for each each trial. Inside this function
# - Select a value for each of the hyperparameters to be tunned using for that the optuna functions that allow to suggest values from int, float and categorical variables.
# - Then, create the keras regressor instance using the get_reg function previoulsy created and using the values of the hyperparameters previoulsy selected.
# - Introduce this regressor in a pipeline. This pipeline includes the tranformer to do the scaling so it is fit (get mean and SD) to the the training sets and not evaluation/validation sets. When using grid search ([link](https://stackoverflow.com/questions/51459406/how-to-apply-standardscaler-in-pipeline-in-scikit-learn-sklearn)) and cross_val_score ([link](https://stackoverflow.com/questions/44446501/how-to-standardize-data-with-sklearns-cross-val-score)), it uses the mean and SD of the training set of a given iteration and apply the transformation to the test set, then in the next iteration does the same but with the new training and test set.
# - The calculate the R2 of several CV partitions using cross_val_score.

# In[72]:


import optuna 
from scikeras.wrappers import KerasRegressor
from sklearn.pipeline import Pipeline
from sklearn.compose import TransformedTargetRegressor
from sklearn.model_selection import cross_val_score
from sklearn import preprocessing

def objective(trial):
    """return the r2-score"""
    
    #search space including integers, float and categorical
    batch_sizesT =  trial.suggest_int('regressor__regressor__batch_size', 
                                      low=np.min(batch_sizes), 
                                      high=np.max(batch_sizes), 
                                      step=2) 
        #Step: A step of discretization. 
        #Note that high is modified if the range is not divisible by step. 
        #So I understand that 1,10,5 means that only 1,5 and 10 will be sampled
        #you should select a step size that is relevant for the parameter
        #e.g., a batch of 1 or 4 will be very similar, so we use a step of 10
        #in the case of batch, we are making smaller steps in the second run
        #because it seems, the best trails are around 10-20, we need 
        #more resolution there (see above).
        #Log: You can also apply log, according to the manual
        #this makes more likely to sample smaller values
        #this maybe useful for learning rate
            #https://optuna.readthedocs.io/en/stable/reference/generated/optuna.trial.Trial.html#optuna.trial.Trial.suggest_int
    epoch_sizesT =  trial.suggest_int('regressor__regressor__epochs', 
                                      low=np.min(epoch_sizes), 
                                      high=np.max(epoch_sizes), 
                                      step=20)
    optimizersT = trial.suggest_categorical('regressor__regressor__optimizer', 
                                            optimizers)
    alphasT = trial.suggest_float('regressor__regressor__optimizer__learning_rate', 
                                  low=np.min(alphas),
                                  high=np.max(alphas), 
                                  step=None, 
                                  log=True) #for log see alphas section in this notebook
    init_modeT = trial.suggest_categorical('regressor__regressor__model__init_mode', 
                                           init_mode)
    activationsT = trial.suggest_categorical('regressor__regressor__model__activation', 
                                           activations)
    dropout_ratesT = trial.suggest_float('regressor__regressor__model__dropout_rate', 
                                            low=np.min(dropout_rates), 
                                            high=np.max(dropout_rates), 
                                            step=None)
    weight_constraintsT = trial.suggest_float('regressor__regressor__model__weight_constraint', 
                                            low=np.min(weight_constraints), 
                                            high=np.max(weight_constraints), 
                                            step=None)
    regu_L1_valuesT = trial.suggest_float('regressor__regressor__model__regu_L1', 
                                            low=np.min(regu_L1_values), 
                                            high=np.max(regu_L1_values), 
                                            step=None, log=True)
    regu_L2_valuesT = trial.suggest_float('regressor__regressor__model__regu_L2', 
                                            low=np.min(regu_L2_values), 
                                            high=np.max(regu_L2_values), 
                                            step=None, log=True)
        #The most common type of regularization is L2, also called simply 
        #“weight decay,” with values often on a logarithmic scale between 
        #0 and 0.1, such as 0.1, 0.001, 0.0001, etc.        
        #https://machinelearningmastery.com/how-to-reduce-overfitting-in-deep-learning-with-weight-regularization/
    n_layersT =  trial.suggest_int('regressor__regressor__model__n_layers', 
                                      low=np.min(n_layers), 
                                      high=np.max(n_layers), 
                                      step=1)
    n_unitsT =  trial.suggest_int('regressor__regressor__model__n_units', 
                                      low=np.min(n_units), 
                                      high=np.max(n_units), 
                                      step=20)
    lossesT = trial.suggest_categorical('regressor__regressor__loss', 
                                            losses)

    #make an instance of keras regressors using the previous hyperparameters
    #and get reg function for our DNN 
    keras_regressor = KerasRegressor(
        model=get_reg,
        batch_size=batch_sizesT,
        epochs=epoch_sizesT,
        optimizer=optimizersT,
        optimizer__learning_rate=alphasT,
        model__init_mode=init_modeT,
        model__activation=activationsT,
        model__dropout_rate=dropout_ratesT, 
        model__weight_constraint=weight_constraintsT,
        model__regu_L1=regu_L1_valuesT,
        model__regu_L2=regu_L2_valuesT,
        model__n_layers=n_layersT,
        model__n_units=n_unitsT,
        loss=lossesT,
        random_state=1, #the same architecture will give the same result
        warm_start=False, 
            #If False, subsequent fit calls will reset the entire model.
            #This has no impact on partial_fit, which always trains
            #for a single epoch starting from the current epoch.
        verbose=0) 

    #include the model in the pipeline
    final_estimator = TransformedTargetRegressor(
        regressor=Pipeline([
            ('scale', preprocessing.StandardScaler()),
            ('regressor',keras_regressor)]),
        transformer=preprocessing.StandardScaler())
    
    #cross validation
    score = cross_val_score(estimator=final_estimator, 
        X=X_train, 
        y=y_train, 
        scoring="r2",
        cv=shuffle_split,
        n_jobs=number_jobs).mean() 
        #This gives nan if any of the folds gives nan. We want this because
        #if an average R2 score is calculated with only 6 folds (the rest were nan)
        #this is is not comparable with an R2 calculated with 10 folds
        #this is not a fair comparison. So we remove any combination that gives
        #nan for any fold.

    #return the R2 score obtained from cross_val_score
    return score


# **Note**
# 
# Note that Optuna support multiobjective optimization! i.e., optimize more than scores, like R2 and MSE! A genetic algorithm is used in this case.

# Then, we have to create the optuna study. Here, we have to indicate: 
# - the sampler dedicated to select different hyperparameter combinations. We will use the default sampler (TPESampler), which works well for integers, float and categorical variables ([link](https://optuna.readthedocs.io/en/stable/reference/samplers/index.html#module-optuna.samplers)).
# - the direction, meaining if you want to select hyperparameter combinations with higher or lower scores. In our case, we need maximize because we want higher R2 scores.

# We have also to indicate where the study will be stored, i.e., the different combinations of hyperparameters in eac trial, the score.... This is an important decision because its implications for parallelizing. 
# 
# If we are going to use the old good `n_jobs` in `study.optimize()`, we can just save the study's information in memory. You will use several cores to run several trials at the same time. This, however, has two problems. First, you can encounter threating problems with python's GIL if your calculations require high CPU time (see `study.optimize()` help), and this could be our case. In addition, we are running the study one time with one seed, so if you want to perform more searched from different starting points to explore from different positions the parametric space, you need to run different scripts but they WILL NOT be connected between them. According to this [guy](https://github.com/optuna/optuna/issues/2351#issuecomment-782972977), he got the better models running different processes with less steps than just one with more steps.
# 
# There is an interesting alternative, which is the use of relational databases (RDB; [link](https://optuna.readthedocs.io/en/stable/tutorial/20_recipes/001_rdb.html#rdb)). When creating a study, you can set its name and the create a SQLite database and set load_if_exists as True. This works as follows:
# 
# You run a first python script with the objective function, which will run each trial, and the optuna study. As it is the first time, the SQLite database will be created and the information of the different trials developed by this script will be saved there. 
# 
# You can run another python script with the code but changing the seed of the TPESampler, i.e., we are starting another search from a different starting point and different suggested parameters, but this belongs to the same study! because load_if_exists=True! So the information of this new process is saved there, and it has available the information of the previous process run. Note that if you set n_trials in study.optimize to 10, and run two different processes, you get a total of 20 trials!
# 
# This means that you can stop of study and resume it later. You can also run multiple processes at the same time, and you can see in real time how the best trial is selected considering all the processes ([link](https://optuna.readthedocs.io/en/stable/tutorial/10_key_features/004_distributed.html#sphx-glr-download-tutorial-10-key-features-004-distributed-py)). IMPORTANTLY, they get parameter suggestions based on shared trials’ history, so we are taking advtange of the full exploration of the parametric space we are doing!
# 
# We can easily use SQLlite using `sqlite:///data_base.db`. This create a SQLite database that you can access easily both from linux and from optuna (using create_study with load_if_exists=True and look into `study.trials_dataframe`). SQLite has cool points like the easy access and speed but it works worse if there is a lot of data (because everything is in a single file) and multiple processes need to access to the database. In those cases it is better to create a MySQL database. This could be done with the following lines:
# 
# `mysql -u root -e "CREATE DATABASE IF NOT EXISTS example"`
# `optuna create-study --study-name "distributed-example" --storage "mysql://root@localhost/example"`
# 
# We are going to stick to SQLite for now. We can try MySQL in the future if required.
# - Most of the work will be done within each trail with several nets fitted to 10 CV folds. Every 10 of these fits only gives a mean score that will be the output for optuna. So, if we run just 5-10 different processes using each one 10 cores for cv_score, then we would use 50-100 cores. In summary, this is not a big deal for optuna, only 5-10 processes are saving info. If using 5-10 processes with different seeds, we get a model with R2 above 0.7 in test, I do not think that any reviewer would have problems with possible unexplored parametric space, because we already got good models and we have explored different seeds, different searches.
# - In addition, I am not sure if I could make work a MySQL database in singularity.
# 
# You can check this thread for an example script to be run optuna in parallele from command line ([link](https://stackoverflow.com/questions/73569369/optuna-hyperparameter-search-repeats-hyperparameters-across-studies-with-paralle)).

# In[73]:

# define a function to run the optuna study taking as argument the seed and 
#the number of trials

from optuna.samplers import TPESampler
from sklearn.pipeline import Pipeline
from sklearn.compose import TransformedTargetRegressor
from sklearn.metrics import r2_score
import joblib

#optuna_seed=1
#n_trials=2
def main_optuna(optuna_seed, n_trials):
    """first argument is optuna seed, while second is the number of optuna trial"""
    
    # create a study
    study = optuna.create_study(study_name = "yoruba_flex_sweep_closest_window_center_optimization", 
        direction="maximize",
        storage='sqlite:///results/optuna_optimization/yoruba_flex_sweep_closest_window_center_optimization.db', 
        sampler=optuna.samplers.TPESampler(seed=optuna_seed),
        load_if_exists=True)

    # Run the search. We are not using prunning becuase I have not found a way to get metrics while the cross-validation is performed in scikit-learn. So we need to finish a trial to get a score and decide if it is good or not.
    study.optimize(func=objective, 
                   n_trials=n_trials,
                   n_jobs=1, #we will parallelize using RDB, not n_jobs
                   gc_after_trial=True, #for recollectin of garbage
                   show_progress_bar=False)

    #we do not do anything else because even if the current process is finished, others could be still working and the best model is still not present. Of course, it would be printed at the end of that final process, but we would have several best models in the output. Better to do it in a separate script when the whole study is finished.

# Note that if you have an error with a hyperparameter suggestion, like indicating the wrong name of an optimizer, you will keep getting the error even if you change unless you remove the original database. I guess that this is caused because the new process uses information from the database, so it can still read the wrong optimizer name, which was saved during the process with the wrong name.

# I have checked that trials of different processes are all saved into the same SQLite database. I have also checked that the parameters combinations are completely dependent of the seed. So if you call again the same study and start a new process with different seed, you get different hypeterparameters suggested. If we repeat again with other process but using the seed of the first process, we get the exact same combination than in the first process. This, along with the fact that the whole database is used for all process in order to suggest new parameters, makes this approach really convinient.

#run main_optuna parsing the required arugments
#if __name__ == "__main__": #used in our source script but I think not needed right no
    #https://www.geeksforgeeks.org/what-does-the-if-__name__-__main__-do/

#define input arugments to be passed by name
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--optuna_seed", help="Seed for Optuna process")
parser.add_argument("--n_trials", help="Number of trials in Optuna process")
args=parser.parse_args()

#get the arguments of the function that have been passed through command line
optuna_seed = np.int64(args.optuna_seed)
    #convert to int64 to avoid problems with optuna.study
n_trials = np.int64(args.n_trials)
    #https://stackoverflow.com/questions/40001892/reading-named-command-arguments
    #alternative with unnamed arguments
        #https://www.saltycrane.com/blog/2007/12/how-to-pass-command-line-arguments-to/

#optuna_seed=1
#n_trials=2
main_optuna(optuna_seed=optuna_seed, n_trials=n_trials) 
    #https://stackoverflow.com/questions/73569369/optuna-hyperparameter-search-repeats-hyperparameters-across-studies-with-paralle

#we are going to use 5 kfold for now to seep up things and explore more hyperparameter combinations, then when we get good architecutres could narrow the search

'''