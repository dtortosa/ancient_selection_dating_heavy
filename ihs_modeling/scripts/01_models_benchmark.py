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

#This script will perform model comparison for analysing iHS.
#We used as reference the script for model comparison in flex-sweep:
    #"/home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating_heavy_analyses/dating_climate_adaptation/flex_sweep_modeling/scripts/01_models_benchmark.py" 

#The idea is to get a sense of the best model for iHS using Yoruba as reference, then deploy in the rest of populations. Ideally, if we can do the modelling and exploration of result fast and automatic enough, we could apply across the 26 pops for which we have average iHS across gene windows. Then check increases of selection in cold regions for BAT. Finally, check whether CVD markers in combat_genes show more signals of association in BAT closing the circle (we are not going to look for the sense of the selection and associaton of alleles).




##############################
#region INITIAL STEPS ########
##############################


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



#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--model_name", type=str, default="elastic_net", help="Selected model. Integer string, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
model_name = args.model_name



############################
# starting with the script #
############################
print_text("starting with model " + str(model_name), header=1)



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
import tensorflow as tf # type: ignore
tf.random.set_seed(seed_value)


print_text("Configure a new global `tensorflow` session: NOT DOING IT", header=2)
#we are not running this because it gets me an error, we will see if we can run tensorflow without this
if False:
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
        ./results/model_comparison/individual_results; \
    ls -l ./results")

# endregion






##################################
# region data preparation ########
##################################
print_text("data preparation", header=1)

print_text("define window size and response", header=4)
gene_window_size = "1000kb"
    #see MDR paper about window selection and why we focused in 1000kb gene windows
print(f"The selected window size is {gene_window_size}")
response="mean_ihs_" + gene_window_size
print(f"The response is {response}")

print_text("load the data", header=4)
import pandas as pd
final_data_yoruba = pd.read_csv( \
    "./data/YRID_modeling_dataset_v1.tsv.gz", \
    sep="\t", \
    header=0 \
)


print_text("clean the data", header=3)
print_text("exclude some columns we are not interested in", header=4)
columns_to_exclude = ["gene_id"]
final_data_yoruba_subset = final_data_yoruba[[column for column in final_data_yoruba.columns if column not in columns_to_exclude]]
print(final_data_yoruba_subset)
print(f"Columns excluded: {columns_to_exclude}")

print_text("make deep copy of the data to do further operations", header=4)
modeling_data = final_data_yoruba_subset.copy(deep=True)
    #deep=True: 
        #a new object will be created with a copy of the calling object's data and indices. Modifications to the data or indices of the copy will not be reflected in the original object (see notes below).
print(modeling_data)


print_text("Apply log transformation to the target variable using the original DF as source", header=4)
import numpy as np
modeling_data[response] = final_data_yoruba_subset[response].apply(lambda x: np.log(x))
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

# endregion






########################################
# region prepare models and HPs ########
########################################

print_text("define function to run across folds and model classes within a nested CV schema", header=1)
print_text("prepare models and HPs", header=2)
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
import xgboost as xgb # type: ignore
from scikeras.wrappers import KerasRegressor # type: ignore
dict_models = { \
    "elastic_net": {"instance": "ElasticNet()"}, \
    "random_forest": {"instance": "RandomForestRegressor()"}, \
    "xgboost": {"instance": "xgb.XGBRegressor()"}, \
    "neural_nets": {"instance": "KerasRegressor()"}}
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
#set the HPs
dict_models["elastic_net"]["HPs"] = { \
    "regressor__l1_ratio": np.arange(0, 1.01, 0.01), #1.01 is not included, it ends at 1, which is like Lasso
    "regressor__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]}
print(dict_models["elastic_net"])
    #HPs
        #l1_ratio
            #alpha in elastic net theory (different from alpha in scikit learn) determines how much weight is given to each of the L1 and L2 penalties (see above). Alpha is a value between 0 and 1 and is used to weight the contribution of the L1 penalty and one minus the alpha value is used to weight the L2 penalty.
            #For example, an alpha of 0.5 would provide a 50 percent contribution of each penalty to the loss function. An alpha value of 0 gives all weight to the L2 penalty and a value of 1 gives all weight to the L1 penalty.
            #the alpha hyperparameter can be set via the “l1_ratio” argument that controls the contribution of the L1 and L2 penalties.
            #Specifically, l1_ratio = 1 is the lasso (L1) penalty. Currently, l1_ratio <= 0.01 is not reliable, unless you supply your own sequence of alpha.
                #No problem because we supply the alpha values.
            #In contrast, l1_ratio=0 is Ridge (L2) penalty.
        #alpha
            #Another hyperparameter is provided called “lambda” (alpha in scikit learn) that controls the weighting of the sum of both penalties to the loss function. A default value of 1.0 is used to use the fully weighted penalty; a value of 0 excludes the penalty. Very small values of lambada, such as 1e-3 or smaller, are common.
            #In scikit learn, “alpha” argument controls the contribution of the sum of both penalties to the loss function
            #``alpha = 0`` is equivalent to an ordinary least square, solved by the :class:`LinearRegression` object
                #Therefore, typical linear model
            #For numerical reasons, using ``alpha = 0`` with the ``Lasso`` object (l1_ratio=1) is not advised. Given this, you should use the :class:`LinearRegression` object.
                #No problem because we are not including alpha=0.
        #Confusingly, as you can see the notation is different in scikit learn, as indicated in the help:
            #The parameter l1_ratio corresponds to alpha in the glmnet R package while alpha corresponds to the lambda parameter in glmnet.
        #fit_intercept (default=True): 
            #Whether the intercept should be estimated or not. If ``False``, the data is assumed to be already centered.
            #our data is scaled (so I guess centered), but I see no reason to avoid the intercept in a linear model.
        #max_iter (default=1000):
            #The maximum number of iterations.
            #Using the default because increasing to 2000 does not improve convergence. Convergence problems could be caused by some combinations of alpha and l1_ratio. For example, alpha=0 (no penalty) gives convergence problems.
        #copy_X (default=True)
            #If ``True``, X will be copied; else, it may be overwritten.
            #Ensure we copy X always? Not sure, just use default.
        #tol (default=1e-4)
            #The tolerance for the optimization: if the updates are smaller than ``tol``, the optimization code checks the dual gap for optimality and continues until it is smaller than ``tol``, see Notes below.
            #It seems this is the maximum level of change in the last weight modified that is required in order to stop the model and finish. But I do not fully understand this argument, so I am going to use the default.
                #https://stats.stackexchange.com/questions/445831/how-is-tol-used-in-scikit-learns-lasso-and-elasticnet
        #warm_start (default=False)
            #When set to ``True``, reuse the solution of the previous call to fit as initialization, otherwise, just erase the previous solution. See :term:`the Glossary <warm_start>`.
            #Avoid reusing any previous call, just fresh start
        #positive (default=False)
            #When set to ``True``, forces the coefficients to be positive.
            #No makes sense for us. Recombination rate should have negative coefficients...
        #There are many other arguments used to speed up calculations, but we do not need them as we can run everything for this model in less than 1 hour.
    #Convergence problems
        #We have convergence problems for some combinations of hyperparameters
        #I have increased the number of iterations from 1000 to 50000 without changes. The tolerance (change in the last weight modified) decreases but it gets stuck after 10000 iterations.
            #https://stats.stackexchange.com/questions/445831/how-is-tol-used-in-scikit-learns-lasso-and-elasticnet
        #The other solution that gives the help is to increase the regularization and we do that as we explore a wide range of penalization levels.
    #Decision about using Elastic net over Ridge, Lasso and LinearRegression
        #alpha=0
            #with alpha=0 we exactly the same R2 scores for ElasticNet as Lasso (l1_ratio=1), ElasticNet as Ridge (l1_ratio=0) and LinearRegression.
            #with alpha=1e-5 we get almost the same results with the three approaches.
            #This is despite the warning of convergence in ElasticNet
        #l1_ratio=1
            #with l1_ratio=1 we get exactly the same results with ElasticNet as Lasso (l1_ratio=1) and Lasso().
        #l1_ratio=0
            #with l1_ratio=0 we get results that are close but not the same with ElasticNet as Lasso (l1_ratio=1) and Ridge (~0.4 of difference).
            #if we set alpha=alpha/X_train.shape[0] in ElasticNet, we get exactly the same results in ElasticNet and Ridge. I have also check what happens in the lowest alpha value I am going to use (1e-5). Both approaches gives almost the same even if we do not divide alpha by X_train.shape[0] for ElasticNet. So it seems we are doing the same here.
                #https://stackoverflow.com/questions/47365978/scikit-learn-elastic-net-approaching-ridge
            #This is despite the warning of convergence in ElasticNet
        #It is important to note that I am ONLY going to use this for selecting the best model class. I am not going to extract coefficients or make inferences about the data. I just need comparable score results in order to select the model class that consistently perform better across folds. These results indicate that ElasticNet can produce similar results than the other functions in terms of score. Therefore, I think we can just use ElasticNet with a wide range of parameters in order to cover Ridge, Lasso and LinearRegression.
    #grid
        #Janizek 2023
            #the ‘alpha’ parameter was tuned over values ranging from 0.1 to 100, whereas the ‘l1_ratio’ parameter was tuned from 0.25 to 0.75.
        #ML mastery
            #One approach would be to gird search l1_ratio values between 0 and 1 with a 0.1 or 0.01 separation and alpha values from perhaps 1e-5 to 100 on a log-10 scale and discover what works best for a dataset.
        #My opinion:
            #The grid of ML Mastery is wider but it takes less than 1 hour to run across all folds without parallelization!
            #It checks 
                #from very low penalty (alpha close to zero; regular linear model) to very high penalty (alpha at 100)
                #at each level of penalty, it consider more impact for L2 penalty (lambda close to zero), equal impact (lambda=0.5) and more impact for L1 penalty (lambda close to 1) 
            #It explores a lot of scenarios that includes pure Ridge, pure Lasso, pure LinearRegression and many intermediate steps in less than 1 hour!
    #we get warning about not converging. ML Mastery says that is ok.


print_text("random forest", header=4)
#general notes about random forest
    #Random forest is based on Decision trees:
        #In a decision tree, you have several predictors (columns, X) used to predict a response (y).
        #A decision tree is made by looking what value of a given predictor can be used to effectively split data points that are different (e.g., low and high flex-sweep probability).
        #For example, using as decision criteria a recombination rate of 4 will be very likely useful to separate genes with very low flex sweep probability because above 4, it is very difficult to detect any sweep due to the eroding effect of recombination. Of course, we will have low and high flew-sweep probabilities in both groups, but we have created to groups that are more homogeneous. In other words, we have increased the impurity of the groups.
        #To decide that criteria use in each decision tree we calculate the variance reduction, i.e., impurity.
            #We calculate the variance of flew-sweep probability in the whole dataset (sum of difference between y minus average divided by sample size). 
            #Then calculate the variance in the two child nodes using criteria one
            #Calculate the variance in the two child nodes using criteria two.
            #Then calculate the variance reducing in both cases with respect to the full dataset. We select the criteria that reduce more variance, i.e., decrease impurity.
            #Note that the model calculates the variance reduction for every possible split!!
        #the process continue for each branch, using other values of the same predictor or other predictor to split each of the two groups in smaller groups that are more homogeneous, i.e., have less impurity.
        #At the end, what we have done is to divided the parametric space in the different regions based on the flex-sweep probabiltiy values. For example, very low recombination and very low vip distance will probably have genes with very high flex-sweep probability. Therefore, even though we are doing regressions, we are effectively classifying the parametric space and the samples.
        #this process will continue until we have reached our desired depth. 
        #At the end, we have different groups with (hopefully) low impurity. For example, a group with very low flex-sweep probs, a group with low but a bit higher flex-sweep probs... and so on.
        #Then to predict a new gene, 
            #you check what conditions the gene satisfies for each decision node (i.e., each threshold of a given predictor), advancing through the tree until reaching the terminal node.
            #once the gene is classified in a terminal node, its prediction is then obtained by calculating the average flex-sweep probability of genes classified in that group during training. For example, if our gene under prediction is classified as very-low flex-sweep prob and we have there 10 training genes, we take their flex-sweep probabilties and calculate the average, that is the prediction for the new gene.
        #very good and simple videos about decision trees classifier and regressors
            #https://youtu.be/ZVR2Way4nwQ
            #https://www.youtube.com/watch?v=UhY5vPfQIrA&ab_channel=NormalizedNerd
    #Random forest
        #What is the problem with a decision tree? that is very dependent on the data that have seen and has low very capacity for generalizying. If you just change a few values of the predictors for a sample, you will get a completely different decision tree.
        #The solution is to use a random set of different decision trees generating a random forest
        #First you do bootstrap
            #start randomly selecting samples (rows) with replacement (the same sample can be selected two times) until you reach a given sample size (usually the size of the original sample size). 
            #for each of this random set of sample you will fit a decision tree. As the decision tree is influenced by small changes in training data, the different decision trees will be different between them.
        #Another random step
            #you randomly select the predictors considered to make decision in each decision tree.
            #if all trees uses the same predictors, it is likely that they will tend to be similar.
            #therefore, each decision tree uses different samples and different predictors as input. 
            #this greatly reduces the correlation between trees and errors.
        #Once you have all your decision trees fitted to the corresponding random set of samples and using a random set of predictors, you can predict.
        #you have a new row to predict flex-sweep probability. Do the prediction across all trees and then obtain the average (if classification it would be a vote, class with more votes across trees wins).
        #This last steps is called aggregation and gives name in combination to bootstrap to bagging (bootstrap+aggregation)
        #very good video about random forest
            #https://youtu.be/v6VJ2RO66Ag
#set the HPs
dict_models["random_forest"]["HPs"] = { \
    "regressor__n_estimators": np.arange(10, 1500, 400), \
    "regressor__max_features": np.arange(2, 16, 4), \
    "regressor__max_depth": [i for i in range(1,8,3)] + [None], \
    "regressor__max_samples": [i for i in np.arange(0.1,0.8,0.3)] + [1]}
        #the last number (second argument) in np.arrange is not included. So we end at 0.7 in max_samples, then we and None to have 100% of samples.
print(dict_models["random_forest"])
    #The most important HPs in general
        #The following HPs are those that are being consistently recommended as those with the greatest impact on predictive power.
            #https://machinelearningmastery.com/random-forest-ensemble-in-python/
            #https://stackoverflow.com/questions/36107820/how-to-tune-parameters-in-random-forest-using-scikit-learn
            #https://www.analyticsvidhya.com/blog/2020/03/beginners-guide-random-forest-hyperparameter-tuning/
        #max_features
            #Perhaps the most important hyperparameter to tune for the random forest is the number of random features to consider at each split point. A good heuristic for regression is to set this hyperparameter to 1/3 the number of input features.
            #It is set via the max_features argument. In this case, for our test dataset, the heuristic would be 19/3 or about 6 features.
            #We will explore the effect of the number of features randomly selected at each split point on predictive power. We will look for around 6 but including more values above 6 because I have seen in CV that higher number of features works better, but there is not great improvement after 12-14...
        #n_estimators
            #the number of decision trees in the ensemble can be set. Often, this is increased until no further improvement is seen. When fitting a final model, it may be desirable to either increase the number of trees until the variance of the model is reduced across repeated evaluations, or to fit multiple final models and average their predictions.
            #Typically, the number of trees is increased until the model performance stabilizes. Intuition might suggest that more trees will lead to overfitting, although this is not the case. Both bagging and random forest algorithms appear to be somewhat immune to overfitting the training dataset given the stochastic nature of the learning algorithm.
            #The number of trees can be set via the “n_estimators” argument and defaults to 100.
            #We will explore the effect of the number of trees with values up to 1000. If the cross-validation performance profiles are still improving at 1,000 trees, then incorporate more trees until performance levels off.
        #max_depth
            #Another important hyperparameter to tune is the depth of the decision trees. Deeper trees are often more overfit to the training data, but also less correlated, which in turn may improve the performance of the ensemble. Depths from 1 to 10 levels may be effective.
            #By default, trees are constructed to an arbitrary depth and are not pruned. This is a sensible default, although we can also explore fitting trees with different fixed depths.
            #The maximum tree depth can be specified via the max_depth argument and is set to None (no maximum depth) by default.
            #We will explore "from 1 to 7 and None=full". Instead of doing 1,2,3... we are going to select a few values between 1 and 7.
        #max_samples
            #The “max_samples” argument can be set to a float between 0 and 1 to control the percentage of the size of the training dataset to make the bootstrap sample used to train each decision tree. For example, if the training dataset has 100 rows, the max_samples argument could be set to 0.5 and each decision tree will be fit on a bootstrap sample with (100 * 0.5) or 50 rows of data. 
            #A smaller sample size will make trees more different, and a larger sample size will make the trees more similar. Setting max_samples to “None” will make the sample size the same size as the training dataset and this is the default.
            #In general, It is good practice to make the bootstrap sample as large as the original dataset size. We will include this option for sure.
    #Other HPs that could be used in a more detailed search if this model class is selected
        #min_samples_split and min_samples_leaf
            #These two were recommended during TDI in case more tunning was required.
            #From documentation "The main difference between the two is that min_samples_leaf guarantees a minimum number of samples in a leaf, while min_samples_split can create arbitrary small leaves, though min_samples_split is more common in the literature."
            #we should make the distinction between a leaf (also called external node) and an internal node. An internal node will have further splits (also called children), while a leaf is by definition a node without any children (without any further splits).
            #min_samples_split specifies the minimum number of samples required to split an internal node, while min_samples_leaf specifies the minimum number of samples required to be at a leaf node.
            #For instance, if min_samples_split = 5, and there are 7 samples at an internal node, then the split is allowed. But let's say the split results in two leaves, one with 1 sample, and another with 6 samples. If min_samples_leaf = 2, then the split won't be allowed (even if the internal node has 7 samples) because one of the leaves resulted will have less then the minimum number of samples required to be at a leaf node.
            #As the documentation referenced above mentions, min_samples_leaf guarantees a minimum number of samples in every leaf, no matter the value of min_samples_split.
            #In other words,
                #The min_samples_split parameter will evaluate the number of samples in the node, and if the number is less than the minimum the split will be avoided and the node will be a leaf.
                #The min_samples_leaf parameter checks before the node is generated, that is, if the possible split results in a child with fewer samples, the split will be avoided (since the minimum number of samples for the child to be a leaf has not been reached) and the node will be replaced by a leaf.
            #https://stackoverflow.com/a/46488222/12772630
            #https://towardsdatascience.com/hyperparameter-tuning-the-random-forest-in-python-using-scikit-learn-28d2aa77dd74
        #max_leaf_nodes
            #This hyperparameter sets a condition on the splitting of the nodes in the tree and hence restricts the growth of the tree. If after splitting we have more terminal nodes than the specified number of terminal nodes, it will stop the splitting and the tree will not grow further.
            #Let’s say we set the maximum terminal nodes as 2 in this case. As there is only one node when starting, it will allow the tree to grow further. Now, after the first split, you can see that there are 2 nodes here and we have set the maximum terminal nodes as 2. Hence, the tree will terminate here and will not grow further. This is how setting the maximum terminal nodes or max_leaf_nodes can help us in preventing overfitting.
            #When the parameter value is very small, the tree is underfitting and as the parameter value increases, the performance of the tree over both test and train increases, but beyond a point, the tree starts to overfit as the parameter value goes beyond 25.
            #This makes sense to me because as the tree grows and you have more terminal leaves, i.e., more groups of samples. Therefore you are splitting the space in more groups, creating small groups defined by very specific combinations of predictors based on the training data. As more grow, more specific conditions. If we reduce the growth, we should limit the overfitting to the specific conditions of the training data and increase generalization.
                #https://www.analyticsvidhya.com/blog/2020/03/beginners-guide-random-forest-hyperparameter-tuning/
            #This parameter is complementary to max_depth, because two trees with the same number of terminal nodes can have very different depths
                #https://stats.stackexchange.com/questions/544111/max-depth-vs-max-leaf-nodes-in-scikit-learns-randomforestclassifier
                #https://www.geeksforgeeks.org/random-forest-hyperparameter-tuning-in-python/
        #bootstrap: NOT RECOMMENDED CHANGING THE DEFAULT
            #Each decision tree in the ensemble is fit on a bootstrap sample drawn from the training dataset. This can be turned off by setting the “bootstrap” argument to False, if you desire. In that case, the whole training dataset will be used to train each decision tree. This is not recommended.
        #there are more HPs that you can explore when working with the selected model class
            #criterion...

print_text("XGboost", header=4)
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
#set the HPs
dict_models["xgboost"]["HPs"] = { \
    "regressor__max_depth": [i for i in range(1,15,4)] + [None], \
    "regressor__min_child_weight": [1, 10, 100], \
    "regressor__subsample":  [i for i in np.arange(0.1,0.8,0.3)]+[1], \
    "regressor__colsample_bytree": [i for i in np.arange(0.1,0.8,0.3)]+[1], \
    "regressor__eta": [0.001, 0.01, 0.1, 0.2, 0.3], \
    "regressor__n_estimators": [10]+[i for i in np.arange(200, 2200, 500)]}
print(dict_models["xgboost"])
    #About the HPs. To fully understand their impact see general notes and, in particular, regularization methods used in order to avoid overfitting.
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
        #there are more HPs that can be useful to combat overfitting and improve performance like lamba or alpha (for regularization) or gamma, but we will not explore them here. We can do it if this class is finally selected.
            #https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
        #You can also narrow the hyperparametric space by using a sequential approach: All default, just tune max_depth and min_child_weight. See what range of values give the best and select it. Repeat now with subsample and colsample leaving all default except max_depth and min_child_weight. Select the best range for subsample and colsample. Then tune eta and n_estimators.
            #https://datascience.stackexchange.com/a/108242


print_text("Neural Networks", header=4)
print("Neural Networks - General Notes")
#general notes about Neural Networks (you can also look at your DNN course in physalia)
    #Artificial neural networks can be considered as function approximation algorithms. 
    #In a supervised learning setting, when presented with many input observations representing the problem of interest, together with their corresponding target outputs, the artificial neural network will seek to approximate the mapping that exists between the two. 
    #The human brain consists of a massive network of interconnected neurons (around one hundred billion of them), with each comprising a cell body, a set of fibres called dendrites, and an axon:
    #The dendrites act as the input channels to a neuron, whereas the axon acts as the output channel. Therefore, a neuron would receive input signals through its dendrites, which in turn would be connected to the (output) axons of other neighbouring neurons. In this manner, a sufficiently strong electrical pulse (also called an action potential) can be transmitted along the axon of one neuron, to all the other neurons that are connected to it. This permits signals to be propagated along the structure of the human brain. 
    #An artificial neural network is analogous to the structure of the human brain, because (1) it is similarly composed of a large number of interconnected neurons that, (2) seek to propagate information across the network by, (3) receiving sets of stimuli from neighbouring neurons and mapping these to outputs, to be fed to the next layer of neurons. 
    #The structure of an artificial neural network is typically organised into layers of neurons (recall the depiction of a tree diagram). For example, the following diagram illustrates a fully-connected  neural network, where all the neurons in one layer are connected to all the neurons in the next layer:
    #The inputs are presented on the left hand side of the network, and the information propagates  (or flows) rightward towards the outputs at the opposite end. Since the information is, hereby, propagating in the forward direction through the network, then we would also refer to such a network as a feedforward neural network. 
    #The layers of neurons in between the input and output layers are called hidden layers, because they are not directly accessible. 
    #Each connection (represented by an arrow in the diagram) between two neurons is attributed a weight, which acts on the data flowing through the network, as we will see shortly. 
    #Math behind
        #More specifically, let’s say that a particular artificial neuron (or a perceptron, as Frank Rosenblatt had initially named it) receives n inputs, [x1, …, xn], where each connection is attributed a corresponding weight, [w1, …, wn]. 
        #The first operation that is carried out multiplies the input values by their corresponding weight, and adds a bias term, b, to their sum, producing an output, z: z = ((x1 × w1) + (x2 × w2) + … + (xn × wn)) + b
        #This weighted sum calculation that we have performed so far is a linear operation. If every neuron had to implement this particular calculation alone, then the neural network would be restricted to learning only linear input-output mappings. 
        #Hence, a second operation is performed by each neuron that transforms the weighted sum by the application of a nonlinear activation function, a(.).
            #The activation determines whether the particular combination of weights and inputs + bias is enough to activate that neuron and send output to a neuron in the next layer.
            #that output will be then multiply for the weight of the next neuron + bias and then the activation of that second neuron will determine if that new neuron is fired up.
            #the previous neuron will do this for each of the neurons in the next layer if the NN is fully connected.
        #Therefore, each neuron can be considered to implement a nonlinear function that maps a set of inputs to an output activation. 
    #training the network
        #Training an artificial neural network involves the process of searching for the set of weights that model best the patterns in the data. It is a process that employs the backpropagation and gradient descent algorithms in tandem. Both of these algorithms make extensive use of calculus. 
        #Each time that the network is traversed in the forward (or rightward) direction, the error of the network can be calculated as the difference between the output produced by the network and the expected ground truth, by means of a loss function (such as the sum of squared errors (SSE)). The backpropagation algorithm, then, calculates the gradient (or the rate of change) of this error to changes in the weights. In order to do so, it requires the use of the chain rule and partial derivatives. 
        #For simplicity, consider a network made up of two neurons connected by a single path of activation. If we had to break them open, we would find that the neurons perform the following operations in cascade:
            #The neuro receives input, which is multiplied by the weight and then bias is sum. Activation is applied and the output is sent to the new neuron and multiply by its weight.
            #This connects the errors...
            #This is suppose to be connected to the fact that deirvatives are used, but I do not fully understand this. What I got in this point that the network does small changes in the weights and this produce an overall error when generating the output and comparing the observed value.
            #The loss (difference between predicted and observed) is minimized by doing small changes in the weights and check whether we have less error now.
    #link
        #https://machinelearningmastery.com/calculus-in-action-neural-networks/

print("define a function to generate DNNs to be used in scikeras. This will be an input argument for KerasRegressor()")
#We are going to use scikeras, a wrapper for using keras on scipy ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). You can use it to run keras neural networks in pipelines and then use gridsearch to perform parameter optimization ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). As we will see below, we can optimize any parameter of the network. Not only the loss function or the learning rate, but also the number of layers and units, the dropout rate... just by using loops and ifelse within the function to create the model. Therefore, we can select the best combination of hyperparamenters to maximize predictive prower in the different evaluation datasets. Of course, as we can use pipelines, we can apply the scaling in response and predictors separately in each training-evaluation set.
#First, we are going to create a function to get the keras models. This function will be in turn used as input in scikeras.KerasRegressor().
    #https://coderzcolumn.com/tutorials/artificial-intelligence/scikeras-give-scikit-learn-like-api-to-your-keras-networks#1
    #https://www.adriangb.com/scikeras/stable/quickstart.html
from tensorflow.keras.layers import Dropout # type: ignore
from tensorflow.keras.constraints import MaxNorm # type: ignore
from tensorflow.keras import regularizers # type: ignore
def get_neural_reg(meta, n_layers, n_units, activation, init_mode, dropout_rate, weight_constraint, regu_L1, regu_L2):
    #note that meta is a dict with input metadata. I think this generated by KerasRegressor automatically once the input data is included
    
    #we import keras inside the function to avoid problems with parallelization
    from tensorflow import keras # type: ignore
        #https://stackoverflow.com/questions/42504669/keras-tensorflow-and-multiprocessing-in-python/42506478#42506478

    #start the sequential model
    model = keras.Sequential()
    
    #add the input layer
    model.add(keras.layers.Input(shape=(meta["n_features_in_"],)))
        #input shape obtained from meta, a dict containing input metadata
        #The error you're encountering indicates that the shape provided to the Input layer is not in the correct format. The shape parameter should be a tuple, even if it contains only one dimension.

    #for each of the layer sizes we have
    for i in range(n_layers):
        
        #add layer
        model.add(keras.layers.Dense(
            units=n_units, \
            activation=activation, \
            kernel_initializer=init_mode, #method to add initial weight values \
            kernel_constraint=MaxNorm(weight_constraint), \
            kernel_regularizer=regularizers.L1L2(
                l1=regu_L1, \
                l2=regu_L2)))
        model.add(Dropout(dropout_rate))        
            #You can add a parameter for the dropout rate and add dropout layers after some layers within the loop (using ifelse). You can also stop the network at certain level of accuracy, i.e., callbacks.
                #https://machinelearningmastery.com/weight-initialization-for-deep-learning-neural-networks/
                #https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#6.-Callbacks

    #add the output layer with just one node (predicting 1 value)
    model.add(keras.layers.Dense(1))
        #No activation function is used for the output layer because it is a regression problem, and you are interested in predicting  numerical values directly without transformation.
            #https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/
    
    #compile the model using the KerasRegressor arguments for
    #loss and optimization
    #model.compile(loss=compile_kwargs["loss"], optimizer=compile_kwargs["optimizer"])
        #to run this line you need to add to the function an argument called "compile_kwargs", but we are going to do this outside of this function
        #According to the author of scikeras, it is easier and safer to allow SciKeras to compile your model for you by passing the loss to KerasRegressor directly (KerasRegressor(loss="")) instead of compiling it inside get_reg. It is less flexible because you can not use custom logic, e.g., if else to select a different learning rate according to the optimizer, but he says it is safer, so we are going to select this approach. In case we need more flexibility, we can do the other approach.
            #https://www.adriangb.com/scikeras/stable/notebooks/MLPClassifier_MLPRegressor.html#2.4-Losses-and-optimizer
            #https://www.adriangb.com/scikeras/stable/advanced.html
    return model
#Now we use KerasRegressor. This takes as input a callable function (get_reg) that returns a keras model (model argument). KerasRegressor also takes some arguments needed to prepare a keras model like the activation function or the optimizer. We must also pass all of the arguments to get_reg as keyword arguments to KerasRegressor. Note that if you do not pass an argument to KerasRegressor, it will not be available for hyperparameter tuning.
dummy_dnn = KerasRegressor(model=get_neural_reg, optimizer="Adam", optimizer__learning_rate=0.001, loss="mse", model__n_layers=2, model__n_units=10, model__activation="relu", model__init_mode="uniform", model__dropout_rate=0, model__weight_constraint=0, model__regu_L1=0, model__regu_L2=0, batch_size=200, epochs=50, verbose=1)
    #https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#2.3-Defining-and-training-the-neural-net-classifier
    #Note about keyword arguments: Keyword arguments (or named arguments) are values that, when passed into a function, are identifiable by specific parameter names.
print(dummy_dnn)
#You can apply the regressor to the data and obtain predictions (user verbose=1 to see losses in each epoch) ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#2.3-Defining-and-training-the-neural-net-classifier)).
dummy_dnn.fit(modeling_data.iloc[:,1:], modeling_data.iloc[:,0])
dummy_dnn.model_.summary()
    #we get exactly the architecture we wanted
#We are going to perform a gridsearch in order to optimize the hyperparameters. SciKeras allows to direct access to all parameters passed to the wrapper constructors. This allows tunning of parameters like hidden_layer_sizes as well as optimizer__learning_rate. 
#The model__ prefix can be used to specify that a paramter is destined only for get_reg. In addition, to differentiate paramters like callbacks which are accepted by both tf.keras.Model.fit and tf.keras.Model.predict you can add a fit__ or predict__ routing suffix respectively ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)).
# We can access the history saved after fitting.
loss_pipe = dummy_dnn.history_
print("print loss")
print(loss_pipe)
# We can see how there is only one loss, meaning that the loss is calculated in the data used as input (X_train), considering it as a whole, not partioning. The splitting in evaluation and training will be done later with gridsearch.
print("set the HPs")
dict_models["neural_nets"]["HPs"] = { \
    "regressor__verbose": [0], \
    "regressor__model": [get_neural_reg], \
    "regressor__optimizer__learning_rate": [0.0006, 0.002], \
    "regressor__optimizer": ["Adamax"], \
    "regressor__model__weight_constraint": [4], \
    "regressor__model__regu_L1": [3.0e-08], \
    "regressor__model__regu_L2": [9.0e-08], \
    "regressor__model__n_units": [1000], \
    "regressor__model__n_layers": [5], \
    "regressor__model__init_mode": ["orthogonal"], \
    "regressor__model__dropout_rate": [0.07], \
    "regressor__model__activation": ["tanh"], \
    "regressor__loss": ["mse"], \
    "regressor__epochs": [630], \
    "regressor__batch_size": [10, 100, 200]} #
print(dict_models["neural_nets"])
#explanation about the specific numbers of the HPs
    #we are going to take advantage of the two optuna attempts we did (runs 01-08 and runs 09-11). 
        #We have invested a looot of computational resources in these analyses, getting good curves of performance. We have done several independent optuna processes, lots of iterations...
        #If after investing all this effort, the resulting network cannot defeat XGBoost solving the same problems (i.e., same folds), then this is not the best class. Remember that we have done much less tunning on XGBoost.
            #Also note that DNNs were tuned using CV with the whole dataset, not a nested CV. This means that the data used for test in the nested CV (this script) has been already used for training for the CV of optuna. Therefore, the scores of the DNN can be too optimistic, there is more risk of overfitting to Yoruba. If, despite this, XGboost can fit better the test data that has not been NEVER seen before and from which has not received any tuning information, then this is clearly the best class.
            #This is ok for the poster. We can just say that we used optuna for DNN only because of the great number of HPs that are important to be considered. A typical gridsearch could take too long and random gridsearch will be likely unable to fully explore the parametric space, so we select an approach that is random but consider information of previous steps. We have to think how to deal with this when analyzing the 25 pops...
        #we are going to take the optimum vale for each HP according to the two optuna attempts. We will further tune just those parameters that differ a lot between the first and the second attempt or when there is no clear peak.
            #The first attempt has more weight because it has much more datapoints and a higher R2.
            #Although we will try to get a compromise in general and check in GS if the results are too different.
    #learning rate
        #first optuna attempt
            #We have a clear peak around 0.002. This is the MOST clear peak.
            #There is a small peak at 0.0006 but it is much smaller and only has R2=0.531 instead of 0.526.
        #second optuna attempt
            #we have to peaks around 0.0002 and 0.005.
        #Decision
            #test 0.0006 and 0.002.
    #optimizer
        #first optuna attempt
            #Adamax, followed by Adam and Nadam.
        #second optuna attempt
            #Adamax, followed by RSMprop and Adam.
        #solution
            #Adamax.
    #weight constrain 
        #first optuna attempt
            #there is a complete plateau, maybe a mild increase at 7 
        #second optuna attempt
            #clear peak in 2
        #solution
            #as this parameter does not seem to influence in the first attempt, we could select any value within the range explored. 
            #we could select 2 and 7, but given we already have too many fits and this parameter is less important than the other two parameters we are tunning here (learning rate and batch size), we are just going to select a value within the plateau and the range of values recommended for this HP (3-4).
            #we will take 4 which is more or less in the middle of the plateau and included in the recommended range.
    #L1 regularization
        #first optuna attempt
            #plateau with the middle being at 3e-8 and then decrease after 1.3e-5
        #second optuna attempt
            #there are two peaks flanking 3e-8
        #solution
            #3e-8
    #L2 regularization
        #first optuna attempt
            #plateau with the middle being at 9e-8 and then decrease after 1.3e-5
        #second optuna attempt
            #there is a peak at 6.1e-5, which is close to the end of the plateau in the first attempt
        #solution
            #9e-8. The peak of the second attempt is already downhill of the plateau in the first attempt. Note that the second attempt did not explore fully the space as the first one, so it is possible that with more iterations it could also get high R2 under lower L2 values.
    #number of units
        #first optuna attempt
            #increase up to 700, then plateau until the limit of the epxlored parametric space (900)
        #second optuna attempt
            #increase up to 1000
        #solution
            #Note that, in the first attempt, we did not explore over 900 units, but we did it in the second attempt, getting the best values over 1000.
            #select 1000 which is the peak of the second attempt and still close to the peak of the first attempt.
    #number layers
        #first optuna attempt
            #peak at 4 and 6
        #second optuna attempt
            #peak at 5
        #solution
            #5. It is not the best in the first attempt, but in the second does very good and the numbers at both sides do also well. If 4 and 6 works, it would be strange that 5 does not work.
    #init mode
        #first optuna attempt
            #orthogonal followed by identity, normal, and truncated normal
        #second optuna attempt
            #orthogonal, lecun uniform and normal
        #solution
            #orthogonal
    #dropout rate
        #first optuna attempt
            #peak at 0.07 and a bit in 0.2, then goes dooown
        #second optuna attempt
            #peak at 0.03
        #solution
            #0.07 which top in the first and very close to the top in the second attempt.
    #model activation
        #first optuna attempt
            #tanh followed by elu and softsign
        #second optuna attempt
            #selu followed by hard sigmoid
        #solution
            #tanh
    #loss function
        #first optuna attempt
            #mse followed by huber and log cosh
        #second optuna attempt
            #mse followed by log cosh and huber
        #solution
            #mse
    #epochs
        #first optuna attempt
            #wide peak around 630, then decrease
        #second optuna attempt
            #peak around 610
        #solution
            #630 so we are close to both peaks
    #batch size
        #first optuna attempt
            #two peaks around 80-100 and 200
        #second optuna attempt
            #peak at 5
        #solution
            #test 10, 100 and 200 to cover all peaks
#about HPs:
    #link 
        #https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/
    #Batch vs. Epoch 
        #https://machinelearningmastery.com/difference-between-a-batch-and-an-epoch/
        #Stochastic Gradient descent: 
            #The job of the algorithm is to find a set of internal model parameters (predictors) that perform well against some performance measure such as logarithmic loss or mean squared error. 
            #Optimization is a type of searching process and you can think of this search as learning. The optimization algorithm is called “gradient descent“, where “gradient” refers to the calculation of an error gradient or slope of error and “descent” refers to the moving down along that slope towards some minimum level of error.
            #The algorithm is iterative. This means that the search process occurs over multiple discrete steps, each step hopefully slightly improving the model parameters.
            #Each step involves using the model with the current set of internal parameters to make predictions on some samples, comparing the predictions to the real expected outcomes, calculating the error, and using the error to update the internal model parameters.
            #This update procedure is different for different algorithms, but in the case of artificial neural networks, the backpropagation update algorithm is used.
            #Stochastic gradient descent is an optimization algorithm that estimates the error gradient for the current state of the model using examples from the training dataset, then updates the weights of the model using the back-propagation of errors algorithm, referred to as simply backpropagation.
            #The amount that the weights are updated during training is referred to as the step size or the “learning rate.”
            #Specifically, the learning rate is a configurable hyperparameter used in the training of neural networks that has a small positive value, often in the range between 0.0 and 1.0.
            #The learning rate controls how quickly the model is adapted to the problem. Smaller learning rates require more training epochs given the smaller changes made to the weights each update, whereas larger learning rates result in rapid changes and require fewer training epochs.
            #A learning rate that is too large can cause the model to converge too quickly to a suboptimal solution, whereas a learning rate that is too small can cause the process to get stuck.
            #The challenge of training deep learning neural networks involves carefully selecting the learning rate. It may be the most important hyperparameter for the model.
        #Sample
            #A sample is a single row of data.
            #It contains inputs that are fed into the algorithm and an output that is used to compare to the prediction and calculate an error.
            #A training dataset is comprised of many rows of data, e.g. many samples. A sample may also be called an instance, an observation, an input vector, or a feature vector.
        #Batch
            #The batch size is a hyperparameter that defines the number of samples (rows) to work through before updating the internal model parameters.
            #Think of a batch as a for-loop iterating over one or more samples and making predictions. In each iteration, a portion of the data is considered and predictions are made considering the current combination of parameters. At the end of the batch, the predictions are compared to the expected output variables and an error is calculated. From this error, the update algorithm is used to improve the model, e.g. move down along the error gradient.
            #A training dataset can be divided into one or more batches.
            #When all training samples are used to create one batch, the learning algorithm is called batch gradient descent. When the batch is the size of one sample, the learning algorithm is called stochastic gradient descent. When the batch size is more than one sample and less than the size of the training dataset, the learning algorithm is called mini-batch gradient descent.
                #Batch Gradient Descent. Batch Size = Size of Training Set
                #Stochastic Gradient Descent. Batch Size = 1
                #Mini-Batch Gradient Descent. 1 < Batch Size < Size of Training Set
            #In the case of mini-batch gradient descent, popular batch sizes include 32, 64, and 128 samples. You may see these values used in models in the literature and in tutorials.
        #Epochs
            #The number of epochs is a hyperparameter that defines the number times that the learning algorithm will work through the entire training dataset.
            #One epoch means that each sample in the training dataset has had an opportunity to update the internal model parameters. An epoch is comprised of one or more batches. For example, as above, an epoch that has one batch is called the batch gradient descent learning algorithm.
            #You can think of a for-loop over the number of epochs where each loop proceeds over the whole training dataset. Within this for-loop is another nested for-loop that iterates over each batch of samples, where one batch has the specified “batch size” number of samples.
            #The number of epochs is traditionally large, often hundreds or thousands, allowing the learning algorithm to run until the error from the model has been sufficiently minimized. You may see examples of the number of epochs in the literature and in tutorials set to 10, 100, 500, 1000, and larger.
            #It is common to create line plots that show epochs along the x-axis as time and the error or skill of the model on the y-axis. These plots are sometimes called learning curves. These plots can help to diagnose whether the model has over learned, under learned, or is suitably fit to the training dataset.  
        #Batch vs. Epoch
            #The batch size is a number of samples processed before the model is updated.
            #The number of epochs is the number of complete passes through the whole training dataset.
            #The size of a batch must be more than or equal to one and less than or equal to the number of samples in the training dataset.
            #The number of epochs can be set to an integer value between one and infinity. You can run the algorithm for as long as you like and even stop it using other criteria besides a fixed number of epochs, such as a change (or lack of change) in model error over time.
            #They are both integer values and they are both hyperparameters for the learning algorithm, e.g. parameters for the learning process, not internal model parameters found by the learning process.
            #You must specify the batch size and number of epochs for a learning algorithm.
            #There are no magic rules for how to configure these parameters. You must try different values and see what works best for your problem.
        #Worked Example
            #Assume you have a dataset with 200 samples (rows of data) and you choose a batch size of 5 and 1,000 epochs.
            #This means that the dataset will be divided into 40 batches (5*40=200), each with five samples. The model weights will be updated after each batch of five samples, i.e., in each learning step, only 5 random rows will be considered.
            #This also means that one epoch will involve 40 batches or 40 updates to the model.
            #With 1,000 epochs, the model will be exposed to or pass through the whole dataset 1,000 times. That is a total of 40,000 batches during the entire training process.
            #If the sample size is not divisible by the batch size without remainder, the last batch would have a smaller number of samples (https://machinelearningmastery.com/difference-between-a-batch-and-an-epoch/). I guess this is ok because we have already fitted the model multiple times with the previous batches and we will repeate this several epochs. In addition, it is usual to tune the batch size testing different sizes, so people usually have this situation.
        #Summary
            #Stochastic gradient descent is an iterative learning algorithm that uses a training dataset to update a model.
            #The batch size is a hyperparameter of gradient descent that controls the number of training samples to work through before the model’s internal parameters are updated.
            #The number of epochs is a hyperparameter of gradient descent that controls the number of complete passes through the training dataset. 
            #The idea behind is that most observations in a large data set will have many neighbors that impart almost the same information. We don't need all of them in our learning process. 
            #BUT CAREFUL IN YOUR CASE, because you do not have millions of samples, although using multiple epochs you will go through the data many times using maany batches.
        #Epochs and batches (number of rows in each step) are not a replacement of cross-validation. It is an iterative process that works better in a subset of the data if you have a looot of data, because you will probably have repeated data. You can increase the batch size until the time of each step starts to increase. By default, tensorflow uses no replacement by default for batches. 
        #Notes about the search I initially did in optuna: 
            # We are only using mini-batches, i.e., batch size>1 and < sample size. Using the whole sample size could make a lot of use of memory and I think batch size=1 could have problems to converge. I have used still a wider range than people use for minibatches.
            #Modification for the second run: The best trials have a batch size between 10-20, so we are in the lower limit of this parameter. I am going to extend the space reaching Stochastic Gradient Descent, i.e., batch size=1. Note that the previous space was 10 to 220, but best models appears very close to 10. Maybe, we do not have so much data, so there is no so high redundancy. I am also using smaller steps to have more resolution.
                #Indeed, best values tend to be around a batch size of 5.
            # The epoch sizes usually go from 1 to 1000, I am not using the whole range to reduce computation time. In the second run I am extending the upper limit because best trials were not far from that 
    #Optimization method 
        # There is a lot of literature and discussion about the selection of the optimization. Adam seems to be an approach that combine the strengths of other methods (i.e., Adagrad and RMSProp; [link](https://towardsdatascience.com/a-visual-explanation-of-gradient-descent-methods-momentum-adagrad-rmsprop-adam-f898b102325c); [link](https://datascience.stackexchange.com/questions/10523/guidelines-for-selecting-an-optimizer-for-training-neural-networks)). There is controversy about it with some results showing poor performance, but then other results show good performance (similar to Stochastic gradient descent + momentum) when doing a well parameter optimization ([link](https://www.fast.ai/posts/2018-07-02-adam-weight-decay.html), [link](http://ruder.io/optimizing-gradient-descent/index.html#adam)).
        # Advantages of Adam ([link](https://machinelearningmastery.com/adam-optimization-algorithm-for-deep-learning/)):
            #Straightforward to implement.
            #Computationally efficient.
            #Little memory requirements.
            #Invariant to diagonal rescale of the gradients.
            #Well suited for problems that are large in terms of data and/or parameters.
            #Appropriate for non-stationary objectives.
            #Appropriate for problems with very noisy/or sparse gradients.
            #Hyper-parameters have intuitive interpretation and typically require little tuning.
            #The authors describe Adam as combining the advantages of two other extensions of stochastic gradient descent. Specifically:
                #Adaptive Gradient Algorithm (AdaGrad) that maintains a per-parameter learning rate that improves performance on problems with sparse gradients (e.g. natural language and computer vision problems).
                    #Our selective pressure variables has many zeros or values close to zero, being the most important values, so an approach like this that try not to underwegight this predictors could be useful.
                #Root Mean Square Propagation (RMSProp) that also maintains per-parameter learning rates that are adapted based on the average of recent magnitudes of the gradients for the weight (e.g. how quickly it is changing). This means the algorithm does well on online and non-stationary problems (e.g. noisy).
                    #We could consider our problem as noisy becasue we have multiple factors influencing selection.
        # In a review of optimizers ([link](https://arxiv.org/abs/1609.04747)), they say "Insofar, RMSprop, Adadelta, and Adam are very similar algorithms that do well in similar circumstances. […] its bias-correction helps Adam slightly outperform RMSprop towards the end of optimization as gradients become sparser. Insofar, Adam might be the best overall choice."
        # It seems this approach work fast and good for non-shallow problems, being usually used as first option currently.
        # We are going to consider different optimizers just in case this makes any difference. 
        #Note about the search done with optuna: 
            #There are two optimizers (AdamW and Adafactor) that cannot be used as input strings in kerasregressor. We can use a class but that cannot be optimized in the tunning process. For doing that, I should compile the model inside get_reg() and add an argument for optimizer. Then include an Ifelse, if optimizer="AdamW", uses class AdamW, which has been previoulsy used. However, scikeras author recommend to compile the model outside, i.e., let kerasregressor to do it.
            #Given that these are only two optimizers and we already tested a battery of optimizers, including several Adam versions, we are going to skip this for now.
    #Adam parameters ([link](https://machinelearningmastery.com/adam-optimization-algorithm-for-deep-learning/))
        #alpha. Also referred to as the learning rate or step size. The proportion that weights are updated (e.g. 0.001). Larger values (e.g. 0.3) results in faster initial learning before the rate is updated. Smaller values (e.g. 1.0E-5) slow learning right down during training
        #beta1. The exponential decay rate for the first moment estimates (e.g. 0.9).
        #beta2. The exponential decay rate for the second-moment estimates (e.g. 0.999). This value should be set close to 1.0 on problems with a sparse gradient (e.g. NLP and computer vision problems).
        #epsilon. Is a very small number to prevent any division by zero in the implementation (e.g. 10E-8).
        # The Adam paper suggests:
            #Good default settings for the tested machine learning problems are alpha=0.001, beta1=0.9, beta2=0.999 and epsilon=10−8
            #The TensorFlow documentation suggests some tuning of epsilon:
                #The default value of 1e-8 for epsilon might not be a good default in general. For example, when training an Inception network on ImageNet a current good choice is 1.0 or 0.1.
            #Keras uses these values as default except epsilon, being 10-7.
        #They say that the default configuration parameters do well on most problems, so we are just going to tune the learning rate alpha. It seems that beta should not be change except to change specific reasons for your data ([link](https://stats.stackexchange.com/questions/499013/adam-adaptive-optimizers-learning-rate-tuning?rq=1)). Maybe in the future we can tune the betas usng bayesian optimization ([link](https://stats.stackexchange.com/questions/265400/deep-learning-how-does-beta-1-and-beta-2-in-the-adam-optimizer-affect-its-lear)). We could even not tune learning rate because it is adaptive in Adam (it changes along the process), but we are going to do it just in case, because it is an important parameter, and it could influence the initila learning rate before updating the parameters. In addition, we are using other optimizers in which learning rate could be more relevant, so it is important to tune this paramater.
        #More general info about tunning learning rate ([link](https://www.bdhammel.com/learning-rates/)).
            # Generally, it is a good idea to also include the number of epochs in an optimization like this as there is a dependency between the amount of learning per batch (learning rate), the number of updates per epoch (batch size), and the number of epochs ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)).
            # We do that as we are also tunning the number of batches and epochs.
        #Note after optuna: 
            # We are not optimizing some specific parameters of the optimizers, like epsilon, etc... In order to do that, we should compile the model inside get_reg() and add ifelse for the value of a new "optimizer" argument. For example, if optimizer="Adam", use Adam(epsilon_...) or something like that. So the function only take epsilon for Adam, but not sure if this would work if you are not using adam and epsilon is still a parameter in keras regressor. 
            # We are going to stick to default parameters (except learning rate) for optimizers for now. Maybe after the big search, having already an optimizer selected, we can optimize its own parameters.
            #I have seen that people using optuna (bayesian) sample learning rate values
            #considering the algorithm scale, I think to prioritize small values. Our reference tutorial goes in the same line with grid searchs typically consisting in picking numbers between 10^−5 and 0.3 on a logaritmic scale ([link](https://machinelearningmastery.com/learning-rate-for-deep-learning-neural-networks/)). The wider range would be 1 - 10^−6: "Typical values for a neural network with standardized inputs (or inputs mapped to the (0,1) interval) are less than 1 and greater than 10^−6". 
            #Given we are going to use the narrower (but still wide) range of 10^−5 to 0.3. If we see that the best models are close to the limits, we can run a more detailed search around these values.
            #the best trials are not veeery close but still close to the lower limit, so we are extending the lower limit in the second run (1e-6).
    #Weight initiallization ([link](https://machinelearningmastery.com/weight-initialization-for-deep-learning-neural-networks/))
        # The optimization algorithm requires a starting point in the space of possible weight values from which to begin the optimization process. Weight initialization is a procedure to set the weights of a neural network to small random values that define the starting point for the optimization (learning or training) of the neural network model.
        # " training deep models is a sufficiently difficult task that most algorithms are strongly affected by the choice of initialization. The initial point can determine whether the algorithm converges at all, with some initial points being so unstable that the algorithm encounters numerical difficulties and fails altogether"
        # Neural network weight initialization used to be simple: use small random values ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)). Now there is a suite of different techniques to choose from. Keras provides a laundry list.
        # Ideally, it may be good to use different weight initialization schemes according to the activation function used on each layer. In our case, however, we will use the same activation across layers because we are working with regression, not classification. Maybe if we use a initialization method that is not good for tahn, that combination will not work and will have low R2, but we will also run other DNNs with tahn and other initialization methods, and other with relu...
    #Activation function ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/))
        # The activation function controls the non-linearity of individual neurons and when to fire ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)). 
    #Dropout rate ([link](https://machinelearningmastery.com/dropout-for-regularizing-deep-neural-networks/))
        # We are going to add dropout so some neurons are shutdown. In this way, we will try to get more generalization. 
        # During training, some number of layer outputs are randomly ignored or “dropped out.” This has the effect of making the layer look-like and be treated-like a layer with a different number of nodes and connectivity to the prior layer. In effect, each update to a layer during training is performed with a different “view” of the configured layer.
        # Dropout has the effect of making the training process noisy, forcing nodes within a layer to probabilistically take on more or less responsibility for the inputs. 
        # The idea behind this is to keep individual neurons from becoming too specialized.  Because each neuron may or may not be in any given run, it cannot be the only neuron detecting a particular feature.  If that feature is important, responsibility for its detection needs to be spread out among several neurons.  When we make a prediction, we will keep all of the neurons, in order to make the best prediction possible.
        # This conceptualization suggests that perhaps dropout breaks-up situations where network layers co-adapt to correct mistakes from prior layers, in turn making the model more robust.
        # Dropout may be implemented on any or all hidden layers in the network as well as the visible or input layer. It is not used on the output layer.
        # A common value is a probability of 0.5 for retaining the output of each node in a hidden layer and a value close to 1.0, such as 0.8, for retaining inputs from the visible layer.
        # When using dropout regularization, it is possible to use larger networks with less risk of overfitting. In fact, a large network (more nodes per layer) may be required as dropout will probabilistically reduce the capacity of the network.
        # A good rule of thumb is to divide the number of nodes in the layer before dropout by the proposed dropout rate and use that as the number of nodes in the new network that uses dropout. For example, a network with 100 nodes and a proposed dropout rate of 0.5 will require 200 nodes (100 / 0.5) when using dropout.
        # Dropotu can be applied to the input layer, so some samples are removed ([link](https://machinelearningmastery.com/dropout-regularization-deep-learning-models-keras/)), but I have not seens this before in my previous courses so I will pass.
    #weight constrain
        # You also need weight constrain:
        # Network weights will increase in size in response to the probabilistic removal of layer activations.Large weight size can be a sign of an unstable network.
        # To counter this effect a weight constraint can be imposed to force the norm (magnitude) of all weights in a layer to be below a specified value. For example, the maximum norm constraint is recommended with a value between 3-4.
        # Constraining the weight matrix directly is another kind of regularization. If you use a simple L2 regularization term (see below) you penalize high weights with your loss function. With this constraint, you regularize directly ([link](https://www.kdnuggets.com/2015/04/preventing-overfitting-neural-networks.html/2), [link](https://stackoverflow.com/questions/45970888/what-does-kernel-constraint-max-norm3-do)). Therefore, this is another layer of regularization that is considered. Our search could combine it with the rest of thechiques or not, depending on the combination.
        # We will try dropout percentages between 0.0 and 0.9 (1.0 does not make sense) and maxnorm weight constraint values between 0 and 5.    
        # **Note**
            # We are going to apply the same dropout rate for all inner layers. Once we have an optimized architecture, we can try to improve it by setting the dropout only in specific layers.
    #Early Stop
        # The idea behind early stop is to stop the fitting process if a given metric does not improve after several steps. 
        # Ideally, you would monitor the loss/accuracy of the validation dataset, so you can if the validation dataset is not improving and hence we would have a high risk of overfitting (i.e., getting better predictions in training than in validation). 
        # This can be implemented in KerasRegressor with the argument callbacks, for example calling the class EarlyStopping, and then add routed parameters like callbacks__patience (how many steps we have to wait to stop if there is no improvement).
        # There is, however, a **problem** if we use this in combination with cross_val_score, grid search.... There is no validaition set identified in KerasRegressor (the validation sets are created are created by grid search), so we can only monitor the loss/accuracy in the training dataset.
        # There are also criticisms to combine CV and early stopping ([link](https://stackoverflow.com/questions/48127550/early-stopping-with-keras-and-sklearn-gridsearchcv-cross-validation)), because you could loss your best model because stopping early. I am not sure about that, but there is an more important point:
        # You are already comparing models with different number of epocs, so if there is overfitting after some epocs, the model with more epocs will get a lower R2 in the validation dataset. 
        # In summary, we are not going to use early stopping to avoid overfitting. We stick to the selection of the best number of epochs according the validation set.
    #Regularization
        # Overfitting is caused by the model having too much flexibility, and the main source of flexibility in a neural network is the weights in the neurons. Specifically, non-zero weights indicate some relationship between the input and the output. Regularization penalizes this flexibility by penalizing non-zero weights. Thus, the network will only have a weight be non-zero if the benefit to the loss function is greater than the penalty applied for the weight.
        # There are two main types of regularization:  𝐿2 -regularization adds a penalty proportional to the sum of the squares of the weights, while  𝐿1 -regularization uses the sum of the absolute values of the weights. (The biases are generally not regularized.). Therefore, making less likely to have non-positive weights and hence associations between predictors and target. In other words, it penalizes “big” weights.
        # The hyperparameter alpha is the regularization parameter. Its size needs to be set to provide the right amount of flexibility that the net avoids both overfitting and its converse, underfitting. In general, you will need to do a bit of a search to find the appropriate value for your problem. Info from Optizimation notebook of TDI.
        # If the hyperparameter is zero, then the L1/2 sum goes to zero, and there is no regulization ([link](https://towardsdatascience.com/intuitions-on-l1-and-l2-regularisation-235f2db4c261)).
        # A linear regression model that implements L1 norm for regularisation is called lasso regression, and one that implements (squared) L2 norm for regularisation is called ridge regression ([link](https://towardsdatascience.com/intuitions-on-l1-and-l2-regularisation-235f2db4c261)). We are going to use L1L2 function so you can have both or one of the two most frequent types of regularization setting the corresponding parameter l1 and l2 for each regularization type ([link](https://keras.io/api/layers/regularizers/)). The default is 0.01 in both cases. 
        # Three different regularizer instances are provided; they are ([link](https://machinelearningmastery.com/how-to-reduce-overfitting-in-deep-learning-with-weight-regularization/)):
            #L1: Sum of the absolute weights.
            #L2: Sum of the squared weights.
            #L1L2: Sum of the absolute and the squared weights.
        # This is achieved by setting the kernel_regularizer argument on each layer. A separate regularizer can also be used for the bias via the bias_regularizer argument, although this is less often used.
        # The most common type of regularization is L2, also called simply “weight decay,” with values often on a **logarithmic scale** between 0 and 0.1, such as 0.1, 0.001, 0.0001, etc. It is a good practice to first grid search through some orders of magnitude between 0.0 and 0.1, then once a level is found, to grid search on that level.
        # Note that in the paper introducing dropout, this techinque was combined with L2 regularization ([link](https://stats.stackexchange.com/questions/241001/deep-learning-use-l2-and-dropout-regularization-simultaneously)). In any case, we have explored in optuna the possibilty of having both or one of them by tunning their corresponding hyperparameters.    
        #Note during optuna optimization
            # In the second search, we can make a finer search around the values of the best models.
            #We have extended the lower limit for the second search because the best trials were close to that limit
    #Batch normalization
        # Another [recently-developed](https://arxiv.org/pdf/1502.03167v3.pdf) tool for deep networks is **batch normalization**.  Although it can help with overfitting, it was originally developed to deal with the vanishing gradient problem.  Recall that activation functions have flat regions, where their gradients are small.  When the input is in these regions, gradient descent will only move the weights small amounts, leaving them stuck in the low-gradient regions.  Intelligent choices for initializations and activation functions try to avoid this as much as possible.
        # Batch normalization takes a more proactive approach, scaling and shifting the inputs so that the average input, over the whole batch, has a target mean and standard deviation.  These target values become parameters of the model, tuned during training.
        # By keeping gradients from vanishing, batch normalization reduces the importance of the weight initialization and the activation function.  Larger learning rates can be used.  Although the initial steps may proceed more slowly, as the correct normalizations must be learned, learning should proceed much faster overall than without.  Batch normalization can also have a regularization effect, reducing the propensity towards overfitting!
        # From TDI optimization notebook. 
            # **Don’t Use With Dropout** ([link](https://machinelearningmastery.com/batch-normalization-for-training-of-deep-neural-networks/), [link](https://stackoverflow.com/a/62806906/12772630))
            # Batch normalization offers some regularization effect, reducing generalization error, perhaps no longer requiring the use of dropout for regularization.
            # Removing Dropout from Modified BN-Inception speeds up training, without increasing overfitting.
            # — Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift, 2015.
            # Further, it may not be a good idea to use batch normalization and dropout in the same network.
            # The reason is that the statistics used to normalize the activations of the prior layer may become noisy given the random dropping out of nodes during the dropout procedure.
            # Batch normalization also sometimes reduces generalization error and allows dropout to be omitted, due to the noise in the estimate of the statistics used to normalize each variable.
                #Page 425, Deep Learning, 2016.
            #**Decision:** We are not going to use batch normalization for now.
    #Number of inner layers and neurons ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/))
        # The number of neurons in a layer is an important parameter to tune. Generally the number of neurons in a layer controls the representational capacity of the network, at least at that point in the topology.
        # Also, generally, a large enough single layer network can approximate any other neural network, at least in theory.
        # A larger network requires more training and at least the batch size and number of epochs should ideally be optimized with the number of neurons (we are doing that).
        # It is recommended in general to start coarse and then fine tune the model, so we are going to compare very different number of layers and neurons.
        #note optuna
            #For the second run, we are extending the maximum number of nodes
            #the best number until now is between 300-400, but some trials
            #are also high close to 500, so maybe the optimum
            #is bigger.
            #in the 10th run we have a lot of good runs around 1000, so we are extending the upper limit a lot.
            #in the 11st run we have a lot of good runs around 1200, so we are extending the upper limit a lot.
        # **Note optuna**
            # We are not looking to networks deeper than 20 layers because computational time. If we see that the number of layers or neurons of the best models is close to the limit, we can extedn this in the finer search.
            # Similarly, in the second search we can try to change a bit the number of units between layers.
    #Losses ([link](https://machinelearningmastery.com/how-to-choose-loss-functions-when-training-deep-learning-neural-networks/))
        # Inside the network, the difference between observed and predicted will be calculated using the loss function, and this information will be used to make the next step while considering the learning rate, etc... 
        # These metrics can be used also as accuracy metrics in the CV, but for the cross-validation, we will use r2 (see below).
        #The Mean Squared Error, or MSE,
            #default loss to use for regression problems. 
            #Mathematically, it is the preferred loss function under the inference framework of maximum likelihood if the distribution of the target variable is Gaussian. It is the loss function to be evaluated first and only changed if you have a good reason.
            #Mean squared error is calculated as the average of the squared differences between the predicted and actual values. The result is always positive regardless of the sign of the predicted and actual values and a perfect value is 0.0. The squaring means that larger mistakes result in more error than smaller mistakes, meaning that the model is punished for making larger mistakes.
            #Mean Squared Logarithmic Error Loss
                #There may be regression problems in which the target value has a spread of values and when predicting a large value, you may not want to punish a model as heavily as mean squared error.
                #Instead, you can first calculate the natural logarithm of each of the predicted values, then calculate the mean squared error. This is called the Mean Squared Logarithmic Error loss, or MSLE for short.
                #It has the effect of relaxing the punishing effect of large differences in large predicted values.
                #As a loss measure, it may be more appropriate when the model is predicting unscaled quantities directly. Nevertheless, we can demonstrate this loss function using our simple regression problem.
                #The model can be updated to use the ‘mean_squared_logarithmic_error‘ loss function and keep the same configuration for the output layer. We will also track the mean squared error as a metric when fitting the model so that we can use it as a measure of performance and plot the learning curve.
            #Mean Absolute Error
                #On some regression problems, the distribution of the target variable may be mostly Gaussian, but may have outliers, e.g. large or small values far from the mean value.
                #The Mean Absolute Error, or MAE, loss is an appropriate loss function in this case as it is more robust to outliers. It is calculated as the average of the absolute difference between the actual and predicted values.
                #The model can be updated to use the ‘mean_absolute_error‘ loss function and keep the same configuration for the output layer.
        
            # There are more losses for regression like Mean Absolute Percentage Error (Mape) or Mean Squared Logarithmic Error (square(log(y_true + 1.) - log(y_pred + 1.))). I used all that are not dedicated to categorical responses in optuna ([link](https://www.tensorflow.org/api_docs/python/tf/keras/losses), [link](https://towardsdatascience.com/understanding-loss-functions-the-smart-way-904266e9393)).
                #"huber", 
                    #In statistics, the Huber loss is a loss function 
                        #used in robust regression, that is less sensitive to
                        #outliers in data than the squared error loss.
                        #it seems to get the strengths of MAE and MSE without 
                        #their weaknesess
                        #https://towardsdatascience.com/understanding-loss-functions-the-smart-way-904266e9393
                #"log_cosh"]
                    #log(cosh(x)) is approximately equal to (x ** 2) / 2 for small x and 
                    #to abs(x) - log(2) for large x. This means that 'logcosh' 
                    #works mostly like the mean squared error, but will not be so 
                    #strongly affected by the occasional wildly incorrect prediction.
    #Scorer for tuning
        # We are going to select ONLY ONE scorer in order to perform the hyperparameter tuning.
        # In regression problems you can not use the same accuracy metrics as in classification problems (e.g. error rate, confusion matrix, etc.): in stead, other metrics are used like:
            # - Pearson linear correlation
            # - Spearman rank correlation
            # - RMSE (root mean squared error)
            # - MAE (mean absolute error)
            # - etc. (there are many more)
        # The R2 is the proportion of the variation in the dependent variable that is predictable from the independent variable(s) ([link](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html#sklearn.metrics.r2_score), [link](https://en.wikipedia.org/wiki/Coefficient_of_determination)). 
        # I have used this metric multiple times including in this very project and I saw a great correlation between R2 and the fit to the distribution. For example, in cases of extreme overfitting, i.e., the distribution of predicted is exactly the same than observed, the R2 calculated with scikitlearn is 1. It is also much easier to interpet and widely use both in machine learning and science in general.
        # We are going to use this metric for now.
    #Future:
        #do like in PRS deep learning paper, you create networks with increasing or decreasing number of neurons in each layer
            #see paper
        #tune parameters for the optimizer (Adamax)
        #dropout in the input layer? 
            #see notes
        # We are going to apply the same dropout rate for all inner layers. Once we have an optimized architecture, we can try to improve it by setting the dropout only in specific layers.
        #remove Dropout and add batch normalization?

# endregion






###################################
# region prepare CV scheme ########
###################################

print_text("prepare CV scheme", header=2)
#It is VERY important that we avoid overfitting given the use we are going to make of the models.
    #In the future, we will use the final model to obtain a probability of selection considering genomic factors and then use it to select genes with the same expected probability of selection (according to these factors) than our genes of interest. The idea is that the interest genes should have the same probability of selection based on genomic features (predicted probability), but if they are target of a selective pressure, their observed probability of selection should be higher. If predicted and observed are exactly the same, there is no room for enrichment of selection in interest genes after controling for confounding factors, and this would be an methodological artifact due to a model that just fit the observed data without any generalization. This can be a problem with algorithms like RF or in deep learning. In other words, more overfitting, less power to detect the impact of selective pressures.

#Held out part of the dataset in each CV round (i.e., test set) 
    #This part will not be used for parameter optimization, but for the final validation after the final model has been optimized. In this way, we avoid potential overfiting in the evaluation sets ([see link](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation)). If you use train in a set of the data and then evaluate in other, you can see if the model trained is not overfitting to the training data and is flexible enough to predict the evaluation data. The problem is that in parameter optimization, we select the best parameters based on the evaluation metrics in the evaluation sets, thus we could get a model that fit too much the evaluation dataset, loosing generalization and thus making the evaluation metrics no longer metrics of generalization. To avoid this, we leave out a set of the data for final evaluation, i.e., the test set. This set will be NOT used in parameter optimization. Once we have selected the best parameters to train a model in the training dataset and predict well in the evaluation datasets of the CV, we use these parameter to create a final model, fit to the whole training data and then predict in the final evaluation dataset, which was not used for anything before. If the model works well, it means it is generalizable and it is not overfitting the data, so, in our case, we can say that it is explaining the variance in selection that it is really explained by the genomic factors and the rest would be variance that could be explained by selective pressures. If there is overfitting, the model fit too much the data, there is not non-explained variance.

#we are going to use nested Cross-Validation, so we can tune hyperparameters and select the best model class without having too optismistic results about model performance. If the same dataset is used to tune hyperparameters and then evaluate the tuned models to compare between model classes, there is a risk of overfitting.
    #https://machinelearningmastery.com/nested-cross-validation-for-machine-learning-with-python/

#in a nested schema, 
    #we create an outer CV
        #create several folds, for example 10
        #You use 9 for training and 1 for evaluatioon
        #using the 9 folds, you perform another CV, this time with less folds
            #from this 90% of the data, we now split in 3 folds
            #take two for training under different hyperparameters
            #evaluate in the third
            #select the best combination of parameter
            #and train using the whole 90% of the data (i.e., the three parts)
        #Use the trained model to predict in the remaining fold from the outer CV, i.e., 10%.
        #This is repeated for each of the 10 folds, i.e, we obtain a metric of performance for each 10% set of the data.
        #we can use these metrics to compare model classes and the hyperparameter optimization was done with a different part of the data.
        #Under this procedure, hyperparameter search does not have an opportunity to overfit the dataset as it is only exposed to a subset of the dataset provided by the outer cross-validation procedure. This reduces, if not eliminates, the risk of the search procedure overfitting the original dataset and should provide a less biased estimate of a tuned model’s performance on the dataset.

#we will use k=10 for now
    #it has been seen in many datasets that 10 folds works relatively well to estimate model performance ([link](https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/#:~:text=configure%20the%20procedure.-,Sensitivity%20Analysis%20for%20k,evaluate%20models%20is%20k%3D10.)). Although there is some debate about it. We will stick to 10 folds for now
    #This is very well explained here (https://stackoverflow.com/questions/45969390/difference-between-stratifiedkfold-and-stratifiedshufflesplit-in-sklearn).

#In the future you could compare different K values and select the best 
    #you could do a sensitivity analysis with elastic net
        #use default parameters of elastic net
        #then run several nested CVs with different number of folds for the inner and outer CVs
        #calculate the average R2 and min-max in the outer folds for each combination of inner-outer
        #select the number of inner and outer folds that give the greatest R2
    #then select that combination and check what happens with the other models (also in default)
    #This will automatically select the best number of inner and outer folds, thus also selecting the selecting the size of the training/evaluation/test sets. For example, 10 outer folds and 3 inner folds means that the test set will be 10% of the data while 90% will be used for training and evaluation. Within that 90%, 2/3 will be used for training and 1/3 for evaluation.
    #link
        #https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/

#Respect to the repetition because of the stochastic nature of the machine learning models
    #If a large number of machine learning algorithms and algorithm configurations are being evaluated systematically on a predictive modeling task, it can be a good idea to fix the random seed of the evaluation procedure. Any value will do.
    #The idea is that each candidate solution (each algorithm or configuration) will be evaluated in an identical manner. This ensures an apples-to-apples comparison. It also allows for the use of paired statistical hypothesis tests later, if needed, to check if differences between algorithms are statistically significant.
    #This what we are gonna do, as we just want to select a model class. We can consider stochasiticity when working with the selected class, running different models with different seeds and get the average prediction. We will likely have several good models, not only one.
        #https://machinelearningmastery.com/different-results-each-time-in-machine-learning/



print_text("define the outer CV", header=3)
print_text("create cross-validator", header=4)
from sklearn.model_selection import KFold
cv_outer = KFold( \
    n_splits=10,  \
    shuffle=True)
    #If you select 10 splits, this means you have 10 test sets. In each split, you have 9 folds for training and 1 for test.
    #random_state is not required as we have already ensured reproducibility by setting the seeds of python, numpy and tensorflow
    #KFold vs. ShuffleSplit
        #We will use KFold because we want to have training and validation datasets that are not overlapped so we can use all the data. KFold with shuffle=True only shuflle one time at the beginning and then make the folds, while ShuffleSplit shuffles every time making possible to select as evaluation two times the same sample. In the case of 5 splits with KFold, one time 1/5 is the validation and the rest is for training, then another 1/5 and so on... until the 5 partitions have been used for validation, i.e., the whole dataset have been used for validation.
        #It is important to us to use the whole training dataset for CV because we need to obtain models that are generalizable enough so we do not overestimate the influence of genomic confounding factors on selection due to overfitting (see above).
print(cv_outer)

print_text("extract the indices for each train-test set", header=4)
indexes_cv_outer = []
#split, train_index, test_index = [(split, train_index, test_index) for split, (train_index, test_index) in enumerate(cv_outer.split(X=modeling_data.iloc[:,1:])) if split==0][0]
    #split Generate indices to split data into training and test set.
        #we get row indexes from the dataset (i.e., genes). Given we have the same genes for y and X, we can just use X (predictors) to get the indexes of rows for training and test
    #for each split
        #get the training/test indexes along with the number of the split (using enumerate for that)
        #select only the first split for debugging
for split, (train_index, test_index) in enumerate(cv_outer.split(X=modeling_data.iloc[:,1:])):
    for model_class in dict_models.keys():
        #save the same train/test indexes for all model class
        indexes_cv_outer.append((split, train_index, test_index, model_class))
print(indexes_cv_outer)
print("Do we have the correct number of splits*models combinations?")
print(len(indexes_cv_outer) == cv_outer.n_splits * len(dict_models.keys()))
    #that number should be the number of splits times the number of model classes
    
#endregion






#################################################################
# region Define function and run it only for elastic net ########
#################################################################

print_text("Define function and run it only for elastic net", header=2)
from sklearn.pipeline import Pipeline
from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score
#split, train_index, test_index, model_class = indexes_cv_outer[0]
#split, train_index, test_index, model_class = indexes_cv_outer[3]
#IMPORTANT:
    #In the case of DNNs, I have the Hyperparameter search has been narrowed based on the best values for predicting in Yoruba with optuna. DNNs are too slow to be train, so running a GridSearch with all possible combinations of HPs is not possible. For the rest of model classes is not cases, we have used a relatively wide range of values, without narrowing the search.
    #This means that the search for DNNs in each population should be different! We are not going to tune manually each pop, so we need to find a way to do this even for DNNs
        #reduce list of HPs
            #we could select a reduced list of DNN HPs only for model class comparison, using more HPs later if DNNs are selected. See for example (https://academic.oup.com/bioinformatics/article/34/9/1538/4747884).
            #Indeed, this is what we do for RF and XGBoost. We have selected the most important parameters. In the final selected class, we can do a more detailed search (with optuna?). 
        #Make model class selection with optuna to deal with the high number of HPs in DNNs
            #Run optuna search for each model and fold during the model class selection. 
            #We know that the three first model classes can be run fast, so we should not have any problem using optuna with them.
            #Maybe we could increase the number of iterations/processes for DNNs if we use more HPs.
        #Reduce the list of HPs in DNNs and also use optuna for all model classes
def model_evaluation(split, train_index, test_index, model_class):

    print_text(f"Starting with split: {split}, model: {model_class}", header=3)
    print_text("split the data", header=4)
    X_train, y_train = modeling_data.iloc[train_index, 1:].to_numpy(), modeling_data.iloc[train_index, 0].to_numpy()
    X_test, y_test = modeling_data.iloc[test_index, 1:].to_numpy(), modeling_data.iloc[test_index, 0].to_numpy()


    print_text("check correct shape training arrays", header=4)
    print((X_train.shape[0] == len(train_index)) & (X_train.shape[1] == modeling_data.shape[1]-1))
    print((y_train.shape[0] == len(train_index)) & (len(y_train.shape) == 1))
        #the shape of the response is (XXXX,), so we only have 1 number (row number) and not column number.

    print_text("check correct shape test arrays", header=4)
    print((X_test.shape[0] == len(test_index)) & (X_test.shape[1] == modeling_data.shape[1]-1))
    print((y_test.shape[0] == len(test_index)) & (len(y_test.shape) == 1))


    print_text("set the inner CV", header=4)
    cv_inner = KFold( \
        n_splits=3,  \
        shuffle=True)
        #random_state is not required as we have already ensured reproducibility by setting the seeds of python, numpy and tensorflow
        #I have run two times elastic net across all folds and I get exactly the same R2 for the best model. This indicates that both inner and outer CV schemas are the same across runs.
        #I have also run several times the training of random forest with this CV inner schema for the same outer folds, and I get the same results always.
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
    
    print_text("set the number of jobs", header=4)
    n_jobs=10
    pre_dispatch_value="1*n_jobs" 
        #pre_dispatch_value="1*n_jobs" to avoid memory explosion, see below function help
        #When running optuna on Yoruba for first time I had to reduce the number of jobs
            #Using more jobs we get "The exit codes of the workers are {SIGABRT(-6)}"
            #A lot of people is having the same problem even RAM seems to be ok. Some solve it by decreasing n_jobs and others by increase RAM usage.    
            #"Turns out allocating all CPUs can be unstable, specially when there are other independent programs running that can suddenly have an uncontrolled spike in memory usage."
            #it seems a problem of memory usage with scikit
                #https://github.com/scikit-learn-contrib/skope-rules/issues/18
            #In optuna-Yoruba, I detected that the decreasing the number of cores reduces the usage of memory even having the same number_jobs and optuna processes. Using for example 20 cores, I can run 10 optuna processes with number_jobs=2 using just 200GB per core. If I increase the number of cores to 40, things seems to seed up, but I get the workers error. So we are keeping things conservative with just 20 cores.


    print_text("define the GridSearch", header=4)
    search = GridSearchCV( \
        estimator=pipeline, \
        param_grid=space, \
        scoring="r2", \
        n_jobs=n_jobs, \
        cv=cv_inner, \
        verbose=0,
        refit=True, \
        pre_dispatch=pre_dispatch_value) #if n_jobs>1, then pre_dispatch should be "1*n_jobs", if not, the dataset will be copied many times increasing a lot memory usage
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
    print(search)


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


    print_text("create tuple for results", header=4)
    tuple_results = (split, model_class, best_params, score)
    print(tuple_results)

    print_text("convert the tuple of results to DF", header=3)
    tuple_results_df = pd.DataFrame([tuple_results], columns=["split", "model_class", "best_params", "score"])
        #the tuple has to be in a list in order to be transformed to DF
    print(tuple_results_df)

    print_text("save the table for the selected model and split", header=3)
    tuple_results_df.to_csv( \
        "./results/model_comparison/individual_results/selection_yoruba_hg19_" + model_class + "_" + str(split) + ".tsv", \
        sep='\t', \
        header=True, \
        index=False)
            #header=False because we are appending in a existing file
            #mode="a"
                #'a' for appending (which on some Unix systems, means that all writes append to the end of the file regardless of the current seek position).

    #return the tuple so it can be combined with the results of the other models
    return tuple_results

#run it across 2 splits for elastic net only
#[model_evaluation(*i) for i in indexes_cv_outer if (i[3]=="elastic_net") & (i[0] in [0, 1])]
    #"*" is used to unpack a tuple and use its elements as arguments
    #in our case, each tuple has as elements the split, indexes and name of the model class, which the arguments of model_evaluation
        #https://stackoverflow.com/a/1993732/12772630
#run it in 1 split for random forest
#[model_evaluation(*i) for i in indexes_cv_outer if (i[3]=="xgboost") & (i[0]==0)]

# endregion






##########################################
# region parallelize the function ########
##########################################

print_text("parallelize the function", header=2)
print_text("open pool with as many cores as splits*model class combinations we have", header=3)
import multiprocessing as mp
#parallelize across CV outer splits we have
pool = mp.Pool(cv_outer.n_splits)
print(pool)
    #40 splits*model classes combinations multiplied by 10 jobs/combination = 400 cores 


print_text("run the function using starmap, which is useful to apply function across iterable whose elements are in turn also iterables storing arguments (tuples in our case)", header=3)
#results = pool.starmap(model_evaluation, indexes_cv_outer)
indexes_cv_outer_subset = [i for i in indexes_cv_outer if i[3]==model_name]
results = pool.starmap(model_evaluation, indexes_cv_outer_subset)
    #map: Apply `func` to each element in `iterable`, collecting the results in a list that is returned
    #starmap: Like `map()` method but the elements of the `iterable` are expected to be iterables as well and will be unpacked as arguments. Hence `func` and (a, b) becomes func(a, b).
        #this is our case, as we have a list of tuples. Therefore, each element of the iterable (list) is in turn another iterable (tuple).
        #elements of the tuple are unpacked as argument which is exactly what we want.
        #https://stackoverflow.com/a/47506842/12772630
        #https://discuss.python.org/t/differences-between-pool-map-pool-apply-and-pool-apply-async/6575
print(results)

#close the pool
pool.close()


print_text("convert the list of results to DF", header=3)
results_df = pd.DataFrame(results, columns=["split", "model_class", "best_params", "score"])
print(results_df)



print_text("save the table", header=3)
import os
path_final_results = "./results/model_comparison/model_class_selection_yoruba_hg19.tsv"
if(not os.path.exists(path_final_results)):
    results_df.to_csv( \
        path_final_results, \
        sep='\t', \
        header=True, \
        index=False)
else:
    results_df.to_csv( \
        path_final_results, \
        mode="a", \
        sep='\t', \
        header=False, \
        index=False)
    #To append the DataFrame to an existing CSV file, you can use the mode='a' and header=False parameters in the to_csv method. We avoid the header because it is already there.


#THIS IS TOO SLOW DUE TO THE NEURAL NETWORKS, ALL MODELS FINISH EXCEPT NEURAL NETS. FROM THE MODELS FINSHED, XGBOOS IS THE BEST, AND ALSO HAS MUCH HIGHE R2 THAN WHAT I SAW WITH OPTUNA AND DEEP NETS BEFORE. IT IS AROUND 0.7!!
#MAYBE IT WOULD BE A GOOD IDEA TO DO THE BENCHMARK USING RMSE INSTEAD R2, you used RMSE for the fine tuning of the selected model class (i.e., XGBOOST; 03_explore_selected_model_class.py).


#WHEN YOU ARE DONE HERE, you can re-run the best model on the subset of non-overlapped 1000Kb windows used in the MDR revision to confirm non-independence of the genes is not affecting our results, although probably not worth it.
    #if you would wanted to move foward you could even use the simulation data of the MDR revision and BAT distance as a factor and model with xgboost to check that neutral simulations do not get the same results.

# endregion


##########
# FINISH #
##########
print_text("finish", header=1)
#chmod +x ./scripts/01_models_benchmark.py; ./scripts/01_models_benchmark.py --model_name="elastic_net" > ./scripts/01_models_benchmark.out 2>&1

