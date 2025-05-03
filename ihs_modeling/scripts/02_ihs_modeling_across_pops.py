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



##########################################
######## IHS MODELING ACROSS POPS ########
##########################################

#In the previous script, we found the best algorithm and hyperparameters for modeling iHS in Yoruba, one of the population where selection should the most visible. We are using this approach to model in other populations.

#The best approach will be used to train the model in the 75% of the data of each pop and then predict in the remaining 25%. Permutation importance will be done in the test set while ALE plots will be done in the whole dataset.



##############################
# region INITIAL STEPS #######
##############################

###########
# imports #
###########

from scipy.stats import spearmanr
import random
import tensorflow as tf # type: ignore
from scikeras.wrappers import KerasRegressor #type: ignore
from tensorflow.keras.layers import Dropout # type: ignore
from tensorflow.keras.constraints import MaxNorm # type: ignore
from tensorflow.keras import regularizers # type: ignore
from scikeras.wrappers import KerasRegressor # type: ignore
from sklearn.pipeline import Pipeline
from sklearn import preprocessing
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from alibi.explainers import PermutationImportance # type: ignore
import pickle
import matplotlib.pyplot as plt
from alibi.explainers import plot_permutation_importance # type: ignore
from alibi.explainers import ALE # type: ignore
from alibi.explainers import plot_ale # type: ignore


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
os.environ["PYTHONHASHSEED"]=str(seed_value)
print(os.environ["PYTHONHASHSEED"])


print_text("Set the `python` built-in pseudo-random generator at a fixed value", header=2)
random.seed(seed_value)


print_text("Set the `numpy` pseudo-random generator at a fixed value", header=2)
np.random.seed(seed_value)


print_text("Set the `tensorflow` pseudo-random generator at a fixed value", header=2)
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



#######################################
# Passing arguments of python program #
#######################################

#define input arguments to be passed when running this script in bash
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--pop_name", type=str, default="CEUD", help="Selected population. String, None does not work!")
parser.add_argument("--n_iterations", type=int, default=1, help="The number of training-evaluation set to run. Int, None does not work!")
parser.add_argument("--energy_type", type=str, default="all_thermogenic", help="The variable related to energy metabolism to be included. String, None does not work!")
    #type=str to use the input as string
    #type=int converts to integer
    #default is the default value when the argument is not passed
args=parser.parse_args()
    #https://docs.python.org/3/library/argparse.html

#get the arguments of the function that have been passed through command line
pop_name = args.pop_name
n_iterations = args.n_iterations
energy_type = args.energy_type

# endregion






###################################################################
# region PREPARE THE BEST REGRESSOR WITH THE BEST PARAMTERS #######
###################################################################
print_text("prepare the best regressor with the best parameters", header=1)

#The best model we found in Yoruba was a DNN with 5 layers, 1000 units, tanh activation function, orthogonal initialization, dropout rate of 0.07, weight constraint of 4, L1 regularization of 3e-08 and L2 regularization of 9e-08. 

#We are going to use this model to predict iHS in other populations. My intuition is that, given positive selection (in the form of iHS signals) is more visible in Yoruba, the selection of models and hyperparameters could discover a good approach able to detect true associations between positive selection and genomic factors than then could be found in other populations with less visibility.

#We are going to train the selected approach in the 75% of the data of the new populations and then test in the 25% just to check the approach is still good.

print_text("initial preparations", header=2)
print_text("define a function to generate DNNs to be used in scikeras. This will be an input argument for KerasRegressor()", header=3)
#We are going to use scikeras, a wrapper for using keras on scipy ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). You can use it to run keras neural networks in pipelines and then use gridsearch to perform parameter optimization ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). As we will see below, we can optimize any parameter of the network. Not only the loss function or the learning rate, but also the number of layers and units, the dropout rate... just by using loops and ifelse within the function to create the model. Therefore, we can select the best combination of hyperparamenters to maximize predictive prower in the different evaluation datasets. Of course, as we can use pipelines, we can apply the scaling in response and predictors separately in each training-evaluation set.
#First, we are going to create a function to get the keras models. This function will be in turn used as input in scikeras.KerasRegressor().
    #https://coderzcolumn.com/tutorials/artificial-intelligence/scikeras-give-scikit-learn-like-api-to-your-keras-networks#1
    #https://www.adriangb.com/scikeras/stable/quickstart.html
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

print_text("define pipeline with scaling and regressor", header=3)
pipeline = Pipeline( \
    steps=[ \
    ('scale', preprocessing.StandardScaler()), \
    ('regressor', KerasRegressor(model=get_neural_reg, optimizer="Adamax", optimizer__learning_rate=0.0006, loss="mse", model__n_layers=5, model__n_units=1000, model__activation="tanh", model__init_mode="orthogonal", model__dropout_rate=0.07, model__weight_constraint=4, model__regu_L1=3e-08, model__regu_L2=9e-08, batch_size=100, epochs=630, verbose=0))])
    #Pipeline:
        #Sequentially apply a list of transforms and a final estimator. Intermediate steps of the pipeline must be 'transforms', that is, they must implement `fit` and `transform` methods. The final estimator only needs to implement `fit`.
    #The hyperparamters of the models are obtained from cross validation in the Yoruba population.
print(pipeline)

print_text("get a list of all populations for which we have average iHS across gene windows", header=3)
print_text("list of populations", header=4)
list_ihs_files = os.listdir("./data/mean_ihs_gene_windows/")
    #this is the folder where we have the mean iHS files for all populations. These were generated for the MDR paper
    #copied to our working folder
        #run_bash("cp --recursive ../../../method_deep_heavy_analyses/ihs_deep_learning/ihs_calculation/results/mean_ihs_gene_windows/ ./data/")

print_text("extract the POP name for each file and then get the unique names", header=4)
#i=list_ihs_files[0]
list_pops = np.unique([i.split("_")[0] for i in list_ihs_files])

# endregion






##########################################
# region TRAIN MODELS AND EVALUATE #######
##########################################
print_text("train models and evaluate", header=1)
print_text("prepare the data for the selected population", header=2)
print_text("load the IHS data for the selected population along with the predictors from the yourba dataset", header=3)
mean_ihs = pd.read_csv( \
    f"./data/mean_ihs_gene_windows/{pop_name}_mean_ihs_gene_windows_final_v1.txt.gz", \
    sep = "\t", \
    header = 0 \
)
    #response variable
n_ihs = pd.read_csv( \
    f"./data/mean_ihs_gene_windows/{pop_name}_n_ihs_gene_windows_final_v1.txt.gz", \
    sep = "\t", \
    header = 0 \
)
    #the number of iHS datapoints for each gene as predictor
yoruba_modeling_data = pd.read_csv( \
    "./data/YRID_modeling_dataset_v1.tsv.gz", \
    sep = "\t", \
    header = 0 \
)
    #just took this from the previous step, 00c_data_preparation.py. There I took predictors from elise´s flex-sweep dataset (she combined my predictors) and also BAT....

print_text("decide what energy predictors to be included", header=3)
if(energy_type=="thermogenic"):
    energy_predictors_to_remove = ["bat_distance_percentile_1", "smt_distance_percentile_1"]
elif(energy_type=="bat"):
    energy_predictors_to_remove = ["thermogenic_distance", "smt_distance_percentile_1"]
elif(energy_type=="smt"):
    energy_predictors_to_remove = ["thermogenic_distance", "bat_distance_percentile_1"]
elif(energy_type=="all_thermogenic"):
    energy_predictors_to_remove = []

print_text("merge the three DFs using 1000kb window data", header=3)
merged_ihs_data = pd.merge( \
    mean_ihs.loc[:, ["gene_id", "mean_ihs_1000kb"]], \
    n_ihs.loc[:, ["gene_id", "n_ihs_1000kb"]], \
    on="gene_id", \
    how="inner" \
)
modeling_data = pd.merge( \
    merged_ihs_data, \
    yoruba_modeling_data.loc[:,~yoruba_modeling_data.columns.isin(["mean_ihs_1000kb", "n_ihs_1000kb"]+energy_predictors_to_remove)], \
    on="gene_id", \
    how="inner" \
)
    #we only want genes that are present in the three datasets, so inner merge in all cases.

print_text("check the correlation between the thermogenic variables", header=3)
if(energy_type=="all_thermogenic"):
    
    print_text("calculate Spearman's rank correlation and overlapping", header=4)
    #thermo_pair=["thermogenic_distance", "bat_distance_percentile_1"]
    for thermo_pair in [["thermogenic_distance", "bat_distance_percentile_1"], ["thermogenic_distance", "smt_distance_percentile_1"], ["bat_distance_percentile_1", "smt_distance_percentile_1"]]:
        
        #correlation
        correlation, p_value = spearmanr(modeling_data[thermo_pair[0]], modeling_data[thermo_pair[1]])
        print(f"Spearman's correlation between {thermo_pair[0]} and {thermo_pair[1]}: Rho: {correlation} and P-value: {p_value}")

print_text("Apply log transformation to the target variable using the original DF as source", header=3)
modeling_data["mean_ihs_1000kb"] = modeling_data["mean_ihs_1000kb"].apply(lambda x: np.log(x))
    #It is should be ok to apply the log before splitting the dataset. There is a problem if you use a transformation that requires learn something from the rest of the data. For example, if you scale the whole dataset, you are using the mean and sd of the whole dataset, influencing data that will be used for test. In other words, there is room for a data leak. In this case, however, log(1.5) is always 0.4, independently of the rest of the data, so I think no data leak is possible. You could do a pipeline with log but it is a little bit more complicated (see [link](https://stats.stackexchange.com/questions/402470/how-can-i-use-scaling-and-log-transforming-together)), so we leave it for now.
        #Indeed I have found people in stack exchange saying this: However, yours (i.e. np.log1p) is a simple transformation that doesn't use any learnable parameters, and it won't matter if you do it before or after the split. It's like dividing a feature by 1000. 
            #https://stats.stackexchange.com/a/456056
    #From all these follow that if you use other transformations like preprocessing.PowerTransformer or QuantileTransformer ([link](https://yashowardhanshinde.medium.com/what-is-skewness-in-data-how-to-fix-skewed-data-in-python-a792e98c0fa6)), it is possible to have data leaks, so be careful.
    #In previous versions I was not using scaling or log transform for deep learning, because I assumed that the DNNs can deal with that, but maybe that was too much and in any case, we are going to use here also more simpler models that can be helped by scaling
    #Update: If I apply the log transformation within the pipeline I get a much lower R2 both in the training and test datasets! Not sure what is going on, but given this transformation does not summarize anything from the whole dataset, I can use it before splitting in training and evaluation. If this transformation was helping the training model to learn from the test set, the R2 in the test would be higher, but we have the opposite scenario.
    #in case you want to apply the log transformation within the pipeline, but this make the model MUCH WORSE compared to just apply the log to the original DF. Do not know why.
        #you can do it with func=np.log and inverse_func=np.exp in TransformedTargetRegressor


print_text("train, eval and visualize predictor's effect across iterations", header=2)
print_text("empty list to save results", header=3)
results = []

print_text("loop across iterations", header=3)
#iteration = [iteration for iteration in range(0, n_iterations)][0]
for iteration in range(0, n_iterations):

    print_text(f"starting iteration {iteration}", header=4)
    print_text("split the data into training and test sets", header=4)
    train_data, test_data = train_test_split(
        modeling_data.loc[:, ~modeling_data.columns.isin(["gene_id"])], 
        test_size=0.25,  #25% for testing
        random_state=iteration  # Ensures reproducibility
    )
        #We are only doing it one time, the hyperparameter selection was done in Yoruba data.

    print_text("get numpy arrays with the reponse and predictors ofr training and test", header=4)
    X_train = train_data.loc[:, ~train_data.columns.isin(["mean_ihs_1000kb"])].to_numpy()
    y_train = train_data.loc[:, train_data.columns.isin(["mean_ihs_1000kb"])].to_numpy()
    X_test = test_data.loc[:, ~test_data.columns.isin(["mean_ihs_1000kb"])].to_numpy()
    y_test = test_data.loc[:, test_data.columns.isin(["mean_ihs_1000kb"])].to_numpy()

    print_text("train the model", header=4)
    training_model = pipeline.fit(X_train, y_train)
    print(training_model)

    print_text("predict on the test set", header=4)
    y_pred = training_model.predict(X_test)

    print_text("calculate evaluation score in the test set", header=4)
    score = r2_score(y_test, y_pred)
    print(f"The R2 is {score}")

    print_text("save the results into the list", header=4)
    results.append({
        "pop_name": pop_name,
        "iteration": iteration,
        "r2_score": score
    })

    print_text("run permutation importance", header=3)
    print_text("define function for using rmse as loss function", header=3)
    #This is the same metric we have used for CV in the training set when modeling in Yoruba. R2 orders the features in the same way but strangely, we get a lot of variability (large error bar) for recombination rate. Therefore, we are using rmse, which does not have this problem and produce similar results for the rest. 
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
        "n_ihs_1000kb": "Number of iHS datapoints (1000kb)",
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
        "smt_distance_percentile_1": "Distance Skeletal Muscle Tissue (SMT) genes"
    }

    print_text("remove the non-selected energy predictors from the dict", header=4)
    #predictor_to_remove=energy_predictors_to_remove[0]
    for predictor_to_remove in energy_predictors_to_remove:
        del dict_nice_feature_names[predictor_to_remove]

    print_text("check that the keys in the dict are the features in the same order than in the data", header=4)
    predictors=[i for i in modeling_data.columns if i not in ["mean_ihs_1000kb", "gene_id"]]
    predictors_nice=[dict_nice_feature_names[i] for i in modeling_data.columns if i not in ["mean_ihs_1000kb", "gene_id"]]
    print(predictors_nice)
    if([i for i in dict_nice_feature_names.keys()] != predictors):
        raise ValueError("ERROR: the keys in the dict are not the same as the features in the data")

    print_text("Initialize the permutation feature importance", header=4)
    #Implementation of the permutation feature importance for tabular datasets. The method measure the importance of a feature as the relative increase/decrease in the loss/score function when the feature values are permuted. Supports black-box models.
        #https://docs.seldon.io/projects/alibi/en/stable/examples/permutation_importance_classification_leave.html
        #https://docs.seldon.io/projects/alibi/en/stable/methods/PermutationImportance.html
        #https://docs.seldon.io/projects/alibi/en/stable/api/alibi.explainers.html#alibi.explainers.PermutationImportance
    explainer_perm = PermutationImportance( \
        predictor=training_model.predict, \
        loss_fns={"rmse": loss_rmse}, \
        score_fns=None, \
        feature_names=predictors_nice, \
        verbose=True \
    )
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

    print_text("compute permutation importance", header=4)
    perm_exp = explainer_perm.explain( \
        X=X_test,  \
        y=y_test,  \
        method="estimate",  \
        kind="ratio", \
        n_repeats=50 \
    )
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

    print_text("save permutations with pickle", header=4)
    run_bash(f" \
        mkdir \
            -p \
            ./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/permutation_importance/ \
    ")
    pickle.dump( \
        perm_exp, \
        open( \
            f"./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/permutation_importance/{pop_name}_{iteration}_{energy_type}_hg19_dnn_permutation_importance_test_set.sav", \
            'wb' \
        ) \
    )
    #import pickle; perm_exp = pickle.load(open(f"./results/ihs_modeling_across_pops/{pop_name}/permutation_importance/{pop_name}_{iteration}_hg19_dnn_permutation_importance_test_set.sav", 'rb'))

    print_text("plot results with Recombination", header=3)
    ax=plot_permutation_importance( \
        exp=perm_exp, \
        features="all", \
        metric_names="all", \
        n_cols=1, \
        sort=True, \
        top_k=None, \
        ax=None, \
        bar_kw=None, \
        fig_kw={'figwidth': 14, 'figheight': 6} \
    )[0][0]
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
    max_x_values = max([i["mean"] for i in perm_exp["feature_importance"][0]])
    ax.update({"xlim": (1, max_x_values+0.2), "xlabel": "Permutation importance (permuted RMSE / original RMSE)"})
        #update the axis to change the xlim and focus on values above 1
        #https://www.geeksforgeeks.org/matplotlib-axes-axes-update-in-python/
    ax.set_title("", fontsize=18)
    ax.xaxis.label.set_size(13.5)
    ax.tick_params(axis='y', labelsize=11)
    plt.savefig( \
        fname=f"./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/permutation_importance/{pop_name}_{iteration}_{energy_type}_permutation_importance_all.png", dpi=300)
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

    #calculate ALE effects and plots only in the first iteration, as we are going to use the wholse dataset to calculate the ALE plots. See below for explanation about why it is ok to do this instead of using the test set
    if iteration == 0:

        print_text("ALE plots", header=3)
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

        print_text("define dict with nice axes names for ALE", header=4)
        dict_nice_feature_names_ale={
            "n_ihs_1000kb": "Number of iHS datapoints (1000kb)",
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
            "smt_distance_percentile_1": "Distance SMT genes (pb)"
        }
        
        print_text("remove the non-selected energy predictors from the dict", header=4)
        #predictor_to_remove=energy_predictors_to_remove[0]
        for predictor_to_remove in energy_predictors_to_remove:
            del dict_nice_feature_names_ale[predictor_to_remove]

        print_text("check that the keys in the dict are the features in the same order than in the data", header=4)
        predictors_ale=[i for i in modeling_data.columns if i not in ["gene_id", "mean_ihs_1000kb"]]
        predictors_nice_ale=[dict_nice_feature_names_ale[i] for i in modeling_data.columns if i not in ["gene_id", "mean_ihs_1000kb"]]
        print(predictors_nice_ale)
        if([i for i in dict_nice_feature_names_ale.keys()] != predictors_ale):
            raise ValueError("ERROR: the keys in the dict are not the same as the features in the data")

        print_text("train the model with the whole dataset", header=4)
        X_full_dataset = modeling_data.loc[:, ~modeling_data.columns.isin(["gene_id", "mean_ihs_1000kb"])].to_numpy()
        y_full_dataset = modeling_data.loc[:, modeling_data.columns.isin(["mean_ihs_1000kb"])].to_numpy()
        training_model_full_dataset = pipeline.fit(X_full_dataset, y_full_dataset)

        print_text("initialize ALE plots using alibi", header=4)
        #alibi is much more maintained than aleplot and it is much easier to install in container
        #we have already used this package for permutation importance
            #https://docs.seldon.io/projects/alibi/en/stable/index.html
        #we have a warning when importing alibi: "NumbaDeprecationWarning: The 'nopython' keyword argument was not supplied to the 'numba.jit' decorator."
            #It seems that this is ok and it is related to shap. This is going to be solved in new versions of shap and we are not using shap anyways.
                #https://github.com/slundberg/shap/issues/2909
        lr_ale = ALE( \
            predictor=training_model_full_dataset.predict, \
            feature_names=predictors_ale, \
            target_names=None, \
            check_feature_resolution=True, \
            low_resolution_threshold=10, \
            extrapolate_constant=True, \
            extrapolate_constant_perc=10.0, \
            extrapolate_constant_min=0.1 \
        )
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

        print_text("calculate ALE curves", header=4)
        lr_exp = lr_ale.explain( \
            X=modeling_data[predictors_ale].to_numpy(), \
            features=None, \
            min_bin_points=4, \
            grid_points=None \
        ) 
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

        print_text("save the explanations into pickle", header=4)
        run_bash(f" \
            mkdir \
                -p \
                ./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/aleplots/ \
        ")
        pickle.dump( \
            lr_exp, \
            open( \
                f"./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/aleplots/{pop_name}_{energy_type}_hg19_dnn_fit_to_full_dataset_ale_alibi_explanations.sav", \
                'wb' \
            ) \
        )
        #import pickle; lr_exp = pickle.load(open(f"./results/ihs_modeling_across_pops/{pop_name}/aleplots/{pop_name}_hg19_dnn_fit_to_full_dataset_ale_alibi_explanations.sav", 'rb'))
            #https://machinelearningmastery.com/save-load-machine-learning-models-python-scikit-learn/

        print_text("make ale plots", header=4)
        #Plotting ale_values against feature_values recovers the ALE curves. For convenience we include a plotting function plot_ale which automatically produces ALE plots using matplotlib:
        #Note that the ALE is estimated for each interval edge and linearly interpolated in between, for real applications it is important to have a sufficiently fine grid but also one that has enough points into each interval for accurate estimates. 
        #Important
            #The x-axis also shows feature deciles of the feature to help judge in which parts of the feature space the ALE plot is interpolating more and the estimate might be less trustworthy.
            #Therefore, first bar is percentile 10%, second bar is percentile 20%, and so on...
        #feature=[feature for (pos,feature) in enumerate(dict_nice_feature_names_ale.keys()) if pos==16][0]
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
                fname=f"./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/aleplots/{pop_name}_{energy_type}_aleplot_" + feature + ".png", dpi=300)
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

print_text("process the results and save", header=3)
print_text("convert to DF", header=4)
results_df = pd.DataFrame(results)

print_text("save the results", header=4)
results_df.to_csv( \
    f"./results/ihs_modeling_across_pops/{pop_name}/{energy_type}/{pop_name}_{energy_type}_model_eval.tsv", \
    sep="\t", \
    header=True, \
    index=False \
)

# endregion





print_text("FINISH", header=1)
#to run the script:
#cd /home/dftortosa/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating_heavy_analyses/dating_climate_adaptation/ihs_modeling/
#chmod +x ./scripts/02_ihs_modeling_across_pops.py
#singularity exec ./containers/03_explore_selected_model_class.sif ./scripts/02_ihs_modeling_across_pops.py --pop_name="FIND" --n_iterations=1 --energy_type="thermogenic" > ./02_ihs_modeling_across_pops_FIND_thermogenic.out 2>&1
    #we use the container of previous steps because it has alibi installed. I have been unable to install alibi in the container of this step.
#grep -Ei 'error|false|fail' ./02_ihs_modeling_across_pops_FIND_thermogenic.out
    #grep: The command used to search for patterns in files.
    #-E: Enables extended regular expressions.
    #-i: Makes the search case-insensitive.
    #'error|false|fail': The pattern to search for. The | character acts as an OR operator, so it matches any line containing "error", "false", or "fail".
    #03a_phenotype_prep.out: The file to search in.
