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




print_text("prepare split data", header=1)
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
        #Now, we only have 1 model class that we have already tuned and we just need to combine several instances of the model to get the final estimator. We need to know how good is this final estimator. As we get a larger test set, we should get a more robust measurement about the generalization ability of the model.





print_text("check performance of the final model", header=1)
print_text("define the final model", header=2)
import xgboost as xgb
from sklearn.pipeline import Pipeline

final_model = Pipeline( \
    steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', xgb.XGBRegressor(
                learning_rate=0.01,
                n_estimators=5000,
                max_depth=14,
                min_child_weight=23,
                gamma=0.01,
                subsample=1.0,
                colsample_bytree=0.9,
                colsample_bylevel=0.65,
                colsample_bynode=0.65,
                max_delta_step=14,
                reg_lambda=0.2,
                reg_alpha=1e-6,
                objective="reg:squarederror",
                nthread=10, 
                eval_metric="rmse", 
                seed=0))])


#combine models with different seeds for improving performance?
#Agree but partially. Some thoughts: 1. Though the standard deviations are high, as the mean comes down, their individual values should also come down (though theoretically not necessary). Actually the point is that some basic tuning helps but as we go deeper, the gains are just marginal. If you think practically, the gains might not be significant. But when you in a competition, these can have an impact because people are close and many times the difference between winning and loosing is 0.001 or even smaller. 2. As we tune our models, it becomes more robust. Even is the CV increases just marginally, the impact on test set may be higher. I've seen Kaggle master's taking AWS instances for hyper-parameter tuning to test out very small differences in values. 3. I actually look at both mean and std of CV. There are instances where the mean is almost the same but std is lower. You can prefer those models at times. 4. As I mentioned in the end, techniques like feature engineering and blending have a much greater impact than parameter tuning. For instance, I generally do some parameter tuning and then run 10 different models on same parameters but different seeds. Averaging their results generally gives a good boost to the performance of the model. Hope this helps. Please share your thoughts.


from sklearn.ensemble import VotingRegressor

#The fundamental difference between voting and stacking is how the final aggregation is done. In voting, user-specified weights are used to combine the classifiers whereas stacking performs this aggregation by using a blender/meta classifier
    #same wegiht ot all models?

#people do gridsearch on the weights of voting, so you can select the best importance for each model. In stack, you do CV to se how the meta-learne does selecting the best wegihts for each base learner


'''

final_model = Pipeline( \
    steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ('regressor', VotingRegressor(estimators=[("est_"+str(seed), xgb.XGBRegressor( learning_rate=0.01, n_estimators=826, max_depth=14, min_child_weight=23, gamma=0.01, subsample=1.0, colsample_bytree=0.9, colsample_bylevel=0.65, colsample_bynode=0.65, max_delta_step=14, reg_lambda=0.2, reg_alpha=1e-6, objective="reg:squarederror", nthread=10,  eval_metric="rmse",  seed=seed)) for seed in range(0,5,1)]))])
                #each model in the list has to be in the form ("estimator", instance)
                    #https://stackoverflow.com/questions/74461779/sklearn-votingclassifier-is-throwing-an-issue-about-argument-not-being-iterable

predictors = [x for x in train.columns if x not in ["prob(sweep)"]]
final_model.fit(train[predictors], train["prob(sweep)"])
y_pred=final_model.predict(test[predictors])
from sklearn.metrics import r2_score
score = r2_score(test["prob(sweep)"], y_pred)
print(score)
'''




from sklearn.ensemble import StackingRegressor
stack_pipeline = Pipeline( \
    steps=[ \
            ('scale', preprocessing.StandardScaler()), \
            ("regressor", StackingRegressor( \
                estimators=[("est_"+str(seed), xgb.XGBRegressor( learning_rate=0.01, n_estimators=5000, max_depth=14, min_child_weight=23, gamma=0.01, subsample=1.0, colsample_bytree=0.9, colsample_bylevel=0.65, colsample_bynode=0.65, max_delta_step=14, reg_lambda=0.2, reg_alpha=1e-6, objective="reg:squarederror", nthread=10,  eval_metric="rmse",  seed=seed)) for seed in range(0,5,1)], \
                final_estimator=Ridge(), \
                cv=cv_scheme))])

#I think this goes sequential, first base model, second base model....

from sklearn.linear_model import Ridge
    #Often simple ensembling models (e.g. simple average of predicted probabilities, weighted average of predicted probabilities or logistic regression of logits - possibly with regularization towards a simple average) perform a lot better than trying to fit fancier models on the second level
        #https://stats.stackexchange.com/questions/561584/why-is-my-stacking-meta-learning-not-outperforming-the-best-base-model
from sklearn.model_selection import KFold
cv_scheme = KFold( \
    n_splits=5,  \
    shuffle=True)
#reg = StackingRegressor(estimators=estimators, final_estimator=Ridge(), cv=cv_scheme)
stack_pipeline.fit(train[predictors], train["prob(sweep)"]).score(test[predictors], test["prob(sweep)"])
    #Note that estimators_ are fitted on the full X while final_estimator_ is trained using cross-validated predictions of the base estimators using cross_val_predict.

    #getting 0.67 with 5000! we have gained 1%!!
    

    #bigger test size implies more robust/stable measure of generalizibilty?
        #https://www.google.com/search?q=size+of+the+test+set&oq=size+of+the+test+set&gs_lcrp=EgZjaHJvbWUyBggAEEUYOTIMCAEQIRgPGBYYHRge0gEJMTM3NTlqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8

predictors = [x for x in train.columns if x not in ["prob(sweep)"]]
final_model.fit(train[predictors], train["prob(sweep)"])
y_pred=final_model.predict(test[predictors])
from sklearn.metrics import r2_score
score = r2_score(test["prob(sweep)"], y_pred)
print(score)







#####PIENSA SI USAR TODO EL SET O SOLO TRAINING, PORQUE SALE DIFERENTE Y PERDEMOS DATOS. SI YA TENEMOS VALDIADO EL MODELO EN EL HELD-OUT...

#Christopher says that "The partial dependence plot shows how the model output changes based on changes of the feature and does not rely on the generalization error. It does not matter whether the PDP is computed with training or test data.". I think this also applies for ALE plots as they do not use the performence, i.e., an error metric from the model obtained by comparing observed and predicted, but they show changes in the prediction as the feature changes in certain data points.
#In addition, I saw the follwing comment "Partial dependence plots can be performed over either the training or validation set, and examples of both cases can be found. For instance, fastbook uses the validation set, whereas Interpretable Machine Learning: A Guide for Making Black Box Models Explainable, in chapter 8.1 4, uses the training set. Frankly, in my experience, it ultimately does not matter which strategy you choose, and there are only a couple of important considerations. First, if the training set is large, PDP may take excessively long, in which case you can use a small chunk of it or resort to the validation set. Second, if the validation set is too small, PDP’s results would understandably be not very reliable, so the training set might be the wiser option."
#Therefore, we could apply this approach to data that has been already seen by the model, i.e., training or full dataset (training+test).
    #https://christophm.github.io/interpretable-ml-book/feature-importance.html#feature-importance-data


#https://datascience.stackexchange.com/questions/62913/shap-explanations-in-case-of-repeated-train-test-split
#https://datascience.stackexchange.com/questions/61395/shap-value-analysis-gives-different-feature-importance-on-train-and-test-set




#you can have slightly different results with XGBoost even setting the seeds
    #Changing subsample and colsample_bytree  to '1' and increasing early_stopping_rounds to '1000' (or whatever n_estimators is set to) should do the trick - let me know if this solves your problem or not. – 
        #https://stackoverflow.com/questions/61764057/how-to-get-reproducible-results-from-xgboostregressor-random-state-has-no-effec
    #I have checked that increasing the number of rounds in the final model makes things more stable.



print_text("ALE plots", header=1)

run_bash(" \
    cd ./results/selected_model_class; \
    mkdir \
        --parents \
        ./ale_plots")

#ale plots avoids problems with correlated features and we can use it with any model!!!
    #see jupyter notebooks for details
    #https://christophm.github.io/interpretable-ml-book/ale.html


#The ALE value can be interpreted as the main effect of the feature at a certain value compared to the average prediction of the data. This interpretation is made possible because ALE plots are centered at zero so each point of the ALE curve represents the difference with mean prediction.[4] For instance, an ALE value of 5 for a feature value 12.2 would mean that if the feature of interest equals 12.2, the prediction it would yield is higher by 5 than the average prediction ([link](https://towardsdatascience.com/explainable-ai-xai-methods-part-3-accumulated-local-effects-ale-cf6ba3387fde))

#alibi is much more maintained than aleplot and it is much easier to install in container
from alibi.explainers import ALE, plot_ale
    #warning numba
    #I think it is ok but check
import matplotlib.pyplot as plt
lr_ale = ALE(final_model.predict, feature_names=modeling_data[predictors].columns, target_names=["log(Flex-Sweep probability)"])

lr_exp = lr_ale.explain(modeling_data[predictors].to_numpy(), features=[16,17,18])
    #DECIDE NUMBER OF INTERVALS!!
        #ALE plots can become a bit shaky (many small ups and downs) with a high number of intervals. In this case, reducing the number of intervals makes the estimates more stable, but also smoothes out and hides some of the true complexity of the prediction model. There is no perfect solution for setting the number of intervals. If the number is too small, the ALE plots might not be very accurate. If the number is too high, the curve can become shaky.
plot_ale(lr_exp, n_cols=3, fig_kw={'figwidth':20, 'figheight': 15});
plt.savefig( \
    fname="./results/selected_model_class/ale_plots/aleplots_all_features_train.png")
plt.close()
    #The ALE on the y-axes of the plot above is in the units of the prediction variable.


#I understand that the bins are calculated as percentiles, the 50% is the percentile 50, i.e., the number that separate the data in two parts with the same size. 10% would be the first 10% of the data.
#Therefore, each interval has the same number of datapoints, so longer intervals means less data density
#You can see this for recombination. 50% is the median, while the last interval 90-100 is percentile 90 to 100. You can see this interval is much longer,from 2.6 (percentile 0.9) to 6.2 (percentile 1 or max value), because of this the bin is larger. The last 10% of the data is spread across a larger range. We have more sparse data there.


    #https://docs.seldon.io/projects/alibi/en/stable/examples/ale_regression_california.html

#This also is unable to detect the positive influence of GC content
#check prevous results with aleplot package (jupyter notebook)




#####tree shap

#tree shap also avoids problems with correlated features, but it is only for tree-based methods, this is not the case for the other SHAP methods

#we can use treeshap that does not include unlikely instances!!!
#Feature dependencies can skew the approximations made by KernelSHAP. The algorithm estimates SHAP values by randomly sampling feature values. The issue is that, when features are correlated, the sampled values can be unlikely. This means that when using SHAP values we can put too much weight on unlikely observations.
#TreeSHAP does not have the problem. However, the algorithm has a different issue due to feature dependencies. That is a feature that has no impact on the prediction can get a SHAP value that is non-zero. This can happen when the feature is correlated with another feature that does impact the prediction. In this case, we can make the incorrect conclusion that a feature has contributed to a prediction.
    #that is ok, if a feature is correlated with recombination, I expect that features to show a signal, but in the case of GC, the positive effect would remain?
    #https://christophm.github.io/interpretable-ml-book/shap.html
    #https://towardsdatascience.com/kernelshap-vs-treeshap-e00f3b3a27db


#check if shaps gives effect after controling for the rest of factors.
    #Note that GC content is also negatively associated with selection accoridng to ALE plots!!

import shap

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

    
##for the future
    #we could do like in the DL paper, they ran 100 independent XGBoost models with the same HPs but changing the random seed
        #After selecting XGBoost as the best-performing model class for the prediction of anti-AML drug synergy, we then wanted to account for the full diversity of possible good XGBoost models fit to the highly cor- related AML gene expression data. We therefore trained 100 models and explained the ensemble model. Each individual model had both row and column subsampling turned on for each additional tree fit, and the difference between the models in the ensemble was the random seed given to generate the subsampling.
            #you need to add a DIFFERENT SEED FOR EACH MODEL!!
    