#!/usr/bin/env python3
# coding: utf-8

# # Modeling flex-sweep probability in Yoruba, massive optuna optimization

# Production script obtained from notebook that models the probability sweeps detected by Flex-sweep in Yoruba using multiple genomic factors. 
# 
# Note that this notebook can be converted to python script by using ```jupyter nbconvert --to script [YOUR_NOTEBOOK].ipynb```. See this [link](https://stackoverflow.com/questions/17077494/how-do-i-convert-a-ipython-notebook-into-a-python-file-via-commandline) for further details.


# Set the working directory:

# In[1]:

'''
import os

os.chdir("/media/dftortosa/Windows/Users/dftor/Documents/diego_docs/science/postdoc_enard_lab/projects/ancient_selection_dating/")
os.getcwd()
'''


# ## Set the window size

# Set the window size selected for those variables that are calculated within windows centered around genes.

# In[2]:


gene_window_size = "1000kb"


# ## Load the input data

# Load the input data and split response and predictors

# In[3]:


import pandas as pd

final_data_yoruba = pd.read_csv("data/flex_sweep_average_prob.txt.gz", 
                                compression="gzip")
final_data_yoruba


# In case you want to apply filters by recombination or disease variants

# In[4]:


#final_data_yoruba = final_data_yoruba.loc[
#    (final_data_yoruba["recombination_final_"+gene_window_size] >
#     np.quantile(final_data_yoruba["recombination_final_"+gene_window_size], 1/2)) & 
#    (final_data_yoruba["n_dis_"+gene_window_size]<
#     np.quantile(final_data_yoruba["n_dis_"+gene_window_size], 1/2)), :]


# ## Load input data

# Indicate the response

# In[5]:


response_name = "average_sweep_prob_" + gene_window_size
response_name


# Extract prob(sweep) and the all the predictors

# In[104]:


columns_to_exclude = ["gene_id", 
                     "predicted_class", 
                     "prob(sweep)",
                     "number_vips_"+gene_window_size,
                     "number_thermogenic_"+gene_window_size]


# In[105]:


modeling_data = final_data_yoruba[[column for column in final_data_yoruba.columns if column not in columns_to_exclude]]
modeling_data


# Put the response as first column:

# In[7]:


first_column = modeling_data.pop(response_name)
modeling_data.insert(0, response_name, first_column)
modeling_data
    #https://www.geeksforgeeks.org/how-to-move-a-column-to-first-position-in-pandas-dataframe/


# Convert to numpy array

# In[8]:


modeling_data_array = modeling_data.values
modeling_data_array


# Apply log transformation to the target variable.

# In[9]:


import numpy as np

modeling_data_array[:, 0] = np.log(modeling_data_array[:, 0])


# It is should be ok to apply the log before splitting the dataset. There is a problem if you use a transformation that requires learn someting from the rest of the data. For example, if you scale the whole dataset, you are using the mean and sd of the whole dataset, influencing data that will be used for test. In other words, there is room for a data leak. In this case, however, log(1.5) is always 0.4, independetly of the rest of the data, so I think no data leak is possible. You could do a pipeline with log but it is a little bit more compliated (see [link](https://stats.stackexchange.com/questions/402470/how-can-i-use-scaling-and-log-transforming-together)), so we leave it for now.
# 
# From all these follow that if you use other transformations like preprocessing.PowerTransformer or QuantileTransformer ([link](https://yashowardhanshinde.medium.com/what-is-skewness-in-data-how-to-fix-skewed-data-in-python-a792e98c0fa6)), it is possible to have data leaks, so be careful.

# ## Split the data

# Separate the response and predictors

# In[10]:


y = modeling_data_array[:, 0]
y.shape


# In[11]:


X = modeling_data_array[:, 1:]
X.shape


# In[12]:


f'Do we have the correct number of predictors? {X.shape[1] + 1 == modeling_data.shape[1]}'


# Held out part of the dataset. This part will not be used for parameter optimization, but for the final validation after the final model has been optimized. In this way, we avoid potential overfittin in the evaluation sets ([see link](https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation)). If you use train in a set of the data and then evaluate in other, you can see if the model trained is not overfitting to the training data and flexible enough to predict the evaluation data. The problem is that in parameter optimization, we select the best parameters based on the evaluation metrics in the evaluation sets, thus we could get a model that fit too much the evaluation dataset, loosing generalization and thus making the evaluation metrics no longer metrics of generalization. To avoid this, we leave out a set of the data for final evaluation. This set will be NOT used in parameter optimization. Once we have selected the best parameters to train a model in the training dataset and predict well in the evaluation dataset, we use these parameter to create a final model, fit to the whole training data and then predict in the final evaluation dataset, which was not used for anything before. If the model works well, it means it is generalizable and it is not overfitting the data, so, in our case, we can say that it is explaining the variance in selection that it is really explained by the genomic factors and the rest would be variance that could be explained by selective pressures. If there is overfitting, the model fit too much the data, there is not non-explained variance.

# In the future, we may want to use the final model to obtain a probability of selection considering genomic factors and then use it to select genes with the same expected probability of selection (according to these factors) than our genes of interest. The idea is that the interest genes should have the same probability of selection based on genomic features (predicted probability), but if they are target of a selective pressure, their observed probability of selection should be higher. If predicted and observed are exactly the same, there is no room for enrichment of selection in interest genes after controling for confounding factors, and this would be an methodological artifact due to a model that just fit the observed data without any generalization. This can be a problem with algorithms like RF or in deep learning.

# It is usually recommended a 80-20 ratio for training-test when the dataset is large, but for small datasets is better 70-30 ([link](https://www.researchgate.net/post/Is-there-an-ideal-ratio-between-a-training-set-and-validation-set-Which-trade-off-would-you-suggest)). 
# 
# Note that the test dataset gives information about the generalization ability of the model, and this is important to us as a measure of overfitting. We need to leave enough amount of data to check this but without leaving the training dataset so small that we underfit.
# 
# Because of this, we are going to use 70-30 even though we do not have a specially small dataset. We want to ensure our validation set covers well the variability found in the genome, it should be a difficult problem to the model so we can have more confidence we are not overfitting.

# In[13]:


from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, 
    y,
    test_size=0.30, #train size is automatically calculated using this value 
    random_state=54, 
    shuffle=True)
    #this function uses internally ShuffleSplit
    #so it can randomly split a dataset in training and evaluation
    #but instead of getting an interable (like in ShuffleSplit), 
    #you directly get the datasets splitted
        #https://stackoverflow.com/questions/66757902/differnce-between-train-test-split-and-stratifiedshufflesplit


# In[14]:


print(X_train.shape, y_train.shape)
print(X_test.shape, y_test.shape)


# ## Scaling

# We could apply the preprocessing to the initial dataset and then split into training and test but this is problematic. As previously explained, the scaling is done using the mean and sd of the sample, so if you do it in the whole dataset, the part will be used for test will be also influenced but training data. If you included your test data in the scaling, that means that your new data is treated differently from the training set, which defeats the purpose of the training set. In practice, this is unlikely to have a large impact, as computing mean and standard deviation is relatively stable on well-behaved datasets. However, I recommend to adhere to best practices, and split off the test set before doing any processing ([more info](https://amueller.github.io/aml/01-ml-workflow/03-preprocessing.html)).

# In other words, when you fit the standard scaler on the whole dataset, information from the test set is used to normalize the training set. This is a common case of "data leakage", which means that information from the test set is used while training the model. This often results in overestimates of the model's performance ([link](https://stackoverflow.com/questions/63037248/is-it-correct-to-use-a-single-standardscaler-before-splitting-data?noredirect=1&lq=1)).

# standard.scaler() is similar to preprocessing.scale(). In both cases, scaling means standardizing by removing the mean and scaling to unit variance. The standard score of a sample x is calculated as: z = (x - u) / s, where u is the mean of the training samples or zero if with_mean=False, and s is the standard deviation of the training samples or one if with_std=False.

# Instead of using preprocessing.scale, we will use preprocessing.StandardScaler(), which is a transformer. You can open an instance of this transformer, call the method fit to learn the mean and sd from the data and then apply the trasfomation

# In[15]:

'''
from sklearn import preprocessing

dummy_sample = np.array([1,1,2,2])

scaler = preprocessing.StandardScaler()
scaler.fit(dummy_sample.reshape(-1,1))
scaler.transform(dummy_sample.reshape(-1,1))
'''

# We get the exactly the same if we use preprocessing.scale

# In[16]:

'''
preprocessing.scale(dummy_sample)
'''

# The great adventage of StandardScaler() is that we can call an instance and use it within a pipeline that will be feed into a gridsearch. Internally, gridsearch scale the training dataset using the training mean-sd, fit the data and then predict in evaluation dataset but after scaling that evaluation dataset with the mean sd of the training set. In the next iteration (next training-eval set) the processes is repeated ([link](https://stackoverflow.com/questions/51459406/how-to-apply-standardscaler-in-pipeline-in-scikit-learn-sklearn)).


# ## Run multiple models

# Comparison of mulitple models recommended for regression problems with less than 100K observations ([scikit learn flowchart](https://scikit-learn.org/stable/tutorial/machine_learning_map/index.html)). In each case, fit the model, predict, see the prediction over the observed data and calculated R2 with cross validation.
# 
# For cross validation, we will use ShuffleSplit in order to create random splits for training and test sets. These splits will be used in cross_val_score to calculate the R2 in each test set.

# #### Fix random seeds for reproducibility

# In[17]:


from tensorflow.random import set_seed as tf_set_seed

np_seed = 7
tf_seed = 42
np.random.seed(np_seed)
tf_set_seed(tf_seed)


# #### Create create random splits for training and test sets.
# 
# We are going to make 10 splits. Therefore, we will obtain the average score of 10 CV splits, increasing the probability to catch models that are overfitting. If we exposed the model to model different datasets, there are more possibilities to detect if the model has problems to generalize. 1/10 of 14000 is 1400 which I think it is ok.
# 
# It is very important for us to avoid overfitting becuase we want to quantify the true impact of genomic factors on selection probability. If the model overfits, it will assume that genomic factors like recombination or GC-content explain much selection variability across the genome than they actually explain, reducing in this way the remaining variability that could be explained by other factors. 
# 
# This is specially relevant in the context of botstrapping, where we compare genes related to a selective pressures with control genes that have a similar probability of selection according to confounding factors. The goal is to test whether interest genes have more selection than predicted, but if selection is already fully explained by confounding factors, there is no room for interest genes to be enriched or depleted in selection. They have the exact selection the confounding factors mark.
# 
# In addition, it has been seen in many datasets that 10 folds works relatively well to estimate model performance ([link](https://machinelearningmastery.com/how-to-configure-k-fold-cross-validation/#:~:text=configure%20the%20procedure.-,Sensitivity%20Analysis%20for%20k,evaluate%20models%20is%20k%3D10.)). Although there is some debate about it. We will stick to 10 folds for now, if we have problems with the predictive power in the test set, we can re-evaluate.

# We will use KFold instead of ShuffleSplit because we want to have training and validation datasets that are not overlapped so we can use all the data. KFold with shuffle=True only shuflle one time at the beginning and then make the folds, while ShuffleSplit shuffles every time making possible to select as evaluation two times the same sample. In the case of 5 splits with KFold, one time 1/5 is the validation and the rest is for training, then another 1/5 and so on... until the 5 partitions have been used for validation, i.e., the whole dataset have been used for validation. 
# 
# It is important to us to use the whole training dataset for CV because we need to obtain models that are generalizable enough so we do not overestimate the influence of genomic confounding factors on selection due to overfitting (see above).
# 
# This is very well explained in this [threat](https://stackoverflow.com/questions/45969390/difference-between-stratifiedkfold-and-stratifiedshufflesplit-in-sklearn).

# In[18]:


from sklearn.model_selection import KFold

shuffle_split = KFold(n_splits=10, shuffle=True, random_state=23)

#set also the number of jobs as only 2
#we will use 100GB per core
number_jobs = 2
    #Using as many jobs as folds, we get
    #"The exit codes of the workers are {SIGABRT(-6)}"
    #A lot of people is having the same problem even RAM seems to be ok. Some solve it by decreasing n_jobs and others by increase RAM usage.    
    #"Turns out allocating all CPUs can be unstable, specially when there are other independent programs running that can suddenly have an uncontrolled spike in memory usage."
    #it seems a problem of memory usage with scikit
        #https://github.com/scikit-learn-contrib/skope-rules/issues/18


# #### Note about p-values

# Scikit learn does not give p-values, we have to study the impact of predictors using permutation plots. In case you need p-values, you can get them for linear models using statsmodel ([link](https://www.statsmodels.org/stable/index.html), [link](https://stackoverflow.com/questions/41045752/using-statsmodel-estimations-with-scikit-learn-cross-validation-is-it-possible)).
# 
# You could also run a bootstrap analysis with the final model adding noise to a given predictor and obtaining a distribution of r2, checking if the actual R2 is above of these cases (see stackoverflow posts about this).


## skipping ridge example

'''
# ### Prepare function to run each model

# #### Example with Ridge

# We indicate the model name using the exact name of the function and the required parameters as a string

# In[19]:


model_name = "Ridge"


# In[20]:


parameters = "alpha=0.5"


# Call the model with the corresponding parameters

# In[21]:


from sklearn.linear_model import Ridge

model = eval("".join([model_name, "(", parameters, ")"]))
model


# Create a pipe that includes the tranformer to do the scaling so it is fit (get mean and SD) to the training sets and not evaluation/validation sets. When using grid search ([link](https://stackoverflow.com/questions/51459406/how-to-apply-standardscaler-in-pipeline-in-scikit-learn-sklearn)) and cross_val_score ([link](https://stackoverflow.com/questions/44446501/how-to-standardize-data-with-sklearns-cross-val-score)), it uses the mean and SD of the training set of a given iteration and apply the transformation to the test set, then in the next iteration does the same but with the new training and test set.

# We need to first create a pipeline ([link](https://stackoverflow.com/questions/33091376/what-is-exactly-sklearn-pipeline-pipeline)) including scaling of predictors and the regressor. Then add it to the transformer of the response using TransformedTargetRegressor. We cannot just use StandardScaling on the response ([link](https://stackoverflow.com/questions/67824676/how-to-scale-both-x-and-y-data-in-sklearn-pipeline)).

# In[22]:


from sklearn.pipeline import Pipeline
from sklearn.compose import TransformedTargetRegressor

estim = TransformedTargetRegressor(regressor=Pipeline([('scale', 
                                                            preprocessing.StandardScaler()), 
                                                       ('regressor', model)]), 
                                       transformer=preprocessing.StandardScaler())
estim


# See the keys, you can see how you need to write two times regressor to reach alpha parameter of Ridge ("regressor__regressor__alpha"). Ridge is a regressor of the pipeline, but the pipeline is in turn a regressor of TransformedTargetRegressor.

# Fit and predict the model using the input data

# In[23]:


estim.fit(X_train, y_train)
prediction = estim.predict(X)


# See R2 in the whole dataset and in CV-subsets

# In[24]:


from sklearn.metrics import r2_score
from sklearn.model_selection import cross_val_score

print("R2 in the whole dataset:")
print(r2_score(y, prediction))

print("R2 across CV folds:")
print(np.mean(cross_val_score(estimator=estim, 
                        X=X_train, 
                        y=y_train, 
                        cv=shuffle_split, 
                        scoring="r2",
                        n_jobs=shuffle_split.n_splits)))


# Plot the observed sweep probability and the prediction in the whole dataset

# In[25]:


import matplotlib.pyplot as plt

plt.hist(y, bins=50, color="green", alpha=0.4, label="Observed sweep probability")
plt.hist(prediction, bins=50, color="blue", alpha=0.4, label="Prediction")
plt.annotate("R2 whole dataset: " + str(np.round(r2_score(y, prediction), 4)), 
             xy=(0.05, 0.7),
             xycoords='axes fraction')
plt.legend(loc='upper left')


# See the coefficients/feature importance

# Access to the regressor

# In[26]:


estim.regressor_.named_steps["regressor"].coef_


# In[27]:


#if the model is not VotingRegressor
if model_name not in ["VotingRegressor", "SVR", "KerasRegressor"]:

    #get the list of attributes of the regressor
    list_attributes = dir(estim.regressor_.named_steps["regressor"])    
    
    #select the attribute of interest (coefficient or feature importance)
    attribute_interest = [element for element in list_attributes if element in ["feature_importances_", "coef_"]][0]
    
    #get the attribute and the name of the predictor
    attribute_interest_pd = pd.DataFrame({
        "Predictors": modeling_data.columns[1:], 
        "Coefficient/feature importance": eval("".join(["estim.regressor_.named_steps['regressor'].", attribute_interest]))})
    print("Predictor importance in the training dataset")
    print(attribute_interest_pd)
'''


# #### Example with Deep Learning

# We are going to use scikeras, a wrapper for using keras on scipy ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). You can use it to run keras neural networks in pipelines and then use gridsearch to perform parameter optimization ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)). As we will see below, we can optimize any parameter of the network. Not only the loss function or the learning rate, but also the number of layers and units, the dropout rate... just by using loops and ifelse within the function to create the model. Therefore, we can select the best combination of hyperparamenters to maximize predictive prower in the different evaluation datasets. Of course, as we can use pipelines, we can apply the scaling in response and predictors separately in each training-evaluation set.

# First, we are going to create a function to get the keras models ([link](https://www.adriangb.com/scikeras/stable/quickstart.html)). This function will be in turn used as input in scikeras.KerasRegressor().
# 
# See this [link](https://coderzcolumn.com/tutorials/artificial-intelligence/scikeras-give-scikit-learn-like-api-to-your-keras-networks#1) for an example using scikeras without a get_model function.

# In[28]:


from tensorflow.keras.layers import Dropout
from tensorflow.keras.constraints import MaxNorm
from tensorflow.keras import regularizers

def get_reg(meta, n_layers, n_units, activation, init_mode, dropout_rate, weight_constraint, regu_L1, regu_L2):
    #note that meta is a dict with input metadata
    #i think this generated by KerasRegressor automatically once the input 
    #data is included
    
    #we import keras inside the function to avoid problems with parallelizinf
    from tensorflow import keras
        #https://stackoverflow.com/questions/42504669/keras-tensorflow-and-multiprocessing-in-python/42506478#42506478

    #start the sequential model
    model = keras.Sequential()
    
    #add the input layer
    model.add(keras.layers.Input(shape=(meta["n_features_in_"])))
        #input shape obtained from meta, a dict containing input metadata

    #for each of the layer sizes we have
    for i in range(n_layers):
        
        #add a Dense layer with the corresponding number of units
        #the selected activation
        #the method to add initial weight values
            #https://machinelearningmastery.com/weight-initialization-for-deep-learning-neural-networks/
        #dropout
        #regularization L1L2
        model.add(keras.layers.Dense(n_units, 
                                     activation=activation, 
                                     kernel_initializer=init_mode, 
                                     kernel_constraint=MaxNorm(weight_constraint), 
                                     kernel_regularizer=regularizers.L1L2(l1=regu_L1, 
                                                                           l2=regu_L2)))
        model.add(Dropout(dropout_rate))        
        
    #add the output layer with just one node (predicting 1 value)
    model.add(keras.layers.Dense(1))
        #No activation function is used for the output layer because 
        #it is a regression problem, and you are interested in predicting 
        #numerical values directly without transformation.
        #https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/
    
    #compile the model using the KerasRegressor arguments for
    #loss and optimization
    #model.compile(loss=compile_kwargs["loss"], optimizer=compile_kwargs["optimizer"])
        #to run this line you need to add to the function 
        #an argument called "compile_kwargs"
    return model


# According to the author of scikeras ([link](https://www.adriangb.com/scikeras/stable/notebooks/MLPClassifier_MLPRegressor.html#2.4-Losses-and-optimizer)), it is easier and safer to allow SciKeras to compile your model for you by passing the loss to KerasRegressor directly (KerasRegressor(loss="")) instead of compiling it inside get_reg. It is lessflexible because you can not use custom logic, e.g., if else to select a different learning rate according to the optimizer ([link](https://www.adriangb.com/scikeras/stable/advanced.html)), but he says it is safer, so we are going to select this approach. In case we need more flexibility, we can do the other approach.

# You can add a parameter for the dropout rate and add dropout layers after some layers within the loop (using ifelse). You can also stop the network at certain level of accuracy, i.e., callbacks ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#6.-Callbacks)).

# Now we use KerasRegressor. This takes as input a callable function (get_reg) that returns a keras model (model argument). KerasRegressor also takes some arugments needed to prepare a keras model like the activation function or the optimizer. We must also pass all of the arguments to get_reg as keyword arguments to KerasRegressor. Note that if you do not pass an argument to KerasRegressor, it will not be avilable for hyperparameter tuning ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#2.3-Defining-and-training-the-neural-net-classifier)).
# 
# Note about keyword arguments: Keyword arguments (or named arguments) are values that, when passed into a function, are identifiable by specific parameter names.

# We are going to perform a gridsearch in order to optimize the hyperparameters. SciKeras allows to direct access to all parameters passed to the wrapper constructors. This allows tunning of parameters like hidden_layer_sizes as well as optimizer__learning_rate.
# 
# The model__ prefix can be used to specify that a paramter is destined only for get_reg. In addition, to differentiate paramters like callbacks which are accepted by both tf.keras.Model.fit and tf.keras.Model.predict you can add a fit__ or predict__ routing suffix respectively ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)).

# In[29]:


from scikeras.wrappers import KerasRegressor
from tensorflow.keras import regularizers

dnn_regressor = KerasRegressor(
    model=get_reg,
    optimizer="Adam",
    optimizer__learning_rate=0.001,
    loss="mse",
    model__n_layers=2,
    model__n_units=10,
    model__activation="relu",
    model__init_mode="uniform",
    model__dropout_rate=0,
    model__weight_constraint=0,
    model__regu_L1=0,
    model__regu_L2=0,
    batch_size=200,
    epochs=50,
    verbose=0
)


# You can apply the regressor to the data and obtain predictions (user verbose=1 to see losses in each epoch) ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#2.3-Defining-and-training-the-neural-net-classifier)).

# In[30]:


dnn_regressor.fit(X_train,y_train)
dnn_regressor.predict(X_train)


# Once you have fitted the keras regressor, get_reg has been efectively called upon. So, you can check the network created inside the keras regressor using get_reg, which should have a hidden layer with 100 units (as indicated in KerasRegressor args) and a output layer 1 unit (as indicated in get_reg)

# In[31]:


dnn_regressor.model_.summary()


# You need to have run kerasregressor.fit in order to see the summary of the model. I guess this is caused because compilation is done outside of get_reg, so there is no compiled model until keras regressor does something.

# Now we can include the dnn_regressor in a pipeline just after the scaling of the predictors and then include everything in TransformedTargetRegressor to also transform the response, see above.

# In[32]:


from sklearn.pipeline import Pipeline
from sklearn.compose import TransformedTargetRegressor
from sklearn import preprocessing

full_dnn_regressor = TransformedTargetRegressor(regressor=Pipeline([('scale', preprocessing.StandardScaler()), 
                                                                    ('regressor', dnn_regressor)]),  
                                                transformer=preprocessing.StandardScaler())
full_dnn_regressor


# You can fit the whole pipeline to the data

# In[33]:


full_dnn_regressor.fit(X_train, y_train)


# We can also predict and obtain R2

# In[34]:


full_dnn_regressor.predict(X_train)


# In[35]:


from sklearn.metrics import r2_score

r2_score(y_train, full_dnn_regressor.predict(X_train))


# We can access the different steps of the pipeline and finally reach the Kerasregressor, which has the history saved after fitting.

# In[36]:


loss_pipe = full_dnn_regressor.regressor_.named_steps["regressor"].history_
loss_pipe


# We can see how there is only one loss, meaning that the loss is calculated in the data used as input (X_train), considering it as a whole, not partioning. The splitting in evaluation and training will be done later with gridsearch.

# In[37]:

'''
import matplotlib.pyplot as plt

plt.plot(loss_pipe["loss"])
'''


# Scikeras can work within cross_val_score, splitting the data in training and evaluation several times.

# In[38]:


from sklearn.model_selection import cross_val_score

print(np.mean(cross_val_score(estimator=full_dnn_regressor, 
                        X=X_train, 
                        y=y_train, 
                        cv=shuffle_split, 
                        scoring="r2",
                        n_jobs=number_jobs)))


##skip example gridsearc
'''

# We can also use it as input in gridsearch ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#7.-Usage-with-sklearn-GridSearchCV)).

# We first set the parameters to be optimized. We need to indicate the different steps to reach the keras regressor. First, we have a reponse transformer in which the pipe is the regressor, then the regressor in the pipe is the keras regressor. There, we have arguments that are part of keras regressor natively (e.g., loss) and also arguments we use in the get_reg function (e.g., size of the hidden layer), which is the input of the keras regressor and last level.

# In[39]:


params = {
    "regressor__regressor__model__n_layers": [2, 3],
    "regressor__regressor__model__n_units": [10, 20],
    "regressor__regressor__model__activation": ["relu", "tanh"],
    "regressor__regressor__loss": ["mse", "mae"]
} 


# Define the gridsearch

# In[40]:


from sklearn.model_selection import GridSearchCV

gs = GridSearchCV(estimator=full_dnn_regressor, 
                  param_grid=params, 
                  refit=True, #to fit a model with the best combination of parameters
                  cv=shuffle_split, 
                  scoring='r2', 
                  n_jobs=shuffle_split.n_splits, 
                  verbose=0)


# Perform the gridsearch

# In[41]:


gs.fit(X_train, y_train)


# There is a warning about memory leak when using more than 1 core, but the results are the same using 1 or several cores (n_jobs argument), so no problem.

# See the best score and params

# In[42]:


print(gs.best_score_, gs.best_params_)


# Get the estimator (keras regressor network) with the best parameters fitted to the whole data used as input

# In[43]:


gs.best_estimator_


# In[44]:


gs.best_estimator_.regressor_.named_steps["regressor"]


# See the loss of the best keras regressor

# In[45]:


plt.plot(gs.best_estimator_.regressor_.named_steps["regressor"].history_["loss"])


# In[46]:


params_best_dnn = gs.best_estimator_.regressor_.named_steps["regressor"].get_params()
    #get a dict with the parameters of a fitted model
    #https://stackoverflow.com/questions/54572468/python-sklearn-get-list-of-available-hyper-parameters-for-model


# See the architecture of the best keras regressor

# In[47]:


gs.best_estimator_.regressor_.named_steps["regressor"].model_.summary()


# Apply to the whole data

# In[48]:


prediction_deep = gs.best_estimator_.predict(X)


# In[49]:


print("R2 in the whole dataset:")
print(r2_score(y, prediction_deep))


# In[50]:


import matplotlib.pyplot as plt

plt.hist(y, bins=50, color="green", alpha=0.4, label="Observed sweep probability")
plt.hist(prediction_deep, bins=50, color="blue", alpha=0.4, label="Prediction")
plt.annotate("R2 whole dataset: " + str(np.round(r2_score(y, prediction_deep), 4)), 
             xy=(0.05, 0.7),
             xycoords='axes fraction')
plt.legend(loc='upper left')


# Save the model. You can use pickle or the keras method.
# 
# Pickle is not recommended if you want a high performance saving or share your data with others, becuase pickle is unsafe and tend to save everything ([link](https://www.adriangb.com/scikeras/stable/notebooks/Basic_Usage.html#4.-Saving-and-loading-a-model)). However, we not only have a kreas regressor, but also the whole pipeline, so we need to use pickle.

# In[51]:


import joblib

joblib.dump(gs.best_estimator_, 'results/dump_objects/initial_deep_learning_gridsearch.pkl')
    #https://stackoverflow.com/questions/34143829/sklearn-how-to-save-a-model-created-from-a-pipeline-and-gridsearchcv-using-jobli
#joblib.load("results/dump_objects/initial_deep_learning_gridsearch.pkl")
    #to load the pipeline you need all the required modules loaded
'''


#skip the function and plots of multiple models
'''


# #### Create the function

# In[52]:


import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.model_selection import cross_val_score

from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import VotingRegressor

def model_comparison(model_name, parameters):
    
    #open an instance of the selected model with the selected parameters
    model = eval("".join([model_name, "(", parameters, ")"]))
        
    #make estimator
    estim = TransformedTargetRegressor(regressor=Pipeline([('scale', 
                                                            preprocessing.StandardScaler()), 
                                                           ('regressor',
                                                            model)]), 
                                       transformer=preprocessing.StandardScaler())
    
    #fit in the training dataset
    estim.fit(X_train, y_train)
    
    #predict in the whole dataset
    prediction = estim.predict(X)    

    #plot
    plt.hist(y, bins=50, color="green", alpha=0.4, label="Observed sweep probability")
    plt.hist(prediction, bins=50, color="blue", alpha=0.4, label="Prediction")
    plt.annotate("R2 whole dataset: " + str(np.round(r2_score(y, prediction), 4)), 
                 xy=(0.05, 0.7),
                 xycoords='axes fraction')
    plt.legend(loc='upper left')
        
    #print R2 in test set
    print("R2 in the whole dataset:")
    print(r2_score(y, prediction))

    #print R2 in CV splits
    print("R2 across CV folds:")
    print(np.mean(cross_val_score(estimator=estim, 
                            X=X_train, 
                            y=y_train, 
                            cv=shuffle_split, 
                            scoring="r2",
                            n_jobs=shuffle_split.n_splits)))
    
    #if the model is not VotingRegressor
    if model_name not in ["VotingRegressor", "SVR", "KerasRegressor"]:
    
        #list attributes
        list_attributes = dir(estim.regressor_.named_steps["regressor"])
        
        #select the attribute of interest (coefficient or feature importance)
        attribute_interest = [element for element in list_attributes if element in ["feature_importances_", "coef_"]][0]

        #get the attribute and the name of the predictor
        attribute_interest_pd = pd.DataFrame({
            "Predictors": modeling_data.columns[1:], 
            "Coefficient/feature importance": eval("".join(["estim.regressor_.named_steps['regressor'].", attribute_interest]))})
        print("Predictor importance in the training dataset")
        print(attribute_interest_pd)


# ### Apply the function

# Apply the function to several models

# In[53]:


model_comparison(model_name="Lasso", parameters="")


# In[54]:


model_comparison(model_name="ElasticNet", parameters="")


# In[55]:


model_comparison(model_name="Ridge", parameters="")


# In[56]:


model_comparison(model_name="SVR", parameters="kernel='rbf'")


# In[57]:


get_ipython().run_cell_magic('time', '', 'model_comparison(model_name="RandomForestRegressor", parameters="random_state=1")\n')


# In[58]:


get_ipython().run_cell_magic('time', '', 'model_comparison(model_name="ExtraTreesRegressor", parameters="random_state=1")\n')


# In[59]:


model_comparison(model_name="GradientBoostingRegressor", parameters="random_state=1")


# In[60]:


model_comparison(model_name="VotingRegressor", 
                 parameters="estimators=[('gb', GradientBoostingRegressor()), ('et', ExtraTreesRegressor()), ('rdg', Ridge())]")


# In[61]:


initial_keras_params = """model=get_reg, 
    optimizer="adam",
    optimizer__learning_rate=0.001,
    loss="mean_absolute_error",
    model__n_layers=1,
    model__n_units=5,
    model__activation="relu",
    model__init_mode="uniform",
    model__dropout_rate=0,
    model__weight_constraint=0,
    model__regu_L1=0,
    model__regu_L2=0,
    batch_size=200,
    epochs=50,
    verbose=0
"""

model_comparison(model_name="KerasRegressor", parameters=initial_keras_params)


# Random forest and extra trees show the best R2, but with clear signs of overfitting. We are going to make an attemp with neural netoworks. Although they do not show the best performance with, we have to bearin mind that there is no default architecture, and we just used one inner layer or increase the epochs. If we do small changes like just adding more layers and units, the performance improves a lot. In addition, the structure of the data suggest it can work, as we have multiple predictors with complex relationships. 
# 
# If deep learning does not work, we can just go to RF, and try to reduce overffiting. RF partial plots show a nice decrease of selection as we go far away from thermogenic genes and vips.

'''


# ### Hyperparameter tunning

# #### Select hyperparameters to be optimized

# Set a list of parameters to be optimized. I am following this [tutorial](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/) about DL optimization. They perform the grid search for each parameter separately, and then agregate, but as they say, some parameters can interact (i.e., some specific combinations can be better), thus we will do the search considering all of them but within a bayesian framework, so we do not exhaustively look evey hyperparameter combination.

# ##### Batch vs. Epoch ([link](https://machinelearningmastery.com/difference-between-a-batch-and-an-epoch/)):
# 
# 
# - Stochastic Gradient descent: 
#     - The job of the algorithm is to find a set of internal model parameters (predictors) that perform well against some performance measure such as logarithmic loss or mean squared error. 
#     - Optimization is a type of searching process and you can think of this search as learning. The optimization algorithm is called “gradient descent“, where “gradient” refers to the calculation of an error gradient or slope of error and “descent” refers to the moving down along that slope towards some minimum level of error.
#     - The algorithm is iterative. This means that the search process occurs over multiple discrete steps, each step hopefully slightly improving the model parameters.
#     - Each step involves using the model with the current set of internal parameters to make predictions on some samples, comparing the predictions to the real expected outcomes, calculating the error, and using the error to update the internal model parameters.
#     - This update procedure is different for different algorithms, but in the case of artificial neural networks, the backpropagation update algorithm is used.
#    
# - More details about Stochastic gradient descent ([link](https://machinelearningmastery.com/understand-the-dynamics-of-learning-rate-on-deep-learning-neural-networks/))
#     - Stochastic gradient descent is an optimization algorithm that estimates the error gradient for the current state of the model using examples from the training dataset, then updates the weights of the model using the back-propagation of errors algorithm, referred to as simply backpropagation.
#     - The amount that the weights are updated during training is referred to as the step size or the “learning rate.”
#     - Specifically, the learning rate is a configurable hyperparameter used in the training of neural networks that has a small positive value, often in the range between 0.0 and 1.0.
#     - The learning rate controls how quickly the model is adapted to the problem. Smaller learning rates require more training epochs given the smaller changes made to the weights each update, whereas larger learning rates result in rapid changes and require fewer training epochs.
#     - A learning rate that is too large can cause the model to converge too quickly to a suboptimal solution, whereas a learning rate that is too small can cause the process to get stuck.
#     - The challenge of training deep learning neural networks involves carefully selecting the learning rate. It may be the most important hyperparameter for the model.
#    
# - Sample
#     - A sample is a single row of data.
#     - It contains inputs that are fed into the algorithm and an output that is used to compare to the prediction and calculate an error.
#     - A training dataset is comprised of many rows of data, e.g. many samples. A sample may also be called an instance, an observation, an input vector, or a feature vector.
# 
# - Batch
#     - The batch size is a hyperparameter that defines the number of samples (rows) to work through before updating the internal model parameters.
#     - Think of a batch as a for-loop iterating over one or more samples and making predictions. In each iteration, a portion of the data is considered and predictions are made considering the current combination of parameters. At the end of the batch, the predictions are compared to the expected output variables and an error is calculated. From this error, the update algorithm is used to improve the model, e.g. move down along the error gradient.
#     - A training dataset can be divided into one or more batches.
#     - When all training samples are used to create one batch, the learning algorithm is called batch gradient descent. When the batch is the size of one sample, the learning algorithm is called stochastic gradient descent. When the batch size is more than one sample and less than the size of the training dataset, the learning algorithm is called mini-batch gradient descent.
#         - Batch Gradient Descent. Batch Size = Size of Training Set
#         - Stochastic Gradient Descent. Batch Size = 1
#         - Mini-Batch Gradient Descent. 1 < Batch Size < Size of Training Set
#     - In the case of mini-batch gradient descent, popular batch sizes include 32, 64, and 128 samples. You may see these values used in models in the literature and in tutorials.
# 
# - Epochs
#     - The number of epochs is a hyperparameter that defines the number times that the learning algorithm will work through the entire training dataset.
#     - One epoch means that each sample in the training dataset has had an opportunity to update the internal model parameters. An epoch is comprised of one or more batches. For example, as above, an epoch that has one batch is called the batch gradient descent learning algorithm.
#     - You can think of a for-loop over the number of epochs where each loop proceeds over the whole training dataset. Within this for-loop is another nested for-loop that iterates over each batch of samples, where one batch has the specified “batch size” number of samples.
#     - The number of epochs is traditionally large, often hundreds or thousands, allowing the learning algorithm to run until the error from the model has been sufficiently minimized. You may see examples of the number of epochs in the literature and in tutorials set to 10, 100, 500, 1000, and larger.
#     - It is common to create line plots that show epochs along the x-axis as time and the error or skill of the model on the y-axis. These plots are sometimes called learning curves. These plots can help to diagnose whether the model has over learned, under learned, or is suitably fit to the training dataset.
#     
# - Batch vs. Epoch
#     - The batch size is a number of samples processed before the model is updated.
#     - The number of epochs is the number of complete passes through the whole training dataset.
#     - The size of a batch must be more than or equal to one and less than or equal to the number of samples in the training dataset.
#     - The number of epochs can be set to an integer value between one and infinity. You can run the algorithm for as long as you like and even stop it using other criteria besides a fixed number of epochs, such as a change (or lack of change) in model error over time.
#     - They are both integer values and they are both hyperparameters for the learning algorithm, e.g. parameters for the learning process, not internal model parameters found by the learning process.
#     - You must specify the batch size and number of epochs for a learning algorithm.
#     - There are no magic rules for how to configure these parameters. You must try different values and see what works best for your problem.
# 
# - Worked Example
#     - Finally, let’s make this concrete with a small example.
#     - Assume you have a dataset with 200 samples (rows of data) and you choose a batch size of 5 and 1,000 epochs.
#     - This means that the dataset will be divided into 40 batches (5*40=200), each with five samples. The model weights will be updated after each batch of five samples, i.e., in each learning step, only 5 random rows will be considered.
#     - This also means that one epoch will involve 40 batches or 40 updates to the model.
#     - With 1,000 epochs, the model will be exposed to or pass through the whole dataset 1,000 times. That is a total of 40,000 batches during the entire training process.
#     - If the sample size is not divisible by the batch size without remainder, the last batch would have a smaller number of samples ([link](https://machinelearningmastery.com/difference-between-a-batch-and-an-epoch/)). I guess this is ok because we have already fitted the model multiple times with the previous batches and we will repeate this several epochs. In addition, it is usual to tune the batch size testing different sizes, so people usually have this situation.
# 
# - Summary
#     - Stochastic gradient descent is an iterative learning algorithm that uses a training dataset to update a model.
#     - The batch size is a hyperparameter of gradient descent that controls the number of training samples to work through before the model’s internal parameters are updated.
#     - The number of epochs is a hyperparameter of gradient descent that controls the number of complete passes through the training dataset. 
#     - The idea behind is that most observations in a large data set will have many neighbors that impart almost the same information. We don't need all of them in our learning process. 
#     - BUT CAREFUL IN YOUR CASE, because you do not have millions of samples, although using multiple epochs you will go through the data many times using maany batches.
#     
# Epochs and batches (number of rows in each step) are not a replacement of cross-validation. It is an iterative process that works better in a subset of the data if you have a looot of data, because you will probably have repeated data. You can increase the batch size until the time of each step starts to increase. By default, tensorflow uses no replacement by default for batches. 

# In[62]:


batch_sizes = [10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220]
epoch_sizes = [10, 50, 100, 200, 400, 610]


# **Note**:
# 
# We are only using mini-batches, i.e., batch size>1 and < sample size. Using the whole sample size could make a lot of use of memory and I think batch size=1 could have problems to converge. I have used still a wider range than people use for minibatches. 
# 
# The epoch sizes usually go from 1 to 1000, I am not using the whole range to reduce computation time.
# 
# If we see the best models are around the limits, we can do a finer search around these values.

# ##### Optimization method 
# 
# There is a lot of literature and discussion about the selection of the optimization. Adam seems to be an approach that combine the strengths of other methods (i.e., Adagrad and RMSProp; [link](https://towardsdatascience.com/a-visual-explanation-of-gradient-descent-methods-momentum-adagrad-rmsprop-adam-f898b102325c); [link](https://datascience.stackexchange.com/questions/10523/guidelines-for-selecting-an-optimizer-for-training-neural-networks)). There is controversy about it with some results showing poor performance, but then other results show good performance (similar to Stochastic gradient descent + momentum) when doing a well parameter optimization ([link](https://www.fast.ai/posts/2018-07-02-adam-weight-decay.html), [link](http://ruder.io/optimizing-gradient-descent/index.html#adam)).
#  
# Advantages of Adam ([link](https://machinelearningmastery.com/adam-optimization-algorithm-for-deep-learning/)):
# - Straightforward to implement.
# - Computationally efficient.
# - Little memory requirements.
# - Invariant to diagonal rescale of the gradients.
# - Well suited for problems that are large in terms of data and/or parameters.
# - Appropriate for non-stationary objectives.
# - Appropriate for problems with very noisy/or sparse gradients.
# - Hyper-parameters have intuitive interpretation and typically require little tuning.
# 
# The authors describe Adam as combining the advantages of two other extensions of stochastic gradient descent. Specifically:
# 
# - Adaptive Gradient Algorithm (AdaGrad) that maintains a per-parameter learning rate that improves performance on problems with sparse gradients (e.g. natural language and computer vision problems).
#     - Our selective pressure variables has many zeros or values close to zero, being the most important values, so an approach like this that try not to underwegight this predictors could be useful.
# - Root Mean Square Propagation (RMSProp) that also maintains per-parameter learning rates that are adapted based on the average of recent magnitudes of the gradients for the weight (e.g. how quickly it is changing). This means the algorithm does well on online and non-stationary problems (e.g. noisy).
#     - We could consider our problem as noisy becasue we have multiple factors influencing selection.
#  
# In a review of optimizers ([link](https://arxiv.org/abs/1609.04747)), they say "Insofar, RMSprop, Adadelta, and Adam are very similar algorithms that do well in similar circumstances. […] its bias-correction helps Adam slightly outperform RMSprop towards the end of optimization as gradients become sparser. Insofar, Adam might be the best overall choice."
# 
# It seems this approach work fast and good for non-shallow problems, being usually used as first option currently. 
# 
# We are going to consider different optimizers just in case this makes any difference. 
# 
# **Note**: 
# 
# There are two optimizers (AdamW and Adafactor) that cannot be used as input strings in kerasregressor. We can use a class but that cannot be optimized in the tunning process. For doing that, I should compile the model inside get_reg() and add an argument for optimizer. Then include an Ifelse, if optimizer="AdamW", uses class AdamW, which has been previoulsy used. However, scikeras author recommend to compile the model outside, i.e., let kerasregressor to do it.
# 
# Given that these are only two optimizers and we are already including a battery of optimizers, including several Adam versions, we are going to skip this for now. 

# In[63]:


optimizers = ["SGD", 
              "RMSprop", 
              "Adagrad", 
              "Adadelta", 
              "Adam", 
              "Adamax", 
              "Nadam", 
              "Ftrl"]


# ##### Adam parameters ([link](https://machinelearningmastery.com/adam-optimization-algorithm-for-deep-learning/))
# 
# - alpha. Also referred to as the learning rate or step size. The proportion that weights are updated (e.g. 0.001). Larger values (e.g. 0.3) results in faster initial learning before the rate is updated. Smaller values (e.g. 1.0E-5) slow learning right down during training
# - beta1. The exponential decay rate for the first moment estimates (e.g. 0.9).
# - beta2. The exponential decay rate for the second-moment estimates (e.g. 0.999). This value should be set close to 1.0 on problems with a sparse gradient (e.g. NLP and computer vision problems).
# - epsilon. Is a very small number to prevent any division by zero in the implementation (e.g. 10E-8).
# 
# 
# The Adam paper suggests:
# - Good default settings for the tested machine learning problems are alpha=0.001, beta1=0.9, beta2=0.999 and epsilon=10−8
# - The TensorFlow documentation suggests some tuning of epsilon:
#     - The default value of 1e-8 for epsilon might not be a good default in general. For example, when training an Inception network on ImageNet a current good choice is 1.0 or 0.1.
# - Keras uses these values as default except epsilon, being 10-7.
# 
# They say that the default configuration parameters do well on most problems, so we are just going to tune the learning rate alpha. It seems that beta should not be change except to change specific reasons for your data ([link](https://stats.stackexchange.com/questions/499013/adam-adaptive-optimizers-learning-rate-tuning?rq=1)). Maybe in the future we can tune the betas usng bayesian optimization ([link](https://stats.stackexchange.com/questions/265400/deep-learning-how-does-beta-1-and-beta-2-in-the-adam-optimizer-affect-its-lear)). We could even not tune learning rate because it is adaptive in Adam (it changes along the process), but we are going to do it just in case, because it is an important parameter, and it could influence the initila learning rate before updating the parameters. In addition, we are using other optimizers in which learning rate could be more relevant, so it is important to tune this paramater.
# 
# More general info about tunning learning rate ([link](https://www.bdhammel.com/learning-rates/)).
# 
# Generally, it is a good idea to also include the number of epochs in an optimization like this as there is a dependency between the amount of learning per batch (learning rate), the number of updates per epoch (batch size), and the number of epochs ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)).
# 
# We do that as we are also tunning the number of batches and epochs.
# 
# **Note**: 
# 
# We are not optimizing some specific parameters of the optimizers, like epsilon, etc... In order to do that, we should compile the model inside get_reg() and add ifelse for the value of a new "optimizer" argument. For example, if optimizer="Adam", use Adam(epsilon_...) or something like that. So the function only take epsilon for Adam, but not sure if this would work if you are not using adam and epsilon is still a parameter in keras regressor. 
# 
# We are going to stick to default parameters (except learning rate) for optimizers for now. Maybe after the big search, having already an optimizer selected, we can optimize its own parameters.

# In[64]:


alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.3]


# I have seen that people using optuna (bayesian) sample learning rate values
# considering the algorithm scale, I think to prioritize small values. Our reference tutorial goes in the same line with grid searchs typically consisting in picking numbers between 10^−5 and 0.3 on a logaritmic scale ([link](https://machinelearningmastery.com/learning-rate-for-deep-learning-neural-networks/)). The wider range would be 1 - 10^−6: "Typical values for a neural network with standardized inputs (or inputs mapped to the (0,1) interval) are less than 1 and greater than 10^−6". 
# 
# Given we are going to use the narrower (but still wide) range of 10^−6 to 0.3. If we see that the best models are close to the limits, we can run a more detailed search around these values.

# ##### Weight initiallization ([link](https://machinelearningmastery.com/weight-initialization-for-deep-learning-neural-networks/))
# 
# The optimization algorithm requires a starting point in the space of possible weight values from which to begin the optimization process. Weight initialization is a procedure to set the weights of a neural network to small random values that define the starting point for the optimization (learning or training) of the neural network model.
# 
# " training deep models is a sufficiently difficult task that most algorithms are strongly affected by the choice of initialization. The initial point can determine whether the algorithm converges at all, with some initial points being so unstable that the algorithm encounters numerical difficulties and fails altogether"
# 
# Neural network weight initialization used to be simple: use small random values ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)). Now there is a suite of different techniques to choose from. Keras provides a laundry list.
# 
# Ideally, it may be good to use different weight initialization schemes according to the activation function used on each layer. In our case, however, we will use the same activation across layers because we are working with regression, not classification. Maybe if we use a initialization method that is not good for tahn, that combination will not work and will have low R2, but we will also run other DNNs with tahn and other initialization methods, and other with relu...

# In[65]:


init_mode = ["uniform",
             "lecun_uniform",
             "normal",
             "zero",
             "glorot_normal",
             "glorot_uniform",
             "he_normal",
             "he_uniform",
             "truncated_normal",
             "ones",
             "identity",
             "orthogonal",
             "constant",
             "variance_scaling"]


# ##### Activation function ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/))
# 
# The activation function controls the non-linearity of individual neurons and when to fire ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/)). 
# 
# We are going to check a big battery of functions

# In[66]:


activations = ['softmax', 
               'softplus', 
               'softsign', 
               'relu', 
               'LeakyReLU',
               'selu',
               'elu',
               'exponential',
               'tanh', 
               'sigmoid', 
               'hard_sigmoid', 
               'linear']


# ##### Dropout rate ([link](https://machinelearningmastery.com/dropout-for-regularizing-deep-neural-networks/))
# 
# We are going to add dropout so some neurons are shutdown. In this way, we will try to get more generalization. 
# 
# During training, some number of layer outputs are randomly ignored or “dropped out.” This has the effect of making the layer look-like and be treated-like a layer with a different number of nodes and connectivity to the prior layer. In effect, each update to a layer during training is performed with a different “view” of the configured layer.
# 
# Dropout has the effect of making the training process noisy, forcing nodes within a layer to probabilistically take on more or less responsibility for the inputs. 
# 
# The idea behind this is to keep individual neurons from becoming too specialized.  Because each neuron may or may not be in any given run, it cannot be the only neuron detecting a particular feature.  If that feature is important, responsibility for its detection needs to be spread out among several neurons.  When we make a prediction, we will keep all of the neurons, in order to make the best prediction possible.
# 
# 
# This conceptualization suggests that perhaps dropout breaks-up situations where network layers co-adapt to correct mistakes from prior layers, in turn making the model more robust.
# 
# Dropout may be implemented on any or all hidden layers in the network as well as the visible or input layer. It is not used on the output layer.
# 
# A common value is a probability of 0.5 for retaining the output of each node in a hidden layer and a value close to 1.0, such as 0.8, for retaining inputs from the visible layer.
# 
# When using dropout regularization, it is possible to use larger networks with less risk of overfitting. In fact, a large network (more nodes per layer) may be required as dropout will probabilistically reduce the capacity of the network.
# 
# A good rule of thumb is to divide the number of nodes in the layer before dropout by the proposed dropout rate and use that as the number of nodes in the new network that uses dropout. For example, a network with 100 nodes and a proposed dropout rate of 0.5 will require 200 nodes (100 / 0.5) when using dropout.
# 
# Dropotu can be applied to the input layer, so some samples are removed ([link](https://machinelearningmastery.com/dropout-regularization-deep-learning-models-keras/)), but I have not seens this before in my previous courses so I will pass.

# **weight constrain**
# 
# You also need weight constrain:
#     
# Network weights will increase in size in response to the probabilistic removal of layer activations.Large weight size can be a sign of an unstable network.
# 
# To counter this effect a weight constraint can be imposed to force the norm (magnitude) of all weights in a layer to be below a specified value. For example, the maximum norm constraint is recommended with a value between 3-4.
# 
# Constraining the weight matrix directly is another kind of regularization. If you use a simple L2 regularization term (see below) you penalize high weights with your loss function. With this constraint, you regularize directly ([link](https://www.kdnuggets.com/2015/04/preventing-overfitting-neural-networks.html/2), [link](https://stackoverflow.com/questions/45970888/what-does-kernel-constraint-max-norm3-do)). Therefore, this is another layer of regularization that is considered. Our search could combine it with the rest of thechiques or not, depending on the combination.

# We will try dropout percentages between 0.0 and 0.9 (1.0 does not make sense) and maxnorm weight constraint values between 0 and 5.

# In[67]:


dropout_rates = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
weight_constraints = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]


# **Note**
# 
# We are going to apply the same dropout rate for all inner layers. Once we have an optimized architecture, we can try to improve it by setting the dropout only in specific layers.

# ##### Early Stop
# 
# The idea behind early stop is to stop the fitting process if a given metric does not improve after several steps. 
# 
# Ideally, you would monitor the loss/accuracy of the validation dataset, so you can if the validation dataset is not improving and hence we would have a high risk of overfitting (i.e., getting better predictions in training than in validation). 
# 
# This can be implemented in KerasRegressor with the argument callbacks, for example calling the class EarlyStopping, and then add routed parameters like callbacks__patience (how many steps we have to wait to stop if there is no improvement).
# 
# There is, however, a **problem** if we use this in combination with cross_val_score, grid search.... There is no validaition set identified in KerasRegressor (the validation sets are created are created by grid search), so we can only monitor the loss/accuracy in the training dataset.
# 
# There are also criticisms to combine CV and early stopping ([link](https://stackoverflow.com/questions/48127550/early-stopping-with-keras-and-sklearn-gridsearchcv-cross-validation)), because you could loss your best model because stopping early. I am not sure about that, but there is an more important point:
# 
# You are already comparing models with different number of epocs, so if there is overfitting after some epocs, the model with more epocs will get a lower R2 in the validation dataset. 
# 
# In summary, we are not going to use early stopping to avoid overfitting. We stick to the selection of the best number of epochs according the validation set.

# #### Regularization
# 
# Overfitting is caused by the model having too much flexibility, and the main source of flexibility in a neural network is the weights in the neurons. Specifically, non-zero weights indicate some relationship between the input and the output. Regularization penalizes this flexibility by penalizing non-zero weights. Thus, the network will only have a weight be non-zero if the benefit to the loss function is greater than the penalty applied for the weight.
# 
# There are two main types of regularization:  𝐿2 -regularization adds a penalty proportional to the sum of the squares of the weights, while  𝐿1 -regularization uses the sum of the absolute values of the weights. (The biases are generally not regularized.). Therefore, making less likely to have non-positive weights and hence associations between predictors and target. In other words, it penalizes “big” weights.
# 
# The hyperparameter  𝛼  is the regularization parameter. Its size needs to be set to provide the right amount of flexibility that the net avoids both overfitting and its converse, underfitting. In general, you will need to do a bit of a search to find the appropriate value for your problem. Info from Optizimation notebook of TDI.
# 
# If the hyperparameter is zero, then the L1/2 sum goes to zero, and there is no regulization ([link](https://towardsdatascience.com/intuitions-on-l1-and-l2-regularisation-235f2db4c261)).
# 
# A linear regression model that implements L1 norm for regularisation is called lasso regression, and one that implements (squared) L2 norm for regularisation is called ridge regression ([link](https://towardsdatascience.com/intuitions-on-l1-and-l2-regularisation-235f2db4c261)). We are going to use L1L2 function so you can have both or one of the two most frequent types of regularization setting the corresponding parameter l1 and l2 for each regularization type ([link](https://keras.io/api/layers/regularizers/)). The default is 0.01 in both cases. 

# Three different regularizer instances are provided; they are ([link](https://machinelearningmastery.com/how-to-reduce-overfitting-in-deep-learning-with-weight-regularization/)):
# 
# - L1: Sum of the absolute weights.
# - L2: Sum of the squared weights.
# - L1L2: Sum of the absolute and the squared weights.
# 
# This is achieved by setting the kernel_regularizer argument on each layer. A separate regularizer can also be used for the bias via the bias_regularizer argument, although this is less often used.
# 
# The most common type of regularization is L2, also called simply “weight decay,” with values often on a **logarithmic scale** between 0 and 0.1, such as 0.1, 0.001, 0.0001, etc. It is a good practice to first grid search through some orders of magnitude between 0.0 and 0.1, then once a level is found, to grid search on that level.
# 
# Note that in the paper introducing dropout, this techinque was combined with L2 regularization ([link](https://stats.stackexchange.com/questions/241001/deep-learning-use-l2-and-dropout-regularization-simultaneously)). In any case, we are going to explore the possibilty of having both or one of them by tunning their corresponding hyperparameters.

# In[68]:


regu_L1_values=[1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7] #use log
regu_L2_values=[1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7] #use log


# **Note**
# 
# In the second search, we can make a finer search around the values of the best models.

# #### Batch normalization
# 
# Another [recently-developed](https://arxiv.org/pdf/1502.03167v3.pdf) tool for deep networks is **batch normalization**.  Although it can help with overfitting, it was originally developed to deal with the vanishing gradient problem.  Recall that activation functions have flat regions, where their gradients are small.  When the input is in these regions, gradient descent will only move the weights small amounts, leaving them stuck in the low-gradient regions.  Intelligent choices for initializations and activation functions try to avoid this as much as possible.
# 
# Batch normalization takes a more proactive approach, scaling and shifting the inputs so that the average input, over the whole batch, has a target mean and standard deviation.  These target values become parameters of the model, tuned during training.
# 
# By keeping gradients from vanishing, batch normalization reduces the importance of the weight initialization and the activation function.  Larger learning rates can be used.  Although the initial steps may proceed more slowly, as the correct normalizations must be learned, learning should proceed much faster overall than without.  Batch normalization can also have a regularization effect, reducing the propensity towards overfitting!
# 
# From TDI optimization notebook. 
# 
# **Don’t Use With Dropout** ([link](https://machinelearningmastery.com/batch-normalization-for-training-of-deep-neural-networks/), [link](https://stackoverflow.com/a/62806906/12772630))
# 
# Batch normalization offers some regularization effect, reducing generalization error, perhaps no longer requiring the use of dropout for regularization.
# 
# Removing Dropout from Modified BN-Inception speeds up training, without increasing overfitting.
# 
# — Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift, 2015.
# 
# Further, it may not be a good idea to use batch normalization and dropout in the same network.
# 
# The reason is that the statistics used to normalize the activations of the prior layer may become noisy given the random dropping out of nodes during the dropout procedure.
# 
# Batch normalization also sometimes reduces generalization error and allows dropout to be omitted, due to the noise in the estimate of the statistics used to normalize each variable.
# 
# — Page 425, Deep Learning, 2016.
# 
# **Decision:** We are not going to use batch normalization for now.

# ##### Number of inner layers and neurons ([link](https://machinelearningmastery.com/grid-search-hyperparameters-deep-learning-models-python-keras/))
# 
# The number of neurons in a layer is an important parameter to tune. Generally the number of neurons in a layer controls the representational capacity of the network, at least at that point in the topology.
# 
# Also, generally, a large enough single layer network can approximate any other neural network, at least in theory.
# 
# A larger network requires more training and at least the batch size and number of epochs should ideally be optimized with the number of neurons (we are doing that).
# 
# It is recommended in general to start coarse and then fine tune the model, so we are going to compare very different number of layers and neurons.

# In[69]:


n_layers = np.arange(1, 21, 3)
n_units = np.arange(1, 500, 20)


# **Note**
# 
# We are not looking to networks deeper than 20 layers because computational time. If we see that the number of layers or neurons of the best models is close to the limit, we can extedn this in the finer search.
# 
# Similarly, in the second search we can try to change a bit the number of units between layers.

# #### Losses ([link](https://machinelearningmastery.com/how-to-choose-loss-functions-when-training-deep-learning-neural-networks/))
# 
# Inside the network, the difference between observed and predicted will be calculated using the loss function, and this information will be used to make the next step while considering the learning rate, etc... 
# 
# These metrics can be used also as accuracy metrics in the CV, but for the cross-validation, we will use r2 (see below).
# 
# - The Mean Squared Error, or MSE,
#     - loss is the default loss to use for regression problems. 
#     - Mathematically, it is the preferred loss function under the inference framework of maximum likelihood if the distribution of the target variable is Gaussian. It is the loss function to be evaluated first and only changed if you have a good reason.
#     - Mean squared error is calculated as the average of the squared differences between the predicted and actual values. The result is always positive regardless of the sign of the predicted and actual values and a perfect value is 0.0. The squaring means that larger mistakes result in more error than smaller mistakes, meaning that the model is punished for making larger mistakes.
# 
# - Mean Squared Logarithmic Error Loss
#     - There may be regression problems in which the target value has a spread of values and when predicting a large value, you may not want to punish a model as heavily as mean squared error.
#     - Instead, you can first calculate the natural logarithm of each of the predicted values, then calculate the mean squared error. This is called the Mean Squared Logarithmic Error loss, or MSLE for short.
#     - It has the effect of relaxing the punishing effect of large differences in large predicted values.
#     - As a loss measure, it may be more appropriate when the model is predicting unscaled quantities directly. Nevertheless, we can demonstrate this loss function using our simple regression problem.
#     - The model can be updated to use the ‘mean_squared_logarithmic_error‘ loss function and keep the same configuration for the output layer. We will also track the mean squared error as a metric when fitting the model so that we can use it as a measure of performance and plot the learning curve.
# 
# - Mean Absolute Error
#     - On some regression problems, the distribution of the target variable may be mostly Gaussian, but may have outliers, e.g. large or small values far from the mean value.
#     - The Mean Absolute Error, or MAE, loss is an appropriate loss function in this case as it is more robust to outliers. It is calculated as the average of the absolute difference between the actual and predicted values.
#     - The model can be updated to use the ‘mean_absolute_error‘ loss function and keep the same configuration for the output layer.
# 
# There are more losses for regression like Mean Absolute Percentage Error (Mape) or Mean Squared Logarithmic Error (square(log(y_true + 1.) - log(y_pred + 1.))). I have selected all that are not dedicated to categorical responses ([link](https://www.tensorflow.org/api_docs/python/tf/keras/losses), [link](https://towardsdatascience.com/understanding-loss-functions-the-smart-way-904266e9393)).

# In[70]:


losses = ["mse", 
          "msle", 
          "mae", 
          "mape",
          "mean_absolute_error",
          "cosine_similarity", 
          "huber", #In statistics, the Huber loss is a loss function 
                  #used in robust regression, that is less sensitive to
                  #outliers in data than the squared error loss.
                  #it seems to get the strengths of MAE and MSE without 
                  #their weaknesess
                  #https://towardsdatascience.com/understanding-loss-functions-the-smart-way-904266e9393
          "log_cosh"]
            #log(cosh(x)) is approximately equal to (x ** 2) / 2 for small x and 
            #to abs(x) - log(2) for large x. This means that 'logcosh' 
            #works mostly like the mean squared error, but will not be so 
            #strongly affected by the occasional wildly incorrect prediction.


# #### Scorer for tuning
# 
# We are going to select ONLY ONE scorer in order to perform the hyperparameter tuning.
# 
# In regression problems you can not use the same accuracy metrics as in classification problems (e.g. error rate, confusion matrix, etc.): in stead, other metrics are used like:
# 
# - Pearson linear correlation
# - Spearman rank correlation
# - RMSE (root mean squared error)
# - MAE (mean absolute error)
# - etc. (there are many more)
# 
# The R2 is the proportion of the variation in the dependent variable that is predictable from the independent variable(s) ([link](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html#sklearn.metrics.r2_score), [link](https://en.wikipedia.org/wiki/Coefficient_of_determination)). 
# 
# I have used this metric multiple times including in this very project and I saw a great correlation between R2 and the fit to the distribution. For example, in cases of extreme overfitting, i.e., the distribution of predicted is exactly the same than observed, the R2 calculated with scikitlearn is 1. It is also much easier to interpet and widely use both in machine learning and science in general.
# 
# We are going to use this metric for now.

# #### Bayesian optimization ([link](https://towardsdatascience.com/hyperparameter-optimization-with-scikit-learn-scikit-opt-and-keras-f13367f3e796))
# 
# Grid search is ok if you are looking a small parameter space, because you just can test all the possible combinations of the few hyperparameters you are using. If you need a little more parameter space you can use random grid search, which randomly select combinations of parameters (not all of them) and can reach a good combination after some time exploring, altough it is unlikely to get the best. Using randomized search is not too hard, and it works well for many fairly simple problems.
# 
# When training is slow, however, (e.g., for more complex problems with larger datasets), this approach will only explore a tiny portion of the hyperparameter space. You can partially alleviate this problem by assisting the search process manually: first, run a quick random search using wide ranges of hyperparameter values (corse search), then run another search using smaller ranges of values centered on the best ones found during the first run, and so on. This approach will hopefully zoom in on a good set of hyperparameters. However, it’s very time consuming, and probably not the best use of your time.
# 
# Fortunately, there are many techniques to explore a search space much more efficiently than randomly. Their core idea is simple: when a region of the space turns out to be good, it should be explored more. Such techniques take care of the “zooming” process for you and lead to much better solutions in much less time. Just like we wanted, first coarse and then zoom, but automatically,not manually.
# 
# One such technique is called Bayesian Optimization. In this [link](https://towardsdatascience.com/grid-search-vs-random-search-vs-bayesian-optimization-2e68f57c3c46#:~:text=Unlike%20the%20grid%20search%20and,is%20determined%20by%20the%20user.), they compare gridsearch, random gridsearch and bayesian sarch (with Optuna), showing that bayesian reach the same solution that gridsearch but in much less time. Random is the fastest, but it did not find the solution. 
# 
# In our case, we will use Optuna ([link](https://optuna.org/), [link](https://towardsdatascience.com/state-of-the-art-machine-learning-hyperparameter-optimization-with-optuna-a315d8564de1)), but there are alternatives:
# 
# - scikit-optimize ([link](https://scikit-optimize.github.io/stable/)).
# - Keras tuner ([link](https://keras.io/keras_tuner/)) from keras.
# - Genetic algorithms (considering each hyperparameter combination as an individual that have fitness, e.g., R2, and will be selected or not based on that; [link](https://towardsdatascience.com/hyperparameters-tuning-from-grid-search-to-optimization-a09853e4e9b8))
# - BayesianOptimization
# 
# I was initially using scikit-optimize, but it gave errors for some hyperparameter combinations and does not show the progress, i.e., the score of the current combination and what is the currelty best.

# #### Bayesian optimization with Optuna

# We are going to hide Info messages fromtensorflow, but of course, we want warnings and errors ([link](https://stackoverflow.com/questions/35911252/disable-tensorflow-debugging-information/42121886#42121886)).
# 
# 
# I keep getting an info message about the potentialities I could activate in tensorflow ([link](https://stackoverflow.com/questions/65298241/what-does-this-tensorflow-message-mean-any-side-effect-was-the-installation-su)). The rest are warning we get for missing libraries needed for GPU tensorflow, which is ok because we are not using GPUs.

# In[71]:


import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1' 
    #0 = all messages are logged (default behavior)
    #1 = INFO messages are not printed
    #2 = INFO and WARNING messages are not printed
    #3 = INFO, WARNING, and ERROR messages are not printed


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
                                      step=10) 
        #Step: A step of discretization. 
        #Note that high is modified if the range is not divisible by step. 
        #So I understand that 1,10,5 means that only 1,5 and 10 will be sampled
        #you should select a step size that is relevant for the parameter
        #e.g., a batch of 1 or 4 will be very similar, so we use a step of 10
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
                                      step=2)
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
    final_estimator = TransformedTargetRegressor(regressor=Pipeline([('scale', 
                                                                      preprocessing.StandardScaler()), 
                                                                     ('regressor',
                                                                      keras_regressor)]), 
                                                 transformer=preprocessing.StandardScaler())
    
    #cross validation
    score = cross_val_score(estimator=final_estimator, 
                                 X=X_train, 
                                 y=y_train, 
                                 scoring="r2",
                                 cv=shuffle_split,
                                 n_jobs=number_jobs).mean() 
    
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
    study = optuna.create_study(study_name = "yoruba_avg_1000kb_flex_sweep_optimization", 
                                    direction="maximize",
                                    storage='sqlite:///results/optuna_optimization/yoruba_avg_1000kb_flex_sweep_optimization.db', 
                                    sampler=optuna.samplers.TPESampler(seed=optuna_seed),
                                    load_if_exists=True)

    # Run the search. We are not using prunning becuase I have not found a way to get metrics while the cross-validation is performed in scikit-learn. So we need to finish a trial to get a score and decide if it is good or not.
    study.optimize(func=objective, 
                   n_trials=n_trials,
                   n_jobs=1, #we will parallelize using RDB, not n_jobs
                   gc_after_trial=True, #for recollectin of garbage
                   show_progress_bar=False)

    #we do not do anything else because even if the current process is finished, others could be still working and the best model is still not present. Of course, it would be printed at the end of that final process, but we would have several best models in the output. Better to do it in a separate script when the whole study is finished.


    '''
    #get the pruned and complete trials
    pruned_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED]
    complete_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]

    #print the numbers
    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    #see the best trial
    trial = study.best_trial
    f"Best trial is trial number {trial.number}"
    
    #see its parameters
    print("Best Score: ", trial.value)
    print("Best Params: ")
    for key, value in trial.params.items():
        print("  {}: {}".format(key, value))

    #make plots about the parameters
    #optuna.visualization.plot_contour(study)
    #optuna.visualization.plot_param_importances(study)

    # Create a kerasregressor instance with default parameters but including a get_reg function and verbose=0 (not info printed across epochs) and include it in a pipeline.
    best_estimator = TransformedTargetRegressor(regressor=Pipeline([('scale', 
                                                                preprocessing.StandardScaler()),     
                                                               ('regressor',
                                                                KerasRegressor(model=get_reg,
                                                                               verbose=0))]), 
                                           transformer=preprocessing.StandardScaler())
    best_estimator

    # Add now the best parameters of the keras regressor
    best_estimator.set_params(**study.best_trial.params)

    # Fit to train data    
    best_estimator.fit(X_train, y_train)
    
    # Calculate R2 in the training set    
    f"R2 in the training set {r2_score(y_train, best_estimator.predict(X_train))}"
    
    # Now predict in the final validation dataset. I have checked that predicting on validation does not influence the fit of the model, but if you fit the model to test, then the fit of the model changes and hence any prediction it makes
    
    # We are going to calculate the R2. Given this set has not been used for ANYTHING during training and evaluation (not even scaling), we can safely say that the R2 obtained is a generalization measure of our model.
    f"R2 in the test set {r2_score(y_test, best_estimator.predict(X_test))}"
    
    # Save the model:    
    joblib.dump(best_estimator, "results/optuna_optimization/optimized_model_train_dataset_avg_sweep_prob_optuna_"+str(optuna_seed)+".pkl")
       #https://stackoverflow.com/questions/51424312/how-to-save-gridsearchcv-object
    #joblib.load("results/dump_objects/optimized_model_train_dataset_optuna_avg_sweep_prob.pkl")
    '''

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