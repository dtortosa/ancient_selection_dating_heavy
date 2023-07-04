eso = pd.read_csv("./results/model_comparison/model_class_selection_yoruba_hg19.tsv", sep="\t", header=0)


import json
eso["best_params"].apply(lambda x: json.loads(x.replace("\'", "\"")))
eso["best_params"].apply(lambda x: json.loads(x.replace("\'", "\""))["regressor__alpha"])


#extract dict values from pandas
	#https://stackoverflow.com/q/57629435/12772630

#conver first to json 
#JSON only allows enclosing strings with double quotes you can manipulate the string
#https://stackoverflow.com/a/50257217/12772630
	#CHECK THIS TO AVOID ERRORS WITH REPLACE, maybe '' is used in some HP names?


