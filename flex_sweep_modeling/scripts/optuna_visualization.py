import sqlite3
con = sqlite3.connect("./results/optuna_optimization/prob_close_center/07_run/yoruba_flex_sweep_closest_window_center_optimization.db")


import pandas as pd

pd.DataFrame(con.execute("SELECT * FROM trial_values").fetchall())
#pd.read_sql_query("SELECT * FROM trial_values", con)
    #https://stackoverflow.com/questions/37051516/printing-a-properly-formatted-sqlite-table-in-python

con.execute("UPDATE trial_values SET value = -1 WHERE value < -1")
    #https://www.digitalocean.com/community/tutorials/how-to-use-the-sqlite3-module-in-python-3

con.commit()
    #https://stackoverflow.com/questions/4840772/python-sqlite3-update-not-updating

pd.DataFrame(con.execute(" \
    SELECT * \
    FROM trial_values \
    WHERE value < 0").fetchall())



#load study
import optuna
study = optuna.create_study(study_name = "yoruba_flex_sweep_closest_window_center_optimization", 
    direction="maximize",
    storage="sqlite:///results/optuna_optimization/prob_close_center/07_run/yoruba_flex_sweep_closest_window_center_optimization.db", 
    sampler=optuna.samplers.TPESampler(seed=45324),
    load_if_exists=True)

#
opt_history = optuna.visualization.plot_optimization_history(study)
opt_history.show()
    #https://optuna.readthedocs.io/en/stable/reference/visualization/generated/optuna.visualization.plot_optimization_history.html


plot_slice = optuna.visualization.plot_slice(study)
    #use params=[] to select specific parameters
plot_slice.show()
    #https://optuna.readthedocs.io/en/stable/reference/visualization/generated/optuna.visualization.plot_slice.html

#https://neptune.ai/blog/optuna-guide-how-to-monitor-hyper-parameter-optimization-runs

#changes for the new search. 
    #batch size
        #extend up to 300 to be sure that the peak at 200 stops.
    #epochs
        #extend to 900 to ensure is consistently decreasing after 500
    #loss
        #mse is the one with more negatives, but also the one with the highest values, and it is an usual loss, so we are going to keep it along with the rest.
    #activation func
        #tanh is the one with more negatives, but also the one with the highest values, and it is an usual act function, so we are going to keep it along with the rest.
    #init mode
        #identity is the one with more negatives, but also the one with the highest values, so we are going to keep it along with the rest.
    #number layers
        #add 2 and 4, because 3 is the highest
    #number of units
        #1000 because the highest is at the upper limit
    #regu L1/L2
        #100 pico (1e-10), to recreated the decrease from 1 micro to 10 nano. There we sill have high values.
    #weight constrain
        #go to 8 just in case, although very plateau
    #optimizer
        #adamax is the one with more negatives, but also the one with the highest values, and it is an usual act function, so we are going to keep it along with the rest.