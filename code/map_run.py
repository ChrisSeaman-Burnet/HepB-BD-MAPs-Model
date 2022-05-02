#   Set Working Directory
import os
os.chdir("C:/Users/seama/Dropbox/PLOS_MAP Update/Code and Results/code")

# Import Required Packages
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import itertools
import re
from data_extract import extract, ve_calcs, dataprep
from model_prep import vac_weight, vax_dist, model_init
from model_run import dis_model,outcome_calcs

#   Start: Model Run Time (64 bit, Windows 10 AMD7, 16GB ram)
start_time=datetime.now()

#"Yes": LMICs with no birth dose coverage use imputed estimates (WHO Region Averages)
#"No": LMICs with no birth dose coverage excluded
impute_bd = "No"
eff_map=1    #Baseline assumption of equal effectiveness
eff_ctc=1    #Baseline assumption of equal effectiveness

#Import Excel Data
databook="map_inputs_revision.xlsx"
raw=pd.ExcelFile(databook)

#Model Runs and Level
runs=1000#pd.read_excel(raw,0).iloc[1,1] #can be manually overwritten
level="LMICs"#pd.read_excel(raw,0).iloc[2,1] #"Regions" or "LMICs", can be manually overwritten


#Range of MAP price points to be analyzed --check that code is malliable to changes (baseline 0,10.25, 0.25); make a loop for MAP analyses so that price point indexing is done automatically
wastage=0.04 
map_price_unit=np.round(np.arange(0,5.05,0.05),2)
map_price=np.round(np.arange(0,5.05,0.05)*(1/(1-wastage)),2)

# Price now equal to MI4A Uniject via UNICEF-SD (US$1.55), CPI adjusted from 2017 USD to 2020
for idx,val in enumerate(map_price):
    if val==np.round(1.65*(1/(1-wastage)),2):
        cpad=int(idx)
    if val==np.round(3.30*(1/(1-wastage)),2):
        tcpad=int(idx)
    if val==np.round(5*(1/(1-wastage)),2):
        five=int(idx)

if runs ==1:
    print ("Model configuration: "+ level+"; "+str(runs)+" iteration")
else:
    print ("Model configuration: "+ level+"; "+str(runs)+" iterations")

#Extract and prepare the data
settings,model_assumptions,raw_inputs=extract(raw, level, impute_bd)
data_prepped=dataprep(raw_inputs, runs)
vac_eff=ve_calcs(raw, runs, eff_map, eff_ctc)
cov_val=model_assumptions["cov_val"]
rep_val=model_assumptions["rep_val"]

#   Model Preperation (steps, vaccine weights etc)

# Time Steps (24 steps in first year consistent with L-->A transition, then annual time steps)
y1_steps=24

t_steps=np.zeros((99+y1_steps))
t_steps[0:y1_steps]=np.linspace(0,1,24)
t_steps[y1_steps:]=np.arange(2,101,1)

# Annualizing Factor
dt=np.zeros((99+y1_steps))
dt[0:y1_steps]=1/y1_steps
dt[y1_steps:]=1

#Main Analysis
exec(open("main_analysis.py").read())   #includes supplemental analysis for coverage and replacement combinations

#One-way sensitivity analysis
if runs==1:
    exec(open("sensitivity.py").read())

#Supplemental Analyses (data presented only for aggregate all LMICs)
if level=="Regions":
    data_prepped["sett_pars"]=data_prepped["sett_pars"][-1,:,:].reshape(1,runs,41) #all LMICs will always be last
    g_ref=len(settings)-1
    settings=["All LMICs"]
    exec(open("supp_incs.py").read())
    exec(open("supp_CTC.py").read())
    exec(open("supp_rep.py").read())
    exec(open("supp_rd1.py").read())
    exec(open("supp_afc.py").read())

#   End: Model Run Time (64 bit, Windows 10 i7, 16GB ram)
print("Model Run Time ="+ str(datetime.now()-start_time))


