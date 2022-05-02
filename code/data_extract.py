#Get data from excel into python
def extract(raw, level, impute_bd):
    """ Extracts required data from the excel workbook for modelling in Python"""
    
    import numpy as np
    import pandas as pd
    import re
    
    
    #Baseline timing, MAP and CTC coverage assumptions
    f_timing=np.array(pd.read_excel(raw,0).iloc[1,8:13]).astype(float)
    c_timing=np.array(pd.read_excel(raw,0).iloc[2,8:13]).astype(float)
    
    f_ctc_adc=pd.read_excel(raw,0).iloc[4,8]
    f_ctc_adt=np.array(pd.read_excel(raw,0).iloc[5,8:13]).astype(float)
    
    c_ctc_adc=pd.read_excel(raw,0).iloc[13,8]
    c_ctc_adt=np.array(pd.read_excel(raw,0).iloc[14,8:13]).astype(float)
    
    f_map_adc=pd.read_excel(raw,0).iloc[22,8]
    f_map_adt=np.array(pd.read_excel(raw,0).iloc[23,8:13]).astype(float)
    
    c_map_adc=pd.read_excel(raw,0).iloc[31,8]
    c_map_adt=np.array(pd.read_excel(raw,0).iloc[32,8:13]).astype(float)
    
    f_ctc_rep=np.array(pd.read_excel(raw,0).iloc[8:12, 8:12]).astype(float)
    c_ctc_rep=np.array(pd.read_excel(raw,0).iloc[17:21, 8:12]).astype(float)
    f_map_rep=np.array(pd.read_excel(raw,0).iloc[26:30, 8:12]).astype(float)
    c_map_rep=np.array(pd.read_excel(raw,0).iloc[35:39, 8:12]).astype(float)
    
    #Assumptions to iterate through for scenarios (coverage and timing); can be manually overwritten in code
    cov_val=pd.read_excel(raw,0).iloc[40,8]
    rep_val=pd.read_excel(raw,0).iloc[41,8]
    
    cov_val=list(np.float_(re.split(",", cov_val)))
    rep_val=list(np.float_(re.split(",", rep_val)))
    
    #Global Variables (excluding vaccine effectiveness)
    glob_input=pd.read_excel(raw,2).iloc[0:24,:]
    
    #Setting Specific Variables (remember to include/exclude imputed coverage values but drop variable afterwards)
    
    if level=="LMICs" and impute_bd=="No":
        set_pars=pd.read_excel(raw,3).iloc[0:136, 0:52].sort_values(by=["setting"])
        set_pars_lb=pd.read_excel(raw,3).iloc[138:274, 0:52].sort_values(by=["setting"])
        set_pars_ub=pd.read_excel(raw,3).iloc[276:412, 0:52].sort_values(by=["setting"])
        sens_cost=pd.read_excel(raw,5).sort_values(by=["setting"])
        
        set_pars=set_pars[set_pars["cov_impute"]=="No"]
        set_pars_lb=set_pars_lb[set_pars_lb["cov_impute"]=="No"]
        set_pars_ub=set_pars_ub[set_pars_ub["cov_impute"]=="No"]
        sens_cost=sens_cost[sens_cost["cov_impute"]=="No"]
        
        set_pars=set_pars.iloc[:,0:50].reset_index(drop=True)
        set_pars_lb=set_pars_lb.iloc[:,0:50].reset_index(drop=True)
        set_pars_ub=set_pars_ub.iloc[:,0:50].reset_index(drop=True)
        sens_cost=sens_cost.iloc[:,0:35].reset_index(drop=True)
        settings=list(sens_cost.iloc[:,0])
        
    if level=="LMICs" and impute_bd=="Yes":
        set_pars=pd.read_excel(raw,3).iloc[0:136, 0:50].sort_values(by=["setting"]).reset_index(drop=True)
        set_pars_lb=pd.read_excel(raw,3).iloc[138:274, 0:50].sort_values(by=["setting"]).reset_index(drop=True)
        set_pars_ub=pd.read_excel(raw,3).iloc[276:412, 0:50].sort_values(by=["setting"]).reset_index(drop=True)
        sens_cost=pd.read_excel(raw,5).iloc[0:136, 0:35].sort_values(by=["setting"]).reset_index(drop=True)
        settings=list(sens_cost.iloc[:,0])
    
    if level=="Regions":
        set_pars=pd.read_excel(raw,4).iloc[0:7, 0:50].reset_index(drop=True)
        set_pars_lb=pd.read_excel(raw,4).iloc[9:16,0:50].reset_index(drop=True)
        set_pars_ub=pd.read_excel(raw,4).iloc[18:25,0:50].reset_index(drop=True)
        sens_cost=pd.read_excel(raw,6).iloc[0:7, 0:35].reset_index(drop=True)
        settings=list(sens_cost.iloc[:,0])
   
    #Dictionary 1: Model Assumptions
    model_assumptions={"f_timing": f_timing, "c_timing": c_timing, "f_ctc_adc": f_ctc_adc,
                       "f_ctc_adt": f_ctc_adt, "c_ctc_adc":c_ctc_adc, "c_ctc_adt":c_ctc_adt,
                       "f_map_adc":f_map_adc, "f_map_adt": f_map_adt, "c_map_adc":c_map_adc,
                       "c_map_adt": c_map_adt, "f_ctc_rep":f_ctc_rep, "c_ctc_rep":c_ctc_rep,
                       "f_map_rep":f_map_rep, "c_map_rep": c_map_rep, "cov_val":cov_val,
                       "rep_val": rep_val}
    
    
    #Dictionary 2: Raw Model Inputs
    raw_inputs={"glob_input":glob_input, "set_pars":set_pars, "set_pars_lb":set_pars_lb,
                "set_pars_ub":set_pars_ub, "sens_cost":sens_cost}
    
    return settings,model_assumptions, raw_inputs

#Vaccine Effectiveness:
def ve_calcs(raw, runs, eff_map, eff_ctc):
    
    """Returns arrays of vaccine success/failure at each time strata, under different effectiveness assumptions. 
    Utilized for model initialization"""
    
    import numpy as np
    import pandas as pd
    
    ve_raw=np.array(pd.read_excel(raw,2).iloc[24:29,1:4])
    
    if runs==1:
        ve_scc=ve_raw[:,0].reshape(5,runs)
        vf_scc=1-ve_scc
        
        ve_ctc=ve_raw[:,0].reshape(5,runs)*eff_ctc
        vf_ctc=1-ve_ctc
        
        ve_map=ve_raw[:,0].reshape(5,runs)*eff_map
        vf_map=1-ve_map
       
   
    if runs > 1:
   #Order of operations: Uncertainity draws using ve_raw, calculate success and failure for each iteration.
       ve_pert=np.zeros((5,runs))
       
       for i in range(5):
           if i>=0 and i<=3:
               ve_pert[i,:]=np.random.triangular(ve_raw[i,1],ve_raw[i,0], ve_raw[i,2], size=runs)
           else:
               ve_pert[i,:]=np.random.uniform(ve_raw[i,1], ve_raw[i,2], size=runs)
        
       
       ve_scc=ve_pert[:,:]
       vf_scc=1-ve_pert[:,:]
       
       ve_ctc=ve_pert[:,:]*eff_ctc
       vf_ctc=1-ve_ctc[:,:]
       
       ve_map=ve_pert[:,:]*eff_map
       vf_map=1-ve_map[:,:]
      
    ve={"ve_cold": ve_scc, "ve_ctc": ve_ctc, "ve_map": ve_map,
       "vf_cold": vf_scc, "vf_ctc": vf_ctc, "vf_map": vf_map}
        
    return ve

#Monte Carlo Simulations:
def dataprep(raw_inputs, runs):
    
    """Draws from uncertainity bounds included in the imported spreadsheet. Random values not seeded, so some variations from run to run
    (if > 1 run) is to be expected"""
    
    import numpy as np
    import pandas as pd
    
    glob_input=raw_inputs["glob_input"]
    set_pars=raw_inputs["set_pars"]
    set_pars_lb=raw_inputs["set_pars_lb"]
    set_pars_ub=raw_inputs["set_pars_ub"]
    
    if runs==1:
        glob_pars=np.zeros((len(glob_input), runs))
        sett_pars=np.zeros((len(set_pars),runs,41)) #41 different (non-regional description) variables; excludes total births and population
        
        for run in range(runs):
            glob_pars[:,run]=glob_input.iloc[:,1]
        
        for loc in range(len(set_pars)):
            for run in range(runs):
                for i in range(41):
                    sett_pars[loc,run,i]=set_pars.iloc[loc,(i+9)]
    
    if runs > 1:
        np.random.seed(22072021)
        glob_pars=np.zeros((len(glob_input),runs))
        sett_pars=np.zeros((len(set_pars),runs,41))
        
        ## To Update (04-02-2022)
        ## Triangular Distribution: DALYs, Vaccine Effectiveness, Costs, Prevalence
        ## Uniform Distribution: Disease Progression, Vaccination Coverage, Facility Births 
        
        for i in range(len(glob_input)):
            if i>=20 and i<=23: #DALYs now triangular distribution
                glob_pars[i,:]=np.random.triangular(glob_input.iloc[i,2],glob_input.iloc[i,1], glob_input.iloc[i,3], runs)
            else:
                glob_pars[i,:]=np.random.uniform(glob_input.iloc[i,2], glob_input.iloc[i,3], runs)
    
        for loc in range(len(set_pars)):
            for i in range(41):
                if (i==3 or i==4) or (i>=7 and i<=17):
                    sett_pars[loc,:,i]=np.random.triangular(set_pars_lb.iloc[loc,i+9], set_pars.iloc[loc,i+9], set_pars_ub.iloc[loc,i+9], runs)
                else:
                    sett_pars[loc,:,i]=np.random.uniform(set_pars_lb.iloc[loc,i+9],set_pars_ub.iloc[loc,i+9], runs)
        
    data_prepped={"glob_pars":glob_pars, "sett_pars":sett_pars}
    
    return data_prepped