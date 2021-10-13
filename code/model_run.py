def dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs):

    """ Runs the disease model, over time horizon t. dt parameter used to annualize values where required"""
    
    import numpy as np
    
    model_run=np.zeros((len(settings*runs), len(t_steps), 13))
    
    #Initial Model Conditions
    init_conds= model_inputs["init_state"]
    
    for loc in range(len(settings)):
        for run in range(runs):
            model_run[(loc*runs)+run,0,0]=init_conds[(loc*runs)+run,0,0]    #Immune
            model_run[(loc*runs)+run,0,1]=init_conds[(loc*runs)+run,0,1]    #Susceptible
            model_run[(loc*runs)+run,0,2]=init_conds[(loc*runs)+run,0,2]    #Latent 
            model_run[(loc*runs)+run,0,12]=init_conds[(loc*runs)+run,0,3]   #Cohort size at time(t)
            
    #Age Dependent Acute --> Chronic
    p_AC=data_prepped["glob_pars"][6,:]
    p_ACt=np.zeros((runs, len(t_steps)))
    
    # To only have a transfer at six months of age
    # for idx, yr in enumerate(t_steps):
    #     for run in range(runs):
    #        if yr==0.5:
    #            p_ACt[run,idx]=p_AC[run]
    
    for idx, yr in enumerate(t_steps):
        for run in range(runs):
            if yr < 0.5:
                p_ACt[run,idx]=p_AC[run]
            else:
                p_ACt[run,idx]=np.exp(-0.645*yr**(0.455))
    
    #All Cause Mortality
    acm=data_prepped["sett_pars"][:,:,19:40]
    acm_t=np.zeros((len(settings),runs,len(t_steps)))
    
    for loc in range(len(settings)):
        for run in range(runs):
            for idx,yr in enumerate(t_steps):
                if yr <1:
                    acm_t[loc,run,idx]=acm[loc,run,0]
                elif yr>= 1 and yr< 5:
                    acm_t[loc,run,idx]=acm[loc,run,1]
                elif yr>=5 and yr <95:
                    acm_t[loc,run,idx]=acm[loc,run,1+int(np.floor(yr)/5)]
                else:
                    acm_t[loc,run,idx]=acm[loc,run,-1]
    
    #Disease Progression and HBV mortality
    r_LA=data_prepped["glob_pars"][4,:].reshape(runs,1)
    r_AC=data_prepped["glob_pars"][5,:].reshape(runs,1)
    r_AZ=data_prepped["glob_pars"][7,:].reshape(runs,1)
    r_CZ=data_prepped["glob_pars"][8,:].reshape(runs,1)
    r_CCC=data_prepped["glob_pars"][9,:].reshape(runs,1)
    r_CHCC=data_prepped["glob_pars"][10,:].reshape(runs,1)
    r_CDC=data_prepped["glob_pars"][11,:].reshape(runs,1)
    r_CCHCC=data_prepped["glob_pars"][12,:].reshape(runs,1)
    r_DCHCC=data_prepped["glob_pars"][13,:].reshape(runs,1)
    
    mu_A=data_prepped["glob_pars"][14,:].reshape(runs,1)
    mu_C=data_prepped["glob_pars"][15,:].reshape(runs,1)
    mu_CC=data_prepped["glob_pars"][16,:].reshape(runs,1)
    mu_DC=data_prepped["glob_pars"][17,:].reshape(runs,1)
    mu_HCC=data_prepped["glob_pars"][18,:].reshape(runs,1)
    
    #Run the model:
    for loc in range(len(settings)):
        for run in range(runs):
            for t in range(len(t_steps)-1):
                #Immune
                    model_run[(loc*runs)+run,t+1,0]=model_run[(loc*runs)+run,t,0]+dt[t]*(((1-p_ACt[run,t])*r_AZ[run,0]*model_run[(loc*runs)+run,t,3])+(r_CZ[run,0]*model_run[(loc*runs)+run,t,4])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,0]))
                #Susceptible
                    model_run[(loc*runs)+run,t+1,1]=model_run[(loc*runs)+run,t,1]-dt[t]*acm_t[loc,run,t]*model_run[(loc*runs)+run,t,1]
                #Latent
                    model_run[(loc*runs)+run,t+1,2]=max(model_run[(loc*runs)+run,t,2]-dt[t]*((r_LA[run,0]*model_run[(loc*runs)+run,t,2])+(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,2])),0)
                #Acute
                    model_run[(loc*runs)+run,t+1,3]=max(model_run[(loc*runs)+run,t,3]+dt[t]*((r_LA[run,0]*model_run[(loc*runs)+run,t,2])-((1-p_ACt[run,t])*r_AZ[run,0]*model_run[(loc*runs)+run,t,3])-(p_ACt[run,t]*r_AC[run,0]*model_run[(loc*runs)+run,t,3])-(mu_A[run,0]*model_run[(loc*runs)+run,t,3])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,3])),0)
                #Chronic
                    model_run[(loc*runs)+run,t+1,4]=model_run[(loc*runs)+run,t,4]+dt[t]*((p_ACt[run,t]*r_AC[run,0]*model_run[(loc*runs)+run,t,3])-((r_CZ[run,0]+r_CCC[run,0]+r_CHCC[run,0])*model_run[(loc*runs)+run,t,4])-(mu_C[run,0]*model_run[(loc*runs)+run,t,4])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,4]))
                #Compensated Cirrhosis
                    model_run[(loc*runs)+run,t+1,5]=model_run[(loc*runs)+run,t,5]+dt[t]*((r_CCC[run,0]*model_run[(loc*runs)+run,t,4])-((r_CDC[run,0]+r_CCHCC[run,0])*model_run[(loc*runs)+run,t,5])-(mu_CC[run,0]*model_run[(loc*runs)+run,t,5])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,5]))
                #Decompensated Cirrhosis
                    model_run[(loc*runs)+run,t+1,6]=model_run[(loc*runs)+run,t,6]+dt[t]*((r_CDC[run,0]*model_run[(loc*runs)+run,t,5])-(r_DCHCC[run,0]*model_run[(loc*runs)+run,t,6])-(mu_DC[run,0]*model_run[(loc*runs)+run,t,6])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,6]))
                #HCC
                    model_run[(loc*runs)+run,t+1,7]=model_run[(loc*runs)+run,t,7]+dt[t]*((r_CHCC[run,0]*model_run[(loc*runs)+run,t,4])+(r_CCHCC[run,0]*model_run[(loc*runs)+run,t,5])+(r_DCHCC[run,0]*model_run[(loc*runs)+run,t,6])-(mu_HCC[run,0]*model_run[(loc*runs)+run,t,7])-(acm_t[loc,run,t]*model_run[(loc*runs)+run,t,7]))
                #Deaths (instant)
                    model_run[(loc*runs)+run,t+1,8]=dt[t]*((mu_A[run,0]*model_run[(loc*runs)+run,t,3])+(mu_C[run,0]*model_run[(loc*runs)+run,t,4])+(mu_CC[run,0]*model_run[(loc*runs)+run,t,5])+(mu_DC[run,0]*model_run[(loc*runs)+run,t,6])+(mu_HCC[run,0]*model_run[(loc*runs)+run,t,7]))
                #Deaths (cumulative)
                    model_run[(loc*runs)+run,t+1,9]=model_run[(loc*runs)+run,t,9]+dt[t]*((mu_A[run,0]*model_run[(loc*runs)+run,t,3])+(mu_C[run,0]*model_run[(loc*runs)+run,t,4])+(mu_CC[run,0]*model_run[(loc*runs)+run,t,5])+(mu_DC[run,0]*model_run[(loc*runs)+run,t,6])+(mu_HCC[run,0]*model_run[(loc*runs)+run,t,7]))
                #CHB Incidence (instant)
                    model_run[(loc*runs)+run,t+1,10]=dt[t]*p_ACt[run,t]*r_AC[run,0]*model_run[(loc*runs)+run,t,3]
                #CHB Incidence (cumulative)
                    model_run[(loc*runs)+run,t+1,11]=model_run[(loc*runs)+run,t,11]+dt[t]*p_ACt[run,t]*r_AC[run,0]*model_run[(loc*runs)+run,t,3]
                #Cohort Size 
                    model_run[(loc*runs)+run,t,12]=sum(model_run[(loc*runs)+run,t,0:8]) #starts at 1000, should always reduce - NOTED: Only minor error, will amend down the track.
    
    results={}
    results["model"]=model_run
    return results

def outcome_calcs(data_prepped, model_inputs, results, map_price, t_steps, dt, runs, settings):
    """Returns outcomes from model runs, with uncertainity bounds. Has been updated to include all the iteration for MAP price points"""
    
    import numpy as np
    import pandas as pd
    
    model_run=results["model"]
    fac_cov=model_inputs["fac_cov"]
    com_cov=model_inputs["com_cov"]
    discount=data_prepped["glob_pars"][3,:].reshape(runs,1)
    
    #Discounting Array
    disc=np.zeros((runs,len(t_steps)))
    
    for run in range(runs):
        for t in range(len(t_steps)):
            if t < 24:
                disc[run,t]=np.exp(-discount[run,0]*t*dt[t])
            else:
                disc[run,t]=np.exp(-discount[run,0]*(t-23)*dt[t])
    
    #Deaths and CHB cases
    ce_outcomes=np.zeros((len(settings), runs, 2))
    
    for loc in range(len(settings)):
        for run in range(runs):
            ce_outcomes[loc,run,0]=model_run[(loc*runs)+run,-1,11] #Chronic Hepatitis B Cases
            ce_outcomes[loc,run,1]=model_run[(loc*runs)+run,-1,9] #Hepatitis B Deaths
            
    #Vaccine Costs
    vac_cov_cost=np.zeros((len(settings), runs, 16)) #coverage for each intervention, then costs and total vaccine costs as final column 
    
    cc_fac=data_prepped["sett_pars"][:,:,12]
    cc_com=data_prepped["sett_pars"][:,:,13]
    ctc_fac=data_prepped["sett_pars"][:,:,17]
    ctc_com=data_prepped["sett_pars"][:,:,18]
    map_fac=data_prepped["sett_pars"][:,:,14]
    map_qhw=data_prepped["sett_pars"][:,:,15]
    map_lhw=data_prepped["sett_pars"][:,:,16]
    
    
    for loc in range(len(settings)):
        for run in range(runs):
            vac_cov_cost[loc,run,0]=sum(fac_cov[(loc*runs)+run,0,0:4])              #cold chain in facilities
            vac_cov_cost[loc,run,1]=sum(com_cov[(loc*runs)+run,0,0:4])              #cold chain in the community
            vac_cov_cost[loc,run,2]=sum(fac_cov[(loc*runs)+run,1,0:4])              #CTC in facilities
            vac_cov_cost[loc,run,3]=sum(com_cov[(loc*runs)+run,1,0:4])              #CTC in the community
            vac_cov_cost[loc,run,4]=sum(fac_cov[(loc*runs)+run,2,0:4])              #MAPs in facilities
            vac_cov_cost[loc,run,5]=sum(com_cov[(loc*runs)+run,2,0:4])              #MAPs in community delivered by QHW
            vac_cov_cost[loc,run,6]=sum(com_cov[(loc*runs)+run,3,0:4])              #MAPs in community delivered by LHW
            vac_cov_cost[loc,run,7]=vac_cov_cost[loc,run,0]*cc_fac[loc,run]         #Cost: cold chain in facilities
            vac_cov_cost[loc,run,8]=vac_cov_cost[loc,run,1]*cc_com[loc,run]         #Cost: cold chain in the community
            vac_cov_cost[loc,run,9]=vac_cov_cost[loc,run,2]*ctc_fac[loc,run]        #Cost: CTC in facilities
            vac_cov_cost[loc,run,10]=vac_cov_cost[loc,run,3]*ctc_com[loc,run]       #Cost: CTC in the community
            vac_cov_cost[loc,run,11]=vac_cov_cost[loc,run,4]*map_fac[loc,run]       #Cost: MAPs in facilities
            vac_cov_cost[loc,run,12]=vac_cov_cost[loc,run,5]*map_qhw[loc,run]       #Cost: MAPS in community delivered by QHW
            vac_cov_cost[loc,run,13]=vac_cov_cost[loc,run,6]*map_lhw[loc,run]       #Cost: MAPS in community delivered by LHW
            vac_cov_cost[loc,run,14]=sum(vac_cov_cost[loc,run,4:7])                 #Total number of MAPs delivered in each scenario, needed for price threshold analysis
            vac_cov_cost[loc,run,15]=sum(vac_cov_cost[loc,run,7:14])                #Cost: Total for scenario in each setting
    
    #DALYs
    life_exp=data_prepped["sett_pars"][:,:,40]
    d_A=data_prepped["glob_pars"][19,:].reshape(runs,1)
    d_C=data_prepped["glob_pars"][20,:].reshape(runs,1)
    d_CC=data_prepped["glob_pars"][21,:].reshape(runs,1)
    d_DC=data_prepped["glob_pars"][22,:].reshape(runs,1)
    d_HCC=data_prepped["glob_pars"][23,:].reshape(runs,1)
    
    life_rem=np.zeros((len(settings),runs, len(t_steps)))
    
    for loc in range(len(settings)):
        for run in range(runs):
            for t in range(len(t_steps)):
                life_rem[loc,run,t]=life_exp[loc,run]-t_steps[t]
                #Valuation of life-years after life expectancy as 0
                if life_rem[loc,run,t] <0:
                    life_rem[loc,run,t]=0
    
    yll=np.zeros((len(settings)*runs,len(t_steps),1))
    yld=np.zeros((len(settings)*runs,len(t_steps),5))
    
    for loc in range(len(settings)):
        for run in range(runs):
            for t in range(len(t_steps)):
                yld[(loc*runs)+run,t,0]=model_run[(loc*runs)+run,t,3]*(d_A[run,0]*dt[t])*disc[run,t]
                yld[(loc*runs)+run,t,1]=model_run[(loc*runs)+run,t,4]*(d_C[run,0]*dt[t])*disc[run,t]
                yld[(loc*runs)+run,t,2]=model_run[(loc*runs)+run,t,5]*(d_CC[run,0]*dt[t])*disc[run,t]
                yld[(loc*runs)+run,t,3]=model_run[(loc*runs)+run,t,6]*(d_DC[run,0]*dt[t])*disc[run,t]
                yld[(loc*runs)+run,t,4]=model_run[(loc*runs)+run,t,7]*(d_HCC[run,0]*dt[t])*disc[run,t]
                yll[(loc*runs)+run,t,0]=model_run[(loc*runs)+run,t,8]*life_rem[loc,run,t]*disc[run,t]
    
    dalys=np.zeros((len(settings),runs,3))
    
    for loc in range(len(settings)):
        for run in range(runs):
            dalys[loc,run,0]=sum(yld[(loc*runs)+run,:,0])+sum(yld[(loc*runs)+run,:,1])+sum(yld[(loc*runs)+run,:,2])+sum(yld[(loc*runs)+run,:,3])+sum(yld[(loc*runs)+run,:,4]) #YLD
            dalys[loc,run,1]=sum(yll[(loc*runs)+run,:,0]) #YLL
            dalys[loc,run,2]=sum(yld[(loc*runs)+run,:,0])+sum(yld[(loc*runs)+run,:,1])+sum(yld[(loc*runs)+run,:,2])+sum(yld[(loc*runs)+run,:,3])+sum(yld[(loc*runs)+run,:,4])+sum(yll[(loc*runs)+run,:,0]) #DALYs
    
    #Disease Management Costs (discounted at 3% per annum)
    c_diag=data_prepped["sett_pars"][:,:,7]
    c_C=data_prepped["sett_pars"][:,:,8]
    c_CC=data_prepped["sett_pars"][:,:,9]
    c_DC=data_prepped["sett_pars"][:,:,10]
    c_HCC=data_prepped["sett_pars"][:,:,11]
    
    t_d_cost=np.zeros((len(settings)*runs,len(t_steps),5))
    
    for loc in range(len(settings)):
        for run in range(runs):
            for t in range(len(t_steps)):
                t_d_cost[(loc*runs)+run,t,0]=model_run[(loc*runs)+run,t,10]*c_diag[loc,run]*disc[run,t]
                t_d_cost[(loc*runs)+run,t,1]=model_run[(loc*runs)+run,t,4]*c_C[loc,run]*dt[t]*disc[run,t]
                t_d_cost[(loc*runs)+run,t,2]=model_run[(loc*runs)+run,t,5]*c_CC[loc,run]*dt[t]*disc[run,t]
                t_d_cost[(loc*runs)+run,t,3]=model_run[(loc*runs)+run,t,6]*c_DC[loc,run]*dt[t]*disc[run,t]
                t_d_cost[(loc*runs)+run,t,4]=model_run[(loc*runs)+run,t,7]*c_HCC[loc,run]*dt[t]*disc[run,t]
    
    dis_cost=np.zeros((len(settings),runs,1))
    
    for loc in range(len(settings)):
        for run in range(runs):
            dis_cost[loc,run,0]=sum(t_d_cost[(loc*runs)+run,:,0])+sum(t_d_cost[(loc*runs)+run,:,1])+sum(t_d_cost[(loc*runs)+run,:,2])+sum(t_d_cost[(loc*runs)+run,:,3])+sum(t_d_cost[(loc*runs)+run,:,4])
    
    #MAP Pricepoint Analysis (note ICERs need to be done seperately!)
    map_costs=np.zeros((len(settings)*runs, len(map_price),5))
    
    for loc in range(len(settings)):
        for run in range(runs):
            for c in range(len(map_price)):
                map_costs[(loc*runs)+run,:,0]=dalys[loc,run,2]
                map_costs[(loc*runs)+run,c,1]=map_price[c]
                map_costs[(loc*runs)+run,c,2]=vac_cov_cost[loc,run,14]
                map_costs[(loc*runs)+run,c,3]=vac_cov_cost[loc,run,15]+dis_cost[loc,run,0]
                map_costs[(loc*runs)+run,c,4]=(map_costs[(loc*runs)+run,c,1]*map_costs[(loc*runs)+run,c,2])+map_costs[(loc*runs)+run,c,3]
    
    outcomes={}
    outcomes["cov_cost"]=vac_cov_cost
    outcomes["ce_out"]=ce_outcomes
    outcomes["DALYs"]=dalys
    outcomes["dis_cost"]=dis_cost
    outcomes["MAP_PTA"]=map_costs
    
    return outcomes