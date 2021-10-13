#Vaccine Weighting (Facility vs Community):
def vac_weight(data_prepped,runs):
    """ As per supplementary formula, this function weights birth dose coverage by
    birth location (facility or community)"""
    
    import numpy as np
    
    #Extract the Required Data
    bd_cov=data_prepped["sett_pars"][:,:,1]
    fac_b=data_prepped["sett_pars"][:,:,0]
    weight=data_prepped["glob_pars"][2,:]
    
    #Calculate Weighted Coverage
    w_fac_cov=np.zeros((len(bd_cov), runs, 1))
    w_com_cov=np.zeros((len(bd_cov), runs, 1))
    
    for loc in range(len(bd_cov)):
        for run in range(runs):
            if bd_cov[loc,run] <=0.5:
                w_fac_cov[loc, run, 0]=bd_cov[loc,run]+bd_cov[loc,run]*(((weight[run]-1)/(weight[run]+1))*(2*(1-fac_b[loc,run])))
                w_com_cov[loc, run, 0]=bd_cov[loc,run]-bd_cov[loc,run]*(((weight[run]-1)/(weight[run]+1))*(2*fac_b[loc,run]))
            elif bd_cov[loc,run] > 0.5:
                w_fac_cov[loc, run, 0]=bd_cov[loc,run]+(1-bd_cov[loc,run])*(((weight[run]-1)/(weight[run]+1))*(2*(1-fac_b[loc,run])))
                w_com_cov[loc, run, 0]=bd_cov[loc,run]-(1-bd_cov[loc,run])*(((weight[run]-1)/(weight[run]+1))*(2*fac_b[loc,run]))
                
    #Store in a dictionary 
    model_inputs={"w_fac_cov": w_fac_cov,
                  "w_com_cov": w_com_cov}
    
    return model_inputs

#Vaccine Distribution for Model Scenarios
def vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen):
    """Calculates additional and replacement coverage for each scenario, distributing across birth dose time strata. 
    scen: iterator for running the model scenarios.
    Updated: 20-04-2021; calculations double checked by hand, matches expected for set assumptions"""
    
    import numpy as np
    import pandas as pd
    
    #Facility Births
    pop=data_prepped["glob_pars"][0,:]
    fac_b=data_prepped["sett_pars"][:,:,0]
    
    #Weighted Vaccine Coverage
    fac_w=model_inputs["w_fac_cov"]
    com_w=model_inputs["w_com_cov"]
    
    b_fac_cov=np.zeros((len(settings)*runs,1,5))   
    b_com_cov=np.zeros((len(settings)*runs,1,5))   
    
    fac_cov=np.zeros((len(settings)*runs,3,5))
    com_cov=np.zeros((len(settings)*runs,4,5))
    
    fac_dist=model_assumptions["f_timing"]
    com_dist=model_assumptions["c_timing"]
        
    for loc in range(len(settings)):
        for run in range(runs):
            b_fac_cov[(loc*runs)+run,0,:]=pop[run]*fac_b[loc,run]*fac_w[loc,run,0]*fac_dist
            b_com_cov[(loc*runs)+run,0,:]=pop[run]*(1-fac_b[loc,run])*com_w[loc,run,0]*com_dist
                
            b_fac_cov[(loc*runs)+run,0,4]=(pop[run]*fac_b[loc,run])-sum(b_fac_cov[(loc*runs)+run,0,0:4])
            b_com_cov[(loc*runs)+run,0,4]=(pop[run]*(1-fac_b[loc,run]))-sum(b_com_cov[(loc*runs)+run,0,0:4])

    if scen==0: #cold chain baseline 
        for loc in range(len(settings)):
            for run in range(runs):
                for i in range(5):
                #Facility, Cold Chain
                    fac_cov[(loc*runs)+run,0,i]=b_fac_cov[(loc*runs)+run,0,i]
                #Community Cold Chain
                    com_cov[(loc*runs)+run,0,i]=b_com_cov[(loc*runs)+run,0,i]
                
    if scen==3: #CTC baseline (supplemental analysis)
        #Additional coverage from CTC
        fac_ctc_ac=model_assumptions["f_ctc_adc"]
        com_ctc_ac=model_assumptions["c_ctc_adc"]
        
        #Timing of additional CTC coverage
        fac_ctc_at=model_assumptions["f_ctc_adt"]
        com_ctc_at=model_assumptions["c_ctc_adt"]
        
        #Replacement CTC distribution
        fac_ctc_r=model_assumptions["f_ctc_rep"]
        com_ctc_r=model_assumptions["c_ctc_rep"]
        
        for loc in range(len(settings)):
            for run in range(runs):
                #Facility, Cold Chain
                fac_cov[(loc*runs)+run,0,0]=b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))
                fac_cov[(loc*runs)+run,0,1]=b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))
                fac_cov[(loc*runs)+run,0,2]=b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))
                fac_cov[(loc*runs)+run,0,3]=b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))
                fac_cov[(loc*runs)+run,0,4]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)
                #Community Cold Chain
                com_cov[(loc*runs)+run,0,0]=b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))
                com_cov[(loc*runs)+run,0,1]=b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))
                com_cov[(loc*runs)+run,0,2]=b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))
                com_cov[(loc*runs)+run,0,3]=b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))
                com_cov[(loc*runs)+run,0,4]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)
                #Facility CTC
                fac_cov[(loc*runs)+run,1,0]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[0])+(b_fac_cov[(loc*runs)+run,0,1]*fac_ctc_r[1,0])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,0])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,0])
                fac_cov[(loc*runs)+run,1,1]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[1])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,1])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,1])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,1])
                fac_cov[(loc*runs)+run,1,2]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[2])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,2])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,2])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,2])
                fac_cov[(loc*runs)+run,1,3]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[3])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,3])
                #Community CTC
                com_cov[(loc*runs)+run,1,0]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[0])+(b_com_cov[(loc*runs)+run,0,1]*com_ctc_r[1,0])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,0])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,0])
                com_cov[(loc*runs)+run,1,1]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[1])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,1])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,1])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,1])
                com_cov[(loc*runs)+run,1,2]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[2])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,2])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,2])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,2])
                com_cov[(loc*runs)+run,1,3]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[3])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,3])

    if scen==1: #Additional MAPs only, cold chain baseline
        #Additional Coverage from MAPs (note: additional coverage is ONLY provided by lay-health workers in the community)
        fac_map_ac=model_assumptions["f_map_adc"]
        com_map_ac=model_assumptions["c_map_adc"]
        
        fac_map_at=model_assumptions["f_map_adt"]
        com_map_at=model_assumptions["c_map_adt"]
        
        for loc in range(len(settings)):
            for run in range(runs):
                #Facility Cold Chain
                fac_cov[(loc*runs)+run,0,0]=b_fac_cov[(loc*runs)+run,0,0]
                fac_cov[(loc*runs)+run,0,1]=b_fac_cov[(loc*runs)+run,0,1]
                fac_cov[(loc*runs)+run,0,2]=b_fac_cov[(loc*runs)+run,0,2]
                fac_cov[(loc*runs)+run,0,3]=b_fac_cov[(loc*runs)+run,0,3]
                fac_cov[(loc*runs)+run,0,4]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_map_ac)
                #Community Cold Chain
                com_cov[(loc*runs)+run,0,0]=b_com_cov[(loc*runs)+run,0,0]
                com_cov[(loc*runs)+run,0,1]=b_com_cov[(loc*runs)+run,0,1]
                com_cov[(loc*runs)+run,0,2]=b_com_cov[(loc*runs)+run,0,2]
                com_cov[(loc*runs)+run,0,3]=b_com_cov[(loc*runs)+run,0,3]
                com_cov[(loc*runs)+run,0,4]=b_com_cov[(loc*runs)+run,0,4]*(1-com_map_ac)
                #Facility MAPs (Additional)
                fac_cov[(loc*runs)+run,2,0]=b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[0]
                fac_cov[(loc*runs)+run,2,1]=b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[1]
                fac_cov[(loc*runs)+run,2,2]=b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[2]
                fac_cov[(loc*runs)+run,2,3]=b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[3]
                #Community MAPs (Additional)
                com_cov[(loc*runs)+run,3,0]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[0]
                com_cov[(loc*runs)+run,3,1]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[1]
                com_cov[(loc*runs)+run,3,2]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[2]
                com_cov[(loc*runs)+run,3,3]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[3]
    
    if scen==4: #Additional MAPs only, CTC baseline (supplemental analysis)
        #Additional coverage from CTC
        fac_ctc_ac=model_assumptions["f_ctc_adc"]
        com_ctc_ac=model_assumptions["c_ctc_adc"]
        
        #Timing of additional CTC coverage
        fac_ctc_at=model_assumptions["f_ctc_adt"]
        com_ctc_at=model_assumptions["c_ctc_adt"]
        
        #Replacement CTC distribution
        fac_ctc_r=model_assumptions["f_ctc_rep"]
        com_ctc_r=model_assumptions["c_ctc_rep"]
        
        #Additional Coverage from MAPs (note: additional coverage is ONLY provided by lay-health workers in the community)
        fac_map_ac=model_assumptions["f_map_adc"]
        com_map_ac=model_assumptions["c_map_adc"]
        
        fac_map_at=model_assumptions["f_map_adt"]
        com_map_at=model_assumptions["c_map_adt"]
        
        for loc in range(len(settings)):
            for run in range(runs):
                #Facility Cold Chain
                fac_cov[(loc*runs)+run,0,0]=b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))
                fac_cov[(loc*runs)+run,0,1]=b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))
                fac_cov[(loc*runs)+run,0,2]=b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))
                fac_cov[(loc*runs)+run,0,3]=b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))
                fac_cov[(loc*runs)+run,0,4]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*(1-fac_map_ac) #As we are assuming CTC as a BASELINE, additional coverage needs to be apportioned like this
                #Community Cold Chain
                com_cov[(loc*runs)+run,0,0]=b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))
                com_cov[(loc*runs)+run,0,1]=b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))
                com_cov[(loc*runs)+run,0,2]=b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))
                com_cov[(loc*runs)+run,0,3]=b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))
                com_cov[(loc*runs)+run,0,4]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*(1-com_map_ac)
                #Facility CTC
                fac_cov[(loc*runs)+run,1,0]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[0])+(b_fac_cov[(loc*runs)+run,0,1]*fac_ctc_r[1,0])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,0])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,0])
                fac_cov[(loc*runs)+run,1,1]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[1])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,1])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,1])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,1])
                fac_cov[(loc*runs)+run,1,2]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[2])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,2])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,2])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,2])
                fac_cov[(loc*runs)+run,1,3]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[3])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,3])
                #Community CTC
                com_cov[(loc*runs)+run,1,0]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[0])+(b_com_cov[(loc*runs)+run,0,1]*com_ctc_r[1,0])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,0])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,0])
                com_cov[(loc*runs)+run,1,1]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[1])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,1])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,1])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,1])
                com_cov[(loc*runs)+run,1,2]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[2])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,2])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,2])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,2])
                com_cov[(loc*runs)+run,1,3]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[3])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,3])
                #Facility MAPs (Additional)
                fac_cov[(loc*runs)+run,2,0]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_map_at[0]
                fac_cov[(loc*runs)+run,2,1]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_map_at[1]
                fac_cov[(loc*runs)+run,2,2]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_map_at[2]
                fac_cov[(loc*runs)+run,2,3]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_map_at[3]
                #Community MAPs (Additional)
                com_cov[(loc*runs)+run,3,0]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[0]
                com_cov[(loc*runs)+run,3,1]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[1]
                com_cov[(loc*runs)+run,3,2]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[2]
                com_cov[(loc*runs)+run,3,3]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[3]
    
    if scen==2: #Additional and replacement MAPs, cold chain baseline
        #Additional Coverage from MAPs (note: additional coverage is ONLY provided by lay-health workers in the community)
        fac_map_ac=model_assumptions["f_map_adc"]
        com_map_ac=model_assumptions["c_map_adc"]
        
        fac_map_at=model_assumptions["f_map_adt"]
        com_map_at=model_assumptions["c_map_adt"]
        
        #Replacement Coverage from MAPs (note: replacement coverage is ONLY provided by qualified health workers in the community)
        fac_map_r=model_assumptions["f_map_rep"]
        com_map_r=model_assumptions["c_map_rep"]
        
        for loc in range(len(settings)):
            for run in range(runs):
                #Facility Cold Chain
                fac_cov[(loc*runs)+run,0,0]=b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_map_r[0,1]+fac_map_r[0,2]+fac_map_r[0,3]))
                fac_cov[(loc*runs)+run,0,1]=b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_map_r[1,0]+fac_map_r[1,2]+fac_map_r[1,3]))
                fac_cov[(loc*runs)+run,0,2]=b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_map_r[2,0]+fac_map_r[2,1]+fac_map_r[2,3]))
                fac_cov[(loc*runs)+run,0,3]=b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_map_r[3,0]+fac_map_r[3,1]+fac_map_r[3,2]))
                fac_cov[(loc*runs)+run,0,4]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_map_ac)
                #Community Cold Chain
                com_cov[(loc*runs)+run,0,0]=b_com_cov[(loc*runs)+run,0,0]*(1-(com_map_r[0,1]+com_map_r[0,2]+com_map_r[0,3]))
                com_cov[(loc*runs)+run,0,1]=b_com_cov[(loc*runs)+run,0,1]*(1-(com_map_r[1,0]+com_map_r[1,2]+com_map_r[1,3]))
                com_cov[(loc*runs)+run,0,2]=b_com_cov[(loc*runs)+run,0,2]*(1-(com_map_r[2,0]+com_map_r[2,1]+com_map_r[2,3]))
                com_cov[(loc*runs)+run,0,3]=b_com_cov[(loc*runs)+run,0,3]*(1-(com_map_r[3,0]+com_map_r[3,1]+com_map_r[3,2]))
                com_cov[(loc*runs)+run,0,4]=b_com_cov[(loc*runs)+run,0,4]*(1-com_map_ac)
                #Facility MAPs (additional and replacement)
                fac_cov[(loc*runs)+run,2,0]=(b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[0])+(b_fac_cov[(loc*runs)+run,0,1]*fac_map_r[1,0])+(b_fac_cov[(loc*runs)+run,0,2]*fac_map_r[2,0])+(b_fac_cov[(loc*runs)+run,0,3]*fac_map_r[3,0])
                fac_cov[(loc*runs)+run,2,1]=(b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[1])+(b_fac_cov[(loc*runs)+run,0,0]*fac_map_r[0,1])+(b_fac_cov[(loc*runs)+run,0,2]*fac_map_r[2,1])+(b_fac_cov[(loc*runs)+run,0,3]*fac_map_r[3,1])
                fac_cov[(loc*runs)+run,2,2]=(b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[2])+(b_fac_cov[(loc*runs)+run,0,0]*fac_map_r[0,2])+(b_fac_cov[(loc*runs)+run,0,1]*fac_map_r[1,2])+(b_fac_cov[(loc*runs)+run,0,3]*fac_map_r[3,2])
                fac_cov[(loc*runs)+run,2,3]=(b_fac_cov[(loc*runs)+run,0,4]*fac_map_ac*fac_map_at[3])+(b_fac_cov[(loc*runs)+run,0,0]*fac_map_r[0,3])+(b_fac_cov[(loc*runs)+run,0,1]*fac_map_r[1,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_map_r[2,3])
                #Community MAPs (replacement - QHW)
                com_cov[(loc*runs)+run,2,0]=(b_com_cov[(loc*runs)+run,0,1]*com_map_r[1,0])+(b_com_cov[(loc*runs)+run,0,2]*com_map_r[2,0])+(b_com_cov[(loc*runs)+run,0,3]*com_map_r[3,0])
                com_cov[(loc*runs)+run,2,1]=(b_com_cov[(loc*runs)+run,0,0]*com_map_r[0,1])+(b_com_cov[(loc*runs)+run,0,2]*com_map_r[2,1])+(b_com_cov[(loc*runs)+run,0,3]*com_map_r[3,1])
                com_cov[(loc*runs)+run,2,2]=(b_com_cov[(loc*runs)+run,0,0]*com_map_r[0,2])+(b_com_cov[(loc*runs)+run,0,1]*com_map_r[1,2])+(b_com_cov[(loc*runs)+run,0,3]*com_map_r[3,2])
                com_cov[(loc*runs)+run,2,3]=(b_com_cov[(loc*runs)+run,0,0]*com_map_r[0,3])+(b_com_cov[(loc*runs)+run,0,1]*com_map_r[1,3])+(b_com_cov[(loc*runs)+run,0,2]*com_map_r[2,3])
                #Community MAPs (additional - LHW)
                com_cov[(loc*runs)+run,3,0]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[0]
                com_cov[(loc*runs)+run,3,1]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[1]
                com_cov[(loc*runs)+run,3,2]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[2]
                com_cov[(loc*runs)+run,3,3]=b_com_cov[(loc*runs)+run,0,4]*com_map_ac*com_map_at[3]
    
    if scen==5: #Additional and replacement MAPs, CTC baseline (supplemental analysis)
        #Additional coverage from CTC
        fac_ctc_ac=model_assumptions["f_ctc_adc"]
        com_ctc_ac=model_assumptions["c_ctc_adc"]
        
        fac_ctc_at=model_assumptions["f_ctc_adt"]
        com_ctc_at=model_assumptions["c_ctc_adt"]
        
        #Replacement CTC distribution
        fac_ctc_r=model_assumptions["f_ctc_rep"]
        com_ctc_r=model_assumptions["c_ctc_rep"]
        
        #Additional Coverage from MAPs (note: additional coverage is ONLY provided by lay-health workers in the community)
        fac_map_ac=model_assumptions["f_map_adc"]
        com_map_ac=model_assumptions["c_map_adc"]
        
        fac_map_at=model_assumptions["f_map_adt"]
        com_map_at=model_assumptions["c_map_adt"]
        
        #Replacement Coverage from MAPs (note: replacement coverage is ONLY provided by qualified health workers in the community)
        fac_map_r=model_assumptions["f_map_rep"]
        com_map_r=model_assumptions["c_map_rep"]
        
        for loc in range(len(settings)):
            for run in range(runs):
                #Facility Cold Chain
                fac_cov[(loc*runs)+run,0,0]=b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))*(1-(fac_map_r[0,1]+fac_map_r[0,2]+fac_map_r[0,3]))
                fac_cov[(loc*runs)+run,0,1]=b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))*(1-(fac_map_r[1,0]+fac_map_r[1,2]+fac_map_r[1,3]))
                fac_cov[(loc*runs)+run,0,2]=b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))*(1-(fac_map_r[2,0]+fac_map_r[2,1]+fac_map_r[2,3]))
                fac_cov[(loc*runs)+run,0,3]=b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))*(1-(fac_map_r[3,0]+fac_map_r[3,1]+fac_map_r[3,2]))
                fac_cov[(loc*runs)+run,0,4]=b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*(1-fac_map_ac)
                #Community Cold Chain
                com_cov[(loc*runs)+run,0,0]=b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))*(1-(com_map_r[0,1]+com_map_r[0,2]+com_map_r[0,3]))
                com_cov[(loc*runs)+run,0,1]=b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))*(1-(com_map_r[1,0]+com_map_r[1,2]+com_map_r[1,3]))
                com_cov[(loc*runs)+run,0,2]=b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))*(1-(com_map_r[2,0]+com_map_r[2,1]+com_map_r[2,3]))
                com_cov[(loc*runs)+run,0,3]=b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))*(1-(com_map_r[3,0]+com_map_r[3,1]+com_map_r[3,2]))
                com_cov[(loc*runs)+run,0,4]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*(1-com_map_ac)
                #Facility CTC
                fac_cov[(loc*runs)+run,1,0]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[0])+(b_fac_cov[(loc*runs)+run,0,1]*fac_ctc_r[1,0])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,0])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,0])
                fac_cov[(loc*runs)+run,1,1]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[1])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,1])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,1])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,1])
                fac_cov[(loc*runs)+run,1,2]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[2])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,2])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,2])+(b_fac_cov[(loc*runs)+run,0,3]*fac_ctc_r[3,2])
                fac_cov[(loc*runs)+run,1,3]=(b_fac_cov[(loc*runs)+run,0,4]*fac_ctc_ac*fac_ctc_at[3])+(b_fac_cov[(loc*runs)+run,0,0]*fac_ctc_r[0,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[1,3])+(b_fac_cov[(loc*runs)+run,0,2]*fac_ctc_r[2,3])
                #Community CTC
                com_cov[(loc*runs)+run,1,0]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[0])+(b_com_cov[(loc*runs)+run,0,1]*com_ctc_r[1,0])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,0])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,0])
                com_cov[(loc*runs)+run,1,1]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[1])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,1])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,1])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,1])
                com_cov[(loc*runs)+run,1,2]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[2])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,2])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,2])+(b_com_cov[(loc*runs)+run,0,3]*com_ctc_r[3,2])
                com_cov[(loc*runs)+run,1,3]=(b_com_cov[(loc*runs)+run,0,4]*com_ctc_ac*com_ctc_at[3])+(b_com_cov[(loc*runs)+run,0,0]*com_ctc_r[0,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[1,3])+(b_com_cov[(loc*runs)+run,0,2]*com_ctc_r[2,3])
                #Facility MAPs (additional and replacement coverage)
                fac_cov[(loc*runs)+run,2,0]=(b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_ctc_at[0])+(b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))*fac_map_r[1,0])+(b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))*fac_map_r[2,0])+(b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))*fac_map_r[3,0])
                fac_cov[(loc*runs)+run,2,1]=(b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_ctc_at[1])+(b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))*fac_map_r[0,1])+(b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))*fac_map_r[2,1])+(b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))*fac_map_r[3,1])
                fac_cov[(loc*runs)+run,2,2]=(b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_ctc_at[2])+(b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))*fac_map_r[0,2])+(b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))*fac_map_r[1,2])+(b_fac_cov[(loc*runs)+run,0,3]*(1-(fac_ctc_r[3,0]+fac_ctc_r[3,1]+fac_ctc_r[3,2]))*fac_map_r[3,2])
                fac_cov[(loc*runs)+run,2,3]=(b_fac_cov[(loc*runs)+run,0,4]*(1-fac_ctc_ac)*fac_map_ac*fac_ctc_at[3])+(b_fac_cov[(loc*runs)+run,0,0]*(1-(fac_ctc_r[0,1]+fac_ctc_r[0,2]+fac_ctc_r[0,3]))*fac_map_r[0,3])+(b_fac_cov[(loc*runs)+run,0,1]*(1-(fac_ctc_r[1,0]+fac_ctc_r[1,2]+fac_ctc_r[1,3]))*fac_map_r[1,3])+(b_fac_cov[(loc*runs)+run,0,2]*(1-(fac_ctc_r[2,0]+fac_ctc_r[2,1]+fac_ctc_r[2,3]))*fac_map_r[2,3])
                #Community MAPs (replacement: QHW only)
                com_cov[(loc*runs)+run,2,0]=(b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))*com_map_r[1,0])+(b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))*com_map_r[2,0])+(b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))*com_map_r[3,0])
                com_cov[(loc*runs)+run,2,1]=(b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))*com_map_r[0,1])+(b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))*com_map_r[2,1])+(b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))*com_map_r[3,1])
                com_cov[(loc*runs)+run,2,2]=(b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))*com_map_r[0,2])+(b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))*com_map_r[1,2])+(b_com_cov[(loc*runs)+run,0,3]*(1-(com_ctc_r[3,0]+com_ctc_r[3,1]+com_ctc_r[3,2]))*com_map_r[3,2])
                com_cov[(loc*runs)+run,2,3]=(b_com_cov[(loc*runs)+run,0,0]*(1-(com_ctc_r[0,1]+com_ctc_r[0,2]+com_ctc_r[0,3]))*com_map_r[0,3])+(b_com_cov[(loc*runs)+run,0,1]*(1-(com_ctc_r[1,0]+com_ctc_r[1,2]+com_ctc_r[1,3]))*com_map_r[1,3])+(b_com_cov[(loc*runs)+run,0,2]*(1-(com_ctc_r[2,0]+com_ctc_r[2,1]+com_ctc_r[2,3]))*com_map_r[2,3])
                #Community MAPs (additional: LHW only)
                com_cov[(loc*runs)+run,3,0]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[0]
                com_cov[(loc*runs)+run,3,1]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[1]
                com_cov[(loc*runs)+run,3,2]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[2]
                com_cov[(loc*runs)+run,3,3]=b_com_cov[(loc*runs)+run,0,4]*(1-com_ctc_ac)*com_map_ac*com_map_at[3]
    
    model_inputs["bl_f"]=b_fac_cov
    model_inputs["bl_c"]=b_com_cov
    model_inputs["fac_cov"]=fac_cov
    model_inputs["com_cov"]=com_cov
    
    return model_inputs

#Model Initialization
def model_init(data_prepped, model_inputs, vac_eff, settings, runs):
    """Calculates coverage using vaccine effectiveness, and initializes model into Protected, Susceptible and Latent compartments.
    Checked to be accurate: 21/04/2021"""
    
    import numpy as np
    import pandas as pd
    
    fac_cov=model_inputs["fac_cov"]
    com_cov=model_inputs["com_cov"]
    
    cc_ve,cc_vf=vac_eff["ve_cold"],vac_eff["vf_cold"]
    ctc_ve, ctc_vf=vac_eff["ve_ctc"],vac_eff["vf_ctc"]
    map_ve, map_vf=vac_eff["ve_map"],vac_eff["vf_map"]
    
    bd_yes=np.zeros((len(settings)*runs,1,1))
    bd_no=np.zeros((len(settings)*runs,1,1))
    bd_total=np.zeros((len(settings)*runs,1,1)) #this can be a final check, and should always sum to 1000.
    
    for loc in range(len(settings)):
        for run in range(runs):
            bd_yes[(loc*runs)+run, 0,0]=sum(fac_cov[(loc*runs)+run,0,:]*cc_ve[:,run])+sum(fac_cov[(loc*runs)+run,1,:]*ctc_ve[:,run])+sum(fac_cov[(loc*runs)+run,2,:]*map_ve[:,run])+sum(com_cov[(loc*runs)+run,0,:]*cc_ve[:,run])+sum(com_cov[(loc*runs)+run,1,:]*ctc_ve[:,run])+sum(com_cov[(loc*runs)+run,2,:]*map_ve[:,run])+sum(com_cov[(loc*runs)+run,3,:]*map_ve[:,run])
            bd_no[(loc*runs)+run, 0,0]=sum(fac_cov[(loc*runs)+run,0,:]*cc_vf[:,run])+sum(fac_cov[(loc*runs)+run,1,:]*ctc_vf[:,run])+sum(fac_cov[(loc*runs)+run,2,:]*map_vf[:,run])+sum(com_cov[(loc*runs)+run,0,:]*cc_vf[:,run])+sum(com_cov[(loc*runs)+run,1,:]*ctc_vf[:,run])+sum(com_cov[(loc*runs)+run,2,:]*map_vf[:,run])+sum(com_cov[(loc*runs)+run,3,:]*map_vf[:,run])
            bd_total[(loc*runs)+run,0,0]=bd_yes[(loc*runs)+run,0,0]+bd_no[(loc*runs)+run,0,0]
    
    init_state=np.zeros((len(settings)*runs,1,4)) #Immune, Susceptible and Latent Infections plus a sum to test if still = 1000
    
    hbv3=data_prepped["sett_pars"][:,:,2]
    prev=data_prepped["sett_pars"][:,:,3]
    hbe=data_prepped["sett_pars"][:,:,4]
    epos=data_prepped["sett_pars"][:,:,5]
    eneg=data_prepped["sett_pars"][:,:,6]
    
    for loc in range(len(settings)):
        for run in range(runs):
            init_state[(loc*runs)+run,0,0]=(bd_yes[(loc*runs)+run,0,0]*hbv3[loc,run])+(bd_no[(loc*runs)+run,0,0]*(1-prev[loc,run])*hbv3[loc,run])+(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*hbe[loc,run]*(1-epos[loc,run])*hbv3[loc,run])+(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*(1-hbe[loc,run])*(1-eneg[loc,run])*hbv3[loc,run])
            init_state[(loc*runs)+run,0,1]=(bd_yes[(loc*runs)+run,0,0]*(1-hbv3[loc,run]))+(bd_no[(loc*runs)+run,0,0]*(1-prev[loc,run])*(1-hbv3[loc,run]))+(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*hbe[loc,run]*(1-epos[loc,run])*(1-hbv3[loc,run]))+(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*(1-hbe[loc,run])*(1-eneg[loc,run])*(1-hbv3[loc,run]))
            init_state[(loc*runs)+run,0,2]=(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*hbe[loc,run]*epos[loc,run])+(bd_no[(loc*runs)+run,0,0]*prev[loc,run]*(1-hbe[loc,run])*eneg[loc,run])
            init_state[(loc*runs)+run,0,3]=sum(init_state[(loc*runs)+run,0,0:3])
    
    model_inputs["bd_yes"]=bd_yes
    model_inputs["bd_no"]=bd_no
    model_inputs["init_state"]=init_state
    
    return model_inputs