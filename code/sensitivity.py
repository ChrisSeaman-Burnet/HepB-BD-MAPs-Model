# Sensitivity Analysis
model_assumptions["f_map_adc"]=0    # No additional coverage from MAPs in facilities
model_assumptions["c_map_adc"]=0.01

#Replacement Coverage
f_map_rep=model_assumptions["f_map_rep"]
c_map_rep=model_assumptions["c_map_rep"]
            
f_map_rep[f_map_rep!=0]=0.01
c_map_rep[c_map_rep!=0]=0.01
            
model_assumptions["f_map_rep"]=f_map_rep
model_assumptions["c_map_rep"]=c_map_rep

#Run the analysis
ow_sens={}
sens_names=["prev_lb", "prev_ub", "prev_pe",
      "hbe_lb", "hbe_ub", "hbe_pe",
      "fac_lb", "fac_ub", "fac_pe",
      "bd_lb", "bd_ub", "bd_pe",
      "ve_lb", "ve_ub",
      "supc_lb", "supc_ub", "supc_pe",
      "com_lb", "com_ub", "com_pe",
      "deliv_lb", "deliv_ub", "deliv_pe",
      "outr_lb", "outr_pe",
      "dism_lb", "dism_ub", "dism_pe"]

for scen in range(3):
    for sen in sens_names:
        #HBV Prevalence (using ranges described in Table 2)
        if sen=="prev_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,3]=raw_inputs["set_pars_lb"].iloc[loc,12]
        if sen=="prev_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,3]=raw_inputs["set_pars_ub"].iloc[loc,12]
        if sen=="prev_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,3]=raw_inputs["set_pars"].iloc[loc,12]
        #HBeAg prevalence (using ranges described in Table 2)
        if sen=="hbe_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,4]=raw_inputs["set_pars_lb"].iloc[loc,13]
        if sen=="hbe_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,4]=raw_inputs["set_pars_ub"].iloc[loc,13]
        if sen=="hbe_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,4]=raw_inputs["set_pars"].iloc[loc,13]
        #Facility Births - as +/- 5% in UA; will use really low (5%) and really high (99%) in leiu 
        if sen=="fac_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,0]=0.05
        if sen=="fac_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,0]=0.99
        if sen=="fac_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,0]=raw_inputs["set_pars"].iloc[loc,9]
        #Baseline birth dose coverage - as +/- 5% in UA; will use really low (5%) and really high (99%) in leiu 
        if sen=="bd_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,1]=0.05
        if sen=="bd_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,1]=0.99
        if sen=="bd_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,1]=raw_inputs["set_pars"].iloc[loc,10]
        #Vaccine Effectiveness - excluded for manuscript as when replacement, low VE means an erroneously negative ICER (more DALYs for more cost)
        if sen=="ve_lb":
            eff_map=0.2
            vac_eff=ve_calcs(raw, runs, eff_map, eff_ctc)
        if sen=="ve_ub":
            eff_map=1
            vac_eff=ve_calcs(raw, runs, eff_map, eff_ctc)
        #Supply Chain Isolated to MAPs - maybe assume that MAPs will incur half/double the supply chain costs?
        if sen=="supc_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:17]=raw_inputs["set_pars"].iloc[loc,23:26]-(raw_inputs["sens_cost"].iloc[loc,8]/2)
        if sen=="supc_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:17]=raw_inputs["set_pars"].iloc[loc,23:26]+raw_inputs["sens_cost"].iloc[loc,8]
        if sen=="supc_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:17]=raw_inputs["set_pars"].iloc[loc,23:26]
        #Commodity Costs, by default, limited to cold chain and CTC vaccines (otherwise unknown)
        if sen=="com_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,12:14]=raw_inputs["set_pars"].iloc[loc,21:23]-(raw_inputs["sens_cost"].iloc[loc,10])+(raw_inputs["sens_cost"].iloc[loc,9])
                data_prepped["sett_pars"][loc,0,17:19]=raw_inputs["set_pars"].iloc[loc,25:27]-(raw_inputs["sens_cost"].iloc[loc,13])+(raw_inputs["sens_cost"].iloc[loc,12])
        if sen=="com_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,12:14]=raw_inputs["set_pars"].iloc[loc,21:23]-(raw_inputs["sens_cost"].iloc[loc,10])+(raw_inputs["sens_cost"].iloc[loc,11])
                data_prepped["sett_pars"][loc,0,17:19]=raw_inputs["set_pars"].iloc[loc,25:27]-(raw_inputs["sens_cost"].iloc[loc,13])+(raw_inputs["sens_cost"].iloc[loc,14])
        if sen=="com_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,12:14]=raw_inputs["set_pars"].iloc[loc,21:23]
                data_prepped["sett_pars"][loc,0,17:19]=raw_inputs["set_pars"].iloc[loc,25:27]
        #Human Resource (vaccine administration time) isolated to MAPs - lower bound minimal difference, upper bound a whole two minutes
        if sen=="deliv_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:16]=raw_inputs["set_pars"].iloc[loc,23:25]-(raw_inputs["sens_cost"].iloc[loc,16])+(raw_inputs["sens_cost"].iloc[loc,15]) #delivered by QHW
                data_prepped["sett_pars"][loc,0,16]=raw_inputs["set_pars"].iloc[loc,25]-(raw_inputs["sens_cost"].iloc[loc,19])+(raw_inputs["sens_cost"].iloc[loc,18]) #delivered by LHW
        if sen=="deliv_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:16]=raw_inputs["set_pars"].iloc[loc,23:25]-(raw_inputs["sens_cost"].iloc[loc,16])+(raw_inputs["sens_cost"].iloc[loc,17]) #delivered by QHW
                data_prepped["sett_pars"][loc,0,16]=raw_inputs["set_pars"].iloc[loc,25]-(raw_inputs["sens_cost"].iloc[loc,19])+(raw_inputs["sens_cost"].iloc[loc,20]) #delivered by LHW
        if sen=="deliv_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,14:16]=raw_inputs["set_pars"].iloc[loc,23:25] #delivered by QHW
                data_prepped["sett_pars"][loc,0,16]=raw_inputs["set_pars"].iloc[loc,25] #delivered by LHW
        #We can assume outreach costs for CC and CTC are upper bounds, we can estimate savings by cadre (QHW, LHW) - maybe 50% reduced for LHW, and 20% reduced for QHW?
        if sen=="outr_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,15]=raw_inputs["set_pars"].iloc[loc,24]-(raw_inputs["sens_cost"].iloc[loc,21]*0.2) #QHW in community
                data_prepped["sett_pars"][loc,0,16]=raw_inputs["set_pars"].iloc[loc,25]-(raw_inputs["sens_cost"].iloc[loc,21]*0.5) #LHW in community
        if sen=="outr_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,15]=raw_inputs["set_pars"].iloc[loc,24]
                data_prepped["sett_pars"][loc,0,16]=raw_inputs["set_pars"].iloc[loc,25]
        #Lower Bound - Tordrup, Upper Bound - HBV Calculator estimates
        if sen=="dism_lb":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,8:12]=raw_inputs["sens_cost"].iloc[loc,22:26]
        if sen=="dism_ub":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,8:12]=raw_inputs["sens_cost"].iloc[loc,26:30]
        if sen=="dism_pe":
            for loc in range(len(settings)):
                data_prepped["sett_pars"][loc,0,8:12]=raw_inputs["set_pars"].iloc[loc,17:21]
           
        ow_sens[str("s%.1d_"+sen)%scen]={}
        model_inputs=vac_weight(data_prepped,runs)
        model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
        model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
        ow_sens[str("s%.1d_"+sen)%scen]["fac_cov"]=model_inputs["fac_cov"]
        ow_sens[str("s%.1d_"+sen)%scen]["com_cov"]=model_inputs["com_cov"]
        ow_sens[str("s%.1d_"+sen)%scen]["init_conds"]=model_inputs["init_state"]
        results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
        outcomes=outcome_calcs(data_prepped, model_inputs, results, map_price, t_steps, dt, runs, settings)
        ow_sens[str("s%.1d_"+sen)%scen]["cov_cost"]=outcomes["cov_cost"]
        ow_sens[str("s%.1d_"+sen)%scen]["ce_out"]=outcomes["ce_out"]
        ow_sens[str("s%.1d_"+sen)%scen]["DALYs"]=outcomes["DALYs"]
        ow_sens[str("s%.1d_"+sen)%scen]["dis_cost"]=outcomes["dis_cost"]
        ow_sens[str("s%.1d_"+sen)%scen]["MAP_PTA"]=outcomes["MAP_PTA"]

sens_s1_excel=np.zeros((len(settings),28))
sens_s2_excel=np.zeros((len(settings),28))

for loc in range(len(settings)):
    for idx, name in enumerate(sens_names):
        sens_s1_excel[loc,idx]=-((ow_sens[str("s1_"+name)]["MAP_PTA"][loc,cpad,4]-ow_sens[str("s0_"+name)]["MAP_PTA"][loc,cpad,4]))/(ow_sens[str("s1_"+name)]["MAP_PTA"][loc,cpad,0]-ow_sens[str("s0_"+name)]["MAP_PTA"][loc,cpad,0])
        sens_s2_excel[loc,idx]=-((ow_sens[str("s2_"+name)]["MAP_PTA"][loc,cpad,4]-ow_sens[str("s0_"+name)]["MAP_PTA"][loc,cpad,4]))/(ow_sens[str("s2_"+name)]["MAP_PTA"][loc,cpad,0]-ow_sens[str("s0_"+name)]["MAP_PTA"][loc,cpad,0])

sens_s1_excel=pd.DataFrame(sens_s1_excel)
sens_s1_excel.insert(0, "setting", settings)
sens_s1_excel=sens_s1_excel.rename(columns=dict(zip(sens_s1_excel.iloc[:,1:].columns, sens_names)))
sens_s1_excel.to_csv(str("Additional 1% Only Sensitivity " + level +".csv"))

sens_s2_excel=pd.DataFrame(sens_s2_excel)
sens_s2_excel.insert(0, "setting", settings)
sens_s2_excel=sens_s2_excel.rename(columns=dict(zip(sens_s2_excel.iloc[:,1:].columns, sens_names)))
sens_s2_excel.to_csv(str("Additional and Replacement 1% Sensitivity " + level +".csv"))