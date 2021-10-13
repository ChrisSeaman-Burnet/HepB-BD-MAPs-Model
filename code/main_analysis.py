#Main Analysis
main={}
for scen in range(3):
    if scen==0:
        main["cc_baseline"]={}
        model_inputs=vac_weight(data_prepped,runs)
        model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
        model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
        main["cc_baseline"]["fac_cov"]=model_inputs["fac_cov"]
        main["cc_baseline"]["com_cov"]=model_inputs["com_cov"]
        main["cc_baseline"]["init_conds"]=model_inputs["init_state"]
        results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
        outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
        main["cc_baseline"]["cov_cost"]=outcomes["cov_cost"]
        main["cc_baseline"]["ce_out"]=outcomes["ce_out"]
        main["cc_baseline"]["DALYs"]=outcomes["DALYs"]
        main["cc_baseline"]["dis_cost"]=outcomes["dis_cost"]
        main["cc_baseline"]["MAP_PTA"]=outcomes["MAP_PTA"]
        main["cc_baseline"]["Table"]=np.zeros((len(settings),43))
        
        #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
        for loc in range(len(settings)):
            for run in range(runs):
                for i in range(3):
                    main["cc_baseline"]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[loc,8])/10e2#births
                    main["cc_baseline"]["Table"][loc,i+1]=np.percentile(main["cc_baseline"]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                    main["cc_baseline"]["Table"][loc,i+4]=np.percentile( main["cc_baseline"]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                    main["cc_baseline"]["Table"][loc,i+7]=np.percentile( main["cc_baseline"]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                    main["cc_baseline"]["Table"][loc,i+10]=np.percentile( main["cc_baseline"]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                    main["cc_baseline"]["Table"][loc,i+13]=np.percentile( main["cc_baseline"]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                    main["cc_baseline"]["Table"][loc,i+16]=np.percentile( main["cc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                    main["cc_baseline"]["Table"][loc,i+19]=np.percentile( main["cc_baseline"]["DALYs"][loc,0:runs,2]-main["cc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                    main["cc_baseline"]["Table"][loc,i+22]=np.percentile( main["cc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                    main["cc_baseline"]["Table"][loc,i+25]=np.percentile( main["cc_baseline"]["dis_cost"][loc,0:runs,0]-main["cc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Cost
                    main["cc_baseline"]["Table"][loc,i+28]=np.percentile( main["cc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                    main["cc_baseline"]["Table"][loc,i+31]=np.percentile( main["cc_baseline"]["cov_cost"][loc,0:runs,15]-main["cc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
                 
    if scen==1: #Additional coverage only
        for cov in cov_val:
            
            #Additional Coverage
            model_assumptions["f_map_adc"]=0    #Additional coverage limited to births outside facility settings
            model_assumptions["c_map_adc"]=cov
            
            main["s1_c%.1d" % (cov*100)]={}
            model_inputs=vac_weight(data_prepped,runs)
            model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
            model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
            main["s1_c%.1d" % (cov*100)]["fac_cov"]=model_inputs["fac_cov"]
            main["s1_c%.1d" % (cov*100)]["com_cov"]=model_inputs["com_cov"]
            main["s1_c%.1d" % (cov*100)]["init_conds"]=model_inputs["init_state"]
            results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
            outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
            main["s1_c%.1d" % (cov*100)]["cov_cost"]=outcomes["cov_cost"]
            main["s1_c%.1d" % (cov*100)]["ce_out"]=outcomes["ce_out"]
            main["s1_c%.1d" % (cov*100)]["DALYs"]=outcomes["DALYs"]
            main["s1_c%.1d" % (cov*100)]["dis_cost"]=outcomes["dis_cost"]
            main["s1_c%.1d" % (cov*100)]["MAP_PTA"]=outcomes["MAP_PTA"]
            main["s1_c%.1d" % (cov*100)]["ICER"]=np.zeros((len(settings*runs), len(map_price), 6))
            main["s1_c%.1d" % (cov*100)]["Table"]=np.zeros((len(settings),43))
            
            #ICERs
            for loc in range(len(settings)):
                for run in range(runs):
                    for c in range(len(map_price)):
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,0]=main["s1_c%.1d" % (cov*100)]["MAP_PTA"][(loc*runs)+run,c,4]
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,1]=main["cc_baseline"]["MAP_PTA"][(loc*runs)+run,c,4]
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,2]=main["s1_c%.1d" % (cov*100)]["MAP_PTA"][(loc*runs)+run,c,0]
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,3]=main["cc_baseline"]["MAP_PTA"][(loc*runs)+run,c,0]
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,4]=map_price[c]
                        main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,5]=-((main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,0]- main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,1])/(main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,2]-main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,3]))

        #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
            for loc in range(len(settings)):
                for run in range(runs):
                    for i in range(3):
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[loc,8])/10e2#Births('000)
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+1]=np.percentile(main["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+4]=np.percentile(main["s1_c%.1d" % (cov*100)]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+7]=np.percentile(main["s1_c%.1d" % (cov*100)]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+10]=np.percentile(main["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+13]=np.percentile(main["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+16]=np.percentile(main["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+19]=np.percentile(main["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,2]-main["cc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+22]=np.percentile( main["s1_c%.1d" % (cov*100)]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+25]=np.percentile( main["s1_c%.1d" % (cov*100)]["dis_cost"][loc,0:runs,0]-main["cc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Costs
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+28]=np.percentile( main["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+31]=np.percentile( main["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,15]-main["cc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+34]=np.percentile(main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,cpad,5],25+(25*i))#ICER CPAD
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+37]=np.percentile(main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,tcpad,5],25+(25*i))#ICER 2xCPAD
                       main["s1_c%.1d" % (cov*100)]["Table"][loc,i+40]=np.percentile(main["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,five,5],25+(25*i))#ICER US$5
                        

    if scen==2: #Additional and replacement coverage only
        for cov in cov_val:
            for rep in rep_val:
                
                #Additional Coverage
                model_assumptions["f_map_adc"]=0    #Additional coverage limited to births outside facility settings
                model_assumptions["c_map_adc"]=cov
                
                #Replacement Coverage
                f_map_rep=model_assumptions["f_map_rep"]
                c_map_rep=model_assumptions["c_map_rep"]
            
                f_map_rep[f_map_rep!=0]=rep
                c_map_rep[c_map_rep!=0]=rep
            
                model_assumptions["f_map_rep"]=f_map_rep
                model_assumptions["c_map_rep"]=c_map_rep
                
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]={}
                model_inputs=vac_weight(data_prepped,runs)
                model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
                model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["fac_cov"]=model_inputs["fac_cov"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["com_cov"]=model_inputs["com_cov"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["init_conds"]=model_inputs["init_state"]
                results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
                outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"]=outcomes["cov_cost"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"]=outcomes["ce_out"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"]=outcomes["DALYs"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"]=outcomes["dis_cost"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"]=outcomes["MAP_PTA"]
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"]=np.zeros((len(settings*runs), len(map_price), 6))
                main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"]=np.zeros((len(settings),43))
                
                #ICERs
                for loc in range(len(settings)):
                    for run in range(runs):
                        for c in range(len(map_price)):
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]=main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,4]
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1]=main["cc_baseline"]["MAP_PTA"][(loc*runs)+run,c,4]
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]=main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,0]
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]=main["cc_baseline"]["MAP_PTA"][(loc*runs)+run,c,0]
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,4]=map_price[c]
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,5]=-((main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]- main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1])/(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]-main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]))
            #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
                for loc in range(len(settings)):
                    for run in range(runs):
                        for i in range(3):
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[loc,8])/10e2#Births('000)
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+1]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+4]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+7]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+10]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+13]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+16]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+19]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,2]-main["cc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+22]=np.percentile( main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+25]=np.percentile( main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"][loc,0:runs,0]-main["cc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Costs
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+28]=np.percentile( main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+31]=np.percentile( main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,15]-main["cc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+34]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,cpad,5],25+(25*i))#ICER CPAD
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+37]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,tcpad,5],25+(25*i))#ICER 2xCPAD
                            main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+40]=np.percentile(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,five,5],25+(25*i))#ICER US$5

##### Outputs
#Table 3: Main Analysis (will be regional in the paper, but need LMICs for supplement)
tab_labs=["births ('000)",
          "dMAP_lb", "dMAP_pe", "dMAP_ub",
          "CHB_lb", "CHB_pe", "CHB_ub", 
          "Deaths_lb", "Deaths_pe", "Deaths_ub",
          "YLD_lb", "YLD_pe", "YLD_ub", 
          "YLL_lb", "YLL_pe", "YLL_ub",
          "DALY_lb", "DALY_pe", "DALY_ub",
          "dDALY_lb", "dDALY_pe", "dDALY_ub",
          "DMC_lb", "DMC_pe", "DMC_ub",
          "dDMC_lb", "dDMC_pe", "dDMC_ub",
          "VC_lb", "VC_pe", "VC_ub", 
          "dVC_lb", "dVC_pe", "dVC_ub",
          "CPAD_ICER_lb", "CPAD_ICER_pe", "CPAD_ICER_ub",
          "CPAD2_ICER_lb", "CPAD2_ICER_pe", "CPAD2_ICER_ub", 
          "US5_ICER_lb", "US5_ICER_pe", "US5_ICER_ub"]

main_excel={}
for scen in range(3):
    if scen==0:
        main_excel["baseline"]=pd.DataFrame(main["cc_baseline"]["Table"])
        main_excel["baseline"].insert(0,"setting",settings)
        main_excel["baseline"]=main_excel["baseline"].rename(columns=dict(zip(main_excel["baseline"].iloc[:,1:].columns, tab_labs)))
    if scen==1:
        for cov in cov_val:
            main_excel["s1_c%.1d" % (cov*100)]=pd.DataFrame(main["s1_c%.1d" % (cov*100)]["Table"])
            main_excel["s1_c%.1d" % (cov*100)].insert(0,"setting", settings)
            main_excel["s1_c%.1d" % (cov*100)]=main_excel["s1_c%.1d" % (cov*100)].rename(columns=dict(zip(main_excel["s1_c%.1d" % (cov*100)].iloc[:,1:].columns, tab_labs)))
    if scen==2:
        for cov in cov_val:
            for rep in rep_val:
                main_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)]=pd.DataFrame(main["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"])
                main_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].insert(0,"setting", settings)
                main_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)]=main_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].rename(columns=dict(zip(main_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].iloc[:,1:].columns, tab_labs)))

writer=pd.ExcelWriter(str("main "+ level+" IQR"+ str(runs) + " runs impute "+str(impute_bd)+".xlsx"))
for df_name, df in main_excel.items():
    df.to_excel(writer, sheet_name=df_name)
writer.save()

#Figure 3: WHO Regions Price Threshold Plot  (AFRO, AMRO, EMRO, EURO, SEARO, WPRO, all LMICs)
if level=="Regions":
    main_s1_plot=np.zeros((7,len(map_price),4))
    main_s2_1_plot=np.zeros((7,len(map_price),4))
    main_s2_5_plot=np.zeros((7,len(map_price),4))
    main_s2_10_plot=np.zeros((7,len(map_price),4))
    
    
    for idx,reg in enumerate(settings):
        for run in range(runs):
            for c in range(len(map_price)):
                for i in range(3):
                    main_s1_plot[idx,c,0]=map_price[c]
                    main_s1_plot[idx,c,i+1]=np.percentile(main["s1_c1"]["ICER"][(idx*runs):(idx*runs)+runs, c,5], (i*25)+25)
                    
                    main_s2_1_plot[idx,c,0]=map_price[c]
                    main_s2_1_plot[idx,c,i+1]=np.percentile(main["s1_c1_r1"]["ICER"][(idx*runs):(idx*runs)+runs, c,5], (i*25)+25)
                    
                    main_s2_5_plot[idx,c,0]=map_price[c]
                    main_s2_5_plot[idx,c,i+1]=np.percentile(main["s1_c1_r5"]["ICER"][(idx*runs):(idx*runs)+runs, c,5], (i*25)+25)
                    
                    main_s2_10_plot[idx,c,0]=map_price[c]
                    main_s2_10_plot[idx,c,i+1]=np.percentile(main["s1_c1_r10"]["ICER"][(idx*runs):(idx*runs)+runs, c,5], (i*25)+25)
    
    plt.figure(figsize=(20, 10))
    for idx,reg in enumerate(settings):
        plt.subplot(int(np.ceil(len(settings)/4)),int(np.ceil(len(settings)/2)), idx+1)
        p1=plt.plot(map_price_unit[:], main_s1_plot[idx,:,2], color='red', marker="^", alpha=0.5)                                        #main_s1_plot[idx,:,0]
        p5=plt.fill_between(map_price_unit[:],main_s1_plot[idx,:,1], main_s1_plot[idx,:,3], alpha=0.3, color="red" )                                      #main_s1_plot[idx,:,0]
        p2=plt.plot(map_price_unit[:], main_s2_1_plot[idx,:,2], color='blue', marker="s", alpha=0.5)                                                       #main_s2_1_plot[idx,:,0]
        p6=plt.fill_between(map_price_unit[:], main_s2_1_plot[idx,:,1], main_s2_1_plot[idx,:,3], alpha=0.3, color="blue" )                                #main_s2_1_plot[idx,:,0]
        p3=plt.plot(map_price_unit[:], main_s2_5_plot[idx,:,2],  color='green', marker="o", alpha=0.5)                                                           #main_s2_5_plot[idx,:,0]
        p7=plt.fill_between(map_price_unit[:], main_s2_5_plot[idx,:,1], main_s2_5_plot[idx,:,3], alpha=0.3, color="green" )                                      #main_s2_5_plot[idx,:,0]
        p4=plt.plot(map_price_unit[:], main_s2_10_plot[idx,:,2], color='orange',marker="h", alpha=0.5)                                                           #main_s2_10_plot[idx,:,0]
        p8=plt.fill_between(map_price_unit[:], main_s2_10_plot[idx,:,1], main_s2_10_plot[idx,:,3], alpha=0.3, color="orange" )                                   #main_s2_10_plot[idx,:,0]
        plt.xlim(0,5)
        plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
        plt.title(str(settings[idx]))
        plt.grid()
        
        if idx==0 or idx==4:
            plt.ylabel("ICER ($US per DALY averted)")
        
        if idx==3 or idx==4 or idx==5 or idx==6:
            plt.xlabel("MAP Price per dose ($US)")
    
    plt.subplots_adjust(top=0.96, bottom=0.08, left=0.055, right=0.959, hspace=0.35, wspace=0.35)
    plt.legend(("Scenario 1: Additional MAPs (1%)","Scenario 1 + 1% Replacement MAPs","Scenario 1 + 5% Replacement MAPs","Scenario 1 + 10% Replacement MAPs"), loc="center left", bbox_to_anchor=(1.1,0.5), shadow=1, fancybox=1, fontsize=12) #put outside the box
    plt.savefig("WHO Regions Fig 2 "+str(runs)+" runs.png", dpi=300) 

# Figure 4: Prepare data for WTP thresholds at MAP price points (turn into bar plot in excel) for median ICER from S1:1%. S2 1% + 1,5,10%
if level=="LMICs":
    wtp_thresholds=np.array(raw_inputs["sens_cost"].iloc[:,31:])
    wtp_yn_daly=np.zeros((len(settings),12))    #Published WTP values
    wtp_yn_hgdp=np.zeros((len(settings),12))    #0.5 GDP per capita
    wtp_yn_ogdp=np.zeros((len(settings),12))    #1 GDP per capita
    wtp_yn_tgdp=np.zeros((len(settings),12))    #3 GDP per capita
   
    for loc in range(len(settings)):
        #WTP per DALY averted (literature)
        if wtp_thresholds[loc,0]>0:
            wtp_yn_daly[loc,0]=main["s1_c1"]["Table"][loc,35]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,1]=main["s1_c1"]["Table"][loc,38]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,2]=main["s1_c1"]["Table"][loc,41]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,3]=main["s1_c1_r1"]["Table"][loc,35]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,4]=main["s1_c1_r1"]["Table"][loc,38]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,5]=main["s1_c1_r1"]["Table"][loc,41]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,6]=main["s1_c1_r5"]["Table"][loc,35]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,7]=main["s1_c1_r5"]["Table"][loc,38]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,8]=main["s1_c1_r5"]["Table"][loc,41]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,9]=main["s1_c1_r10"]["Table"][loc,35]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,10]=main["s1_c1_r10"]["Table"][loc,38]-wtp_thresholds[loc,0]
            wtp_yn_daly[loc,11]=main["s1_c1_r10"]["Table"][loc,41]-wtp_thresholds[loc,0]
        #0.5 GDP threshold
        wtp_yn_hgdp[loc,0]=main["s1_c1"]["Table"][loc,35]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,1]=main["s1_c1"]["Table"][loc,38]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,2]=main["s1_c1"]["Table"][loc,41]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,3]=main["s1_c1_r1"]["Table"][loc,35]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,4]=main["s1_c1_r1"]["Table"][loc,38]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,5]=main["s1_c1_r1"]["Table"][loc,41]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,6]=main["s1_c1_r5"]["Table"][loc,35]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,7]=main["s1_c1_r5"]["Table"][loc,38]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,8]=main["s1_c1_r5"]["Table"][loc,41]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,9]=main["s1_c1_r10"]["Table"][loc,35]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,10]=main["s1_c1_r10"]["Table"][loc,38]-wtp_thresholds[loc,1]
        wtp_yn_hgdp[loc,11]=main["s1_c1_r10"]["Table"][loc,41]-wtp_thresholds[loc,1]
        #1 GDP threshold
        wtp_yn_ogdp[loc,0]=main["s1_c1"]["Table"][loc,35]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,1]=main["s1_c1"]["Table"][loc,38]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,2]=main["s1_c1"]["Table"][loc,41]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,3]=main["s1_c1_r1"]["Table"][loc,35]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,4]=main["s1_c1_r1"]["Table"][loc,38]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,5]=main["s1_c1_r1"]["Table"][loc,41]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,6]=main["s1_c1_r5"]["Table"][loc,35]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,7]=main["s1_c1_r5"]["Table"][loc,38]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,8]=main["s1_c1_r5"]["Table"][loc,41]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,9]=main["s1_c1_r10"]["Table"][loc,35]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,10]=main["s1_c1_r10"]["Table"][loc,38]-wtp_thresholds[loc,2]
        wtp_yn_ogdp[loc,11]=main["s1_c1_r10"]["Table"][loc,41]-wtp_thresholds[loc,2]
        #3 GDP threshold
        wtp_yn_tgdp[loc,0]=main["s1_c1"]["Table"][loc,35]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,1]=main["s1_c1"]["Table"][loc,38]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,2]=main["s1_c1"]["Table"][loc,41]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,3]=main["s1_c1_r1"]["Table"][loc,35]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,4]=main["s1_c1_r1"]["Table"][loc,38]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,5]=main["s1_c1_r1"]["Table"][loc,41]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,6]=main["s1_c1_r5"]["Table"][loc,35]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,7]=main["s1_c1_r5"]["Table"][loc,38]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,8]=main["s1_c1_r5"]["Table"][loc,41]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,9]=main["s1_c1_r10"]["Table"][loc,35]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,10]=main["s1_c1_r10"]["Table"][loc,38]-wtp_thresholds[loc,3]
        wtp_yn_tgdp[loc,11]=main["s1_c1_r10"]["Table"][loc,41]-wtp_thresholds[loc,3]

    wtp_labs=["s1_CPAD", "s1_2CPAD", "s1_5USD",
              "s2_1%_CP", "s2_1%_2CP", "s2_1%_5USD",
              "s2_5%_CP", "s2_5%_2CP", "s2_5%_5USD",
              "s2_10%_CP", "s2_10%_2CP", "s2_10%_5USD"]

    wtp_yn_daly=pd.DataFrame(wtp_yn_daly)
    wtp_yn_daly.insert(0,"setting", settings)
    wtp_yn_daly=wtp_yn_daly.rename(columns=dict(zip(wtp_yn_daly.iloc[:,1:].columns,wtp_labs)))
    wtp_yn_daly.to_csv(str("WTP per DALY averted_"+str(level)+"_"+str(impute_bd)+"_"+str(runs)+".csv"))
    
    wtp_yn_hgdp=pd.DataFrame(wtp_yn_hgdp)
    wtp_yn_hgdp.insert(0,"setting", settings)
    wtp_yn_hgdp=wtp_yn_hgdp.rename(columns=dict(zip(wtp_yn_hgdp.iloc[:,1:].columns,wtp_labs)))
    wtp_yn_hgdp.to_csv(str("0.5 per capita GDP_"+str(level)+"_"+str(impute_bd)+"_"+str(runs)+".csv"))
    
    wtp_yn_ogdp=pd.DataFrame(wtp_yn_ogdp)
    wtp_yn_ogdp.insert(0,"setting", settings)
    wtp_yn_ogdp=wtp_yn_ogdp.rename(columns=dict(zip(wtp_yn_ogdp.iloc[:,1:].columns,wtp_labs)))
    wtp_yn_ogdp.to_csv(str("per capita GDP_"+str(level)+"_"+str(impute_bd)+"_"+str(runs)+".csv"))
    
    wtp_yn_tgdp=pd.DataFrame(wtp_yn_tgdp)
    wtp_yn_tgdp.insert(0,"setting", settings)
    wtp_yn_tgdp=wtp_yn_tgdp.rename(columns=dict(zip(wtp_yn_tgdp.iloc[:,1:].columns,wtp_labs)))
    wtp_yn_tgdp.to_csv(str("3 per capita GDP_"+str(level)+"_"+str(impute_bd)+"_"+str(runs)+".csv"))
    

# Scatterplot for regions (cpad, tcpad, five)
if level == "Regions":
    cpad_plot=np.zeros((7,12,6))
    tcpad_plot=np.zeros((7,12,6))
    fusd_plot=np.zeros((7,12,6))
    
    map_scen=["s1_c1", "s1_c5", "s1_c10",
              "s1_c1_r1", "s1_c1_r5", "s1_c1_r10",
              "s1_c5_r1", "s1_c5_r5", "s1_c5_r10",
              "s1_c10_r1", "s1_c10_r5", "s1_c10_r10",]
    
    # DALYs Averted (0 LB, 1 PE, 2 UB)
    for loc in range(7):
        for idx,val in enumerate(map_scen):
            for pct in range(3):
                cpad_plot[loc,idx,pct]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,cpad,2]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,cpad,3], 25+(pct*25))
                tcpad_plot[loc,idx,pct]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,tcpad,2]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,tcpad,3], 25+(pct*25))
                fusd_plot[loc,idx,pct]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,five,2]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,five,3], 25+(pct*25))
    # Additional Costs (3 LB, 4 PE,  5 UB)
    for loc in range(7):
        for idx,val in enumerate(map_scen):
            for pct in range(3):
                cpad_plot[loc,idx,pct+3]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,cpad,0]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,cpad,1], 25+(pct*25))
                tcpad_plot[loc,idx,pct+3]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,tcpad,0]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,tcpad,1], 25+(pct*25))
                fusd_plot[loc,idx,pct+3]=np.percentile(main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,five,0]-main[str(val)]["ICER"][loc*runs:(loc*runs)+runs,five,1], 25+(pct*25))

    fig_s=plt.figure(figsize=(19.2, 9.08))
    
    markers=["s","D","X","s","s","s","D","D","D","X","X","X"] #s = 1% additional, D=5% additional, *=10% additional
    shades=["red","red","red","limegreen","royalblue", "darkviolet","limegreen","royalblue", "darkviolet","limegreen","royalblue", "darkviolet"] #bisque=no rep, lg=1% rep, rb=5% rep, dv=10% rep
    
    #### CPAD Price
    
    #AFRO CPAD
    ax1=fig_s.add_subplot(3,7,1)
    for i in range(12):
        ax1.scatter(cpad_plot[0,i,4], -cpad_plot[0,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax1.errorbar(cpad_plot[0,i,4], -cpad_plot[0,i,1], xerr=[[(cpad_plot[0,i,4]-cpad_plot[0,i,3])],[(cpad_plot[0,i,5]-cpad_plot[0,i,4])]], yerr=[[(cpad_plot[0,i,1]-cpad_plot[0,i,0])],[(cpad_plot[0,i,2]-cpad_plot[0,i,1])]], color=shades[i], elinewidth=0.5)
        ax1.axvline(0, color='k', linestyle=":")
        ax1.set(title="AFRO")#, ylabel="DALYs Averted (per 1000 births)")
        ax1.set_xlim([-10,150])
    
    #AMRO CPAD
    ax2=fig_s.add_subplot(3,7,2)
    for i in range(12):
        ax2.scatter(cpad_plot[1,i,4], -cpad_plot[1,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax2.errorbar(cpad_plot[1,i,4], -cpad_plot[1,i,1], xerr=[[(cpad_plot[1,i,4]-cpad_plot[1,i,3])],[(cpad_plot[1,i,5]-cpad_plot[1,i,4])]], yerr=[[(cpad_plot[1,i,1]-cpad_plot[1,i,0])],[(cpad_plot[1,i,2]-cpad_plot[1,i,1])]], color=shades[i], elinewidth=0.5)
        ax2.axvline(0, color='k', linestyle=":")
        ax2.set(title="AMRO")
        ax2.set_xlim([-50,80])
    
   #EMRO CPAD
    ax3=fig_s.add_subplot(3,7,3)
    for i in range(12):
        ax3.scatter(cpad_plot[2,i,4], -cpad_plot[2,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax3.errorbar(cpad_plot[2,i,4], -cpad_plot[2,i,1], xerr=[[(cpad_plot[2,i,4]-cpad_plot[2,i,3])],[(cpad_plot[2,i,5]-cpad_plot[2,i,4])]], yerr=[[(cpad_plot[2,i,1]-cpad_plot[2,i,0])],[(cpad_plot[2,i,2]-cpad_plot[2,i,1])]], color=shades[i], elinewidth=0.5)
        ax3.axvline(0, color='k', linestyle=":")
        ax3.set(title="EMRO")
        ax3.set_xlim([-10,150])
        
    #EURO CPAD
    ax4=fig_s.add_subplot(3,7,4)
    for i in range(12):
        ax4.scatter(cpad_plot[3,i,4], -cpad_plot[3,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax4.errorbar(cpad_plot[3,i,4], -cpad_plot[3,i,1], xerr=[[(cpad_plot[3,i,4]-cpad_plot[3,i,3])],[(cpad_plot[3,i,5]-cpad_plot[3,i,4])]], yerr=[[(cpad_plot[3,i,1]-cpad_plot[3,i,0])],[(cpad_plot[3,i,2]-cpad_plot[3,i,1])]], color=shades[i], elinewidth=0.5)
        ax4.axvline(0, color='k', linestyle=":")
        ax4.set(title="EURO")
        ax4.set_xlim([-10,150])
        
   #SEARO CPAD
    ax5=fig_s.add_subplot(3,7,5)
    for i in range(12):
        ax5.scatter(cpad_plot[4,i,4], -cpad_plot[4,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax5.errorbar(cpad_plot[4,i,4], -cpad_plot[4,i,1], xerr=[[(cpad_plot[4,i,4]-cpad_plot[4,i,3])],[(cpad_plot[4,i,5]-cpad_plot[4,i,4])]], yerr=[[(cpad_plot[4,i,1]-cpad_plot[4,i,0])],[(cpad_plot[4,i,2]-cpad_plot[4,i,1])]], color=shades[i], elinewidth=0.5)
        ax5.axvline(0, color='k', linestyle=":")
        ax5.set(title="SEARO")
        ax5.set_xlim([-100,150])
        
    #WPRO CPAD
    ax6=fig_s.add_subplot(3,7,6)
    for i in range(12):
        ax6.scatter(cpad_plot[5,i,4], -cpad_plot[5,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax6.errorbar(cpad_plot[5,i,4], -cpad_plot[5,i,1], xerr=[[(cpad_plot[5,i,4]-cpad_plot[5,i,3])],[(cpad_plot[5,i,5]-cpad_plot[5,i,4])]], yerr=[[(cpad_plot[5,i,1]-cpad_plot[5,i,0])],[(cpad_plot[5,i,2]-cpad_plot[5,i,1])]], color=shades[i], elinewidth=0.5)
        ax6.axvline(0, color='k', linestyle=":")
        ax6.set(title="WPRO")
        ax6.set_xlim([-150,50])
  
    #All LMICs CPAD
    ax7=fig_s.add_subplot(3,7,7)
    for i in range(12):
        ax7.scatter(cpad_plot[6,i,4], -cpad_plot[6,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax7.errorbar(cpad_plot[6,i,4], -cpad_plot[6,i,1], xerr=[[(cpad_plot[6,i,4]-cpad_plot[6,i,3])],[(cpad_plot[6,i,5]-cpad_plot[6,i,4])]], yerr=[[(cpad_plot[6,i,1]-cpad_plot[6,i,0])],[(cpad_plot[6,i,2]-cpad_plot[6,i,1])]], color=shades[i], elinewidth=0.5)
        ax7.axvline(0, color='k', linestyle=":")
        ax7.set(title="All LMICs")
        ax7.set_xlim([-100,100])

    #### 2xCPAD price

    #AFRO TCPAD
    ax8=fig_s.add_subplot(3,7,8)
    for i in range(12):
        ax8.scatter(tcpad_plot[0,i,4], -tcpad_plot[0,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax8.errorbar(tcpad_plot[0,i,4], -tcpad_plot[0,i,1], xerr=[[(tcpad_plot[0,i,4]-tcpad_plot[0,i,3])],[(tcpad_plot[0,i,5]-tcpad_plot[0,i,4])]], yerr=[[(tcpad_plot[0,i,1]-tcpad_plot[0,i,0])],[(tcpad_plot[0,i,2]-tcpad_plot[0,i,1])]], color=shades[i], elinewidth=0.5)
        ax8.axvline(0, color='k', linestyle=":")
        ax8.set(ylabel="DALYs Averted (per 1000 births)")
        ax8.set_xlim([-10,150])
    
    #AMRO TCPAD
    ax9=fig_s.add_subplot(3,7,9)
    for i in range(12):
        ax9.scatter(tcpad_plot[1,i,4], -tcpad_plot[1,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax9.errorbar(tcpad_plot[1,i,4], -tcpad_plot[1,i,1], xerr=[[(tcpad_plot[1,i,4]-tcpad_plot[1,i,3])],[(tcpad_plot[1,i,5]-tcpad_plot[1,i,4])]], yerr=[[(tcpad_plot[1,i,1]-tcpad_plot[1,i,0])],[(tcpad_plot[1,i,2]-tcpad_plot[1,i,1])]], color=shades[i], elinewidth=0.5)
        ax9.axvline(0, color='k', linestyle=":")
        ax9.set_xlim([-50,80])
    
   #EMRO TCPAD
    ax10=fig_s.add_subplot(3,7,10)
    for i in range(12):
        ax10.scatter(tcpad_plot[2,i,4], -tcpad_plot[2,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax10.errorbar(tcpad_plot[2,i,4], -tcpad_plot[2,i,1], xerr=[[(tcpad_plot[2,i,4]-tcpad_plot[2,i,3])],[(tcpad_plot[2,i,5]-tcpad_plot[2,i,4])]], yerr=[[(tcpad_plot[2,i,1]-tcpad_plot[2,i,0])],[(tcpad_plot[2,i,2]-tcpad_plot[2,i,1])]], color=shades[i], elinewidth=0.5)
        ax10.axvline(0, color='k', linestyle=":")
        ax10.set_xlim([-10,150])
        
    #EURO TCPAD
    ax11=fig_s.add_subplot(3,7,11)
    for i in range(12):
        ax11.scatter(tcpad_plot[3,i,4], -tcpad_plot[3,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax11.errorbar(tcpad_plot[3,i,4], -tcpad_plot[3,i,1], xerr=[[(tcpad_plot[3,i,4]-tcpad_plot[3,i,3])],[(tcpad_plot[3,i,5]-tcpad_plot[3,i,4])]], yerr=[[(tcpad_plot[3,i,1]-tcpad_plot[3,i,0])],[(tcpad_plot[3,i,2]-tcpad_plot[3,i,1])]], color=shades[i], elinewidth=0.5)
        ax11.axvline(0, color='k', linestyle=":")
        ax11.set_xlim([-10,150])
        
   #SEARO TCPAD
    ax12=fig_s.add_subplot(3,7,12)
    for i in range(12):
        ax12.scatter(tcpad_plot[4,i,4], -tcpad_plot[4,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax12.errorbar(tcpad_plot[4,i,4], -tcpad_plot[4,i,1], xerr=[[(tcpad_plot[4,i,4]-tcpad_plot[4,i,3])],[(tcpad_plot[4,i,5]-tcpad_plot[4,i,4])]], yerr=[[(tcpad_plot[4,i,1]-tcpad_plot[4,i,0])],[(tcpad_plot[4,i,2]-tcpad_plot[4,i,1])]], color=shades[i], elinewidth=0.5)
        ax12.axvline(0, color='k', linestyle=":")
        ax12.set_xlim([-100,150])
        
    #WPRO TCPAD
    ax13=fig_s.add_subplot(3,7,13)
    for i in range(12):
        ax13.scatter(tcpad_plot[5,i,4], -tcpad_plot[5,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax13.errorbar(tcpad_plot[5,i,4], -tcpad_plot[5,i,1], xerr=[[(tcpad_plot[5,i,4]-tcpad_plot[5,i,3])],[(tcpad_plot[5,i,5]-tcpad_plot[5,i,4])]], yerr=[[(tcpad_plot[5,i,1]-tcpad_plot[5,i,0])],[(tcpad_plot[5,i,2]-tcpad_plot[5,i,1])]], color=shades[i], elinewidth=0.5)
        ax13.axvline(0, color='k', linestyle=":")
        ax13.set_xlim([-150,50])
        
    #All LMICs TCPAD
    ax14=fig_s.add_subplot(3,7,14)
    for i in range(12):
        ax14.scatter(tcpad_plot[6,i,4], -tcpad_plot[6,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax14.errorbar(tcpad_plot[6,i,4], -tcpad_plot[6,i,1], xerr=[[(tcpad_plot[6,i,4]-tcpad_plot[6,i,3])],[(tcpad_plot[6,i,5]-tcpad_plot[6,i,4])]], yerr=[[(tcpad_plot[6,i,1]-tcpad_plot[6,i,0])],[(tcpad_plot[6,i,2]-tcpad_plot[6,i,1])]], color=shades[i], elinewidth=0.5)
        ax14.axvline(0, color='k', linestyle=":")
        ax14.set_xlim([-100,100])
        
    
    ### 5USD per dose
    
    #AFRO FIVE
    ax15=fig_s.add_subplot(3,7,15)
    for i in range(12):
        ax15.scatter(fusd_plot[0,i,4], -fusd_plot[0,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax15.errorbar(fusd_plot[0,i,4], -fusd_plot[0,i,1], xerr=[[(fusd_plot[0,i,4]-fusd_plot[0,i,3])],[(fusd_plot[0,i,5]-fusd_plot[0,i,4])]], yerr=[[(fusd_plot[0,i,1]-fusd_plot[0,i,0])],[(fusd_plot[0,i,2]-fusd_plot[0,i,1])]], color=shades[i], elinewidth=0.5)
        ax15.axvline(0, color='k', linestyle=":")
        #ax15.set(ylabel="DALYs Averted (per 1000 births)", xlabel="Additional Costs (USD per 1000 briths)")
        ax15.set_xlim([-10,150])
    
    #AMRO FIVE
    ax16=fig_s.add_subplot(3,7,16)
    for i in range(12):
        ax16.scatter(fusd_plot[1,i,4], -fusd_plot[1,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax16.errorbar(fusd_plot[1,i,4], -fusd_plot[1,i,1], xerr=[[(fusd_plot[1,i,4]-fusd_plot[1,i,3])],[(fusd_plot[1,i,5]-fusd_plot[1,i,4])]], yerr=[[(fusd_plot[1,i,1]-fusd_plot[1,i,0])],[(fusd_plot[1,i,2]-fusd_plot[1,i,1])]], color=shades[i], elinewidth=0.5)
        ax16.axvline(0, color='k', linestyle=":")
        #ax16.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax16.set_xlim([-50,80])
    
   #EMRO FIVE
    ax17=fig_s.add_subplot(3,7,17)
    for i in range(12):
        ax17.scatter(fusd_plot[2,i,4], -fusd_plot[2,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax17.errorbar(fusd_plot[2,i,4], -fusd_plot[2,i,1], xerr=[[(fusd_plot[2,i,4]-fusd_plot[2,i,3])],[(fusd_plot[2,i,5]-fusd_plot[2,i,4])]], yerr=[[(fusd_plot[2,i,1]-fusd_plot[2,i,0])],[(fusd_plot[2,i,2]-fusd_plot[2,i,1])]], color=shades[i], elinewidth=0.5)
        ax17.axvline(0, color='k', linestyle=":")
        #ax17.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax17.set_xlim([-10,150])
        
    #EURO FIVE
    ax18=fig_s.add_subplot(3,7,18)
    for i in range(12):
        ax18.scatter(fusd_plot[3,i,4], -fusd_plot[3,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax18.errorbar(fusd_plot[3,i,4], -fusd_plot[3,i,1], xerr=[[(fusd_plot[3,i,4]-fusd_plot[3,i,3])],[(fusd_plot[3,i,5]-fusd_plot[3,i,4])]], yerr=[[(fusd_plot[3,i,1]-fusd_plot[3,i,0])],[(fusd_plot[3,i,2]-fusd_plot[3,i,1])]], color=shades[i], elinewidth=0.5)
        ax18.axvline(0, color='k', linestyle=":")
        ax18.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax18.set_xlim([-10,150])
        
   #SEARO FIVE
    ax19=fig_s.add_subplot(3,7,19)
    for i in range(12):
        ax19.scatter(fusd_plot[4,i,4], -fusd_plot[4,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax19.errorbar(fusd_plot[4,i,4], -fusd_plot[4,i,1], xerr=[[(fusd_plot[4,i,4]-fusd_plot[4,i,3])],[(fusd_plot[4,i,5]-fusd_plot[4,i,4])]], yerr=[[(fusd_plot[4,i,1]-fusd_plot[4,i,0])],[(fusd_plot[4,i,2]-fusd_plot[4,i,1])]], color=shades[i], elinewidth=0.5)
        ax19.axvline(0, color='k', linestyle=":")
        #ax19.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax19.set_xlim([-100,150])
        
    #WPRO FIVE
    ax20=fig_s.add_subplot(3,7,20)
    for i in range(12):
        ax20.scatter(fusd_plot[5,i,4], -fusd_plot[5,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax20.errorbar(fusd_plot[5,i,4], -fusd_plot[5,i,1], xerr=[[(fusd_plot[5,i,4]-fusd_plot[5,i,3])],[(fusd_plot[5,i,5]-fusd_plot[5,i,4])]], yerr=[[(fusd_plot[5,i,1]-fusd_plot[5,i,0])],[(fusd_plot[5,i,2]-fusd_plot[5,i,1])]], color=shades[i], elinewidth=0.5)
        ax20.axvline(0, color='k', linestyle=":")
        #ax20.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax20.set_xlim([-150,50])
        
    #All LMICs FIVE
    ax21=fig_s.add_subplot(3,7,21)
    for i in range(12):
        ax21.scatter(fusd_plot[6,i,4], -fusd_plot[6,i,1], marker=markers[i],color=shades[i], alpha=0.5, s=40)
        ax21.errorbar(fusd_plot[6,i,4], -fusd_plot[6,i,1], xerr=[[(fusd_plot[6,i,4]-fusd_plot[6,i,3])],[(fusd_plot[6,i,5]-fusd_plot[6,i,4])]], yerr=[[(fusd_plot[6,i,1]-fusd_plot[6,i,0])],[(fusd_plot[6,i,2]-fusd_plot[6,i,1])]], color=shades[i], elinewidth=0.5)
        ax21.axvline(0, color='k', linestyle=":")
        #ax21.set(xlabel="Additional Costs (USD per 1000 briths)")
        ax21.set_xlim([-100,100])
        
    plt.subplots_adjust(top=0.96, bottom=0.08, left=0.050, right=0.959, hspace=0.35, wspace=0.35)
    plt.savefig("MAP Scatterplot_per 1000 births.png", dpi=300)
    
