# Supplement: CTC
import numpy as np

ctc={}
for scen in range(3,6):
    if scen==3:
        ctc["ctc_baseline"]={}
        model_inputs=vac_weight(data_prepped,runs)
        model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
        model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
        ctc["ctc_baseline"]["fac_cov"]=model_inputs["fac_cov"]
        ctc["ctc_baseline"]["com_cov"]=model_inputs["com_cov"]
        ctc["ctc_baseline"]["init_conds"]=model_inputs["init_state"]
        results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
        outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
        ctc["ctc_baseline"]["cov_cost"]=outcomes["cov_cost"]
        ctc["ctc_baseline"]["ce_out"]=outcomes["ce_out"]
        ctc["ctc_baseline"]["DALYs"]=outcomes["DALYs"]
        ctc["ctc_baseline"]["dis_cost"]=outcomes["dis_cost"]
        ctc["ctc_baseline"]["MAP_PTA"]=outcomes["MAP_PTA"]
        ctc["ctc_baseline"]["Table"]=np.zeros((len(settings),43))

       #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
        for loc in range(len(settings)):
            for run in range(runs):
                for i in range(3):
                    ctc["ctc_baseline"]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[-1,8])/10e2#births
                    ctc["ctc_baseline"]["Table"][loc,i+1]=np.percentile(ctc["ctc_baseline"]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                    ctc["ctc_baseline"]["Table"][loc,i+4]=np.percentile( ctc["ctc_baseline"]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                    ctc["ctc_baseline"]["Table"][loc,i+7]=np.percentile( ctc["ctc_baseline"]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                    ctc["ctc_baseline"]["Table"][loc,i+10]=np.percentile( ctc["ctc_baseline"]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                    ctc["ctc_baseline"]["Table"][loc,i+13]=np.percentile( ctc["ctc_baseline"]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                    ctc["ctc_baseline"]["Table"][loc,i+16]=np.percentile( ctc["ctc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                    ctc["ctc_baseline"]["Table"][loc,i+19]=np.percentile( ctc["ctc_baseline"]["DALYs"][loc,0:runs,2]-ctc["ctc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                    ctc["ctc_baseline"]["Table"][loc,i+22]=np.percentile( ctc["ctc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                    ctc["ctc_baseline"]["Table"][loc,i+25]=np.percentile( ctc["ctc_baseline"]["dis_cost"][loc,0:runs,0]-ctc["ctc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Cost
                    ctc["ctc_baseline"]["Table"][loc,i+28]=np.percentile( ctc["ctc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                    ctc["ctc_baseline"]["Table"][loc,i+31]=np.percentile( ctc["ctc_baseline"]["cov_cost"][loc,0:runs,15]-ctc["ctc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
  
    if scen==4:
        for cov in cov_val:
            
            #Additional Coverage
            model_assumptions["f_map_adc"]=0    #Additional coverage limited to births outside facility settings
            model_assumptions["c_map_adc"]=cov
            
            ctc["s1_c%.1d" % (cov*100)]={}
            model_inputs=vac_weight(data_prepped,runs)
            model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
            model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
            ctc["s1_c%.1d" % (cov*100)]["fac_cov"]=model_inputs["fac_cov"]
            ctc["s1_c%.1d" % (cov*100)]["com_cov"]=model_inputs["com_cov"]
            ctc["s1_c%.1d" % (cov*100)]["init_conds"]=model_inputs["init_state"]
            results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
            outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
            ctc["s1_c%.1d" % (cov*100)]["cov_cost"]=outcomes["cov_cost"]
            ctc["s1_c%.1d" % (cov*100)]["ce_out"]=outcomes["ce_out"]
            ctc["s1_c%.1d" % (cov*100)]["DALYs"]=outcomes["DALYs"]
            ctc["s1_c%.1d" % (cov*100)]["dis_cost"]=outcomes["dis_cost"]
            ctc["s1_c%.1d" % (cov*100)]["MAP_PTA"]=outcomes["MAP_PTA"]
            
            ctc["s1_c%.1d" % (cov*100)]["ICER"]=np.zeros((len(settings*runs), len(map_price), 6))
            ctc["s1_c%.1d" % (cov*100)]["Table"]=np.zeros((len(settings),43))   
               
            #ICERs
            for loc in range(len(settings)):
                for run in range(runs):
                    for c in range(len(map_price)):
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,0]=ctc["s1_c%.1d" % (cov*100)]["MAP_PTA"][(loc*runs)+run,c,4]
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,1]=ctc["ctc_baseline"]["MAP_PTA"][(loc*runs)+run,c,4]
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,2]=ctc["s1_c%.1d" % (cov*100)]["MAP_PTA"][(loc*runs)+run,c,0]
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,3]=ctc["ctc_baseline"]["MAP_PTA"][(loc*runs)+run,c,0]
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,4]=map_price[c]
                        ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,5]=-((ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,0]- ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,1])/(ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,2]-ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs)+run,c,3]))
            
            #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
            for loc in range(len(settings)):
                for run in range(runs):
                    for i in range(3):
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[-1,8])/10e2#Births('000)
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+1]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+4]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+7]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+10]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+13]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+16]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+19]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["DALYs"][loc,0:runs,2]-ctc["ctc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+22]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+25]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["dis_cost"][loc,0:runs,0]-ctc["ctc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Costs
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+28]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+31]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["cov_cost"][loc,0:runs,15]-ctc["ctc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+34]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,cpad,5],25+(25*i))#ICER CPAD
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+37]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,tcpad,5],25+(25*i))#ICER 2xCPAD
                       ctc["s1_c%.1d" % (cov*100)]["Table"][loc,i+40]=np.percentile(ctc["s1_c%.1d" % (cov*100)]["ICER"][(loc*runs):(loc*runs)+runs,five,5],25+(25*i))#ICER US$5
    if scen==5:
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
                
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]={}
                model_inputs=vac_weight(data_prepped,runs)
                model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
                model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["fac_cov"]=model_inputs["fac_cov"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["com_cov"]=model_inputs["com_cov"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["init_conds"]=model_inputs["init_state"]
                results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
                outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"]=outcomes["cov_cost"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"]=outcomes["ce_out"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"]=outcomes["DALYs"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"]=outcomes["dis_cost"]
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"]=outcomes["MAP_PTA"]
                
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"]=np.zeros((len(settings*runs), len(map_price), 6))
                ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"]=np.zeros((len(settings),43))
                
                #ICERs
                for loc in range(len(settings)):
                    for run in range(runs):
                        for c in range(len(map_price)):
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]=ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,4]
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1]=ctc["ctc_baseline"]["MAP_PTA"][(loc*runs)+run,c,4]
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]=ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,0]
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]=ctc["ctc_baseline"]["MAP_PTA"][(loc*runs)+run,c,0]
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,4]=map_price[c]
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,5]=-((ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]- ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1])/(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]-ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]))
                 
                #Collate Outputs for export and tabulation (with IQRs); as baseline no ICER needed
                for loc in range(len(settings)):
                    for run in range(runs):
                        for i in range(3):
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,0]=(raw_inputs["set_pars"].iloc[-1,8])/10e2#Births('000)
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+1]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,14],25+(25*i))#MAP Doses
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+4]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"][loc,0:runs,0],25+(25*i))#CHB
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+7]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"][loc,0:runs,1],25+(25*i))#Deaths
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+10]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,0],25+(25*i))#YLL
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+13]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,1],25+(25*i))#YLD
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+16]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,2],25+(25*i))#DALYs
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+19]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"][loc,0:runs,2]-ctc["ctc_baseline"]["DALYs"][loc,0:runs,2],25+(25*i))#Diff DALYs
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+22]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"][loc,0:runs,0],25+(25*i))#Dis Man Cost
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+25]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"][loc,0:runs,0]-ctc["ctc_baseline"]["dis_cost"][loc,0:runs,0],25+(25*i))#Diff Dis Man Costs
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+28]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,15],25+(25*i))#Vax Cost
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+31]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"][loc,0:runs,15]-ctc["ctc_baseline"]["cov_cost"][loc,0:runs,15],25+(25*i))#Diff Vax Cost
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+34]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,cpad,5],25+(25*i))#ICER CPAD
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+37]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,tcpad,5],25+(25*i))#ICER 2xCPAD
                            ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"][loc,i+40]=np.percentile(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs):(loc*runs)+runs,five,5],25+(25*i))#ICER US$5            
# Outputs:
ctc_s1_plot=np.zeros((1,len(map_price),4))
ctc_s2_1_plot=np.zeros((1,len(map_price),4))
ctc_s2_5_plot=np.zeros((1,len(map_price),4))
ctc_s2_10_plot=np.zeros((1,len(map_price),4))

for idx,reg in enumerate(settings):
    for run in range(runs):
            for c in range(len(map_price)):
                for i in range(3):
                    ctc_s1_plot[idx,c,0]=map_price[c]
                    ctc_s1_plot[idx,c,i+1]=np.percentile(ctc["s1_c1"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                    
                    ctc_s2_1_plot[idx,c,0]=map_price[c]
                    ctc_s2_1_plot[idx,c,i+1]=np.percentile(ctc["s1_c1_r1"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                    
                    ctc_s2_5_plot[idx,c,0]=map_price[c]
                    ctc_s2_5_plot[idx,c,i+1]=np.percentile(ctc["s1_c1_r5"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                     
                    ctc_s2_10_plot[idx,c,0]=map_price[c]
                    ctc_s2_10_plot[idx,c,i+1]=np.percentile(ctc["s1_c1_r10"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                    
# Plot
s2_fig=plt.figure(figsize=(20, 10))

ax0=s2_fig.add_subplot(121)
ax0.plot(map_price_unit[:], main_s1_plot[-1,:,2], color='red', marker="^", alpha=0.5)
ax0.fill_between(map_price_unit[:],main_s1_plot[-1,:,1], main_s1_plot[-1,:,3], alpha=0.3, color="red" )
ax0.plot(map_price_unit[:], main_s2_1_plot[-1,:,2], color='blue', marker="s", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_1_plot[-1,:,1], main_s2_1_plot[-1,:,3], alpha=0.3, color="blue" )
ax0.plot(map_price_unit[:], main_s2_5_plot[-1,:,2],  color='green', marker="o", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_5_plot[-1,:,1], main_s2_5_plot[-1,:,3], alpha=0.3, color="green" )
ax0.plot(map_price_unit[:], main_s2_10_plot[-1,:,2], color='orange',marker="h", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_10_plot[-1,:,1], main_s2_10_plot[-1,:,3], alpha=0.3, color="orange" )
plt.ylim(-170, 250)
plt.xlim(0,5)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.xlabel("MAP Price ($US)")
plt.ylabel("ICER (USD/DALY averted)")
plt.title("All LMICs - Main Analysis")

ax1=s2_fig.add_subplot(122)
ax1.plot(map_price_unit[:], ctc_s1_plot[0,:,2], color='red', marker="^", alpha=0.5)
ax1.fill_between(map_price_unit[:],ctc_s1_plot[0,:,1], ctc_s1_plot[0,:,3], alpha=0.3, color="red" )
ax1.plot(map_price_unit[:], ctc_s2_1_plot[0,:,2], color='blue', marker="s", alpha=0.5)
ax1.fill_between(map_price_unit[:], ctc_s2_1_plot[0,:,1], ctc_s2_1_plot[0,:,3], alpha=0.3, color="blue" )
ax1.plot(map_price_unit[:], ctc_s2_5_plot[0,:,2],  color='green', marker="o", alpha=0.5)
ax1.fill_between(map_price_unit[:], ctc_s2_5_plot[0,:,1], ctc_s2_5_plot[0,:,3], alpha=0.3, color="green" )
ax1.plot(map_price_unit[:], ctc_s2_10_plot[0,:,2], color='orange',marker="h", alpha=0.5)
ax1.fill_between(map_price_unit[:], ctc_s2_10_plot[0,:,1], ctc_s2_10_plot[0,:,3], alpha=0.3, color="orange" )
plt.ylim(-170, 250)
plt.xlim(0,5)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.xlabel("MAP Price ($US)")
plt.title("All LMICs - CTC Baseline")
ax1.legend(("Scenario 1: Additional MAPs (1%)","Scenario 1 + 1% Replacement MAPs","Scenario 1 + 5% Replacement MAPs","Scenario 1 + 10% Replacement MAPs"),loc="best")
plt.savefig("Supp CTC "+str(runs)+" runs.png", dpi=300)     

# Table 
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

ctc_excel={}
for scen in range(3,6):
    if scen==3:
        ctc_excel["baseline"]=pd.DataFrame(ctc["ctc_baseline"]["Table"])
        ctc_excel["baseline"].insert(0,"setting",settings)
        ctc_excel["baseline"]=ctc_excel["baseline"].rename(columns=dict(zip(ctc_excel["baseline"].iloc[:,1:].columns, tab_labs)))
    if scen==4:
        for cov in cov_val:
            ctc_excel["s1_c%.1d" % (cov*100)]=pd.DataFrame(ctc["s1_c%.1d" % (cov*100)]["Table"])
            ctc_excel["s1_c%.1d" % (cov*100)].insert(0,"setting", settings)
            ctc_excel["s1_c%.1d" % (cov*100)]=ctc_excel["s1_c%.1d" % (cov*100)].rename(columns=dict(zip(ctc_excel["s1_c%.1d" % (cov*100)].iloc[:,1:].columns, tab_labs)))
    if scen==5:
        for cov in cov_val:
            for rep in rep_val:
                ctc_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)]=pd.DataFrame(ctc["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["Table"])
                ctc_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].insert(0,"setting", settings)
                ctc_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)]=ctc_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].rename(columns=dict(zip(ctc_excel["s1_c%.1d_r%.1d" % (cov*100, rep*100)].iloc[:,1:].columns, tab_labs)))

writer=pd.ExcelWriter(str("CTC "+ level+" IQR"+ str(runs) + " runs impute "+str(impute_bd)+".xlsx"))
for df_name, df in ctc_excel.items():
    df.to_excel(writer, sheet_name=df_name)
writer.save()
                