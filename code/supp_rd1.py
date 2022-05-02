# Supplement: Replacement on day 1
d_one={}
#Remake the array so coverage is all day 1
model_assumptions["f_map_rep"]=np.array([0,0,0,0,0.01,0,0,0,0.01,0,0,0,0.01,0,0,0]).reshape(4,4)
model_assumptions["c_map_rep"]=np.array([0,0,0,0,0.01,0,0,0,0.01,0,0,0,0.01,0,0,0]).reshape(4,4)

for scen in range(2,3):
    if scen==2:
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
                
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]={}
                model_inputs=vac_weight(data_prepped,runs)
                model_inputs=vax_dist(model_assumptions, data_prepped, settings, model_inputs, runs, scen)
                model_inputs=model_init(data_prepped, model_inputs, vac_eff, settings, runs)
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["fac_cov"]=model_inputs["fac_cov"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["com_cov"]=model_inputs["com_cov"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["init_conds"]=model_inputs["init_state"]
                results=dis_model(data_prepped, model_inputs, settings, t_steps, dt, runs)
                outcomes=outcome_calcs(data_prepped, model_inputs, results,map_price, t_steps, dt, runs, settings)
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["cov_cost"]=outcomes["cov_cost"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ce_out"]=outcomes["ce_out"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["DALYs"]=outcomes["DALYs"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["dis_cost"]=outcomes["dis_cost"]
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"]=outcomes["MAP_PTA"] 
                
                d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"]=np.zeros((len(settings*runs), len(map_price), 6))
                
             #ICERs (update this to factor in only running for all LMICs)
                for loc in range(len(settings)):
                    for run in range(runs):
                        for c in range(len(map_price)):
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]=d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,4]
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1]=main["cc_baseline"]["MAP_PTA"][(g_ref*runs)+run,c,4]
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]=d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["MAP_PTA"][(loc*runs)+run,c,0]
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]=main["cc_baseline"]["MAP_PTA"][(g_ref*runs)+run,c,0]
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,4]=map_price[c]
                            d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,5]=-((d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,0]-d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,1])/(d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,2]-d_one["s1_c%.1d_r%.1d" % (cov*100, rep*100)]["ICER"][(loc*runs)+run,c,3]))
                

#Remake the array so coverage back to baseline assumptions
model_assumptions["f_map_rep"]=np.array([0,0,0,0,0.01,0,0,0,0,0.01,0,0,0,0,0.01,0]).reshape(4,4)
model_assumptions["c_map_rep"]=np.array([0,0,0,0,0.01,0,0,0,0,0.01,0,0,0,0,0.01,0]).reshape(4,4)

# #Outputs:
rd1_s2_1_plot=np.zeros((1,len(map_price),4))
rd1_s2_5_plot=np.zeros((1,len(map_price),4))
rd1_s2_10_plot=np.zeros((1,len(map_price),4))

for idx,reg in enumerate(settings):
    for run in range(runs):
            for c in range(len(map_price)):
                for i in range(3):
                    rd1_s2_1_plot[idx,c,0]=map_price[c]
                    rd1_s2_1_plot[idx,c,i+1]=np.percentile(d_one["s1_c1_r1"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                    
                    rd1_s2_5_plot[idx,c,0]=map_price[c]
                    rd1_s2_5_plot[idx,c,i+1]=np.percentile(d_one["s1_c1_r5"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                     
                    rd1_s2_10_plot[idx,c,0]=map_price[c]
                    rd1_s2_10_plot[idx,c,i+1]=np.percentile(d_one["s1_c1_r10"]["ICER"][(idx*runs):(idx*runs)+runs,c,5], (i*25)+25)
                    
# Plot
s2_fig=plt.figure(figsize=(20, 10))

ax0=s2_fig.add_subplot(121)
ax0.plot(map_price_unit[:], main_s2_1_plot[-1,:,2], color='blue', marker="s", alpha=0.5)
ax0.plot(map_price_unit[:], main_s2_5_plot[-1,:,2],  color='green', marker="o", alpha=0.5)
ax0.plot(map_price_unit[:], main_s2_10_plot[-1,:,2], color='orange',marker="h", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_1_plot[-1,:,1], main_s2_1_plot[-1,:,3], alpha=0.3, color="blue" )
ax0.fill_between(map_price_unit[:], main_s2_5_plot[-1,:,1], main_s2_5_plot[-1,:,3], alpha=0.3, color="green" )
ax0.fill_between(map_price_unit[:], main_s2_10_plot[-1,:,1], main_s2_10_plot[-1,:,3], alpha=0.3, color="orange" )
plt.ylim(-420, 500)
plt.xlim(0,5)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.xlabel("MAP Price ($US)")
plt.ylabel("ICER (USD/DALY averted)")
plt.title("All LMICs - Main Analysis")

ax1=s2_fig.add_subplot(122)
ax1.plot(map_price_unit[:], rd1_s2_1_plot[0,:,2], color='blue', marker="s", alpha=0.5)
ax1.plot(map_price_unit[:], rd1_s2_5_plot[0,:,2],  color='green', marker="o", alpha=0.5)
ax1.plot(map_price_unit[:], rd1_s2_10_plot[0,:,2], color='orange',marker="h", alpha=0.5)
ax1.fill_between(map_price_unit[:], rd1_s2_1_plot[0,:,1], rd1_s2_1_plot[0,:,3], alpha=0.3, color="blue" )
ax1.fill_between(map_price_unit[:], rd1_s2_5_plot[0,:,1], rd1_s2_5_plot[0,:,3], alpha=0.3, color="green" )
ax1.fill_between(map_price_unit[:], rd1_s2_10_plot[0,:,1], rd1_s2_10_plot[0,:,3], alpha=0.3, color="orange" )
plt.ylim(-420, 500)
plt.xlim(0,5)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.xlabel("MAP Price ($US)")
plt.title("All LMICs - Replacement Coverage All Day One")
ax1.legend(("Scenario 1 +1% Replacement MAPs", "Scenario 1+ 5% Replacement MAPs", "Scenario 1+ 10% Replacement MAPs"),loc="upper left")
plt.savefig("Supp Timing "+str(runs)+" runs.png", dpi=300)                     