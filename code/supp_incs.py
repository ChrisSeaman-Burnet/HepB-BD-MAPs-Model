# Supplemental Analysis output
main_s5_plot=np.zeros((1,len(map_price),4))
main_s5_r1_plot=np.zeros((1,len(map_price),4))
main_s5_r5_plot=np.zeros((1,len(map_price),4))
main_s5_r10_plot=np.zeros((1,len(map_price),4))

main_s10_plot=np.zeros((1,len(map_price),4))
main_s10_r1_plot=np.zeros((1,len(map_price),4))
main_s10_r5_plot=np.zeros((1,len(map_price),4))
main_s10_r10_plot=np.zeros((1,len(map_price),4))

#who_reg=[who_reg[-1]]

for idx,reg in enumerate(settings):
    for run in range(runs):
            for c in range(len(map_price)):
                for i in range(3):
                    main_s5_plot[idx,c,0]=map_price[c]
                    main_s5_plot[idx,c,i+1]=np.percentile(main["s1_c5"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                    
                    main_s5_r1_plot[idx,c,0]=map_price[c]
                    main_s5_r1_plot[idx,c,i+1]=np.percentile(main["s1_c5_r1"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                    
                    main_s5_r5_plot[idx,c,0]=map_price[c]
                    main_s5_r5_plot[idx,c,i+1]=np.percentile(main["s1_c5_r5"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                     
                    main_s5_r10_plot[idx,c,0]=map_price[c]
                    main_s5_r10_plot[idx,c,i+1]=np.percentile(main["s1_c5_r10"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                    
                    main_s10_plot[idx,c,0]=map_price[c]
                    main_s10_plot[idx,c,i+1]=np.percentile(main["s1_c10"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                     
                    main_s10_r1_plot[idx,c,0]=map_price[c]
                    main_s10_r1_plot[idx,c,i+1]=np.percentile(main["s1_c10_r1"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                     
                    main_s10_r5_plot[idx,c,0]=map_price[c]
                    main_s10_r5_plot[idx,c,i+1]=np.percentile(main["s1_c10_r5"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)
                     
                    main_s10_r10_plot[idx,c,0]=map_price[c]
                    main_s10_r10_plot[idx,c,i+1]=np.percentile(main["s1_c10_r10"]["ICER"][(g_ref*runs):(g_ref*runs)+runs,c,5], (i*25)+25)

#Produce Plot
s1_fig=plt.figure(figsize=(20,10))

#1% Additional Coverage
ax0=s1_fig.add_subplot(221)
ax0.plot(map_price_unit[:], main_s1_plot[-1,:,2], color="red", marker="^", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s1_plot[-1,:,1], main_s1_plot[-1,:,3], color="red", alpha=0.3)
ax0.plot(map_price_unit[:], main_s2_1_plot[-1,:,2], color="blue", marker="^", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_1_plot[-1,:,1], main_s2_1_plot[-1,:,3], color="blue", alpha=0.3)
ax0.plot(map_price_unit[:], main_s2_5_plot[-1,:,2], color="green", marker="^", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_5_plot[-1,:,1], main_s2_5_plot[-1,:,3], color="green", alpha=0.3)
ax0.plot(map_price_unit[:], main_s2_10_plot[-1,:,2], color="orange", marker="^", alpha=0.5)
ax0.fill_between(map_price_unit[:], main_s2_10_plot[-1,:,1], main_s2_10_plot[-1,:,3], color="orange", alpha=0.3)
plt.xlim(0,5)
plt.ylim(-170, 250)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.ylabel("ICER (USD/DALY averted)")
plt.title("All LMICs - 1% Additional Coverage")

#5% Additional Coverage
ax1=s1_fig.add_subplot(222)
ax1.plot(map_price_unit[:], main_s5_plot[0,:,2], color="red", marker="^", alpha=0.5)
ax1.fill_between(map_price_unit[:], main_s5_plot[0,:,1], main_s5_plot[0,:,3], color="red", alpha=0.3)
ax1.plot(map_price_unit[:], main_s5_r1_plot[0,:,2], color="blue", marker="^", alpha=0.5)
ax1.fill_between(map_price_unit[:], main_s5_r1_plot[0,:,1], main_s5_r1_plot[0,:,3], color="blue", alpha=0.3)
ax1.plot(map_price_unit[:], main_s5_r5_plot[0,:,2], color="green", marker="^", alpha=0.5)
ax1.fill_between(map_price_unit[:], main_s5_r5_plot[0,:,1], main_s5_r5_plot[0,:,3], color="green", alpha=0.3)
ax1.plot(map_price_unit[:], main_s5_r10_plot[0,:,2], color="orange", marker="^", alpha=0.5)
ax1.fill_between(map_price_unit[:], main_s5_r10_plot[0,:,1], main_s5_r10_plot[0,:,3], color="orange", alpha=0.3)
plt.xlim(0,5)
plt.ylim(-170, 250)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.xlabel("MAP Price ($US)")
plt.title("All LMICs - 5% Additional Coverage")
    
ax2=s1_fig.add_subplot(223)
ax2.plot(map_price_unit[:], main_s10_plot[0,:,2], color="red", marker="^", alpha=0.5)
ax2.fill_between(map_price_unit[:], main_s10_plot[0,:,1], main_s10_plot[0,:,3], color="red", alpha=0.3)
ax2.plot(map_price_unit[:], main_s10_r1_plot[0,:,2], color="blue", marker="^", alpha=0.5)
ax2.fill_between(map_price_unit[:], main_s10_r1_plot[0,:,1], main_s10_r1_plot[0,:,3], color="blue", alpha=0.3)
ax2.plot(map_price_unit[:], main_s10_r5_plot[0,:,2], color="green", marker="^", alpha=0.5)
ax2.fill_between(map_price_unit[:], main_s10_r5_plot[0,:,1], main_s10_r5_plot[0,:,3], color="green", alpha=0.3)
ax2.plot(map_price_unit[:], main_s10_r10_plot[0,:,2], color="orange", marker="^", alpha=0.5)
ax2.fill_between(map_price_unit[:], main_s10_r10_plot[0,:,1], main_s10_r10_plot[0,:,3], color="orange", alpha=0.3)
plt.xlim(0,5)
plt.ylim(-170, 250)
plt.hlines(y=0, xmin=0, xmax=10, color="black", linestyle="--")
plt.grid()
plt.ylabel("ICER (USD/DALY averted)")
plt.xlabel("MAP Price ($US)")
plt.title("All LMICs - 10% Additional Coverage")
ax2.legend(("Additional Coverage Only", "1% Replacement", "5% Replacement", "10% Replacement"),loc="center left", bbox_to_anchor=(1.1,0.5), shadow=1, fancybox=1, fontsize=12)
plt.savefig("Supp Increments "+str(runs)+" runs.png", dpi=300) 