import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameters
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res_opt = 400
ab_res_comp = 100
n_plot = 3

Tab_max = 24
T0 = 5

lag_opt  = np.zeros([n_plot, ab_res_opt])
del_opt  = np.zeros([n_plot, ab_res_opt])
lag_comp = np.zeros([n_plot, ab_res_comp])
del_comp = np.zeros([n_plot, ab_res_comp])

index = np.array([10, 14, 18])
T = T0 + index

for i in range(n_plot):
	idx_opt  = int(index[i] * ab_res_opt  / Tab_max)
	idx_comp = int(index[i] * ab_res_comp / Tab_max)
	
	lag_opt[i]  = np.loadtxt('data/high_resolution/optimal_lag-T05')[:,idx_opt]
	del_opt[i]  = np.loadtxt('data/high_resolution/optimal_delta-T05')[:,idx_opt]
	lag_comp[i] = np.loadtxt('data/competition_two_species/competition_lag-T05')[:,idx_comp]
	del_comp[i] = np.loadtxt('data/competition_two_species/competition_delta-T05')[:,idx_comp]
	
color = ['orange', 'red', 'maroon', 'blue']	       
x_opt  = np.linspace(0, 1, ab_res_opt)
x_comp = np.linspace(0, 1, ab_res_comp)

fig, ax = plt.subplots(1,2, figsize=(6.7,2.5))

ax[0].set(xlabel=r"$p$", title=r"$\lambda^* / ~T$")
ax[1].set(xlabel=r"$p$", title=r"$\delta^*$")
# ax[1,0].set(xlabel=r"$p$", ylabel=r"$\lambda^*$")
# ax[1,1].set(xlabel=r"$p$", ylabel=r"$\delta^*$")



for i in range(n_plot):
    ax[0].plot(x_opt, (lag_opt[i] / T[i]), '--', lw=1, c=color[i])
    ax[1].plot(x_opt, del_opt[i], '--', lw=1, c=color[i])

    before = (lag_opt[i] < 1)
    after  = (lag_opt[i] > 0.1)

    ax[0].plot(x_opt[before], (lag_opt[i] / T[i])[before], color=color[i], label=r"$T_{AB}$="+f"{index[i]}")
    ax[0].plot(x_opt[after],  (lag_opt[i] / T[i])[after],  color=color[i])
    
    ax[1].plot(x_opt[before], del_opt[i][before], color=color[i])
    ax[1].plot(x_opt[after],  del_opt[i][after],  color=color[i])
    
    before = (lag_comp[i] < 1)
    after  = (lag_comp[i] > 0.1)
    
    #ax[0].plot(x_comp[before], (lag_comp[i] / T[i])[before], 'o', color=color[i])
    # ax[0].plot(x_comp[after], (lag_comp[i] / T[i])[after], 'o', color=color[i])

    # ax[1].plot(x_comp[before], del_comp[i][before], 'o', color=color[i])
    # ax[1].plot(x_comp[after], del_comp[i][after], 'o', color=color[i])



fig.tight_layout(rect=[0, 0, 0.85, 1])
fig.legend(loc='upper center', bbox_to_anchor=(0.91, 0.75), ncol=1, fancybox=True, shadow=False, handletextpad=0.05)
fig.savefig("compare_competition_T0_5.png")
