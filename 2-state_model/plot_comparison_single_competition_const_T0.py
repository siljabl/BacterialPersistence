import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameters
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res_opt = 400
ab_res_comp = 96
n_plot = 3

Tab_max = 24
T0 = 5

lag_opt  = np.zeros([n_plot, ab_res_opt])
del_opt  = np.zeros([n_plot, ab_res_opt])
lag_comp = np.zeros([n_plot, ab_res_comp])
del_comp = np.zeros([n_plot, ab_res_comp])

Tab_opt = np.array([10, 14, 18])
idx_opt  = (Tab_opt * ab_res_opt  / Tab_max).astype(int)
idx_comp = (Tab_opt * ab_res_comp / Tab_max).astype(int)
T0_comp  = idx_comp * Tab_max / ab_res_comp

T = T0 + Tab_opt

for i in range(n_plot):
	
	lag_opt[i]  = np.loadtxt(f'data/high_resolution/optimal_lag-T0{T0}.txt')[:,idx_opt[i]]
	del_opt[i]  = np.loadtxt(f'data/high_resolution/optimal_delta-T0{T0}.txt')[:,idx_opt[i]]
	lag_comp[i] = np.loadtxt(f'data/competition_two_species/competition_lag-T0{T0}.txt')[:,idx_comp[i]]
	del_comp[i] = np.loadtxt(f'data/competition_two_species/competition_delta-T0{T0}.txt')[:,idx_comp[i]]
	
color = ['orange', 'red', 'maroon', 'blue']	       
x_opt  = np.linspace(0, 1, ab_res_opt)
x_comp = np.linspace(0, 1, ab_res_comp)

fig, ax = plt.subplots(1,2, figsize=(6.7,2.5))
for i in range(n_plot):
    ax[0].plot(x_opt, (lag_opt[i] / T[i]), '--', lw=1, c=color[i])
    ax[1].plot(x_opt, del_opt[i], '--', lw=1, c=color[i])

    before = (lag_opt[i] < 1)
    after  = (lag_opt[i] > 0.1)

    ax[0].plot(x_opt[before], (lag_opt[i] / T[i])[before], color=color[i], label=r"$T_{ab}$="+f"{Tab_opt[i]}")
    ax[0].plot(x_opt[after],  (lag_opt[i] / T[i])[after],  color=color[i])
    
    ax[1].plot(x_opt[before], del_opt[i][before], color=color[i])
    ax[1].plot(x_opt[after],  del_opt[i][after],  color=color[i])
    
    # before = (lag_comp[i] < 1)
    # after  = (lag_comp[i] > 0.1)
    
    # ax[0].plot(x_comp[before], (lag_comp[i] / T[i])[before], 'o', color=color[i])
    # ax[0].plot(x_comp[after],  (lag_comp[i] / T[i])[after], 'o', color=color[i])

    # ax[1].plot(x_comp[before], del_comp[i][before], 'o', color=color[i])
    # ax[1].plot(x_comp[after],  del_comp[i][after],  'o', color=color[i])

ax[0].set(xlabel=r"$p$", title=r"$\lambda^{\star} / ~T$")
ax[1].set(xlabel=r"$p$", title=r"$\delta^{\star}$")


fig.tight_layout(rect=[0, 0, 0.85, 1])
fig.legend(loc='upper center', bbox_to_anchor=(0.91, 0.75), ncol=1, fancybox=True, shadow=False, handletextpad=0.05)
fig.savefig(f"figs/compare_competition_T0_{T0}.png")
