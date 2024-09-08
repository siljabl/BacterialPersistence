import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameters
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res_opt = 400
ab_res_comp = 100
n_plot = 3

T0_max = 12
Tab = 10

lag_opt  = np.zeros([n_plot, ab_res_opt])
del_opt  = np.zeros([n_plot, ab_res_opt])
lag_comp = np.zeros([n_plot, ab_res_comp])
del_comp = np.zeros([n_plot, ab_res_comp])

T0_opt = np.array([0, 3.7, 7])
idx_opt  = (T0_opt * ab_res_opt  / T0_max).astype(int)
idx_comp = (T0_opt * ab_res_comp / T0_max).astype(int) + 6
T0_comp  = idx_comp * T0_max / ab_res_comp

T = Tab + T0_opt

for i in range(n_plot):
	print(idx_opt,  np.linspace(0,12,ab_res_opt)[idx_opt])
	print(idx_comp, np.linspace(0,12,ab_res_comp)[idx_comp])
      
	lag_opt[i]  = np.loadtxt('data/high_resolution/optimal_lag-Tab10')[:,idx_opt[i]]
	del_opt[i]  = np.loadtxt('data/high_resolution/optimal_delta-Tab10')[:,idx_opt[i]]
	lag_comp[i] = np.loadtxt('data/competition_two_species/competition_lag-Tab10')[:,idx_comp[i]]
	del_comp[i] = np.loadtxt('data/competition_two_species/competition_delta-Tab10')[:,idx_comp[i]]
	
color = color = ['dodgerblue', 'blue', 'black']	       
x_opt  = np.linspace(0, 1, ab_res_opt)
x_comp = np.linspace(0, 1, ab_res_comp)

fig, ax = plt.subplots(1,2, figsize=(6.7,2.5))
for i in range(n_plot):
    ax[0].plot(x_opt, (lag_opt[i] / T[i]), '--', lw=1, c=color[i])
    ax[1].plot(x_opt, del_opt[i], '--', lw=1, c=color[i])

    before = (lag_opt[i] < 1)
    middle = (lag_opt[i] > 0.1)*(del_opt[i] > 0.01)
    after  = (lag_opt[i] > 0.1)*(del_opt[i] < 0.01)

    ax[0].plot(x_opt[before], (lag_opt[i] / T[i])[before], c=color[i])
    ax[0].plot(x_opt[middle], (lag_opt[i] / T[i])[middle], c=color[i])
    ax[0].plot(x_opt[after],  (lag_opt[i] / T[i])[after],  c=color[i])
    
    ax[1].plot(x_opt[before], del_opt[i][before], color=color[i], label=r"$T_{0}$="+f"{T0_opt[i]}")
    ax[1].plot(x_opt[middle], del_opt[i][middle], color=color[i])
    ax[1].plot(x_opt[after],  del_opt[i][after],  color=color[i])
    
    before = (lag_comp[i] < 1)
    middle = (lag_comp[i] > 0.1)*(del_comp[i] > 0.01)
    after  = (lag_comp[i] > 0.1)*(del_comp[i] < 0.01)
    
    ax[0].plot(x_comp[before], (lag_comp[i] / T[i])[before], '.', c=color[i], lw=3, alpha=0.5)
    ax[0].plot(x_comp[middle], (lag_comp[i] / T[i])[middle], '.', c=color[i], lw=3, alpha=0.5)
    ax[0].plot(x_comp[after],  (lag_comp[i] / T[i])[after],  '.', c=color[i], lw=3, alpha=0.5)

    ax[1].plot(x_comp[before], del_comp[i][before], '.', c=color[i], alpha=0.5)#, label=r"$T_{0}$="+f"{np.round(T0_comp[i],1)}")
    ax[1].plot(x_comp[middle], del_comp[i][middle], '.', c=color[i], alpha=0.5)
    ax[1].plot(x_comp[after],  del_comp[i][after],  '.', c=color[i], alpha=0.5)


ax[0].set(xlabel=r"$p$", title=r"$\lambda^* / ~T$")
ax[1].set(xlabel=r"$p$", title=r"$\delta^*$")

fig.tight_layout(rect=[0, 0, 0.85, 1])
fig.legend(loc='upper center', bbox_to_anchor=(0.91, 0.75), ncol=1, fancybox=True, shadow=False, handletextpad=0)
fig.savefig("compare_competition_Tab_10.png")

