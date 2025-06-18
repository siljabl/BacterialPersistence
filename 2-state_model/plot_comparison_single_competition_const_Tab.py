import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameters
mpl.rcParams["font.size"]   = "12"

p_res = 96
T_res = 101
n_plot = 3

T0_max = 12
Tab = 10

lag_opt  = np.zeros([n_plot, p_res])
del_opt  = np.zeros([n_plot, p_res])
lag_comp = np.zeros([n_plot, p_res])
del_comp = np.zeros([n_plot, p_res])

T0_opt = np.array([0, 3, 6])
idx_opt  = (T0_opt * T_res / T0_max).astype(int)
idx_comp = (T0_opt * T_res / T0_max).astype(int)

T = Tab + T0_opt

for i in range(n_plot):      
	lag_opt[i]  = np.loadtxt(f'data/low_resolution/optimal_lag-Tab{Tab}.txt')[:p_res,idx_opt[i]]
	del_opt[i]  = np.loadtxt(f'data/low_resolution/optimal_delta-Tab{Tab}.txt')[:p_res,idx_opt[i]]
	lag_comp[i] = np.loadtxt(f'data/competition_two_species/optimal_lag-Tab{Tab}.txt')[:,idx_comp[i]]
	del_comp[i] = np.loadtxt(f'data/competition_two_species/optimal_delta-Tab{Tab}.txt')[:,idx_comp[i]]

color = sns.color_palette("crest", as_cmap=True)([0, 0.5, 1])
#color = ['dodgerblue', 'blue', 'black']	       
x_opt  = np.linspace(0, 1*96/101, p_res)
x_comp = np.linspace(0, 1*96/101, p_res)

sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(1,2, figsize=(6.7,2.5), sharex=True)
for i in range(n_plot):
    ax[0].plot(x_opt, (lag_opt[i] / T[i]), '--', lw=1.5, c=color[i])
    ax[1].plot(x_opt, del_opt[i],          '--', lw=1.5, c=color[i])

    before = (lag_opt[i] < 1)
    middle = (lag_opt[i] > 0.1)*(del_opt[i] > 0.01)
    after  = (lag_opt[i] > 0.1)*(del_opt[i] < 0.01)

    for mask in [before, middle, after]:
        ax[0].plot(x_opt[mask], (lag_opt[i] / T[i])[mask], lw=2, c=color[i])
        ax[1].plot(x_opt[mask], del_opt[i][mask],          lw=2, color=color[i])

    # label  
    ax[1].plot(x_opt[before], del_opt[i][before], color=color[i], label=f"{T0_opt[i]}")
    
    before = (lag_comp[i] < 1)
    middle = (lag_comp[i] > 0.1)*(del_comp[i] > 0.01)
    after  = (lag_comp[i] > 0.1)*(del_comp[i] < 0.01)
    
    for mask in [before, middle, after]:
        ax[0].plot(x_comp[mask], (lag_comp[i] / T[i])[mask], '.', c=color[i], lw=3, alpha=0.5)
        ax[1].plot(x_comp[mask], del_comp[i][mask],          '.', c=color[i], lw=3, alpha=0.5)


ax[0].set(xlabel=r"$p$", title=r"$\lambda^{\star} / ~T$")
ax[1].set(xlabel=r"$p$", title=r"$\delta^{\star}$")

# Ticks
plt.gca().xaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().xaxis.set_ticks([0.25, 0.75], minor=True)

ax[0].yaxis.set_ticks([0, 0.5, 1],  minor=False)
ax[0].yaxis.set_ticks([0.25, 0.75], minor=True)

ax[1].yaxis.set_ticks([0, 0.03, 0.06], minor=False)
ax[1].yaxis.set_ticks([0.015, 0.045],  minor=True)


sns.despine()
fig.tight_layout(rect=[0, 0, 0.85, 1], w_pad=2)
fig.legend(loc='upper center',
           bbox_to_anchor=(0.91, 0.9),
           ncol=1,
           frameon=False,
           handlelength=1, 
           title=r"$T_0$")
fig.savefig(f"figs/compare_competition_Tab_{Tab}.png", dpi=300)

