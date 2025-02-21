import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

# Plot parameters
mpl.rcParams["font.size"]   = "12"

ab_res_opt  = 101
ab_res_comp = 100
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

T = T0 + Tab_opt

for i in range(n_plot):
	
	lag_opt[i]  = np.loadtxt(f'data/low_resolution/optimal_lag-T0{T0}.txt')[:,idx_opt[i]]
	del_opt[i]  = np.loadtxt(f'data/low_resolution/optimal_delta-T0{T0}.txt')[:,idx_opt[i]]
	lag_comp[i] = np.loadtxt(f'data/competition_two_species/optimal_lag-T0{T0}.txt')[:,idx_comp[i]]
	del_comp[i] = np.loadtxt(f'data/competition_two_species/optimal_delta-T0{T0}.txt')[:,idx_comp[i]]


color = sns.color_palette("flare", as_cmap=True)([0, 0.5, 1])
#color = ['orange', 'red', 'maroon', 'blue']	       
x_opt  = np.linspace(0, 1, ab_res_opt)
x_comp = np.linspace(0, 1, ab_res_comp)

sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(1,2, figsize=(6.7,2.5), sharex=True)
for i in range(n_plot):
    ax[0].plot(x_opt, (lag_opt[i] / T[i]), '--', lw=1.5, c=color[i])
    ax[1].plot(x_opt, del_opt[i],          '--', lw=1.5, c=color[i])

    before = (lag_opt[i] < 1)
    after  = (lag_opt[i] > 0.1)

    # legend
    ax[0].plot(x_opt[before], (lag_opt[i] / T[i])[before], color=color[i], label=f"{Tab_opt[i]}")
    for mask in [before, after]:
        ax[0].plot(x_opt[mask],  (lag_opt[i] / T[i])[mask], lw=2, color=color[i])
        ax[1].plot(x_opt[mask],   del_opt[i][mask],         lw=2, color=color[i])
    
    before = (lag_comp[i] < 1)
    after  = (lag_comp[i] > 0.1)

    for mask in [before, after]:
        ax[0].plot(x_comp[mask], (lag_comp[i] / T[i])[mask], '.', c=color[i], lw=3, alpha=0.5)
        ax[1].plot(x_comp[mask],  del_comp[i][mask],         '.', c=color[i], lw=3, alpha=0.5)

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
           title=r"$T_{ab}$")
fig.savefig(f"figs/compare_competition_T0_{T0}_test.png")

