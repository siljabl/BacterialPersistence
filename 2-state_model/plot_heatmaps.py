import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


ab_res = 400
lag_cmap = mpl.colormaps.get_cmap('viridis')
del_cmap = mpl.colormaps.get_cmap('plasma')

T0  = 5
Tab = 10

# importing data 
folder = 'high_resolution'
lag_Tab = np.loadtxt(f'data/{folder}/optimal_lag-T0{T0}.txt')
del_Tab = np.loadtxt(f'data/{folder}/optimal_delta-T0{T0}.txt')
lag_T0  = np.loadtxt(f'data/{folder}/optimal_lag-Tab{Tab}.txt')
del_T0  = np.loadtxt(f'data/{folder}/optimal_delta-Tab{Tab}.txt')


# arrays for scaling imshow
Tab_arr = np.linspace(0, 24, ab_res)
T0_arr  = np.linspace(0, 12, ab_res)
T_Tab = T0  + np.outer(np.ones(ab_res), Tab_arr)
T_T0  = Tab + np.outer(np.ones(ab_res),  T0_arr)


# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(2, 2, figsize=(6.75, 5.5), sharey=True, sharex='col')

# plotting lag
ax[0,0].set(ylabel=r"$p$", title=r"$\lambda^{\star} / ~T$")
ax[1,0].set(xlabel=r'$T_{0}$', ylabel=r"$p$",  title=r"$\delta^{\star}$")
ax[0,1].set(title=r"$\lambda^{\star} / ~T$")
ax[1,1].set(xlabel=r'$T_{ab}$', title=r"$\delta^{\star}$")

im00 = ax[0,0].imshow(lag_T0 / T_T0,   cmap=lag_cmap, origin="lower", vmin=0, vmax=1.01,   aspect="auto", extent=[0, 12, 0, 1]) 
im10 = ax[1,0].imshow(del_T0,          cmap=del_cmap, origin="lower", vmin=0, vmax=0.07, aspect="auto", extent=[0, 12, 0, 1])
im01 = ax[0,1].imshow(lag_Tab / T_Tab, cmap=lag_cmap, origin="lower", vmin=0, vmax=1.01,   aspect="auto", extent=[0, 24, 0, 1])
im11 = ax[1,1].imshow(del_Tab,         cmap=del_cmap, origin="lower", vmin=0, vmax=0.07, aspect="auto", extent=[0, 24, 0, 1])

# Ticks
for j in range(2):
    ax[j,0].xaxis.set_ticks([0, 6, 12],  minor=False)
    ax[j,0].xaxis.set_ticks([3, 9], minor=True)

    ax[j,1].xaxis.set_ticks([0, 12, 24],  minor=False)
    ax[j,1].xaxis.set_ticks([6, 18], minor=True)

plt.gca().yaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().yaxis.set_ticks([0.25, 0.75], minor=True)


# saving
fig.tight_layout(w_pad=2, h_pad=2)

# Colorbars
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((0, 0))  # adjust the limits as needed

cbar0 = fig.colorbar(im01, ax=ax[0,:], aspect=10, anchor=(-.1, 0.5))
cbar1 = fig.colorbar(im11, ax=ax[1,:], aspect=10, anchor=(-.1, 0.5))

cbar0.ax.yaxis.set_ticks([0,0.5, 1],   minor=False)
cbar0.ax.yaxis.set_ticks([0.25, 0.75], minor=True)
cbar1.formatter = cbformat


fig.savefig(f"figs/2state_heatmap_{folder}.png", dpi=300)

