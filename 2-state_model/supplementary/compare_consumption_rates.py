import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


ab_res = 101
lag_cmap = mpl.colormaps.get_cmap('viridis')
del_cmap = mpl.colormaps.get_cmap('plasma')

T0  = 5
Tab = 10

# importing data 
dir_1 = 'low_resolution'
dir_2 = 'correct_consumption_rate'
lag_Tab_paper = np.loadtxt(f'../data/{dir_1}/optimal_lag-T0{T0}.txt')
del_Tab_paper = np.loadtxt(f'../data/{dir_1}/optimal_delta-T0{T0}.txt')
lag_T0_paper  = np.loadtxt(f'../data/{dir_1}/optimal_lag-Tab{Tab}.txt')
del_T0_paper  = np.loadtxt(f'../data/{dir_1}/optimal_delta-Tab{Tab}.txt')

lag_Tab_correct = np.loadtxt(f'../data/{dir_2}/optimal_lag-T0{T0}.txt')
del_Tab_correct = np.loadtxt(f'../data/{dir_2}/optimal_delta-T0{T0}.txt')
lag_T0_correct  = np.loadtxt(f'../data/{dir_2}/optimal_lag-Tab{Tab}.txt')
del_T0_correct  = np.loadtxt(f'../data/{dir_2}/optimal_delta-Tab{Tab}.txt')

lag_Tab = (lag_Tab_paper - lag_Tab_correct) / lag_Tab_paper
del_Tab = (del_Tab_paper - del_Tab_correct) / (del_Tab_paper+1e-2)
lag_T0  = (lag_T0_paper  - lag_T0_correct)  / lag_T0_paper
del_T0  = (del_T0_paper  - del_T0_correct)  / (del_T0_paper+1e-2)


# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(2, 2, figsize=(6.75, 5.5), sharey=True, sharex='col')

# plotting lag
ax[0,0].set(title=r"$(\lambda^{\star}_{eq(1)} - \lambda^{\star}_{eq(2)}) ~/~ \lambda^{\star}_{eq(1)}$", ylabel=r"$p$")
ax[1,0].set(title=r"$(\delta^{\star}_{eq(1)}  - \delta^{\star}_{eq(2)})  ~/~ \delta^{\star}_{eq(1)}$",  xlabel=r'$T_{0}$', ylabel=r"$p$")
ax[0,1].set(title=r"$(\lambda^{\star}_{eq(1)} - \lambda^{\star}_{eq(2)}) ~/~ \lambda^{\star}_{eq(1)}$")
ax[1,1].set(title=r"$(\delta^{\star}_{eq(1)}  - \delta^{\star}_{eq(2)})  ~/~ \delta^{\star}_{eq(1)}$",  xlabel=r'$T_{ab}$')

im00 = ax[0,0].imshow(lag_T0,  cmap=lag_cmap, vmin=-1, vmax=1, origin="lower", aspect="auto", extent=[0, 12, 0, 1]) 
im10 = ax[1,0].imshow(del_T0,  cmap=del_cmap, vmin=-1, vmax=1, origin="lower", aspect="auto", extent=[0, 12, 0, 1])
im01 = ax[0,1].imshow(lag_Tab, cmap=lag_cmap, vmin=-1, vmax=1, origin="lower", aspect="auto", extent=[0, 24, 0, 1])
im11 = ax[1,1].imshow(del_Tab, cmap=del_cmap, vmin=-1, vmax=1, origin="lower", aspect="auto", extent=[0, 24, 0, 1])

# Ticks
for j in range(2):
    ax[j,0].xaxis.set_ticks([0, 6, 12],  minor=False)
    ax[j,0].xaxis.set_ticks([3, 9], minor=True)

    ax[j,1].xaxis.set_ticks([0, 12, 24],  minor=False)
    ax[j,1].xaxis.set_ticks([6, 18], minor=True)

plt.gca().yaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().yaxis.set_ticks([0.25, 0.75], minor=True)


# saving
#fig.suptitle(r"$|x_{eq (1)} - x_{eq (2)}| / x_{eq. (1)}$")
fig.tight_layout(w_pad=2, h_pad=2)

# Colorbars
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((0, 0))  # adjust the limits as needed

cbar0 = fig.colorbar(im01, ax=ax[0,:], aspect=10, anchor=(-.1, 0.5))
cbar1 = fig.colorbar(im11, ax=ax[1,:], aspect=10, anchor=(-.1, 0.5))

# cbar0.ax.yaxis.set_ticks([0,0.5, 1],   minor=False)
# cbar0.ax.yaxis.set_ticks([0.25, 0.75], minor=True)
# cbar1.formatter = cbformat


fig.savefig(f"../figs/supplementary/2state_heatmap_consumption_rates.png") #, dpi=100)

