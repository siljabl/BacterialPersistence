import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

#mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res = 400
lag_cmap = mpl.colormaps.get_cmap('viridis')
del_cmap = mpl.colormaps.get_cmap('plasma')

# importing data 
folder = 'high_resolution'
#lag_Tab = np.loadtxt(f'data/{folder}/competition_lag-T05.txt')
#del_Tab = np.loadtxt(f'data/{folder}/competition_delta-T05.txt')
lag_Tab = np.loadtxt(f'data/{folder}/optimal_lag-T05.txt')
del_Tab = np.loadtxt(f'data/{folder}/optimal_delta-T05.txt')
lag_T0  = np.loadtxt(f'data/{folder}/optimal_lag-Tab10.txt')
del_T0  = np.loadtxt(f'data/{folder}/optimal_delta-Tab10.txt')


# arrays for scaling imshow
Tab_arr = np.linspace(0, 24, ab_res)
T0_arr  = np.linspace(0, 12, ab_res)
T_Tab = 5  + np.outer(np.ones_like(Tab_arr), Tab_arr)
T_T0  = 10 + np.outer(np.ones_like(T0_arr),  T0_arr)


# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(2, 2, figsize=(6.75, 6), sharey=True)

# plotting lag
ax[0,0].set(xlabel=r'$T_{0}$',  title=r"$\lambda^{\star} / ~T$", ylabel=r"$p$")
ax[0,1].set(xlabel=r'$T_{0}$',  title=r"$\delta^{\star}$")
ax[1,0].set(xlabel=r'$T_{ab}$', title=r"$\lambda^{\star} / ~T$", ylabel=r"$p$")
ax[1,1].set(xlabel=r'$T_{ab}$', title=r"$\delta^{\star}$")

im00 = ax[0,0].imshow(lag_T0 / T_T0,   cmap=lag_cmap, origin="lower", vmin=0, vmax=1.1,   aspect="auto", extent=[0, 12, 0, 1]) #, cmap=lag_cmap, vmin=0
im01 = ax[0,1].imshow(del_T0,          cmap=del_cmap, origin="lower", vmin=0, vmax=0.065, aspect="auto", extent=[0, 12, 0, 1])#, cmap=del_cmap, vmin=0, vmax=0.075
im10 = ax[1,0].imshow(lag_Tab / T_Tab, cmap=lag_cmap, origin="lower", vmin=0, vmax=1.1,   aspect="auto", extent=[0, 24, 0, 1]) #, cmap=lag_cmap,vmin=0
im11 = ax[1,1].imshow(del_Tab,         cmap=del_cmap, origin="lower", vmin=0, vmax=0.065, aspect="auto", extent=[0, 24, 0, 1])#, cmap=del_cmap, vmin=0, vmax=0.075

# Colorbars
ims = [im00, im01, im10, im11]
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((-2, 2))  # adjust the limits as needed

for im, axes in zip(ims, ax.flatten()):
    cbar = fig.colorbar(im, ax=axes, aspect=15, anchor=(-.1, 0.5))
    cbar.formatter = cbformat
    if im in [im00, im10]:
        cbar.ax.yaxis.set_ticks([0,0.5, 1],   minor=False)
        cbar.ax.yaxis.set_ticks([0.25, 0.75], minor=True)

    if im in [im01, im11]:
        cbar.ax.yaxis.set_ticks([0, 0.03, 0.06], minor=False)
        cbar.ax.yaxis.set_ticks([0.015, 0.045],  minor=True)


# Ticks
for j in range(2):
    ax[0,j].xaxis.set_ticks([0, 6, 12],  minor=False)
    ax[0,j].xaxis.set_ticks([3, 9], minor=True)

    ax[1,j].xaxis.set_ticks([0, 12, 24],  minor=False)
    ax[1,j].xaxis.set_ticks([6, 18], minor=True)

plt.gca().yaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().yaxis.set_ticks([0.25, 0.75], minor=True)


# saving
fig.tight_layout(w_pad=2)
fig.savefig(f"figs/2state_heatmap.png") #, dpi=100)

