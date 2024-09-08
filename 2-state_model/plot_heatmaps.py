import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res = 400
lag_cmap = mpl.colormaps.get_cmap('viridis')
del_cmap = mpl.colormaps.get_cmap('plasma')

# importing data 
lag_Tab = np.loadtxt(f'data/high_resolution/optimal_lag-T05')
del_Tab = np.loadtxt(f'data/high_resolution/optimal_delta-T05')
lag_T0  = np.loadtxt(f'data/high_resolution/optimal_lag-Tab10')
del_T0  = np.loadtxt(f'data/high_resolution/optimal_delta-Tab10')

# arrays for scaling imshow
Tab_arr = np.linspace(0, 24, ab_res)
T0_arr  = np.linspace(0, 12, ab_res)
T_Tab = 5  + np.outer(np.ones_like(Tab_arr), Tab_arr)
T_T0  = 10 + np.outer(np.ones_like(T0_arr),  T0_arr)


# setting up figure
fig, ax = plt.subplots(2, 2, figsize=(6.75, 6), sharey=True)

# plotting lag
ax[0,0].set(xlabel=r'$T_{AB}$', title=r"$\lambda^* / ~T$", ylabel=r"$p$")
ax[0,1].set(xlabel=r'$T_{AB}$', title=r"$\delta^*$")
ax[1,0].set(xlabel=r'$T_{0}$',  title=r"$\lambda^* / ~T$", ylabel=r"$p$")
ax[1,1].set(xlabel=r'$T_{0}$',  title=r"$\delta^*$")

im00 = ax[0,0].imshow(lag_Tab / T_Tab, cmap=lag_cmap, origin="lower", vmin=0, vmax=1.1,   aspect="auto", extent=[0, 24, 0, 1]) #, cmap=lag_cmap,vmin=0
im01 = ax[0,1].imshow(del_Tab,         cmap=del_cmap, origin="lower", vmin=0, vmax=0.065, aspect="auto", extent=[0, 24, 0, 1])#, cmap=del_cmap, vmin=0, vmax=0.075
im10 = ax[1,0].imshow(lag_T0 / T_T0,   cmap=lag_cmap, origin="lower", vmin=0, vmax=1.1,   aspect="auto", extent=[0, 12, 0, 1]) #, cmap=lag_cmap, vmin=0
im11 = ax[1,1].imshow(del_T0,          cmap=del_cmap, origin="lower", vmin=0, vmax=0.065, aspect="auto", extent=[0, 12, 0, 1])#, cmap=del_cmap, vmin=0, vmax=0.075

fig.colorbar(im00, ax=ax[0,0], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im01, ax=ax[0,1], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im10, ax=ax[1,0], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im11, ax=ax[1,1], aspect=20, anchor=(-.1, 0.5))

# saving
fig.tight_layout()
fig.savefig(f"figs/2state_heatmap.png") #, dpi=100)

