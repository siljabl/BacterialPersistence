import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("src")
from config_functions import read_config

# Plot parameters
font = {'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

d_cmap = mpl.colormaps['viridis']
r_cmap = mpl.colormaps['plasma']

# Import configuration
folder = sys.argv[1]
config = read_config(folder)

Tab    = int(config['Tab'])
T0_min = int(config['T0_min'])
T0_max = int(config['T0_max'])
T_arr  = Tab + np.linspace(T0_min, T0_max, int(config['ab_res']))
T      = np.outer(np.ones_like(T_arr), T_arr)
T[0,:] = 1

p_arr = np.linspace(0, 1, int(config['ab_res']))
p = np.outer(p_arr, np.ones_like(p_arr))

#importing data 
λd_opt = np.loadtxt(f'{folder}/single_optimal_λd-Tab_{Tab}.txt')
λr_opt = np.loadtxt(f'{folder}/single_optimal_λr-Tab_{Tab}.txt')
δ_opt  = np.loadtxt(f'{folder}/single_optimal_δ-Tab_{Tab}.txt')

λr_opt = λr_opt * (δ_opt > 0) 

# setting up figure
fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
# title = rf'$T_{0}$ = {T0}'
# fig.suptitle(title)

# plotting lag
ax[0].set(xlabel=r'$T_0$', ylabel=r"$p$", title=r"$\lambda^*_d / T$")
ax[1].set(xlabel=r'$T_0$', ylabel=r"$p$", title=r"$\lambda^*_r / T$")
ax[2].set(xlabel=r'$T_0$', ylabel=r"$p$", title=r"$\delta^*$")

im0 = ax[0].imshow(λd_opt / T, origin="lower", cmap=d_cmap, aspect="auto", vmin=0, extent=[T0_min, T0_max, 0, 1])
im1 = ax[1].imshow(λr_opt / T, origin="lower", cmap=r_cmap, aspect="auto", vmin=0, extent=[T0_min, T0_max, 0, 1])
im2 = ax[2].imshow(δ_opt,  origin="lower", cmap=r_cmap, aspect="auto", vmin=0, extent=[T0_min, T0_max, 0, 1])
#im2 = ax[2].imshow(δ_opt,  origin="lower", norm=mpl.colors.LogNorm(vmin=10**(-6), vmax=δ_max, clip=False), cmap=δ_cmap, aspect="auto", extent=[T_arr.min(), T_arr.max(), 0, p_max])

# Colorbars
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))


# saving
fig.tight_layout()
fig.savefig(f"figs/single_optimal/optimal_heatmap_Tab_{Tab}_rescaled.png", dpi=100)