import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("src")
from differential_equations import S0, f, p_min, λ_min, δ_max

T_min = 0
T_max = 12

font = {'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

λ_cmap = mpl.cm.get_cmap('viridis')
δ_cmap = mpl.cm.get_cmap('plasma')

#importing data 
data = sys.argv[1]

λd_opt = np.loadtxt(f'data/single_optimal/optimal_λd-T0_{data}')
λr_opt = np.loadtxt(f'data/single_optimal/optimal_λr-T0_{data}')
δ_opt  = np.loadtxt(f'data/single_optimal/optimal_δ-T0_{data}')

# arrays for scaling imshow
#T = np.outer(np.ones_like(Tab), Tab) + np.outer(np.ones_like(T0), T0)
#T[0, 0] = 1


# setting up figure
fig, ax = plt.subplots(1, 3, figsize=(15, 5))#, sharey=True)
#title = r'$T_{0}$'# + ' = ' + str(T_const)
#fig.suptitle(title)

# plotting lag
ax[0].set(xlabel=r'$T_{AB}$', ylabel=r"$p$")#, title=r"$\lambda^*$")
ax[1].set(xlabel=r'$T_{AB}$')#, ylabel=r"$p$")#, title=r"$\lambda^*$")
ax[2].set(xlabel=r'$T_{AB}$')#, ylabel=r"$p$")#, title=r"$\delta^*$")

im0 = ax[0].imshow(λd_opt, origin="lower", cmap=λ_cmap, aspect="auto", vmin=0, extent=[T_min, T_max, 0, 1])
im1 = ax[1].imshow(λr_opt * (δ_opt > 0), origin="lower", cmap=λ_cmap, aspect="auto", vmin=0, extent=[T_min, T_max, 0, 1])
im2 = ax[2].imshow(δ_opt,  origin="lower", cmap=δ_cmap, aspect="auto", extent=[T_min, T_max, 0, 1])
#im2 = ax[2].imshow(δ_opt,  origin="lower", norm=mpl.colors.LogNorm(vmin=10**(-6), vmax=δ_max, clip=False), cmap=δ_cmap, aspect="auto", extent=[T_arr.min(), T_arr.max(), 0, 1])

# Colorbars
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))

#cbar.set_label(r"$\lambda^*$", rotation=0, labelpad=30)
#cbar.set_label(r"$\lambda^*$", rotation=0, labelpad=30)
#cbar.set_label(r"$\delta^*$", rotation=0, labelpad=50)

# saving
fig.tight_layout()
fig.savefig(f"figs/single_optimal/optimal_heatmap_T0_{data}.png", dpi=100)