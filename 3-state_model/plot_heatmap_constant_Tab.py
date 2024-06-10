import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("src")
from differential_equations import n0, S0, f, n_min, λ_min, δ_max, p_max


font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

λ_cmap = mpl.cm.get_cmap('viridis')
δ_cmap = mpl.cm.get_cmap('plasma')

#importing data 
data = sys.argv[1]

λd_opt = np.loadtxt(f'data/optimal_λd-Tab_{data}')
λr_opt = np.loadtxt(f'data/optimal_λr-Tab_{data}')
δ_opt  = np.loadtxt(f'data/optimal_δ-Tab_{data}')

# arrays for scaling imshow
#T = np.outer(np.ones_like(Tab), Tab) + np.outer(np.ones_like(T0), T0)
#T[0, 0] = 1


# setting up figure
fig, ax = plt.subplots(1, 3, figsize=(15, 4))#, sharey=True)
#title = r'$T_{AB}$'# + ' = ' + str(T_const)
#fig.suptitle(title)

# plotting lag
ax[0].set(xlabel=r'$T_0\:$', ylabel=r"$p$")#, title=r"$\lambda^*$")
ax[1].set(xlabel=r'$T_0\:$', ylabel=r"$p$")#, title=r"$\lambda^*$")
ax[2].set(xlabel=r'$T_0\:$', ylabel=r"$p$")#, title=r"$\delta^*$")

im0 = ax[0].imshow(λd_opt, origin="lower", cmap=λ_cmap, aspect="auto", vmin=0) #, extent=[T_arr.min(), T_arr.max(), 0, p_max])
im1 = ax[1].imshow(λr_opt * (δ_opt > 0), origin="lower", cmap=λ_cmap, aspect="auto", vmin=0) #, extent=[T_arr.min(), T_arr.max(), 0, p_max])
im2 = ax[2].imshow(δ_opt,  origin="lower", cmap=δ_cmap, aspect="auto", vmin=0, vmax=δ_max) #, extent=[T_arr.min(), T_arr.max(), 0, p_max])
#im2 = ax[2].imshow(δ_opt,  origin="lower", norm=mpl.colors.LogNorm(vmin=10**(-6), vmax=δ_max, clip=False), cmap=δ_cmap, aspect="auto", extent=[T_arr.min(), T_arr.max(), 0, p_max])

# Colorbars
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))

cbar.set_label(r"$\lambda^*$", rotation=0, labelpad=30)
cbar.set_label(r"$\lambda^*$", rotation=0, labelpad=30)
cbar.set_label(r"$\delta^*$", rotation=0, labelpad=50)

# ax[0].plot(p_arr, λd_opt)
# ax[0].set(xlabel="p", ylabel=r"$λ_d$")
# ax[1].plot(p_arr, λp_opt)
# ax[1].set(xlabel="p", ylabel=r"$λ_p$")
# ax[2].plot(p_arr, δ_opt)
# ax[2].set(xlabel="p", ylabel=r"$δ$")

# # plotting Ts
# ax[2].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$T_S$")
# im2 = ax[2].imshow(fitness, origin="lower", aspect="auto", extent=[T_arr.min(), T_arr.max(), 0, 1])
# cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))
# cbar.set_label(r"$\langle T_S \rangle_p^{-1} ~[h^{-1}]$", rotation=0, labelpad=50)

# saving
fig.tight_layout()
fig.savefig(f"figs/optimal_heatmap_Tab_{data}.png", dpi=100)