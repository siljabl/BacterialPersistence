import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("src")
from config_functions import read_config

# Plot parameters
#mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

ab_res = 41
Tab = 12
T0  = 2

d_cmap = mpl.colormaps['viridis']
r_cmap = mpl.colormaps['plasma']


#importing data
λd_Tab_low  = np.loadtxt(f'data/constant_T0_0_Tab_6/single_optimal_λd-T0_0.txt')
λd_Tab = np.loadtxt(f'data/constant_T0_0/single_optimal_λd-T0_0.txt')
λr_Tab = np.loadtxt(f'data/constant_T0_0/single_optimal_λr-T0_0.txt')
δ_Tab  = np.loadtxt(f'data/constant_T0_0/single_optimal_δ-T0_0.txt')
λd_T0  = np.loadtxt(f'data/constant_Tab_12/single_optimal_λd-Tab_12.txt')
λr_T0  = np.loadtxt(f'data/constant_Tab_12/single_optimal_λr-Tab_12.txt')
δ_T0   = np.loadtxt(f'data/constant_Tab_12/single_optimal_δ-Tab_12.txt')

λd_Tab[:, :11] = λd_Tab_low[:,:11]


# λr_opt = λr_opt * (δ_opt > 0) 

Tab_arr = np.linspace(0, 24, ab_res)
T0_arr  = np.linspace(0, 12, ab_res)
T_Tab = np.outer(np.ones_like(Tab_arr), Tab_arr)
T_T0  = Tab + np.outer(np.ones_like(T0_arr),  T0_arr)
T_Tab[:,0] = 1

p_arr = [0.1, 0.3, 0.5, 0.7, 0.9]

# setting up figure,
fig, ax = plt.subplots(2, 3, figsize=(6.75, 4.5), sharey=True)#, layout='constrained')
# title = rf'$T_{0}$ = {T0}'
# fig.suptitle(title)

# plotting lag
ax[0,0].set(xlabel=r'$T_{ab}$', title=r"$\lambda^{\star} / ~T$", ylabel=r"$p$")
ax[0,1].set(xlabel=r'$T_{ab}$', title=r"$\omega^{\star} / ~T$")
ax[0,2].set(xlabel=r'$T_{ab}$', title=r"$\delta^{\star}$")
ax[1,0].set(xlabel=r'$T_{0}$',  title=r"$\lambda^{\star} / ~T$", ylabel=r"$p$")
ax[1,1].set(xlabel=r'$T_{0}$',  title=r"$\omega^{\star} / T$")
ax[1,2].set(xlabel=r'$T_{0}$',  title=r"$\delta^{\star}$")

im00 = ax[0,0].imshow(λd_Tab / T_Tab, origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 24, 0, 1])
im01 = ax[0,1].imshow(λr_Tab / T_Tab, origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 24, 0, 1])
im02 = ax[0,2].imshow(δ_Tab,          origin="lower", cmap=r_cmap, aspect="auto", vmin=0, vmax=0.06, extent=[0, 24, 0, 1])
im10 = ax[1,0].imshow(λd_T0 / T_T0,   origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 12, 0, 1])
im11 = ax[1,1].imshow(λr_T0 / T_T0,   origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 12, 0, 1])
im12 = ax[1,2].imshow(δ_T0,           origin="lower", cmap=r_cmap, aspect="auto", vmin=0, vmax=0.06, extent=[0, 12, 0, 1])

cmap = sns.color_palette('coolwarm', as_cmap=True)
colors = cmap(np.linspace(0,1,5))
ax[0,0].scatter([Tab, Tab, Tab, Tab, Tab], p_arr, s=20, edgecolors='k', c=colors, marker='o')
ax[0,1].scatter([Tab, Tab, Tab, Tab, Tab], p_arr, s=20, edgecolors='k', c=colors, marker='o')
ax[0,2].scatter([Tab, Tab, Tab, Tab, Tab], p_arr, s=20, edgecolors='k', c=colors, marker='o')
ax[1,0].scatter([T0, T0, T0, T0, T0], p_arr, s=20, edgecolors='k', c='w', marker='o')
ax[1,1].scatter([T0, T0, T0, T0, T0], p_arr, s=20, edgecolors='k', c='w', marker='o')
ax[1,2].scatter([T0, T0, T0, T0, T0], p_arr, s=20, edgecolors='k', c='w', marker='o')

# Colorbars
fig.colorbar(im00, ax=ax[0,0], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im01, ax=ax[0,1], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im02, ax=ax[0,2], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im10, ax=ax[1,0], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im11, ax=ax[1,1], aspect=20, anchor=(-.1, 0.5))
fig.colorbar(im12, ax=ax[1,2], aspect=20, anchor=(-.1, 0.5))

# saving
fig.tight_layout()
fig.savefig(f"figs/3state_heatmaps.png", dpi=100)