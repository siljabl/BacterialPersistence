import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("../src")

# Plot parameters
mpl.rcParams["font.size"]   = "12"

ab_res = 41
Tab = 12
T0  = 2

d_cmap = mpl.colormaps['viridis']
r_cmap = mpl.colormaps['plasma']
p_cmap = sns.color_palette('coolwarm', as_cmap=True)(np.linspace(0,1,5))

#importing data
x_Tab  = np.loadtxt(f'../data/supplementary/bethedging_triggered/single_optimal_x-T0_{T0}.txt')
λr_Tab = np.loadtxt(f'../data/supplementary/bethedging_triggered/single_optimal_λr-T0_{T0}.txt')
δ_Tab  = np.loadtxt(f'../data/supplementary/bethedging_triggered/single_optimal_δ-T0_{T0}.txt')
# λd_T0  = np.loadtxt(f'data/constant_Tab_12/single_optimal_λd-Tab_12.txt')
# λr_T0  = np.loadtxt(f'data/constant_Tab_12/single_optimal_λr-Tab_12.txt')
# δ_T0   = np.loadtxt(f'data/constant_Tab_12/single_optimal_δ-Tab_12.txt')


# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(3, 2, figsize=(6.75, 7.5), sharey=True, sharex='col')

# plotting lag
ax[0,0].set(ylabel=r"$p$", title=r"$x$")
ax[1,0].set(ylabel=r"$p$", title=r"$\omega^{\star}$")
ax[2,0].set(xlabel=r'$T_{ab}$', ylabel=r"$p$", title=r"$\delta^{\star}$")
ax[0,1].set(title=r"$\lambda^{\star} / ~T$")
ax[1,1].set( title=r"$\omega^{\star} / T$")
ax[2,1].set(xlabel=r'$T_{0}$', title=r"$\delta^{\star}$")


im00 = ax[0,0].imshow(x_Tab,  origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 24, 0, 1])
im01 = ax[1,0].imshow(λr_Tab, origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=24, extent=[0, 24, 0, 1])
im02 = ax[2,0].imshow(δ_Tab,          origin="lower", cmap=r_cmap, aspect="auto", vmin=0, vmax=0.06, extent=[0, 24, 0, 1])
# im10 = ax[0,1].imshow(λd_T0 / T_T0,   origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 12, 0, 1])
# im11 = ax[1,1].imshow(λr_T0 / T_T0,   origin="lower", cmap=d_cmap, aspect="auto", vmin=0, vmax=1, extent=[0, 12, 0, 1])
# im12 = ax[2,1].imshow(δ_T0,           origin="lower", cmap=r_cmap, aspect="auto", vmin=0, vmax=0.06, extent=[0, 12, 0, 1])

for j in range(3):
    ax[j,0].set(xticks=[0, 12, 24])
    ax[j,1].set(xticks=[0,  6, 12])

# Ticks
plt.gca().yaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().yaxis.set_ticks([0.25, 0.75], minor=True)

for j in range(3):
    ax[j,0].xaxis.set_ticks([0, 12, 24], minor=False)
    ax[j,0].xaxis.set_ticks([6, 18],     minor=True)
    
    ax[j,1].xaxis.set_ticks([0, 6, 12], minor=False)
    ax[j,1].xaxis.set_ticks([3, 9],     minor=True)

# saving
fig.tight_layout(rect=[0, 0, 1, 1],w_pad=2, h_pad=1.5)

# Colorbars
ims = [im00, im01, im02]
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((-2, 2))

# i = 0
# for im, axes in zip(ims, ax.flatten()):
#     cbar = fig.colorbar(im, ax=ax[i,:], aspect=10, anchor=(-.1, 0.5))
#     cbar.formatter = cbformat
    
#     if im in [im00, im01, im10, im11]:
#         cbar.ax.yaxis.set_ticks([0,0.5, 1],   minor=False)
#         cbar.ax.yaxis.set_ticks([0.25, 0.75], minor=True)

#     i += 1

fig.savefig(f"../figs/supplementary/bethedging_heatmaps.png", dpi=300)