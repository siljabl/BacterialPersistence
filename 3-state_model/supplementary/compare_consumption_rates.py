import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("src")

# Plot parameters
mpl.rcParams["font.size"]   = "12"

ab_res = 41
Tab = 12
T0  = 2

d_cmap = mpl.colormaps['viridis']
r_cmap = mpl.colormaps['plasma']
p_cmap = sns.color_palette('coolwarm', as_cmap=True)(np.linspace(0,1,5))

dir_1 = "high_resolution"
dir_2 = "high_resolution"
dir_3 = "correct_consumption_rate"

λd_Tab_paper = np.loadtxt(f'../data/{dir_1}/single_optimal_λd-T0_0.txt')
λr_Tab_paper = np.loadtxt(f'../data/{dir_1}/single_optimal_λr-T0_0.txt')
δ_Tab_paper  = np.loadtxt(f'../data/{dir_1}/single_optimal_δ-T0_0.txt')
λd_T0_paper  = np.loadtxt(f'../data/{dir_2}/single_optimal_λd-Tab_12.txt')
λr_T0_paper  = np.loadtxt(f'../data/{dir_2}/single_optimal_λr-Tab_12.txt')
δ_T0_paper   = np.loadtxt(f'../data/{dir_2}/single_optimal_δ-Tab_12.txt')

λd_Tab_correct = np.loadtxt(f'../data/{dir_3}/single_optimal_λd-T0_0.txt')
λr_Tab_correct = np.loadtxt(f'../data/{dir_3}/single_optimal_λr-T0_0.txt')
δ_Tab_correct  = np.loadtxt(f'../data/{dir_3}/single_optimal_δ-T0_0.txt')
λd_T0_correct  = np.loadtxt(f'../data/{dir_3}/single_optimal_λd-Tab_12.txt')
λr_T0_correct  = np.loadtxt(f'../data/{dir_3}/single_optimal_λr-Tab_12.txt')
δ_T0_correct   = np.loadtxt(f'../data/{dir_3}/single_optimal_δ-Tab_12.txt')

λd_Tab = (λd_Tab_paper - λd_Tab_correct) / λd_Tab_paper
λr_Tab = (λr_Tab_paper - λr_Tab_correct) / (λr_Tab_paper+1e-2)
δ_Tab  = (δ_Tab_paper  - δ_Tab_correct)  / (δ_Tab_paper+1e-2)
λd_T0  = (λd_T0_paper  - λd_T0_correct)  / λd_T0_paper
λr_T0  = (λr_T0_paper  - λr_T0_correct)  / (λr_T0_paper+1e-2)
δ_T0   = (δ_T0_paper   - δ_T0_correct)   / (δ_T0_paper+1e-2)


# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(3, 2, figsize=(6.75, 7.5), sharey=True, sharex='col')

# plotting lag
ax[0,0].set(title=r"$(\lambda^{\star}_{eq(1)} - \lambda^{\star}_{eq(2)}) ~/~ \lambda^{\star}_{eq(1)}$", ylabel=r"$p$")
ax[1,0].set(title=r"$(\omega^{\star}_{eq(1)}  - \omega^{\star}_{eq(2)})  ~/~ \omega^{\star}_{eq(1)}$",  ylabel=r"$p$")
ax[2,0].set(title=r"$(\delta^{\star}_{eq(1)}  - \delta^{\star}_{eq(2)})  ~/~ \delta^{\star}_{eq(1)}$",  xlabel=r'$T_{ab}$', ylabel=r"$p$")
ax[0,1].set(title=r"$(\lambda^{\star}_{eq(1)} - \lambda^{\star}_{eq(2)}) ~/~ \lambda^{\star}_{eq(1)}$")
ax[1,1].set(title=r"$(\omega^{\star}_{eq(1)}  - \omega^{\star}_{eq(2)})  ~/~ \omega^{\star}_{eq(1)}$")
ax[2,1].set(title=r"$(\delta^{\star}_{eq(1)}  - \delta^{\star}_{eq(2)})  ~/~ \delta^{\star}_{eq(1)}$",  xlabel=r'$T_{0}$')


im00 = ax[0,0].imshow(λd_Tab, origin="lower", cmap=d_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 24, 0, 1])
im01 = ax[1,0].imshow(λr_Tab, origin="lower", cmap=d_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 24, 0, 1])
im02 = ax[2,0].imshow(δ_Tab,  origin="lower", cmap=r_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 24, 0, 1])
im10 = ax[0,1].imshow(λd_T0,  origin="lower", cmap=d_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 12, 0, 1])
im11 = ax[1,1].imshow(λr_T0,  origin="lower", cmap=d_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 12, 0, 1])
im12 = ax[2,1].imshow(δ_T0,   origin="lower", cmap=r_cmap, aspect="auto", vmin=-1, vmax=1, extent=[0, 12, 0, 1])


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
ims = [im10, im11, im12]
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((-2, 2))

i = 0
for im, axes in zip(ims, ax.flatten()):
    cbar = fig.colorbar(im, ax=ax[i,:], aspect=10, anchor=(-.1, 0.5))
    cbar.formatter = cbformat
    i += 1

fig.savefig(f"../figs/supplementary/3state_heatmaps_consumption_rates.png", dpi=300)