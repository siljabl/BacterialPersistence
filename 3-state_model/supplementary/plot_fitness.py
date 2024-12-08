import skimage
from skimage import filters
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import skimage.morphology


def get_edge(heatmap, X, Y):
    mask = (heatmap > 0.01)
    heatmap[mask] = 1
    edge = filters.sobel(heatmap)
    edge[edge < 0.01] = 0
    edge[edge > 0.01] = 1
    edge = skimage.morphology.skeletonize(edge)
    edge = edge.astype(np.int8)
    x = X[edge == 1]
    y = Y[edge == 1]
    idx = np.argsort(x)

    return x[idx], y[idx]

# folder = args.folder
ab_res = 41
x_T0  = np.linspace(0, 12, ab_res)
x_Tab = np.linspace(0, 24, ab_res)
y = np.linspace(0, 1,  ab_res)

X_T0, Y = np.meshgrid(x_T0, y)
X_Tab,_ = np.meshgrid(x_Tab,y)

cmap = mpl.colormaps['viridis']

#importing data
F_T0  = np.loadtxt(f'../data/constant_T0_0/single_optimal_F_max.txt')
F_Tab = np.loadtxt(f'../data/constant_Tab_12/single_optimal_F_max.txt')

heatmap_T0  = np.loadtxt(f'../data/constant_T0_0/single_optimal_λd-T0_0.txt')
heatmap_Tab = np.loadtxt(f'../data/constant_Tab_12/single_optimal_λd-Tab_12.txt')

x_T0,  y_T0  = get_edge(heatmap_T0,  X_T0,  Y)
x_Tab, y_Tab = get_edge(heatmap_Tab, X_Tab, Y)

# setting up figure
sns.set_theme(style='ticks', font_scale=1.2)

fig, ax = plt.subplots(1, 2, figsize=(6.5, 3), sharey=True)

im0=ax[0].imshow(F_T0,  origin="lower", cmap=cmap, aspect="auto", extent=[0, 12, 0, 1])
im1=ax[1].imshow(F_Tab, origin="lower", cmap=cmap, aspect="auto", extent=[0, 24, 0, 1])


ax[0].plot(x_T0,  y_T0,  'k', lw=3)
ax[1].plot(x_Tab, y_Tab, 'k', lw=3)

ax[0].set(xlabel=r'$T_{0}$', ylabel=r"$p$", title=r"$F_{I}(\lambda^{\star}, \omega^{\star}, \delta^{\star};T_{0}=0)$")
ax[1].set(xlabel=r'$T_{ab}$', title=r"$F_{I}(\lambda^{\star}, \omega^{\star}, \delta^{\star};T_{ab}=12)$")

ax[0].set(xticks=[0,  6, 12])
ax[1].set(xticks=[0, 12, 24])

# Colorbars
cbar = fig.colorbar(im0, ax=ax[0], aspect=15)
cbar = fig.colorbar(im1, ax=ax[1], aspect=15)

# Ticks
plt.gca().yaxis.set_ticks([0, 0.5, 1],  minor=False)
plt.gca().yaxis.set_ticks([0.25, 0.75], minor=True)

ax[0].xaxis.set_ticks([0, 6, 12], minor=False)
ax[0].xaxis.set_ticks([3, 9],     minor=True)

ax[1].xaxis.set_ticks([0, 12, 24], minor=False)
ax[1].xaxis.set_ticks([6, 18],     minor=True)


# saving
fig.tight_layout(w_pad=0.45)
fig.savefig(f"3statefitness_optimal_parameters.png", dpi=100)