import sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append("../src")
from config_functions       import read_config
from differential_equations import δ_max
from optimal_from_file import identify_optimal_parameters_const_Tab


################
## Read input ##
################
parser = argparse.ArgumentParser(description='Run and plot fitness only for one set of antibiotic parameters.')
parser.add_argument('p',   type=float, help='probability of antibiotics')
#parser.add_argument('T0_max',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('-bac_res', type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()

p_spon = np.linspace(0,1,11)

###############
## Load data ##
sns.set_theme(style='ticks', font_scale=1.1)
color = sns.color_palette("crest", as_cmap=True)([0, 0.5, 1])
fig, ax = plt.subplots(3, 1, figsize=(5.5,5.5), sharey=False, sharex=True)

T0 = 4.0
δ  = np.loadtxt(f"../data/supplementary/δ_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λd = np.loadtxt(f"../data/supplementary/λd_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λr = np.loadtxt(f"../data/supplementary/λr_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")

T = T0*p_spon + args.Tab
ax[0].plot(p_spon, λd / T, '.-', color=color[0])
ax[1].plot(p_spon, λr / T, '.-', color=color[0], label=4)
ax[2].plot(p_spon, δ, '.-', color=color[0])


T0 = 6.0
δ  = np.loadtxt(f"../data/supplementary/δ_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λd = np.loadtxt(f"../data/supplementary/λd_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λr = np.loadtxt(f"../data/supplementary/λr_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")

T = T0*p_spon + args.Tab
ax[0].plot(p_spon, λd / T, '.-', color=color[1])
ax[1].plot(p_spon, λr / T, '.-', color=color[1], label=6)
ax[2].plot(p_spon, δ, '.-', color=color[1])

T0 = 8.0
δ  = np.loadtxt(f"../data/supplementary/δ_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λd = np.loadtxt(f"../data/supplementary/λd_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")
λr = np.loadtxt(f"../data/supplementary/λr_stochasticT0_p{args.p}_T0max{T0}_Tab{args.Tab}.txt")

T = T0*p_spon + args.Tab
ax[0].plot(p_spon, λd / T, '.-', color=color[2])
ax[1].plot(p_spon, λr / T, '.-', color=color[2], label=8)
ax[2].plot(p_spon, δ, '.-', color=color[2])

ax[0].set(ylabel=r"$\lambda^{\star} / ~T$")
ax[1].set(ylabel=r"$\omega^{\star} / ~T$")
ax[2].set(xlabel=r"Probability of desynchronized antibiotics", ylabel=r"$\delta^{\star}$")

fig.tight_layout(rect=[0, 0, 0.87, 1], h_pad=3)
ax[1].legend(loc='upper center',
               bbox_to_anchor=(1.1, 1.15),
               ncol=1,
               reverse=True,
               frameon=False,
               fancybox=True,
               handlelength=1,
               title=r"$T_0$",
               alignment='center')

sns.despine()
fig.savefig(f"../figs/supplementary/stochastic_T0_p{args.p}.png", dpi=300)
