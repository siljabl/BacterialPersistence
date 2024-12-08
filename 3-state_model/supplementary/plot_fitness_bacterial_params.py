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
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('-bac_res', type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()

T = args.T0 + args.Tab


###############
## Load data ##
###############
F_max = np.loadtxt(f"../data/supplementary/fitness_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt")
δ_max = np.loadtxt(f"../data/supplementary/δ_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt")
λd    = np.loadtxt(f"../data/supplementary/λd_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt")
λr    = np.loadtxt(f"../data/supplementary/λr_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt")


folder = '../data/constant_Tab_12'
config = read_config(folder)    
ab_params = {'p':args.p, 'T0':args.T0, 'Tab':args.Tab, 'T':args.T0+args.Tab}
opt_params = identify_optimal_parameters_const_Tab(ab_params, config, folder)

ab_params = {'p':0.1, 'T0':args.T0, 'Tab':args.Tab, 'T':args.T0+args.Tab}
stoch_params = identify_optimal_parameters_const_Tab(ab_params, config, folder)

ab_params = {'p':args.p, 'T0':0, 'Tab':args.Tab, 'T':args.Tab}
trig_params = identify_optimal_parameters_const_Tab(ab_params, config, folder)


#################
## Plot figure ##
#################
sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True, sharex=True)

im0 = ax[0].imshow(F_max, origin="lower", aspect="auto", extent=[0, T+1,  0, T+1])
im1 = ax[1].imshow(δ_max, origin="lower", aspect="auto", extent=[0, T+1,  0, T+1])
ax[0].set(xlabel=r'$\lambda$',  ylabel=r"$\omega$", title=r"$F(\lambda, \omega, \delta^{\star})$")
ax[1].set(xlabel=r'$\lambda$',  title=r"$\delta^{\star}(\lambda, \omega)$")


# Colorbars
cbformat = mpl.ticker.ScalarFormatter(useMathText=False)
cbformat.set_powerlimits((0, 0))

cbar0 = fig.colorbar(im0, ax=ax[0], aspect=15, anchor=(-.1, 0.5))
cbar1 = fig.colorbar(im1, ax=ax[1], aspect=15, anchor=(-.1, 0.5))
cbar0.formatter = cbformat
cbar1.formatter = cbformat


# Ticks
ax_lim = T
ax[0].xaxis.set_ticks(np.linspace(0, ax_lim, 4, endpoint=True),  minor=False)
ax[0].yaxis.set_ticks(np.linspace(0, ax_lim, 4, endpoint=True),  minor=False)
ax[0].xaxis.set_ticks(np.linspace(0, ax_lim, 7, endpoint=True),  minor=True)
ax[0].yaxis.set_ticks(np.linspace(0, ax_lim, 7, endpoint=True),  minor=True)


# Optimal
for j in range(2):
    #ax[j].scatter(opt_params[0], opt_params[1], s=40, edgecolors='gray', color=color, marker='*')
    ax[j].scatter(stoch_params[0], stoch_params[1], s=50, color=sns.color_palette('colorblind')[2], marker=5, clip_on=False)
    ax[j].scatter(trig_params[0] * T / args.Tab, 0, s=50, color=sns.color_palette('colorblind')[4], marker=6, clip_on=False)


# saving
fig.tight_layout()
fig.savefig(f"../figs/supplementary/fitness_bacterial_params_p{args.p}_T0{args.T0}_Tab{args.Tab}.png", dpi=100)