import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append('src')

#mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

T0  = 0
Tab = 12
T   = T0 + Tab
p   = 0.3



i = 1000
j = 9999
fig, ax = plt.subplots(3,1)

for p in [0.1, 0.5]:
    λd = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    λr = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    δ  = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")



    ax[0].plot(δ[i:j], '.-')
    ax[1].plot(λr[i:j]/12, '.-')
    ax[2].plot(λd[i:j]/12, '.-')

    ax[0].set(yscale='log')
    ax[1].set(yscale='log')
    ax[2].set(yscale='log')

fig.savefig("check.png")
