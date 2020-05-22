import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPAR/kon/'

"""
Functions

"""


def asi_plot(n, name):
    print(name)

    # Import data
    res = np.loadtxt(datadirec + str(n) + '/Res.txt')
    res[res == 1] = np.nan

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(res.T, origin='lower', vmin=2, vmax=6, cmap='viridis')
    # ax.set_xlabel(r'$k_{pos} \; (\mu m^3 \: s^{-1})$')
    # ax.set_ylabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    # ax.tick_params(axis='both', labelsize=7)
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    fig.set_size_inches(3, 3)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def dosage_plot(n, name):
    print(name)

    # Import data
    res0 = (np.loadtxt(datadirec + str(n) + '/Res.txt') == 2).astype(float)
    res1 = (np.loadtxt(datadirec + str(n + 1) + '/Res.txt') == 2).astype(float)
    res2 = (np.loadtxt(datadirec + str(n + 2) + '/Res.txt') == 2).astype(float)
    res3 = (np.loadtxt(datadirec + str(n + 3) + '/Res.txt') == 2).astype(float)
    res4 = (np.loadtxt(datadirec + str(n + 4) + '/Res.txt') == 2).astype(float)
    restotal = res0 + res1 + res2 + res3 + res4
    restotal[restotal == 0] = np.nan

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(restotal.T, origin='lower', vmin=1, vmax=5, cmap='cividis')
    # ax.set_xlabel(r'$k_{pos} \; (\mu m^3 \: s^{-1})$')
    # ax.set_ylabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    # ax.tick_params(axis='both', labelsize=7)
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    fig.set_size_inches(3, 3)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def trigger_plot(n1, n2, name):
    print(name)

    # Import data
    res1 = (np.loadtxt(datadirec + str(n1) + '/Res.txt') == 2).astype(float)
    res2 = (np.loadtxt(datadirec + str(n2) + '/Res.txt') == 2).astype(float)
    restotal = res1 + res2
    restotal[restotal == 0] = np.nan

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(restotal.T, origin='lower', vmin=0, vmax=2, cmap='Greys')
    # ax.set_xlabel('kAP/PA')
    # ax.set_ylabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    # ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(3, 3)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def pattern_plot(vals, cs, name):
    plt.close()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    M = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
            ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

    for i, v in enumerate(vals):
        m = copy.deepcopy(M)
        m.konA = v
        m.konP = v
        m.initiate()
        m.run()
        ax1.plot(m.A)
        ax1.plot(m.P)
        ax2.plot(m.A / max(m.A))
        ax2.plot(m.P / max(m.P))

    fig1.set_size_inches(3, 3)
    fig1.tight_layout()
    fig1.savefig(name + '.png', dpi=300, transparent=True)

    fig2.set_size_inches(3, 3)
    fig2.tight_layout()
    fig2.savefig(name + '_norm.png', dpi=300, transparent=True)

    plt.close()


"""
ASI

"""

# asi_plot(n=5, name='kon_antag_ASI')

"""
Dosage

"""

# dosage_plot(n=0, name='kon_antag_dosage')

"""
Trigger

"""

# trigger_plot(n1=4, n2=6, name='kon_antag_trigger')

"""
Patterns

"""

pattern_plot([10 ** 0, 10 ** -0.5, 10 ** -1, 10 ** -1.5], [], name='kon_patterns')
