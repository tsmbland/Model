import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker
from Funcs import evaluate
import matplotlib as mpl

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPF/xkon0HD/'

"""
Functions

"""


def func(cyt, mem, kon, kpos, koff):
    return kon * cyt - koff * mem + kpos * cyt * mem


def rundown_plot():
    plt.close()
    fig, ax = plt.subplots()

    # Params
    psi = 0.1
    koff = 0.01
    kon0 = 0.1
    xvals = [0, 0.6, 0.9]
    cs = ['0.8', '0.4', '0']
    # alphas = [1, 0.6, 0.2]

    # Conservation of mass
    ax.plot([0, 1], [1 / psi, 0], linestyle='--')

    # Rundown
    for i, x in enumerate(xvals):
        kon = kon0 * (1 - x)
        kpos = x * (psi * kon0 + koff) / 1

        ax.scatter(*evaluate(func, xrange=(0, 0.6), yrange=(0, 10), iterations=5, resolution0=10, resolution_step=5,
                             args=(kon, kpos, koff)), s=0.1, c=cs[i])

    # Legend
    ax.plot([], [], c=cs[0], label='0')
    ax.plot([], [], c=cs[1], label='0.6')
    ax.plot([], [], c=cs[2], label='0.9')
    legend = ax.legend(frameon=False, fontsize=7, title=r'$\epsilon$')
    plt.setp(legend.get_title(), fontsize=7)

    # Finalise figure
    ax.set_xlabel('Cytoplasmic concentration ' + r'$(\mu m^{-3})$')
    ax.set_ylabel('Membrane \n concentration ' + r'$(\mu m^{-2})$')
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    plt.savefig('rundown.png', dpi=300, transparent=True)


def asi_plot(n, n2, name):
    print(name)

    # Import data
    res = np.loadtxt(datadirec + str(n) + '/Res.txt')
    res[res == 1] = np.nan

    # Mask
    res2 = np.loadtxt(datadirec + str(n2) + '/Res.txt')
    res2[res2 != 2] = 1
    res2[res2 == 2] = 0.2

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(res.T, origin='lower', vmin=2, vmax=6, extent=(0, 1, -3, 0), aspect=1 / 3, cmap='viridis')
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$k_{on0} \; (\mu m \: s^{-1})$')
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
            fontsize=6)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def dosage_plot(n, n2, name):
    print(name)

    # Import data
    res0 = (np.loadtxt(datadirec + str(n) + '/Res.txt') == 2).astype(float)
    res1 = (np.loadtxt(datadirec + str(n + 1) + '/Res.txt') == 2).astype(float)
    res2 = (np.loadtxt(datadirec + str(n + 2) + '/Res.txt') == 2).astype(float)
    res3 = (np.loadtxt(datadirec + str(n + 3) + '/Res.txt') == 2).astype(float)
    res4 = (np.loadtxt(datadirec + str(n + 4) + '/Res.txt') == 2).astype(float)
    restotal = res0 + res1 + res2 + res3 + res4
    restotal[restotal == 0] = np.nan

    # Mask
    res2 = np.loadtxt(datadirec + str(n2) + '/Res.txt')
    res2[res2 != 2] = 1
    res2[res2 == 2] = 0.2

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(restotal.T, origin='lower', vmin=1, vmax=5, extent=(0, 1, -3, 0), aspect=1 / 3, cmap='cividis')
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$k_{on0} \; (\mu m \: s^{-1})$')
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
            fontsize=6)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(2.5, 2.5)
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

    # Set up cmap
    colours = {1: 'green', 2: 'grey'}
    cmap = mpl.colors.ListedColormap(list(colours.values()))
    norm = mpl.colors.BoundaryNorm(list(colours.keys()) + [100], cmap.N)

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(restotal.T, origin='lower', vmin=0, vmax=2, extent=(0, 1, -3, 0), aspect=1 / 3, cmap=cmap, norm=norm,
              alpha=0.5)
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$k_{on0} \; (\mu m \: s^{-1})$')
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
            fontsize=6)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def pattern_plot(xs, cs, name):
    print(name)

    plt.close()
    fig, ax = plt.subplots()

    M = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
            ePneg=2, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

    kon0 = 0.1
    alphas = [0.2, 0.6, 1]

    for i, x in enumerate(xs):
        print(i)
        m = copy.deepcopy(M)

        m.konP = kon0 * (1 - x)
        m.konA = kon0 * (1 - x)
        m.kposP = x * (m.psi * kon0 + m.koffP) / 1
        m.kposA = x * (m.psi * kon0 + m.koffA) / 1
        m.initiate()
        m.run()
        ax.plot(np.linspace(0, m.L, m.xsteps), m.A, c='tab:red', alpha=alphas[i])
        ax.plot(np.linspace(0, m.L, m.xsteps), m.P, c='tab:blue', alpha=alphas[i])

    # Legend
    ax.plot([], [], c='k', alpha=alphas[0], label='0')
    ax.plot([], [], c='k', alpha=alphas[1], label='0.6')
    ax.plot([], [], c='k', alpha=alphas[2], label='0.9')
    legend = ax.legend(frameon=False, fontsize=7, title=r'$\epsilon$')
    plt.setp(legend.get_title(), fontsize=7)

    # Finalise figure
    ax.set_xlabel(r'$ x \; (\mu m)$')
    ax.set_ylabel('Concentration ' + r'$(\mu m^{-2})$')
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


# """
# ASI
#
# """
#
# asi_plot(n=27, n2=32, name='xkon0_ASI_0.005')
# asi_plot(n=28, n2=33, name='xkon0_ASI_0.01')
# asi_plot(n=29, n2=34, name='xkon0_ASI_0.1')
#
# """
# Dosage
#
# """
#
# dosage_plot(n=10, n2=32, name='xkon0_dosage_0.005')
# dosage_plot(n=15, n2=33, name='xkon0_dosage_0.01')
dosage_plot(n=20, n2=34, name='xkon0_dosage_0.1')
#
# """
# Trigger
#
# """
#
# trigger_plot(n1=14, n2=32, name='xkon0_trigger_0.005')
# trigger_plot(n1=19, n2=33, name='xkon0_trigger_0.01')
trigger_plot(n1=24, n2=34, name='xkon0_trigger_0.1')

"""
Rundown

"""

# rundown_plot()

"""
Patterns

"""

# pattern_plot(xs=[0, 0.6, 0.9], cs=['0.6', '0.3', '0'], name='patterns')
