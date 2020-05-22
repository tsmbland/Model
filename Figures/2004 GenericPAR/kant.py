import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPAR/kantHD/'

"""
Functions

"""


def asi_plot(n, name):
    print(name)

    plt.close()
    fig, ax = plt.subplots()

    # Import data
    res = np.loadtxt(open(datadirec + str(n) + '/Res.csv', 'rb'), delimiter=',')
    res = res[res[:, 0] > -3.0381282495667246]
    order = np.argsort(res[:, 0])

    # Plot line
    ax.plot(10 ** res[order, 0], res[order, 1])
    ax.axvline(10 ** -3.0381282495667246, c='k', linestyle='--')

    # Finalise figure
    ax.set_xlabel(r'$k_{ant} \; (\mu m^4 \: s^{-1})$')
    ax.set_ylabel('ASI')
    ax.tick_params(axis='both', labelsize=7)
    ax.set_xscale('log')
    ax.set_xlim(0.0001, 10)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    fig.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def dosage_plot(n, name):
    print(name)

    # Import data
    res = (np.loadtxt(datadirec + str(n) + '/Res.txt') == 2).astype(float)

    # Find boundaries
    a = np.nonzero(np.diff(res, axis=0))
    b = np.nonzero(np.diff(res, axis=1))
    unique_points = np.unique(np.c_[np.r_[a[0], b[0]] + 1, np.r_[a[1], b[1]] + 1], axis=0)
    xpoints = unique_points[:, 0]
    ypoints = unique_points[:, 1]

    # Link points
    xdist = np.tile(xpoints, (len(xpoints), 1)) - np.tile(xpoints, (len(xpoints), 1)).T
    ydist = np.tile(ypoints, (len(ypoints), 1)) - np.tile(ypoints, (len(ypoints), 1)).T
    dist = ((xdist ** 2) + (ydist ** 2)) ** 0.5
    dist[dist == 0] = np.nan
    x = np.argmax(ypoints)
    order = [x]
    while np.count_nonzero(~np.isnan(dist)) != 0:
        next_x = np.nanargmin(dist[x, :])
        order.append(next_x)
        dist[:, x] = np.nan
        dist[x, :] = np.nan
        x = next_x

    # Normalise
    p1_range = [-4, 1]
    p2_range = [0, 1]
    xpoints = p1_range[0] + xpoints * (p1_range[1] - p1_range[0]) / res.shape[0]
    ypoints = p2_range[0] + ypoints * (p2_range[1] - p2_range[0]) / res.shape[0]
    print(min(xpoints))

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.plot(10 ** xpoints[order], (1 - ypoints[order]) * 100)
    ax.axvline(10 ** -3.0381282495667246, c='k', linestyle='--')
    ax.set_xlabel(r'$k_{ant} \; (\mu m^4 \: s^{-1})$')
    ax.set_ylabel('Maximum dosage \n imbalance (%)', linespacing=1.5)
    ax.tick_params(axis='both', labelsize=7)
    ax.set_xscale('log')
    ax.set_xlim(0.0001, 10)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    fig.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def trigger_plot(n1, name):
    print(name)
    plt.close()
    fig, ax = plt.subplots()

    # Import data
    res = (np.loadtxt(datadirec + str(n1) + '/Res.txt') == 2).astype(float)

    # Find boundaries
    a = np.nonzero(np.diff(res, axis=0))
    b = np.nonzero(np.diff(res, axis=1))
    unique_points = np.unique(np.c_[np.r_[a[0], b[0]] + 1, np.r_[a[1], b[1]] + 1], axis=0)
    xpoints = unique_points[:, 0]
    ypoints = unique_points[:, 1]

    # Link points
    xdist = np.tile(xpoints, (len(xpoints), 1)) - np.tile(xpoints, (len(xpoints), 1)).T
    ydist = np.tile(ypoints, (len(ypoints), 1)) - np.tile(ypoints, (len(ypoints), 1)).T
    dist = ((xdist ** 2) + (ydist ** 2)) ** 0.5
    dist[dist == 0] = np.nan
    x = np.argmax(ypoints)
    order = [x]
    while np.count_nonzero(~np.isnan(dist)) != 0:
        next_x = np.nanargmin(dist[x, :])
        order.append(next_x)
        dist[:, x] = np.nan
        dist[x, :] = np.nan
        x = next_x

    # Normalise
    p1_range = [-4, 1]
    p2_range = [0, 0.5]
    xpoints = p1_range[0] + xpoints * (p1_range[1] - p1_range[0]) / res.shape[0]
    ypoints = p2_range[0] + ypoints * (p2_range[1] - p2_range[0]) / res.shape[0]

    # Shade regions
    area1 = [10 ** -4, 10 ** -3.0381282495667246]
    area2 = [10 ** -3.0381282495667246, 10 ** -2.22]
    area3 = [10 ** -2.22, 10]
    ax.axvspan(area2[0], area2[1], alpha=0.5, color='grey')
    ax.axvspan(area3[0], area3[1], alpha=0.5, color='green')

    # Plot lines
    ax.plot(10 ** xpoints[order][xpoints[order] > -3.0],
            ypoints[order][xpoints[order] > -3.0])
    ax.axvline(10 ** -3.0381282495667246, c='k', linestyle='--')
    ax.axvline(10 ** -2.22, c='k', linestyle='--')

    # Add text
    ax.text(0.07, 0.5, 'No polarity', transform=ax.transAxes, fontsize=7, verticalalignment='center',
            rotation='vertical')
    ax.text(0.26, 0.5, 'Spontaneous polarity', transform=ax.transAxes, fontsize=7, verticalalignment='center',
            rotation='vertical')
    # ax.text(0.8, 0.5, 'Triggered polarity', transform=ax.transAxes, fontsize=7, verticalalignment='center',
    #         rotation='vertical')
    ax.text(0.72, 0.4, 'Triggered polarity', transform=ax.transAxes, fontsize=7, horizontalalignment='center')

    # Finalise figure
    ax.set_xlabel(r'$k_{ant} \; (\mu m^4 \: s^{-1})$')
    ax.set_ylabel('Minimum trigger \n strength (ASI)', linespacing=1.5)
    ax.tick_params(axis='both', labelsize=7)
    ax.set_xscale('log')
    ax.set_xlim(0.0001, 10)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    fig.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


def pattern_plot(vals, name):
    print(name)

    plt.close()
    fig, ax = plt.subplots()

    M = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.01, kPA=0.01,
            ePneg=2, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

    alphas = [1, 0.6, 0.2]

    for i, v in enumerate(vals):
        print(i)
        m = copy.deepcopy(M)
        m.kAP = v
        m.kPA = v
        m.initiate()
        m.run()
        ax.plot(np.linspace(0, m.L, m.xsteps), m.A, c='tab:red', alpha=alphas[i])
        ax.plot(np.linspace(0, m.L, m.xsteps), m.P, c='tab:blue', alpha=alphas[i], label=str(v))

    ax.set_xlabel(r'$ x \; (\mu m)$')
    ax.set_ylabel('Concentration' + r'$(\mu m^{-2})$')
    legend = ax.legend(frameon=False, fontsize=8, title=r'$k_{ant} \; (\mu m^4 \: s^{-1})$')
    plt.setp(legend.get_title(), fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


"""
ASI

"""

# asi_plot(n='ASI', name='kant_ASI')

"""
Dosage

"""

# dosage_plot(n='Dosage', name='kant_dosage')

"""
Trigger

"""

trigger_plot(n1='Trigger', name='kant_trigger')

"""
Patterns

"""

# pattern_plot([1, 0.1, 0.005], name='kant_patterns')
