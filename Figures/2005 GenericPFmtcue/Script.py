import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PARflow import PAR
import copy
import matplotlib.ticker as ticker

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPFmtcue/'

"""
Functions

"""


def trigger_plot(ns, name):
    print(name)
    plt.close()
    fig, ax = plt.subplots()

    labels = ['0', '0.6', '0.9']
    cs = ['0.8', '0.4', '0']

    for i, n in enumerate(ns):

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
        x = np.argmax(xpoints)
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

        # Plot line
        plt.plot(10 ** xpoints[order], ypoints[order], label=labels[i], c=cs[i])

    # Legend
    legend = ax.legend(frameon=False, fontsize=7, title=r'$\epsilon$')
    plt.setp(legend.get_title(), fontsize=7)

    # Finalise figure
    ax.set_xlabel(r'$k_{ant} \; (\mu m^4 \: s^{-1})$')
    ax.set_ylabel('Critical trigger strength', linespacing=1.5)
    ax.tick_params(axis='both', labelsize=7)
    ax.set_xscale('log')
    ax.set_xlim(0.001, 10)
    ax.set_ylim(top=1)
    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    sns.despine()
    fig.savefig(name + '.png', dpi=300, transparent=True)
    # plt.show()
    plt.close()


# def profiles(name):
#     print(name)
#     plt.close()
#     fig, ax = plt.subplots()
#
#     M = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.00, kPA=0.00,
#             ePneg=2, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01, L=50, psi=0.1, pA=1, pP=0, v=0.04)
#
#     e_vals = [0, 0.6, 0.9]
#     cs = ['0.8', '0.4', '0']
#     kon0 = 0.1
#
#     for i, x in enumerate(e_vals):
#         m = copy.deepcopy(M)
#         m.konP = kon0 * (1 - x)
#         m.konA = kon0 * (1 - x)
#         m.kposP = x * (m.psi * kon0 + m.koffP) / 1
#         m.kposA = x * (m.psi * kon0 + m.koffA) / 1
#
#         m.initiate3(0)
#         m.run()
#         plt.plot(np.linspace(0, 50, 100), m.A, c=cs[i], label=str(x))
#
#     # Legend
#     legend = ax.legend(frameon=False, fontsize=7, title=r'$\epsilon$')
#     plt.setp(legend.get_title(), fontsize=7)
#
#     # Finalise figure
#     ax.set_xlabel(r'$ x \; (\mu m)$')
#     ax.set_ylabel('Concentration ' + r'$(\mu m^{-2})$')
#     ax.set_ylim(bottom=0)
#     ax.tick_params(axis='both', labelsize=7)
#     fig.set_size_inches(3, 2.5)
#     fig.tight_layout()
#     sns.despine()
#     plt.savefig(name + '.png', dpi=300, transparent=True)
#     plt.close()


"""
Trigger plot

"""

trigger_plot([0, 3, 5], name='kant_c')


"""
Example profiles - no antagonism

"""

# profiles(namme='profiles')
