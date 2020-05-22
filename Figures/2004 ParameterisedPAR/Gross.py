import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/ParameterisedPAR/GrossHD/'

"""
Functions

"""


# def find_bounds(array, val1, val2):
#     """
#     To do:
#     - way to automatically specify the best start point
#     - or different way of linking points that doesn't require start point to be specified
#
#     """
#
#     # Find boundaries
#     a = np.nonzero(np.logical_and(array[1:, :] == val1, array[:-1, :] == val2))
#     b = np.nonzero(np.logical_and(array[1:, :] == val2, array[:-1, :] == val1))
#     c = np.nonzero(np.logical_and(array[:, 1:] == val1, array[:, :-1] == val2))
#     d = np.nonzero(np.logical_and(array[1:, :] == val2, array[:-1, :] == val1))
#     unique_points = np.unique(np.c_[np.r_[a[0], b[0], c[0], d[0]] + 1, np.r_[a[1], b[1], c[1], d[1]] + 1], axis=0)
#     xpoints = unique_points[:, 0]
#     ypoints = unique_points[:, 1]
#
#     # Link points
#     xdist = np.tile(xpoints, (len(xpoints), 1)) - np.tile(xpoints, (len(xpoints), 1)).T
#     ydist = np.tile(ypoints, (len(ypoints), 1)) - np.tile(ypoints, (len(ypoints), 1)).T
#     dist = ((xdist ** 2) + (ydist ** 2)) ** 0.5
#     dist[dist == 0] = np.nan
#
#     x = np.argmax(ypoints)
#     order = [x]
#     while np.count_nonzero(~np.isnan(dist)) != 0:
#         next_x = np.nanargmin(dist[x, :])
#         order.append(next_x)
#         dist[:, x] = np.nan
#         dist[x, :] = np.nan
#         x = next_x
#
#     return xpoints[order], ypoints[order]
#
#
# res = np.loadtxt(datadirec + '1/Res.txt')
# xpoints, ypoints = find_bounds(res, 2, 3)


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
    res = (np.loadtxt(datadirec + str(n + 4) + '/Res.txt') == 2).astype(float)
    res += np.logical_and(res, (np.loadtxt(datadirec + str(n + 3) + '/Res.txt') == 2).astype(float))
    res += np.logical_and(res, (np.loadtxt(datadirec + str(n + 2) + '/Res.txt') == 2).astype(float))
    res += np.logical_and(res, (np.loadtxt(datadirec + str(n + 1) + '/Res.txt') == 2).astype(float))
    res += np.logical_and(res, (np.loadtxt(datadirec + str(n + 0) + '/Res.txt') == 2).astype(float))
    res[res == 0] = np.nan

    # Plot
    plt.close()
    fig, ax = plt.subplots()
    ax.imshow(res.T, origin='lower', vmin=1, vmax=5, cmap='cividis')
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
ASI A

"""

asi_plot(n=2, name='Gross_asi_a')

"""
ASI P

"""
asi_plot(n=3, name='Gross_asi_p')

"""
Dosage A

"""

# dosage_plot(n=5, name='Gross_dosage_a')

"""
Dosage P

"""

# dosage_plot(n=10, name='Gross_dosage_p')

"""
Trigger

"""

# trigger_plot(n1=0, n2=4, name='Gross_trigger')
