import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
# from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker
from Funcs import evaluate
import matplotlib as mpl

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPF/xkon0_linearHD/'

"""
Functions

"""


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
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^2 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
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
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^2 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
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
    ax.text(0.05, 0.05, r'$k_{ant} \: = \: %s \; \mu m^2 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
            fontsize=6)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.show()
    # plt.savefig(name + '.png', dpi=300, transparent=True)
    # plt.close()


"""
ASI

"""

# asi_plot(n=27, n2=32, name='linear_xkon0_ASI_0.025')
# asi_plot(n=28, n2=33, name='linear_xkon0_ASI_0.05')
# asi_plot(n=29, n2=34, name='linear_xkon0_ASI_0.5')

"""
Dosage

"""

# dosage_plot(n=10, n2=32, name='linear_xkon0_dosage_0.025')
# dosage_plot(n=15, n2=33, name='linear_xkon0_dosage_0.05')
# dosage_plot(n=20, n2=34, name='linear_xkon0_dosage_0.5')

"""
Trigger

"""

# trigger_plot(n1=14, n2=32, name='linear_xkon0_trigger_0.025')
# trigger_plot(n1=19, n2=33, name='linear_xkon0_trigger_0.05')
trigger_plot(n1=24, n2=34, name='linear_xkon0_trigger_0.5')
