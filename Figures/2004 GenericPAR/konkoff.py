import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR
import copy
import matplotlib.ticker as ticker
import matplotlib as mpl

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/GenericPAR/konkoffHD/'

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
    ax.imshow(res.T, origin='lower', vmin=2, vmax=6, extent=(-3, 0, -3, -1), aspect=3 / 2, cmap='viridis')
    ax.set_xlabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    ax.set_ylabel(r'$k_{off} \; (s^{-1})$')
    ax.text(0.05, 0.9, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
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
    ax.imshow(restotal.T, origin='lower', vmin=1, vmax=5, extent=(-3, 0, -3, -1), aspect=3 / 2, cmap='cividis')
    ax.set_xlabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    ax.set_ylabel(r'$k_{off} \; (s^{-1})$')
    ax.text(0.05, 0.9, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
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
    ax.imshow(restotal.T, origin='lower', vmin=0, vmax=2, extent=(-3, 0, -3, -1), aspect=3 / 2, cmap=cmap, norm=norm,
              alpha=0.5)
    ax.set_xlabel(r'$k_{on} \; (\mu m \: s^{-1})$')
    ax.set_ylabel(r'$k_{off} \; (s^{-1})$')
    ax.text(0.05, 0.9, r'$k_{ant} \: = \: %s \; \mu m^4 \: s^{-1}$' % name.split('_')[-1], transform=ax.transAxes,
            fontsize=6)
    ax.tick_params(axis='both', labelsize=7)
    fig.set_size_inches(2.5, 2.5)
    fig.tight_layout()
    plt.savefig(name + '.png', dpi=300, transparent=True)
    plt.close()


# def pattern_plot(vals, cs, name):
#     plt.close()
#     fig1, ax1 = plt.subplots()
#     fig2, ax2 = plt.subplots()
#
#     M = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
#             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
#
#     for i, v in enumerate(vals):
#         m = copy.deepcopy(M)
#         m.konA = v
#         m.konP = v
#         m.initiate()
#         m.run()
#         ax1.plot(m.A)
#         ax1.plot(m.P)
#         ax2.plot(m.A / max(m.A))
#         ax2.plot(m.P / max(m.P))
#
#     fig1.set_size_inches(3, 3)
#     fig1.tight_layout()
#     fig1.savefig(name + '.png', dpi=300, transparent=True)
#
#     fig2.set_size_inches(3, 3)
#     fig2.tight_layout()
#     fig2.savefig(name + '_norm.png', dpi=300, transparent=True)
#
#     plt.close()


"""
ASI

"""

# asi_plot(n=27, n2=32, name='konkoff_ASI_0.005')
# asi_plot(n=28, n2=33, name='konkoff_ASI_0.01')
# asi_plot(n=29, n2=34, name='konkoff_ASI_0.1')


"""
Dosage

"""

# dosage_plot(n=10, n2=32, name='konkoff_dosage_0.005')
# dosage_plot(n=15, n2=33, name='konkoff_dosage_0.01')
# dosage_plot(n=20, n2=34, name='konkoff_dosage_0.1')

"""
Trigger

"""

# trigger_plot(n1=14, n2=32, name='konkoff_trigger_0.005')
# trigger_plot(n1=19, n2=33, name='konkoff_trigger_0.01')
# trigger_plot(n1=24, n2=34, name='konkoff_trigger_0.1')

# """
# Colourbars - vertical
#
# """
#
# # ASI
# n = 5
# plt.close()
# fig, ax = plt.subplots(figsize=(1.5, 4))
# fig.subplots_adjust(right=0.4)
# cmap = plt.get_cmap('viridis', n)
# norm = matplotlib.colors.BoundaryNorm(np.arange(0, n + 1) - 0.5, n)
# cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=np.arange(0, n))
# cbar.ax.set_yticklabels(['< 0.2', '0.2 - 0.35', '0.35 - 0.45', '0.45 - 0.49', '> 0.49'])
# cbar.ax.tick_params(size=0)
# plt.savefig('ASI_legend.png', dpi=300, transparent=True)
# plt.close()
#
# # Dosage
# n = 5
# plt.close()
# fig, ax = plt.subplots(figsize=(1.5, 4))
# fig.subplots_adjust(right=0.4)
# cmap = plt.get_cmap('cividis', n)
# norm = matplotlib.colors.BoundaryNorm(np.arange(0, n + 1) - 0.5, n)
# cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=np.arange(0, n))
# cbar.ax.set_yticklabels(['< 20%', '20% - 40%', '40% - 60%', '60% - 80%', '> 80%'])
# cbar.ax.tick_params(size=0)
# plt.savefig('Dosage_legend.png', dpi=300, transparent=True)
# plt.close()
#
"""
Colourbars - horizontal

"""

# ASI
n = 5
plt.close()
fig, ax = plt.subplots(figsize=(3, 1.2))
fig.subplots_adjust(bottom=0.6, top=0.8)
cmap = plt.get_cmap('viridis', n)
norm = matplotlib.colors.BoundaryNorm(np.arange(0, n + 1) - 0.7, n)
cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=np.arange(0, n), orientation='horizontal')
cbar.ax.set_xticklabels(['< 0.2', '0.2 - 0.35', '0.35 - 0.45', '0.45 - 0.49', '> 0.49'], rotation=45, ha='right')
cbar.ax.tick_params(size=0)
cbar.set_label('ASI')
ax.xaxis.set_label_position('top')
plt.savefig('ASI_legend_horizontal.png', dpi=300, transparent=True)
plt.close()

# Dosage
n = 5
plt.close()
fig, ax = plt.subplots(figsize=(3, 1.2))
fig.subplots_adjust(bottom=0.6, top=0.8)
cmap = plt.get_cmap('cividis', n)
norm = matplotlib.colors.BoundaryNorm(np.arange(0, n + 1) - 0.7, n)
cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=np.arange(0, n), orientation='horizontal')
cbar.ax.set_xticklabels(['< 20%', '20% - 40%', '40% - 60%', '60% - 80%', '> 80%'], rotation=45, ha='right')
cbar.ax.tick_params(size=0)
cbar.set_label('Maximum dosage imbalance')
ax.xaxis.set_label_position('top')
plt.savefig('Dosage_legend_horizontal.png', dpi=300, transparent=True)
plt.close()

# Trigger
n = 2
plt.close()
fig, ax = plt.subplots(figsize=(3, 1.2))
fig.subplots_adjust(bottom=0.6, top=0.8)
colours = {1: 'green', 2: 'grey'}
cmap = mpl.colors.ListedColormap(list(colours.values()))
norm = matplotlib.colors.BoundaryNorm(np.arange(0, n + 1) - 0.5, n)
cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, alpha=0.5, ticks=np.arange(0, n),
                                        orientation='horizontal')
cbar.ax.set_xticklabels(['Yes', 'No'])
cbar.ax.tick_params(size=0)
cbar.set_label('Requires trigger?')
ax.xaxis.set_label_position('top')
plt.savefig('Trigger_legend_horizontal.png', dpi=300, transparent=True)
plt.close()
