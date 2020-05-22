import matplotlib as mpl

# mpl.use('Agg')

import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit
import scipy.odr as odr
import copy
from Models.ODE.DM import Model

datadirec = '/Volumes/lab-goehringn/working/Tom/ModelData/Dimerisation/'

"""
Functions

"""


def func(map, binary):
    print('i')

    # Import data
    resmap = np.loadtxt(map + '/Res.txt')

    print('j')
    resbin = np.loadtxt(binary + '/Res.txt')

    print('k')

    # # Set up figure
    # plt.close()
    # fig, ax = plt.subplots()
    #
    # # Map plot
    # ax.imshow(resmap.T, origin='lower', vmin=0, vmax=1, extent=(-2, 2, -2, 2), cmap='viridis')

    # # Find boundaries
    # a = np.nonzero(np.diff(resbin, axis=0))
    # b = np.nonzero(np.diff(resbin, axis=1))
    # unique_points = np.unique(np.c_[np.r_[a[0], b[0]] + 1, np.r_[a[1], b[1]] + 1], axis=0)
    # xpoints = unique_points[:, 0]
    # ypoints = unique_points[:, 1]
    #
    # # Link points
    # xdist = np.tile(xpoints, (len(xpoints), 1)) - np.tile(xpoints, (len(xpoints), 1)).T
    # ydist = np.tile(ypoints, (len(ypoints), 1)) - np.tile(ypoints, (len(ypoints), 1)).T
    # dist = ((xdist ** 2) + (ydist ** 2)) ** 0.5
    # dist[dist == 0] = np.nan
    # x = np.argmax(xpoints)
    # order = [x]
    # while np.count_nonzero(~np.isnan(dist)) != 0:
    #     next_x = np.nanargmin(dist[x, :])
    #     order.append(next_x)
    #     dist[:, x] = np.nan
    #     dist[x, :] = np.nan
    #     x = next_x
    #
    # # Normalise
    # p1_range = [-2, 2]
    # p2_range = [-2, 2]
    # xpoints = p1_range[0] + xpoints * (p1_range[1] - p1_range[0]) / resbin.shape[0]
    # ypoints = p2_range[0] + ypoints * (p2_range[1] - p2_range[0]) / resbin.shape[0]
    #
    # # Plot line
    # plt.plot(xpoints[order], ypoints[order])

    # Finalise plot
    # fig.set_size_inches(2.5, 2.5)
    # fig.tight_layout()
    # plt.show()


# func(datadirec + '/e/2', datadirec + '/e_binary/2')
# func(datadirec + '/kon/2', datadirec + '/kon_binary/2')
# plt.show()

# a = (np.loadtxt(datadirec + '/kon/2/Res.txt') > np.log10(22.799483171274492)).astype(int)
# b = (np.loadtxt(datadirec + '/e/2/Res.txt') > 0.9261655550386308).astype(int)
# c = a + b
# plt.imshow(c.T, origin='lower', extent=(-2, 2, -2, 2))
# plt.show()

# c = a + b
# plt.imshow(c.T, origin='lower', extent=(-2, 2, -2, 2))
# plt.show()

# a = np.loadtxt(datadirec + '/e/2/Res.txt')
# plt.imshow(a.T, origin='lower')
# plt.show()

"""
Rundown data best fit (when kon=1, kon2=100)

kd_f 10 ** 0.844
kd_b 10 ** 0.358


"""

BaseModel = Model(kon=1, kon2=100, koff=1, kd_f=1, kd_b=1, psi=0.174, dosage=1)


def pf_profile(x, kon0n, E):
    psi = 0.174
    a = kon0n * (1 - E)
    b = E * (psi * kon0n + 1)
    y = (a * x) / (1 - (b * x))
    return y


def fitting(cyt, mem):
    # Ols fit
    popt, pcov = curve_fit(pf_profile, cyt, mem)
    a_min_0 = popt[0]
    b_min_0 = popt[1]

    # Odr fit
    def perform_odr(x, y):
        quadr = odr.Model(lambda B, x: pf_profile(x, B[0], B[1]))
        mydata = odr.Data(x, y)
        myodr = odr.ODR(mydata, quadr, beta0=[a_min_0, b_min_0])

        output = myodr.run()
        return output

    regression = perform_odr(cyt, mem)
    a_min = regression.beta[0]
    b_min = regression.beta[1]
    return a_min, b_min


def rundown(model):
    # Create results containers
    pm1 = np.zeros([100])
    pm2s = np.zeros([100])
    pm2d = np.zeros([100])
    pc1 = np.zeros([100])
    pc2 = np.zeros([100])

    # Perform simulations
    for i, d in enumerate(np.linspace(0, 1, 100)):
        m = copy.deepcopy(model)
        m.dosage = d

        sol = m.solve()
        pm1[i] = sol[0]
        pm2s[i] = sol[1]
        pm2d[i] = sol[2]
        pc1[i] = sol[3]
        pc2[i] = sol[4]

    # Compile
    cyt = pc1 + 2 * pc2
    mem = pm1 * 2 * pm2s + 2 * pm2d

    # Fit
    a, b = fitting(cyt, mem)

    return cyt, mem


m = copy.deepcopy(BaseModel)
m.kd_f = 10 ** 0.844
m.kd_b = 10 ** 0.358
cd, md = rundown(m)
plt.plot(cd, md)

plt.plot(np.linspace(0, 0.2, 100), pf_profile(np.linspace(0, 0.2, 100), 22.799483171274492, 0.9261655550386308))
plt.show()
