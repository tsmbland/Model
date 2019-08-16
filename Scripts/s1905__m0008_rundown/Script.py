import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '.')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../MembraneQuant')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../ImageAnalysis_Scripts')

import numpy as np
import IA2 as ia2
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import InVivo.s1905__par2_ng as iv
from matplotlib.widgets import Slider
import Models.m0008 as m
import copy

"""
Funcs

"""


def mse(xdata, ydata, xdata_model, ydata_model):
    # Interpolate model data
    xdata_model_itp = ia2.interp_1d_array(xdata_model, 10000)
    ydata_model_itp = ia2.interp_1d_array(ydata_model, 10000)

    # Loop through data points
    errors = np.zeros([len(xdata)])
    for i, x in enumerate(xdata):
        if sum(abs(x - xdata_model_itp) < 0.001) == 0:
            return np.nan
        else:
            try:
                errors[i] = min(ydata[i] - ydata_model_itp[abs(x - xdata_model_itp) < 0.001]) ** 2
            except:
                print('ERROR')
                return np.nan

    return np.mean(errors)


"""
Import data

"""


def import_mem(direcs):
    return ia2.bounded_mean_2d(np.array(ia2.importall(direcs, 'Quantification_g/mem_spa1000.txt')), [0.9, 0.1])


def import_cyt(direcs):
    return ia2.bounded_mean_2d(np.array(ia2.importall(direcs, 'Quantification_g/cyt_spa1000.txt')), [0.9, 0.1])


d = ia2.ddirec + '/Experiments/e1905__par2_ng_rundown_tom9/'
a_mem = import_mem(ia2.direcslist(d, 2))
a_cyt = import_cyt(ia2.direcslist(d, 2))
a_mem_wt = import_mem(ia2.direcslist(d + 'lp637_wt', 1))
a_cyt_wt = import_cyt(ia2.direcslist(d + 'lp637_wt', 1))
a_mem_uniform = import_mem(ia2.direcslist(d + 'nwg201_wt', 1) + ia2.direcslist(d + 'nwg201_rd', 1))
a_cyt_uniform = import_cyt(ia2.direcslist(d + 'nwg201_wt', 1) + ia2.direcslist(d + 'nwg201_rd', 1))
a_mem_polar = import_mem(ia2.direcslist(d + 'lp637_wt', 1) + ia2.direcslist(d + 'lp637_rd', 1))
a_cyt_polar = import_cyt(ia2.direcslist(d + 'lp637_wt', 1) + ia2.direcslist(d + 'lp637_rd', 1))

d2 = ia2.ddirec + '/Experiments/e1901__par2_cai/'
b_mem = import_mem(ia2.direcslist(d2, 2))
b_cyt = import_cyt(ia2.direcslist(d2, 2))
b_mem_wt = import_mem(ia2.direcslist(d2 + 'kk1273_wt', 1))
b_cyt_wt = import_cyt(ia2.direcslist(d2 + 'kk1273_wt', 1))
b_mem_uniform = import_mem(ia2.direcslist(d2 + 'kk1273_par6', 1) + ia2.direcslist(d2 + 'th415_par6', 1))
b_cyt_uniform = import_cyt(ia2.direcslist(d2 + 'kk1273_par6', 1) + ia2.direcslist(d2 + 'th415_par6', 1))
b_mem_polar = import_mem(ia2.direcslist(d2 + 'kk1273_wt', 1) + ia2.direcslist(d2 + 'th415_wt', 1))
b_cyt_polar = import_cyt(ia2.direcslist(d2 + 'kk1273_wt', 1) + ia2.direcslist(d2 + 'th415_wt', 1))

d3 = ia2.ddirec + '/Experiments/e1905__ring_crispr/nwg214/'
ring_mem = import_mem(ia2.direcslist(d3, 1))
ring_cyt = import_cyt(ia2.direcslist(d3, 1))

b_mem *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
b_mem_uniform *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt_uniform *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
b_mem_polar *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))
b_cyt_polar *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
ring_cyt *= (np.mean(a_cyt_wt) / np.mean(b_cyt_wt))
ring_mem *= (np.mean(a_mem_wt) / np.mean(b_mem_wt))

cyt_polar = np.r_[a_cyt_polar, b_cyt_polar] * 0.001 * (0.255 ** 2)
mem_polar = np.r_[a_mem_polar, b_mem_polar] * 0.001 * 0.255 * 1.63672480374
cyt_uniform = np.r_[a_cyt_uniform, b_cyt_uniform] * 0.001 * (0.255 ** 2)
mem_uniform = np.r_[a_mem_uniform, b_mem_uniform] * 0.001 * 0.255 * 1.63672480374
ring_cyt *= 0.001 * (0.255 ** 2)
ring_mem *= 0.001 * 0.255 * 1.63672480374

"""


"""

params_base = {'Da': 0, 'kon_a': 0, 'koff_a': 0, 'ra': 0, 'Dp': 0, 'kon_p': 0.0421342134213, 'koff_p': 0.0073,
               'kon_p_2': 0.1, 'kd_f': 0.1, 'kd_b': 0.1, 'rp': 0, 'xsteps': 0, 'psi': 0.10412, 'Tmax': 1000,
               'deltat': 0.1, 'deltal': 0.01, 'radii': 0, 'am_0': 0, 'ac_0': 0, 'pm1_0': 0,
               'pm2s_0': 0, 'pm2d_0': 0, 'pc1_0': 0, 'pc2_0': 0, 'spatial': False}

"""
Free parameters: 
- kd_f / kd_b
- kon_p_2

"""

# For each parameter set
# For a range of dosages
# Start with all monomeric in cytoplasm
# Simulate for ~ seconds
# For each real point
# Compare with closest point on graph


dosages = np.linspace(0, 1, 10)
kon_p_2_vals = np.linspace(0, 1, 10)
ratio_vals = np.linspace(0, 1, 10)

kon_p_2 = 0.5
ratio = 0.2

res = np.zeros([10, 10])

# for j, kon_p_2 in enumerate(kon_p_2_vals):
#     print(j)
#     for z, ratio in enumerate(ratio_vals):
#         print(z)
cyts = np.zeros([10])
mems = np.zeros([10])
for i, d in enumerate(dosages):
    params = copy.deepcopy(params_base)
    params['pc1_0'] = d
    params['kon_p_2'] = kon_p_2
    params['kd_f'] = ratio * params['kd_b']
    model = m.Model(**params)
    model.run()
    cyts[i] = model.pc1 + 2 * model.pc2
    mems[i] = model.pm1 + 2 * model.pm2s + 2 * model.pm2d
# res[j, z] = mse(cyt_uniform, mem_uniform, cyts, mems)

# plt.imshow(res, cmap='gist_rainbow', vmin=0, vmax=2, origin='lower')
plt.plot(cyts, mems)
plt.scatter(cyt_uniform, mem_uniform)
plt.show()
