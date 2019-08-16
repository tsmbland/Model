import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../..')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../MembraneQuant')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../../../ImageAnalysis_Scripts')

import numpy as np
import IA2 as ia2
import matplotlib.pyplot as plt
import seaborn as sns

"""
Raw data converted to micron units. 
Membrane values converted to cytoplasmic units
Multiplied by 0.001 to keep numbers manageable

"""


def func(direc):
    """
    Calculates r_local

    """
    res = np.zeros([len(ia2.direcslist(direc, 1)), 1000])
    for i, d in enumerate(ia2.direcslist(direc, 1)):
        res[i, :] = abs(ia2.interp_1d_array(np.loadtxt(d + '/ROI_norm.txt')[:, 1], 1000, 'cubic'))
    return np.mean(res, 0)


def fold(array):
    """
    Folds array in half

    """
    return np.flipud((array[:len(array) // 2] + np.flipud(array[len(array) // 2:])) / 2)


"""
PAR-2

"""

d1 = ia2.ddirec + 'Experiments/e1905__par2_ng_rundown_tom9/lp637_wt'
p2_mem = 0.001 * 0.255 * 1.63672480374 * fold(
    np.mean(ia2.importall(ia2.direcslist(d1, 1), 'Quantification_g/mem_spa1000.txt'), 0))
p2_mem_std = 0.001 * 0.255 * 1.63672480374 * fold(
    np.std(ia2.importall(ia2.direcslist(d1, 1), 'Quantification_g/mem_spa1000.txt'), 0))
p2_cyt = 0.001 * (0.255 ** 2) * np.mean(ia2.importall(ia2.direcslist(d1, 1), 'Quantification_g/cym.txt'), 0)
p2_cyt_std = 0.001 * (0.255 ** 2) * np.std(ia2.importall(ia2.direcslist(d1, 1), 'Quantification_g/cym.txt'), 0)
svr = np.mean(ia2.importall(ia2.direcslist(d1, 1), 'Quantification_g/svr.txt'), 0)
spres = np.mean(ia2.importall(ia2.direcslist(d1, 1), 'spres.txt'), 0)
r_local = fold(func(d1))

"""
PAR-6

"""

d2 = ia2.ddirec + 'Experiments/e1810__par6_rundown_stages_nelio/NEBD/kk1248_wt'
p6_mem = 0.001 * 0.255 * 1.63672480374 * fold(
    np.mean(ia2.importall(ia2.direcslist(d2, 1), 'Quantification_g/mem_spa1000.txt'), 0))
p6_mem_std = 0.001 * 0.255 * 1.63672480374 * fold(
    np.std(ia2.importall(ia2.direcslist(d2, 1), 'Quantification_g/mem_spa1000.txt'), 0))
p6_cyt = 0.001 * (0.255 ** 2) * np.mean(ia2.importall(ia2.direcslist(d2, 1), 'Quantification_g/cym.txt'), 0)
p6_cyt_std = 0.001 * (0.255 ** 2) * np.std(ia2.importall(ia2.direcslist(d2, 1), 'Quantification_g/cym.txt'), 0)

"""
Figures

"""

# plt.plot(p2_mem, c='b')
# plt.fill_between(range(500), p2_mem - p2_mem_std, p2_mem + p2_mem_std, color='b', alpha=0.2)
# sns.despine()
# plt.xticks([])
# plt.ylabel('Cortical concentration')
# plt.xlabel('Anterior - Posterior')
# plt.show()

print(svr)