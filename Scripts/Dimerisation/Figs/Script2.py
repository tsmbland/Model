import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

"""
Figs:
1: mem vs cyt full dosage
2: mem vs cyt non-linearity
3: total dimerisation fraction full dosage
4: total dimerisation fraction non-linearity
5: effective off rate full dosage


"""

direc = '/Volumes/lab-goehringn/working/Tom/ModelData/Dimerisation/4/'
p1_vals = np.linspace(-2, 2, 20)
p2_vals = np.linspace(-2, 2, 20)
koff = 0.1
kon = 0.1


def res1():
    """
    Mem vs cyt at full dosage

    """

    res = np.zeros([len(p1_vals), len(p2_vals)])
    for i, p1_val in enumerate(p1_vals):
        for j, p2_val in enumerate(p2_vals):
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
            mems = np.loadtxt(direc + '/' + name + '/pco.txt')
            cyts = np.loadtxt(direc + '/' + name + '/pcy.txt')
            ratios = mems / cyts
            res[i, j] = ratios[-1]
    return res


# def res2():
#     """
#     Mem vs cyt non-linearity
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         print(i)
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             mems = np.loadtxt(direc + '/' + name + '/pco.txt')
#             cyts = np.loadtxt(direc + '/' + name + '/pcy.txt')
#             ratios = mems / cyts
#             res[i, j] = ratios[-1] / ratios[49]
#     return res


def res2b():
    """
    Mem vs cyt non-linearity, log-log fit

    """

    res = np.zeros([len(p1_vals), len(p2_vals)])
    for i, p1_val in enumerate(p1_vals):
        print(i)
        for j, p2_val in enumerate(p2_vals):
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
            mems = np.loadtxt(direc + '/' + name + '/pco.txt')
            cyts = np.loadtxt(direc + '/' + name + '/pcy.txt')
            results = np.polyfit(np.log10(cyts[1:]), np.log10(mems[1:]), 1)
            res[i, j] = results[0]

    return res


def res3():
    """
    Fraction in dimeric form (wt dosage)

    """

    res = np.zeros([len(p1_vals), len(p2_vals)])
    for i, p1_val in enumerate(p1_vals):
        for j, p2_val in enumerate(p2_vals):
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
            pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
            pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
            pm2d = np.loadtxt(direc + '/' + name + '/pm2d.txt')
            pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
            pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')

            pmem_m = pm1
            pmem_d = 2 * pm2s + 2 * pm2d
            pcyt_m = pc1
            pcyt_d = 2 * pc2
            ptot_m = pcyt_m + 0.1 * pmem_m
            ptot_d = pcyt_d + 0.1 * pmem_d

            fractions = ptot_d / (ptot_m + ptot_d)
            res[i, j] = fractions[-1]
    return res


# def res4():
#     """
#     Dimeric fraction non-linearity
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
#             pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
#             pm2d = np.loadtxt(direc + '/' + name + '/pm2d.txt')
#             pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
#             pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')
#
#             pmem_m = pm1
#             pmem_d = 2 * pm2s + 2 * pm2d
#             pcyt_m = pc1
#             pcyt_d = 2 * pc2
#             ptot_m = pcyt_m + 0.1 * pmem_m
#             ptot_d = pcyt_d + 0.1 * pmem_d
#
#             fractions = ptot_d / (ptot_m + ptot_d)
#             res[i, j] = fractions[-1] / fractions[49]
#     return res
#
#
# def res5():
#     """
#     Particle off rate????
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
#             pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
#             pm2d = np.loadtxt(direc + '/' + name + '/pm2d.txt')
#             koffs = koff * (pm1 + pm2s) / (pm1 + pm2s + pm2d)
#             res[i, j] = koffs[-1]
#     return res
#
#
# def res6():
#     """
#     On component
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
#             pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')
#             ons = kon * (pc1 + 2 * pc2)
#             res[i, j] = ons[-1]
#     return res
#
#
# def res7():
#     """
#     Off component
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
#             pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
#             kd_b = 10 ** p2_val
#             offs = koff * (pm1 + 2 * pm2s) + kd_b * pm2s
#             res[i, j] = offs[-1]
#     return res
#
#
# def res8():
#     """
#     Positive feedback component
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
#             pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
#             kd_f = 10 ** p1_val
#             pfs = kd_f * pm1 * pc1
#             res[i, j] = pfs[-1]
#     return res
#
#
# def res9():
#     """
#     Total on (Positive feedback + basal on)
#
#     """
#
#     res = np.zeros([len(p1_vals), len(p2_vals)])
#     for i, p1_val in enumerate(p1_vals):
#         for j, p2_val in enumerate(p2_vals):
#             name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
#             pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
#             pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
#             pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')
#             kd_f = 10 ** p1_val
#             ons = (kd_f * pm1 * pc1) + (kon * (pc1 + 2 * pc2))
#             res[i, j] = ons[-1]
#     return res


def res10():
    """
    Ensemble kon

    """

    res = np.zeros([len(p1_vals), len(p2_vals)])
    for i, p1_val in enumerate(p1_vals):
        for j, p2_val in enumerate(p2_vals):
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
            kd_f = 10 ** p1_val
            kd_b = 10 ** p2_val

            pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
            pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
            pm2d = np.loadtxt(direc + '/' + name + '/pm2d.txt')
            pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
            pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')

            flux_on = kon * (pc1 + 2 * pc2) + (kd_f * pm1 * pc1)
            cyt_conc = pc1 + (2 * pc2)
            ensemble_kon = flux_on / cyt_conc
            res[i, j] = ensemble_kon[-1]
    return res


def res11():
    """
    Ensemble koff

    """

    res = np.zeros([len(p1_vals), len(p2_vals)])
    for i, p1_val in enumerate(p1_vals):
        for j, p2_val in enumerate(p2_vals):
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p2_val))
            kd_f = 10 ** p1_val
            kd_b = 10 ** p2_val

            pm1 = np.loadtxt(direc + '/' + name + '/pm1.txt')
            pm2s = np.loadtxt(direc + '/' + name + '/pm2s.txt')
            pm2d = np.loadtxt(direc + '/' + name + '/pm2d.txt')
            pc1 = np.loadtxt(direc + '/' + name + '/pc1.txt')
            pc2 = np.loadtxt(direc + '/' + name + '/pc2.txt')

            flux_off = koff * (pm1 + 2 * pm2s) + (kd_b * pm2s)
            mem_conc = pm1 + 2 * (pm2s + pm2d)
            ensemble_koff = flux_off / mem_conc
            res[i, j] = ensemble_koff[-1]
    return res


def im_fig(res):
    fig, ax = plt.subplots()
    im = ax.imshow(res.T, origin='lower', cmap='gray',
                   extent=(p1_vals[0], p1_vals[-1], p2_vals[0], p2_vals[-1]))
    ax.set_xlabel('log10(kon / koff)')
    ax.set_ylabel('log10(kon2 / koff)')
    ax.set_xticks(np.arange(p1_vals[0], p1_vals[-1] + 1, 1.0))
    ax.set_yticks(np.arange(p2_vals[0], p2_vals[-1] + 1, 1.0))
    fig.set_size_inches(4, 4)
    cbar = fig.colorbar(im, fraction=0.046, pad=0.04)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Non-linearity score', rotation=270)
    fig.tight_layout()


"""

"""
res = res2b()
im_fig(res)
# plt.scatter(1, 1, c='r', s=70)
# plt.scatter(-1, -1, c='b', s=70)
# plt.scatter(1, -1, c='g', s=70)
# plt.scatter(-1, 1, c='y', s=70)
fig = plt.gcf()
fig.set_size_inches(4, 4)
plt.show()
# plt.savefig('Fig3b.png', dpi=300)

"""
Non-linearity vs exchange rate

"""


def func(direc):
    """
    Mem vs cyt non-linearity, log-log fit

    """

    res = np.zeros([len(p1_vals)])
    for i, p1_val in enumerate(p1_vals):
        print(i)
        try:
            name = '%s_%s' % (float('%.4g' % p1_val), float('%.4g' % p1_val))
            mems = np.loadtxt(direc + '/' + name + '/pco.txt')
            cyts = np.loadtxt(direc + '/' + name + '/pcy.txt')
            results = np.polyfit(np.log10(cyts[1:]), np.log10(mems[1:]), 1)
            # plt.plot(np.log10(cyts[1:]), np.log10(mems[1:]))
            res[i] = results[0]
        except:
            res[i] = np.nan

    plt.plot(p1_vals, res)


#
# # func('/Volumes/lab-goehringn/working/Tom/ModelData/Dimerisation50x50,4oom/0/')
# # func('/Volumes/lab-goehringn/working/Tom/ModelData/Dimerisation/1/')
# # func('/Volumes/lab-goehringn/working/Tom/ModelData/Dimerisation50x50,4oom/2/')
# plt.show()

"""
Single non-linearity fit

"""


# def func(name):
#     mems = np.loadtxt(direc + '/' + name + '/pco.txt')
#     cyts = np.loadtxt(direc + '/' + name + '/pcy.txt')
#     results = np.polyfit(np.log10(cyts[1:]), np.log10(mems[1:]), 1)
#
#     plt.plot(cyts, mems, label='Model data')
#     plt.plot(cyts, (10 ** results[1]) * (cyts ** results[0]), label='Best fit')
#     plt.xlabel('Cytoplasmic concentration')
#     plt.ylabel('Membrane concentration')
#     sns.despine()
#     plt.legend(frameon=False)
#     plt.show()
#
#
# func('1.579_1.579')
