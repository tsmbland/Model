import M as x
import Models.m0008 as m0008
import Models.m0008b as m0008b

x.datadirec = '../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'

# params = ['ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
# ranges = [[0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.001, 1], [0.0001, 1]]
# sims = [544,
#         827,
#         480,
#         532,
#         902,
#         202,
#         88,
#         604,
#         347,
#         348,
#         417,
#         454]
#
#
# x.cam_diagram(7, 0, sims, params, ranges)

# for sim in sims:
#     print(sim)
#     print(vars(x.loaddata(7, 1, sim).params))
#     x.plot_singlesim(7, 1, sim)



# # x.plot_singlesim(7, 0, 323)
# p = x.loaddata(9, 1, 750).params
# print(vars(p))
#
# p.Tmax = 1000
# p.Da = 0.28
# p.Dp = 0.15
# x.alg_singlesim(m0008b.Model(p), compression=0)
# # x.sliderplot()
# x.plot_singlesim()

# params.Tmax = 1000
# x.datadirec = '../ModelData'
# x.alg_singlesim(m0008.Model(params), compression=0)
# x.plot_singlesim()


# params = ['ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
# ranges = [[0.01, 1], [0.001, 0.1], [0.01, 1], [0.01, 1], [0.01, 1], [0.01, 1]]
# simids = [604, 412, 302, 198, 861]
#
# print(vars(x.loaddata(7, 1, 128).params))


# x.cam_diagram(jobid=11, subjobid=0, simids=simids, params=params, ranges=ranges)


import pandas as pd

res1 = pd.read_csv('%s/res.csv' % x.direc_to(11, 0))

# print((res1['domainsize_a'] < 67.3) & (res1['domainsize_p'] < 67.3))

list1 = (res1['simid'][(res1['domainsize_a'] < 67.3) & (res1['domainsize_p'] < 67.3)])

# Need to print results where a and p are polarised in model 1, and p is polarised in model 2


res2 = pd.read_csv('%s/res.csv' % x.direc_to(11, 1))

list2 = (res2['simid'][(res2['domainsize_p'] < 67.3)])

list3 = list(set(list1).intersection(list2))

params = ['ra', 'kon_p', 'kon_p_2', 'kd_f', 'kd_b', 'rp']
ranges = [[0.001, 0.1], [0.0001, 0.01], [0.1, 10], [0.01, 1], [0.01, 1], [0.0001, 0.01]]

# x.cam_diagram(11, 0, list3, params, ranges)

for i in [706]:
    # print(i)
    # res = x.loaddata(11,0,i)
    # print(res.pco[-1, -1])

    ax1 = x.plt.subplot2grid((1, 3), (0, 0))
    a, b = x.parplot_norm(ax1, 11, 0, i)

    ax2 = x.plt.subplot2grid((1, 3), (0, 1))
    x.parplot_norm(ax2, 11, 1, i, a, b)

    ax3 = x.plt.subplot2grid((1, 3), (0, 2))
    x.parplot_norm(ax3, 11, 2, i, a, b)

    x.plt.rcParams['savefig.dpi'] = 600
    x.plt.show()

# sim 512 is best

