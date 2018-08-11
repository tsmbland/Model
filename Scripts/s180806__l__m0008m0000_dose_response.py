"""


"""

import M as x
import Models.m0010 as m0010
import Models.m0008b as m0008b
import numpy as np
import matplotlib.pyplot as plt

x.datadirec = '../../ModelData'

########## m0000

p0 = m0010.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                  eAneg=2, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, starts=[0, 0, 0, 1])

# x.alg_parsim(m=m0010.Model(p0), params=['pc_0'], vals=[np.linspace(0, 1, 8)], jobid=8, subjobid=1)

########## m0008

x.datadirec = '../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'
# p0 = x.loaddata(7, 0, 604).params
p0 = x.loaddata(11, 0, 706).params
p0.Tmax = 1000
p0.eqTmax = 0
p0.deltat = 0.1
p0.ac_0 = 0
x.datadirec = '../../ModelData'

x.alg_parsim(m=m0008b.Model(p0), params=['pc1_0'], vals=[np.linspace(0, 1, 8)], jobid=9, subjobid=0, compression=0)

###### Plot

for sim in range(0, 8):
    res0 = x.loaddata(8, 1, sim)
    plt.scatter(res0.pc, np.mean(res0.pco), c='k', s=10)

    res1 = x.loaddata(9, 0, sim)
    plt.scatter(res1.pcyt, np.mean(res1.pco), c='r', s=10)

x.sns.despine()
plt.xticks([])
plt.yticks([])
# plt.xlabel('[Cytoplasmic PAR-2]')
# plt.ylabel('[Cortical PAR-2]')

plt.show()


# ##############################################
#
# # Control: A dose response
#
# x.datadirec = '../../ModelData'
#
# ########## m0000
#
# p0 = m0010.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
#                   eAneg=2, L=67.3, xsteps=500, psi=0.174, Tmax=100, deltat=0.1, starts=[0, 0, 0, 0])
#
# # x.alg_parsim(m=m0010.Model(p0), params=['ac_0'], vals=[np.linspace(0, 3, 16)], jobid=8, subjobid=1)
#
# ########## m0008
#
# x.datadirec = '../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'
# p0 = x.loaddata(7, 0, 604).params
# p0.Tmax = 100
# p0.deltat = 0.1
# p0.pc1_0 = 0
#
# x.datadirec = '../../ModelData'
#
# # x.alg_parsim(m=m0008b.Model(p0), params=['ac_0'], vals=[np.linspace(0, 3, 16)], jobid=9, subjobid=0, compression=0)
#
# ####### Plot
#
# for sim in range(0, 16):
#     res0 = x.loaddata(8, 1, sim)
#     plt.scatter(res0.ac, np.mean(res0.aco), c='k', s=10)
#
#     res1 = x.loaddata(9, 0, sim)
#     plt.scatter(res1.acyt, np.mean(res1.aco), c='r', s=10)
#
# x.sns.despine()
# # plt.xticks([])
# # plt.yticks([])
# plt.show()
