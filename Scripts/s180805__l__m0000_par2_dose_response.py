"""


"""

import M as x
import Models.m0000 as m0000
import numpy as np
import matplotlib.pyplot as plt

x.datadirec = '../../ModelaData'

########## Polarised

p0 = m0000.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                  eAneg=2,
                  pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5,
                  Peqmin=0.5,
                  Peqmax=1)

# x.alg_parsim(m=m0000.Model(p0), params=['pP'], vals=[np.linspace(0, 2, 16)], jobid=8, subjobid=0)

########## Uniform

p0 = m0000.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                  eAneg=2,
                  pA=0, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.1, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5,
                  Peqmax=1)

x.alg_parsim(m=m0000.Model(p0), params=['pP'], vals=[np.linspace(0, 3, 16)], jobid=8, subjobid=1)

####### Plot

# for sim in x.simidlist(8, 0):
#     res = x.loaddata(8, 0, sim)
#     plt.scatter((res.params.pP - res.params.psi * np.mean(res.pco)), np.mean(res.pco[-1, res.aco[-1, :] < 0.01]), c='r',
#                 s=10)

for sim in x.simidlist(8, 1):
    res = x.loaddata(8, 1, sim)
    plt.scatter((res.params.pP - res.params.psi * np.mean(res.pco)), np.mean(res.pco), c='k', s=10)

x.sns.despine()
plt.xticks([])
plt.yticks([])
plt.show()
