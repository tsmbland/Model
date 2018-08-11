"""


"""

import M as x
import Models.m0008b as m0008b
import numpy as np
import matplotlib.pyplot as plt


# x.datadirec = '../../ModelData'

########## Uniform


# p0 = m0008b.Params(Da=0.28, kon_a=0.0085, koff_a=0.0054, ra=1, Dp=0.15, kon_p=0.001, koff_p=0.0073, kon_p_2=0.1,
#                    kd_f=1, kd_b=1, rp=1, L=67.3,
#                    xsteps=500, psi=0.174, Tmax=100, deltat=0.01, starts=[0, 0, 0, 0, 0, 1, 0])

x.datadirec = '../../../../../../../Volumes/lab-goehringn/working/Tom/ModelData'
p0 = x.loaddata(7, 0, 202).params
p0.Tmax = 100
x.datadirec = '../../ModelData'

# x.alg_parsim(m=m0008b.Model(p0), params=['pc1_0'], vals=[np.linspace(0, 3, 16)], jobid=9, subjobid=0, compression=0)


########## Polarised


p0 = m0008b.Params(Da=0.28, kon_a=0.0085, koff_a=0.0054, ra=1, Dp=0.15, kon_p=0.0474, koff_p=0.0073, kon_p_2=0.0474,
                   kd_f=1, kd_b=1, rp=1, L=67.3,
                   xsteps=500, psi=0.174, Tmax=100, deltat=0.01, starts=[0, 0, 0, 0, 0, 1, 0])

# x.alg_parsim(m=m0008b.Model(p0), params=['pc1_0'], vals=[np.linspace(0, 3, 8)], jobid=9, subjobid=1, compression=0)

####### Plot

for sim in range(0, 16):
    res = x.loaddata(9, 0, sim)
    c = res.pcyt[-1]
    m = np.mean(res.pco[-1, :])
    plt.scatter(c, m, c='k', s=10)


for sim in range(0, 8):
    res = x.loaddata(9, 1, sim)
    c = res.pcyt[-1]
    m = np.mean(res.pco[-1, res.aco[-1, :] < 0.05])
    plt.scatter(c, m, c='r')

x.sns.despine()
plt.yticks([])
plt.xticks([])
plt.show()

