import numpy as np
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR

"""
Functions

"""


def func(cyt, mem, kon, koff, kpos):
    return kon * cyt - koff * mem + kpos * cyt * mem


koff = 0.01
p = 1
kon_base = 0.1
psi = 0.05
y0 = (kon_base * p) / (psi * kon_base + koff)
kposmax = kon_base / y0
# kpos_vals = [0 * kposmax, 0.25 * kposmax, 0.5 * kposmax, 0.75 * kposmax, 0.9 * kposmax]
c = ['0.8', '0.6', '0.4', '0.2', '0']

kpos_vals = [0, 0.005, 0.0075, 0.01, 0.012]
kon_vals = [0.1, 0.065, 0.045, 0.025, 0.015]

# kpos_vals = [0, 0.01, 0.03, 0.05]
# koff_vals = [0.01, 0.015, 0.03, 0.045]

"""
Rundown relationship

"""
for i, kpos in enumerate(kpos_vals):
    kon = kon_base - (kpos * y0)

    plt.scatter(*evaluate(func, xrange=(0, 0.8), yrange=(0, 20), iterations=5, resolution0=10, resolution_step=5,
                          args=(kon, koff, kpos)), s=0.1, c=c[i])

# plt.xlim(left=0, right=0.7)
# plt.ylim(bottom=0)
sns.despine()
plt.xlabel('Cytoplasmic concentration')
plt.ylabel('Cortical concentration')
plt.show()

"""
Model symmetric

"""

# for i, kpos in enumerate(kpos_vals):
#     # kon = kon_vals[i]
#     koff = koff_vals[i]
#     # m = PAR(Da=0.1, Dp=0.1, konA=kon, koffA=0.01, kposA=kpos, konP=kon, koffP=0.01, kposP=kpos, kAP=0.01, kPA=0.01,
#     #         ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
#     m = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=koff, kposA=kpos, konP=0.1, koffP=koff, kposP=kpos, kAP=0.01, kPA=0.01,
#             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
#     m.initiate()
#     m.run(kill_stab=True)
#     plt.plot(m.A, c=c[i])
#     plt.plot(m.P, c=c[i])
#
# plt.show()

"""
Model asymmetric

"""

# for i, kpos in enumerate(kpos_vals):
#     kon = kon_base - (kpos * y0)
#     m = PAR(Da=0.01, Dp=0.01, konA=0.1, koffA=0.01, kposA=0, konP=kon, koffP=0.01, kposP=kpos, kAP=0.01, kPA=0.01,
#             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
#     m.initiate()
#     m.run(kill_stab=True)
#     plt.plot(m.A, c=c[i])
#     plt.plot(m.P, c=c[i])
#     # plt.scatter(m.pP - m.psi * np.mean(m.P), m.P[-1], c=c[i])
#
# plt.show()


"""

"""

# """
# 18: PAR, non-linear feedback, long
#
# """
#
# koffA = 0.0092
# koffP = 0.0073
# psi = 0.10318684114244771
# dosP = 0.294005475663175
# dosA = 1.05143336288 * dosP
#
# m = PAR(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
#         koffP=koffP, kposP=0, kAP=0.01, kPA=0.01, ePneg=1, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50,
#         psi=psi, pA=dosA, pP=dosP)
# m.initiate()
# m.run(kill_stab=True)
# plt.plot(m.A)
# plt.plot(m.P)
# plt.show()
#
# """
# 19: PAR, both (odr), long
#
# """
#
# koffA = 0.0092
# koffP = 0.0073
# psi = 0.10318684114244771
# dosP = 0.294005475663175
# dosA = 1.05143336288 * dosP
#
# m = PAR(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
#         koffP=koffP, kposP=12.7357711156 * koffP, kAP=0.01, kPA=0.01, ePneg=1, eAneg=2, xsteps=100, Tmax=10000,
#         deltat=0.01, L=50, psi=psi, pA=dosA, pP=dosP)
# m.initiate()
# m.run(kill_stab=True)
# plt.plot(m.A)
# plt.plot(m.P)
# plt.show()
