import numpy as np
import matplotlib.pyplot as plt
from Funcs import evaluate
import seaborn as sns
from Models.PDE.PAR import PAR

k_vals = [0.001, 0.01, 0.1]
c = ['0.8', '0.6', '0']

"""
Model symmetric

"""

for i, k in enumerate(k_vals):
    m = PAR(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=k, kPA=k,
            ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)
    m.initiate()
    m.run(kill_stab=True)
    plt.plot(m.A, c=c[i])
    plt.plot(m.P, c=c[i])

plt.show()
