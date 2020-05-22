import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/GenericPFmtcue/'

from Models.PDE.PARmtcue import PAR as PARs
from Models.PDE import States as s
from Funcs import ParamSpaceQual2D
import numpy as np
import itertools
import copy

print(sys.argv[1])

"""

Investigating the effects of positive feedback on polarity triggering
Non-linear model

- kant vs cue strength for different E



"""

# Generic parameter set
BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.001, kPA=0.001,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1, cue=0)

E_vals = [0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99]

"""
Trigger strength

"""

kon0 = 0.1
x = E_vals[int(sys.argv[1])]

BaseS.konP = kon0 * (1 - x)
BaseS.konA = kon0 * (1 - x)
BaseS.kposP = x * (BaseS.psi * kon0 + BaseS.koffP) / 1
BaseS.kposA = x * (BaseS.psi * kon0 + BaseS.koffA) / 1


def func(k, cue):
    m = copy.deepcopy(BaseS)
    m.kAP = 10 ** k
    m.kPA = 10 ** k
    m.initiate3(0)
    m.Tmax = 1000
    m.cue = cue
    m.run()
    m.cue = 0
    m.Tmax = 9000
    m.run()
    return s.par_state_pol2(m.A, m.P)


ParamSpaceQual2D(func, p1_range=[-4, 1], p2_range=[0, 1], cores=32, resolution0=10, resolution_step=2,
                 n_iterations=7, direc=save_direc + sys.argv[1], parallel=True, crange=[1, 6]).run()
