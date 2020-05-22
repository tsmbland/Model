import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/GenericPFflow/'

from Models.PDE.PARflow import PAR as PARs
from Models.PDE import States as s
from Funcs import ParamSpaceQual2D
import numpy as np
import itertools
import copy

# print(sys.argv[1])

"""

Investigating the effects of positive feedback on polarity triggering
Non-linear model

- E vs trigger ASI for different kon0
- E vs trigger ASI for different kant
- E vs flow strength for different kon
- E vs flow strength for different kant

kon values: -2, -1, 0
kant values: 0.005, 0.1



"""

# Generic parameter set
BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.001, kPA=0.001,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1, v=0)

kon_vals = [0.01, 0.1, 1]
kant_vals = [0.005, 0.1]
kon_kant = list(itertools.product(kon_vals, kant_vals))

if int(sys.argv[1]) in range(0, 6):
    kon0 = kon_kant[int(sys.argv[1])][0]
    kant = kon_kant[int(sys.argv[1])][1]
    BaseS.kAP = kant
    BaseS.kPA = kant


    def funcA(x, asi):
        m = copy.deepcopy(BaseS)
        m.konP = kon0 * (1 - x)
        m.konA = kon0 * (1 - x)
        m.kposP = x * (m.psi * kon0 + m.koffP) / 1
        m.kposA = x * (m.psi * kon0 + m.koffA) / 1

        m.initiate3(asi)
        m.run(kill_uni=False, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(funcA, p1_range=[0, 1], p2_range=[0, 0.499], cores=32, resolution0=20, resolution_step=2,
                     n_iterations=1, direc=save_direc + str(sys.argv[1]), parallel=True, crange=[1, 6]).run()

if int(sys.argv[1]) in range(6, 12):
    kon0 = kon_kant[int(sys.argv[1]) - 6][0]
    kant = kon_kant[int(sys.argv[1]) - 6][1]
    BaseS.kAP = kant
    BaseS.kPA = kant


    def funcB(x, flow):
        m = copy.deepcopy(BaseS)
        m.konP = kon0 * (1 - x)
        m.konA = kon0 * (1 - x)
        m.kposP = x * (m.psi * kon0 + m.koffP) / 1
        m.kposA = x * (m.psi * kon0 + m.koffA) / 1

        m.initiate3(0)
        m.Tmax = 1000
        m.v = flow
        m.run()
        m.v = 0
        m.Tmax = 9000
        m.run()
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(funcB, p1_range=[0, 1], p2_range=[0, 0.1], cores=32, resolution0=20, resolution_step=2,
                     n_iterations=1, direc=save_direc + str(sys.argv[1]), parallel=True, crange=[1, 6]).run()
