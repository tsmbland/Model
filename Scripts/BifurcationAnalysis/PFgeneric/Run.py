import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/BifurcationPFgeneric/'

from Models.PDE.PAR import PAR as PARs
from Models.ODE.PAR import PAR as PARu
from Funcs import ParamSpaceQual2D
import numpy as np
import copy

print(sys.argv[1])

"""

OBTAINING DIFFERENT PF SCHEMES

Aim: with different pf values obtain the same domain conc by adjusting either kon or koff
1. Positive feedback (both, lin) vs on rate (both, lin)
2. Positive feedback (both, lin) vs off rate (both, lin)

> 7 parameter sets
a default
b weak pf (via on)
c medium pf (via on)
d high pf (via on)
e weak pf (via off)
f medium pf (via off)
g high pf (via off)



ROBUSTNESS ANALYSIS

With each parameter set perform
3-9. Dosage (P, lin) vs dosage (A, lin)
10-16. Antag (kAP, log) vs antag (kPA, log)
17-23. Dosage (P, lin) vs antag (both, log)



"""

# Generic parameter set

BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

BaseU = PARu(konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01, ePneg=2, eAneg=2,
             psi=0.1, pA=1, pP=1)

"""
OBTAINING DIFFERENT PF SCHEMES

"""

# 1. Positive feedback (both, lin) vs on rate (both, lin)
if int(sys.argv[1]) == 1:
    p1_range = (0, 0.015)
    p2_range = (0, 0.12)


    def spatial(kpos, kon):
        m = copy.deepcopy(BaseS)
        m.kposA = kpos
        m.konA = kon
        m.kposP = kpos
        m.konP = kon
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateC(6.715780138212306)


    def uniform(kpos, kon):
        m = copy.deepcopy(BaseU)
        m.kposA = kpos
        m.konA = kon
        m.kposP = kpos
        m.konP = kon
        return m.bistability_instability()

# 2. Positive feedback (both, lin) vs off rate (both, lin)
if int(sys.argv[1]) == 1:
    p1_range = (0, 0.1)
    p2_range = (0, 0.1)


    def spatial(kpos, koff):
        m = copy.deepcopy(BaseS)
        m.kposA = kpos
        m.koffA = koff + 0.0001
        m.kposP = kpos
        m.koffP = koff
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateC(6.715780138212306)


    def uniform(kpos, koff):
        m = copy.deepcopy(BaseU)
        m.kposA = kpos
        m.koffA = koff + 0.0001
        m.kposP = kpos
        m.koffP = koff
        return m.bistability_instability()

"""
ROBUSTNESS ANALYSIS

"""

# 7 parameter sets (so far these are approximate)
psets = [{},
         {'kposP': 0.0025, 'konP': 0.085, 'kposA': 0.0025, 'konA': 0.085},
         {'kposP': 0.0075, 'konP': 0.048, 'kposA': 0.0075, 'konA': 0.048},
         {'kposP': 0.0125, 'konP': 0.01, 'kposA': 0.0125, 'konA': 0.01},
         {'kposP': 0.02, 'koffP': 0.025, 'kposA': 0.02, 'koffA': 0.025},
         {'kposP': 0.06, 'koffP': 0.05, 'kposA': 0.06, 'koffA': 0.05},
         {'kposP': 0.1, 'koffP': 0.075, 'kposA': 0.1, 'koffA': 0.075}]

# 3-9. Dosage (P, lin) vs dosage (A, lin)
if int(sys.argv[1]) in range(3, 10):
    pset = int(sys.argv[1]) - 3
    for key, value in psets[pset].items():
        setattr(BaseS, key, value)
        setattr(BaseU, key, value)

    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateB()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 10-16. Antag (kAP, log) vs antag (kPA, log)
if int(sys.argv[1]) in range(10, 17):
    pset = int(sys.argv[1]) - 10
    for key, value in psets[pset].items():
        setattr(BaseS, key, value)
        setattr(BaseU, key, value)

    p1_range = (-4, 0)
    p2_range = (-4, 0)


    def spatial(kAP, kPA):
        m = copy.deepcopy(BaseS)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateB()


    def uniform(kAP, kPA):
        m = copy.deepcopy(BaseU)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        return m.bistability_instability()

# 17-23. Dosage (P, lin) vs antag (both, log)
if int(sys.argv[1]) in range(17, 24):
    pset = int(sys.argv[1]) - 17
    for key, value in psets[pset].items():
        setattr(BaseS, key, value)
        setattr(BaseU, key, value)

    p1_range = (0, 2)
    p2_range = (-4, 0)


    def spatial(pP, k):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateB()


    def uniform(pP, k):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        return m.bistability_instability()

###############################################################################################

# Bifurcation2D(uniform, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=50, resolution_step=2,
#               n_iterations=6, direc=save_direc + sys.argv[1] + '/Uniform', parallel=True, crange=[1, 6]).run()

ParamSpaceQual2D(spatial, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=10, resolution_step=2,
                 n_iterations=8, direc=save_direc + sys.argv[1] + '/Spatial', parallel=True, crange=[1, 5]).run()
