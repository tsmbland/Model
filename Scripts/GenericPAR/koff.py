import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/GenericPAR/koff/'

from Models.PDE.PAR import PAR as PARs
from Models.PDE import States as s
from Funcs import ParamSpaceQual2D
import numpy as np
import copy

print(sys.argv[1])

"""
Exploring the effects of off rate variation

"""

# Generic parameter set
BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

dosages = [0.2, 0.4, 0.6, 0.8, 1]

koff_range = (-3, -1)
kant_range = (-4, 0)

"""
Dosage

"""

if int(sys.argv[1]) in range(0, 5):
    dosage = dosages[int(sys.argv[1])]
    BaseS.pP = dosage


    def func(koff, k):
        m = copy.deepcopy(BaseS)
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

"""
ASI

"""

if int(sys.argv[1]) == 5:
    def func(koff, k):
        m = copy.deepcopy(BaseS)
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_asi_p(m.A, m.P)

"""
Spontaneous

"""

if int(sys.argv[1]) == 6:
    def func(koff, k):
        m = copy.deepcopy(BaseS)
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate2()
        m.run(kill_uni=False, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

###############################################################################################


ParamSpaceQual2D(func, p1_range=koff_range, p2_range=kant_range, cores=32, resolution0=10, resolution_step=2,
                 n_iterations=3, direc=save_direc + sys.argv[1], parallel=True, crange=[1, 6]).run()
