import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/GenericPAR/konkoffHD/'

from Models.PDE.PAR import PAR as PARs
from Models.PDE import States as s
from Funcs import ParamSpaceQual2D
import numpy as np
import itertools
import copy

print(sys.argv[1])

"""



"""

# Generic parameter set
BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.001, kPA=0.001,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

dosages = [0.2, 0.4, 0.6, 0.8, 1]
kant_vals = [0.001, 0.002, 0.005, 0.01, 0.1]
kant_dosages = list(itertools.product(kant_vals, dosages))
p1_range = (-3, 0)
p2_range = (-3, -1)

"""
Dosage

"""

if int(sys.argv[1]) in range(0, 25):
    kant = kant_dosages[int(sys.argv[1])][0]
    BaseS.kAP = kant
    BaseS.kPA = kant
    dosage = kant_dosages[int(sys.argv[1])][1]
    BaseS.pP = dosage


    def func(kon, koff):
        m = copy.deepcopy(BaseS)
        m.konP = 10 ** kon
        m.konA = 10 ** kon
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

"""
ASI

"""

if int(sys.argv[1]) in range(25, 30):
    kant = kant_vals[int(sys.argv[1]) - 25]
    BaseS.kAP = kant
    BaseS.kPA = kant


    def func(kon, koff):
        m = copy.deepcopy(BaseS)
        m.konP = 10 ** kon
        m.konA = 10 ** kon
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_asi_p(m.A, m.P)

"""
Spontaneous

"""

if int(sys.argv[1]) in range(30, 35):
    kant = kant_vals[int(sys.argv[1]) - 30]
    BaseS.kAP = kant
    BaseS.kPA = kant


    def func(kon, koff):
        m = copy.deepcopy(BaseS)
        m.konP = 10 ** kon
        m.konA = 10 ** kon
        m.koffP = (10 ** koff) * 1.01
        m.koffA = 10 ** koff
        m.initiate2()
        m.run(kill_uni=False, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

###############################################################################################


ParamSpaceQual2D(func, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=20, resolution_step=2,
                 n_iterations=6, direc=save_direc + sys.argv[1], parallel=True, crange=[1, 6]).run()
