import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/ParameterisedPAR/GrossHD/'

from Models.PDE.PAR import PAR
import Models.PDE.States as s
from Funcs import ParamSpaceQual2D
import numpy as np
import copy

print(sys.argv[1])

"""

For 
alpha=2, beta=1

0. kAP vs kPA, showing polarisible region
1. kAP vs kPA, showing polarisable region and line where domains are equal in size
2. kAP vs kPA, showing ASI of A
3. kAP vs kPA, showing ASI of P
4. kAP vs kPA, showing spontaneous polarity
5-9. kAP vs kPA, showing A dosage dependence
10-14. kAP vs kPA, showing P dosage dependence



"""

# Base parameter set

BaseModel = PAR(Da=0.28, Dp=0.15, konA=0.02115, koffA=0.0092, konP=0.13012, koffP=0.0073, kAP=0, kPA=0, eAneg=2,
                ePneg=1, xsteps=100, psi=0.174, Tmax=10000, deltat=0.01, L=67.3, pA=1.56, pP=1, kposA=0, kposP=0)

dosages = [0.2, 0.4, 0.6, 0.8, 1]

kAP_range = [-3, 0]
kPA_range = [-2, 1]

"""
Polarisible region

"""

if int(sys.argv[1]) == 0:
    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

"""
Equal domains

"""

if int(sys.argv[1]) == 1:
    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_bp(m.A, m.P)

"""
ASI of A

"""

if int(sys.argv[1]) == 2:
    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_asi_a(m.A, m.P)

"""
ASI of P

"""

if int(sys.argv[1]) == 3:
    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_asi_p(m.A, m.P)

"""
Spontaneous polarity

"""

if int(sys.argv[1]) == 4:
    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate2()
        m.run(kill_uni=False, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

"""
A dosage dependence

"""

if int(sys.argv[1]) in range(5, 10):
    dosage = dosages[int(sys.argv[1]) - 5]
    BaseModel.pA *= dosage


    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

"""
P dosage dependence

"""

if int(sys.argv[1]) in range(10, 15):
    dosage = dosages[int(sys.argv[1]) - 10]
    BaseModel.pP *= dosage


    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

####################################################################################


ParamSpaceQual2D(func, p1_range=kAP_range, p2_range=kPA_range, cores=32, resolution0=10, resolution_step=2,
                 n_iterations=7, direc=save_direc + sys.argv[1], parallel=True, crange=[1, 6]).run()
