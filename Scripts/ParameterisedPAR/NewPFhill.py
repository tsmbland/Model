import matplotlib as mpl

mpl.use('Agg')

import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/ParameterisedPAR/NewPFhill/'

from Models.PDE.PARhill import PAR
import Models.PDE.States as s
from Funcs import ParamSpace2D
import numpy as np
import copy

n = int(sys.argv[1])
print(n)

group = n // 20
n_adjusted = n % 20

"""

For 
00-20: eAneg = 1, ePneg = 1
20-40: eAneg = 2, ePneg = 1
40-60: eAneg = 1, ePneg = 2
60-80: eAneg = 2, ePneg = 2



0. kAP vs kPA, showing polarisible region
1. kAP vs kPA, showing polarisable region and line where domains are equal in size
2. kAP vs kPA, showing ASI of A
3. kAP vs kPA, showing ASI of P
4. kAP vs kPA, showing spontaneous polarity
5-9. kAP vs kPA, showing A dosage dependence
10-14. kAP vs kPA, showing P dosage dependence



"""

# Base parameter set

BaseModel = PAR(Da=0.28, Dp=0.15, konA=0.02115, koffA=0.0092, konP=0.00981, koffP=0.0073, kAP=0, kPA=0, eAneg=0,
                ePneg=0, xsteps=100, psi=0.174, Tmax=10000, deltat=0.01, L=67.3, pA=1.56, pP=1, kposA=0, kposP=0.882,
                khillA=0, khillP=22.19)

eAneg_vals = [1, 2, 1, 2]
ePneg_vals = [1, 1, 2, 2]

eAneg = eAneg_vals[group]
ePneg = ePneg_vals[group]
BaseModel.eAneg = eAneg
BaseModel.ePneg = ePneg

dosages = [0.2, 0.4, 0.6, 0.8, 1]

"""
Polarisible region

"""

if n_adjusted == 0:
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

if n_adjusted == 1:
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

if n_adjusted == 2:
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

if n_adjusted == 3:
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

if n_adjusted == 4:
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

if n_adjusted in range(5, 10):
    dosage = dosages[n_adjusted - 5]
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

if n_adjusted in range(10, 15):
    dosage = dosages[n_adjusted - 10]
    BaseModel.pP *= dosage


    def func(kAP, kPA):
        m = copy.deepcopy(BaseModel)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)

####################################################################################


ParamSpace2D(func, p1_range=[-4, 0], p2_range=[-4, 0], cores=32, resolution0=20, resolution_step=2,
             n_iterations=6, direc='%s/%s/%s' % (save_direc, group, n_adjusted), parallel=True, crange=[1, 6]).run()
