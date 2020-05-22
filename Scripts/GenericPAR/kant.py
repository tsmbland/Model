import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/GenericPAR/kantHD/'

from Models.PDE.PAR import PAR as PARs
from Models.PDE.PARflow import PAR as PARflow
from Models.PDE import States as s
from Funcs import ParamSpaceQual2D, ParamSpaceQuant1D
import numpy as np
import copy

print(sys.argv[1])

"""
Exploring the effects of antagonism rate variation

"""

# Generic parameter set
BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

kant_range = (-4, 1)

"""
Dosage

"""

if int(sys.argv[1]) == 0:
    def func(k, pP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(func, p1_range=kant_range, p2_range=[0, 1], cores=32, resolution0=10, resolution_step=2,
                     n_iterations=7, direc=save_direc + 'Dosage', parallel=True, crange=[1, 6]).run()

"""
Dosage wide

"""

if int(sys.argv[1]) == 3:
    def func(k, pP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(func, p1_range=[-4, 2], p2_range=[0, 1], cores=32, resolution0=10, resolution_step=2,
                     n_iterations=4, direc=save_direc + 'DosageWide', parallel=True, crange=[1, 6]).run()

"""
Trigger strength

"""

if int(sys.argv[1]) == 1:
    def func(k, asi):
        m = copy.deepcopy(BaseS)
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate3(asi)
        m.run(kill_uni=False, kill_stab=False)
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(func, p1_range=kant_range, p2_range=[0, 0.499], cores=32, resolution0=10, resolution_step=2,
                     n_iterations=7, direc=save_direc + 'Trigger', parallel=True, crange=[1, 6]).run()

"""
Flow strength

"""

if int(sys.argv[1]) == 4:
    BaseS = PARflow(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.0101, kposP=0, kAP=0.01, kPA=0.01,
                    ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1, v=0)


    def func(k, v):
        m = copy.deepcopy(BaseS)
        m.kAP = 10 ** k
        m.kPA = 10 ** k

        m.initiate3(0)
        m.Tmax = 1000
        m.v = v
        m.run()
        m.v = 0
        m.Tmax = 9000
        m.run()
        return s.par_state_pol2(m.A, m.P)


    ParamSpaceQual2D(func, p1_range=kant_range, p2_range=[0, 0.1], cores=32, resolution0=10, resolution_step=2,
                     n_iterations=2, direc=save_direc + 'Flow', parallel=True, crange=[1, 6]).run()

"""
ASI

"""

if int(sys.argv[1]) == 2:
    kant_vals = np.linspace(-4, 1, 500)


    def func(k):
        m = copy.deepcopy(BaseS)
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        ant = np.mean(m.P[:50])
        post = np.mean(m.P[50:])
        asi = abs((ant - post) / (2 * (ant + post)))
        return asi


    ParamSpaceQuant1D(func=func, p_vals=kant_vals, direc=save_direc + 'ASI', parallel=True, cores=32).run()
