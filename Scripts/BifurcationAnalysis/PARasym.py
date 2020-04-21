import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/BifurcationPARasym/'

from Models.PDE.PAR import PAR as PARs
from Models.ODE.PAR import PAR as PARu
from Funcs import ParamSpaceQual2D
import numpy as np
import copy

print(sys.argv[1])

"""
# Nonlinear antag, positive feedback on P

1. Positive feedback (P, lin) vs on rate (P, lin) - medium antagonism
2. Positive feedback (P, lin) vs off rate (P, lin) - medium antagonism
3. Positive feedback (P, lin) vs kAP (P, lin) - medium antagonism
4. Positive feedback (P, lin) vs kPA (P, lin) - medium antagonism


# 10. Dosage (P, lin) vs antagonism (both, log) - low positive feedback
# 11. Dosage (P, lin) vs antagonism (both, log) - medium positive feedback
# 12. Dosage (P, lin) vs antagonism (both, log) - high positive feedback
# 19. Dosage (P, lin) vs antagonism (both, log) - v. high positive feedback
# 7. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
# 8. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
# 9. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
# 21. Dosage (A, lin) vs positive feedback (P, lin) - low antagonism
# 22. Dosage (A, lin) vs positive feedback (P, lin) - medium antagonism
# 23. Dosage (A, lin) vs positive feedback (P, lin) - high antagonism
# 10. Dosage (P, lin) vs dosage A (lin) - low positive feedback
# 11. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
# 12. Dosage (P, lin) vs dosage A (lin) - high positive feedback
# 19. Dosage (P, lin) vs dosage A (lin) - v. high positive feedback
# 
# 27. Dosage (P, lin) vs antagonism (both, log) - low positive feedback


"""

# Generic parameter set

BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

BaseU = PARu(konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01, ePneg=2, eAneg=2,
             psi=0.1, pA=1, pP=1)

# 1. Positive feedback (P, lin) vs on rate (P, lin) - medium antagonism
if int(sys.argv[1]) == 1:
    p1_range = (0, 0.2)
    p2_range = (0, 0.2)


    def spatial(kposP, konP):
        m = copy.deepcopy(BaseS)
        m.kposP = kposP
        m.konP = konP
        m.initiate()
        m.run(kill_uni=True, kill_stab=True)
        return m.stateB()


    def uniform(kposP, konP):
        m = copy.deepcopy(BaseU)
        m.kposP = kposP
        m.konP = konP
        return m.bistability_instability()

# 2. Positive feedback (P, lin) vs off rate (P, lin) - medium antagonism
if int(sys.argv[1]) == 2:
    p1_range = (0, 0.2)
    p2_range = (0, 0.02)


    def spatial(kposP, koffP):
        m = copy.deepcopy(BaseS)
        m.kposP = kposP
        m.koffP = koffP
        m.initiate()
        m.run(kill_uni=True, kill_stab=True)
        return m.stateB()


    def uniform(kposP, koffP):
        m = copy.deepcopy(BaseU)
        m.kposP = kposP
        m.konP = koffP
        return m.bistability_instability()

# 3. Positive feedback (P, lin) vs kAP (P, lin) - medium antagonism
if int(sys.argv[1]) == 3:
    p1_range = (0, 0.2)
    p2_range = (-4, 0)


    def spatial(kposP, kAP):
        m = copy.deepcopy(BaseS)
        m.kposP = kposP
        m.kAP = 10 ** kAP
        m.initiate()
        m.run(kill_uni=True, kill_stab=True)
        return m.stateB()


    def uniform(kposP, kAP):
        m = copy.deepcopy(BaseU)
        m.kposP = kposP
        m.kAP = 10 ** kAP
        return m.bistability_instability()

# 4. Positive feedback (P, lin) vs kPA (P, lin) - medium antagonism
if int(sys.argv[1]) == 4:
    p1_range = (0, 0.2)
    p2_range = (-4, 0)


    def spatial(kposP, kPA):
        m = copy.deepcopy(BaseS)
        m.kposP = kposP
        m.kPA = 10 ** kPA
        m.initiate()
        m.run(kill_uni=True, kill_stab=True)
        return m.stateB()


    def uniform(kposP, kPA):
        m = copy.deepcopy(BaseU)
        m.kposP = kposP
        m.kPA = 10 ** kPA
        return m.bistability_instability()

###############################################################################################


# Bifurcation2D(uniform, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=50, resolution_step=2,
#               n_iterations=8, direc=save_direc + sys.argv[1] + '/Uniform', parallel=True, crange=[1, 6]).run()

ParamSpaceQual2D(spatial, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=10, resolution_step=2,
                 n_iterations=5, direc=save_direc + sys.argv[1] + '/Spatial', parallel=True, crange=[1, 5]).run()
