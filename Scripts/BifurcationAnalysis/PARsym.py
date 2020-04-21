import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/BifurcationPARsym/'

from Models.PDE.PAR import PAR as PARs
from Models.ODE.PAR import PAR as PARu
from Funcs import ParamSpaceQual2D
import numpy as np
import copy

print(sys.argv[1])

"""

Nonlinear antag

1. Dosage (P, lin) vs antagonism (both, log)
2. Dosage (P, lin) vs antagonism (P, log)
3. Antagonism (P, log) vs antagonism (A, log)
4. Dosage (P, lin) vs dosage (A, lin) - v low antagonism
5. Dosage (P, lin) vs dosage (A, lin) - low antagonism
6. Dosage (P, lin) vs dosage (A, lin) - medium antagonism
7. Dosage (P, lin) vs dosage (A, lin) - high antagonism


Nonlinear antag, positive feedback on both

8. Positive feedback (both, lin) vs on rate (both, lin)
9. Positive feedback (both, lin) vs off rate (both, lin)
# 16. Positive feedback (both, lin) vs antagonism (both, log)


12. Dosage (P, lin) vs positive feedback (both, lin) - v low antagonism
13. Dosage (P, lin) vs positive feedback (both, lin) - low antagonism
14. Dosage (P, lin) vs positive feedback (both, lin) - medium antagonism
15. Dosage (P, lin) vs positive feedback (both, lin) - high antagonism

8. Dosage (P, lin) vs antagonism (both, log) - low positive feedback
9. Dosage (P, lin) vs antagonism (both, log) - medium positive feedback
10. Dosage (P, lin) vs antagonism (both, log) - high positive feedback
11. Dosage (P, lin) vs antagonism (both, log) - v. high positive feedback





"""

# Generic parameter set

BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

BaseU = PARu(konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01, ePneg=2, eAneg=2,
             psi=0.1, pA=1, pP=1)

"""
Nonlinear antag

"""

# 1. Dosage (P, lin) vs antagonism (both, log)
if int(sys.argv[1]) == 1:
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

# 2. Dosage (P, lin) vs antagonism (P, log)
if int(sys.argv[1]) == 2:
    p1_range = (0, 2)
    p2_range = (-4, 0)


    def spatial(pP, k):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        m.initiate()
        m.run(kill_uni=True, kill_stab=False)
        return m.stateB()


    def uniform(pP, k):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.kAP = 10 ** k
        return m.bistability_instability()

# 3. Antagonism (P, log) vs antagonism (A, log)
if int(sys.argv[1]) == 3:
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

# 4. Dosage (P, lin) vs dosage (A, lin) - v low antagonism
if int(sys.argv[1]) == 4:
    BaseS.kAP = 0.0001
    BaseS.kPA = 0.0001
    BaseU.kAP = 0.0001
    BaseU.kPA = 0.0001
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

# 5. Dosage (P, lin) vs dosage (A, lin) - low antagonism
if int(sys.argv[1]) == 5:
    BaseS.kAP = 0.001
    BaseS.kPA = 0.001
    BaseU.kAP = 0.001
    BaseU.kPA = 0.001
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

# 6. Dosage (P, lin) vs dosage (A, lin) - medium antagonism
if int(sys.argv[1]) == 6:
    BaseS.kAP = 0.01
    BaseS.kPA = 0.01
    BaseU.kAP = 0.01
    BaseU.kPA = 0.01
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

# 7. Dosage (P, lin) vs dosage (A, lin) - high antagonism
if int(sys.argv[1]) == 7:
    BaseS.kAP = 0.1
    BaseS.kPA = 0.1
    BaseU.kAP = 0.1
    BaseU.kPA = 0.1
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

"""
Nonlinear antag, positive feedback on both

"""

# 8. Positive feedback (both, lin) vs on rate (both, lin)
if int(sys.argv[1]) == 8:
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

# 9. Positive feedback (both, lin) vs off rate (both, lin)
if int(sys.argv[1]) == 9:
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

# # 12. Dosage (P, lin) vs positive feedback (both, lin) - v low antagonism
# if int(sys.argv[1]) == 7:
#     BaseS.kAP = 0.0001
#     BaseS.kPA = 0.0001
#     BaseU.kAP = 0.0001
#     BaseU.kPA = 0.0001
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kpos):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseS.konP - (kpos * y0)
#         m.konA = BaseS.konA - (kpos * y0)
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kpos):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseU.konP - (kpos * y0)
#         m.konA = BaseU.konA - (kpos * y0)
#         return m.bistability_instability()
#
#
# # 13. Dosage (P, lin) vs positive feedback (both, lin) - low antagonism
# if int(sys.argv[1]) == 7:
#     BaseS.kAP = 0.001
#     BaseS.kPA = 0.001
#     BaseU.kAP = 0.001
#     BaseU.kPA = 0.001
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kpos):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseS.konP - (kpos * y0)
#         m.konA = BaseS.konA - (kpos * y0)
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kpos):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseU.konP - (kpos * y0)
#         m.konA = BaseU.konA - (kpos * y0)
#         return m.bistability_instability()
#
# # 14. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
# if int(sys.argv[1]) == 8:
#     BaseS.kAP = 0.01
#     BaseS.kPA = 0.01
#     BaseU.kAP = 0.01
#     BaseU.kPA = 0.01
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kpos):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseS.konP - (kpos * y0)
#         m.konA = BaseS.konA - (kpos * y0)
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kpos):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseU.konP - (kpos * y0)
#         m.konA = BaseU.konA - (kpos * y0)
#         return m.bistability_instability()
#
# # 15. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
# if int(sys.argv[1]) == 9:
#     BaseS.kAP = 0.1
#     BaseS.kPA = 0.1
#     BaseU.kAP = 0.1
#     BaseU.kPA = 0.1
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kpos):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseS.konP - (kpos * y0)
#         m.konA = BaseS.konA - (kpos * y0)
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kpos):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.kposP = kpos
#         m.kposA = kpos
#         m.konP = BaseU.konP - (kpos * y0)
#         m.konA = BaseU.konA - (kpos * y0)
#         return m.bistability_instability()

# """
# PAR model (generic, nonlinear antag, positive feedback)
#
# """
#
# # 7. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
# if int(sys.argv[1]) == 7:
#     BaseS.kAP = 0.001
#     BaseS.kPA = 0.001
#     BaseU.kAP = 0.001
#     BaseU.kPA = 0.001
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 8. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
# if int(sys.argv[1]) == 8:
#     BaseS.kAP = 0.01
#     BaseS.kPA = 0.01
#     BaseU.kAP = 0.01
#     BaseU.kPA = 0.01
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 9. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
# if int(sys.argv[1]) == 9:
#     BaseS.kAP = 0.1
#     BaseS.kPA = 0.1
#     BaseU.kAP = 0.1
#     BaseU.kPA = 0.1
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pP, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 21. Dosage (A, lin) vs positive feedback (P, lin) - low antagonism
# if int(sys.argv[1]) == 21:
#     BaseS.kAP = 0.001
#     BaseS.kPA = 0.001
#     BaseU.kAP = 0.001
#     BaseU.kPA = 0.001
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pA, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pA, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 22. Dosage (A, lin) vs positive feedback (P, lin) - medium antagonism
# if int(sys.argv[1]) == 22:
#     BaseS.kAP = 0.01
#     BaseS.kPA = 0.01
#     BaseU.kAP = 0.01
#     BaseU.kPA = 0.01
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pA, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pA, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 23. Dosage (A, lin) vs positive feedback (P, lin) - high antagonism
# if int(sys.argv[1]) == 23:
#     BaseS.kAP = 0.1
#     BaseS.kPA = 0.1
#     BaseU.kAP = 0.1
#     BaseU.kPA = 0.1
#     p1_range = (0, 2)
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     # p2_range = (0, BaseS.konP / y0)
#     p2_range = (0, 0.02)
#
#
#     def spatial(pA, kposP):
#         m = copy.deepcopy(BaseS)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pA, kposP):
#         m = copy.deepcopy(BaseU)
#         m.pA = pA
#         m.konP = m.konP - (kposP * y0)
#         m.kposP = kposP
#         return m.bistability_instability()
#
# # 10. Dosage (P, lin) vs dosage A (lin) - low positive feedback
# if int(sys.argv[1]) == 10:
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     BaseS.kposP = 0.005
#     BaseU.kposP = 0.005
#     BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
#     BaseU.konP = BaseU.konP - (BaseU.kposP * y0)
#     p1_range = (0, 2)
#     p2_range = (0, 2)
#
#
#     def spatial(pP, pA):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.pA = pA
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, pA):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.pA = pA
#         return m.bistability_instability()
#
# # 11. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
# if int(sys.argv[1]) == 11:
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     BaseS.kposP = 0.01
#     BaseU.kposP = 0.01
#     BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
#     BaseU.konP = BaseU.konP - (BaseU.kposP * y0)
#     p1_range = (0, 2)
#     p2_range = (0, 2)
#
#
#     def spatial(pP, pA):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.pA = pA
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, pA):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.pA = pA
#         return m.bistability_instability()
#
# # 12. Dosage (P, lin) vs dosage A (lin) - high positive feedback
# if int(sys.argv[1]) == 12:
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     BaseS.kposP = 0.015
#     BaseU.kposP = 0.015
#     BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
#     BaseU.konP = BaseU.konP - (BaseU.kposP * y0)
#     p1_range = (0, 2)
#     p2_range = (0, 2)
#
#
#     def spatial(pP, pA):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.pA = pA
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, pA):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.pA = pA
#         return m.bistability_instability()
#
# # 19. Dosage (P, lin) vs dosage A (lin) - v. high positive feedback
# if int(sys.argv[1]) == 19:
#     y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
#     BaseS.kposP = 0.018
#     BaseU.kposP = 0.018
#     BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
#     BaseU.konP = BaseU.konP - (BaseU.kposP * y0)
#     p1_range = (0, 2)
#     p2_range = (0, 2)
#
#
#     def spatial(pP, pA):
#         m = copy.deepcopy(BaseS)
#         m.pP = pP
#         m.pA = pA
#         m.initiate()
#         m.run(kill_uni=True, kill_stab=False)
#         return m.state()
#
#
#     def uniform(pP, pA):
#         m = copy.deepcopy(BaseU)
#         m.pP = pP
#         m.pA = pA
#         return m.bistability_instability()

###############################################################################################

# Bifurcation2D(uniform, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=50, resolution_step=2,
#               n_iterations=6, direc=save_direc + sys.argv[1] + '/Uniform', parallel=True, crange=[1, 6]).run()

ParamSpaceQual2D(spatial, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=10, resolution_step=2,
                 n_iterations=8, direc=save_direc + sys.argv[1] + '/Spatial', parallel=False, crange=[1, 5]).run()
