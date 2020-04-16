import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/AntagSweeps/'

from Models.PDE.PAR import PAR as PARs
from Models.ODE.PAR import PAR as PARu
from Funcs import Bifurcation2D
import numpy as np
import copy

print(sys.argv[1])

"""

PAR model (generic, nonlinear antag)

1. Dosage (P, lin) vs antagonism (both, log)
2. Dosage (P, lin) vs antagonism (P, log)
3. Antagonism (P, log) vs antagonism (A, log)
4. Dosage (P, lin) vs dosage (A, lin) - low antagonism
5. Dosage (P, lin) vs dosage (A, lin) - medium antagonism
6. Dosage (P, lin) vs dosage (A, lin) - high antagonism


PAR model (generic, nonlinear antag, positive feedback)

7. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
8. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
9. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
10. Dosage (P, lin) vs dosage A (lin) - low positive feedback
11. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
12. Dosage (P, lin) vs dosage A (lin) - high positive feedback


PAR model (generic, linear antag, positive feedback)

13. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
14. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
15. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
16. Dosage (P, lin) vs dosage A (lin) - low positive feedback
17. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
18. Dosage (P, lin) vs dosage A (lin) - high positive feedback


"""

# Generic parameter set

BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

BaseU = PARu(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01, ePneg=2, eAneg=2,
             psi=0.1, pA=1, pP=1)

print(sys.argv[1])

"""
PAR model (generic, nonlinear antag)

"""

# 1. Dosage (P, lin) vs antagonism (both, log)
if int(sys.argv[1]) == 1:
    p1_range = (0, 2)
    p2_range = (-3, 0)


    def spatial(pP, k):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        return m.initiate().run().polarised()


    def uniform(pP, k):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.kAP = 10 ** k
        m.kPA = 10 ** k
        return m.bistability_instability()

# 2. Dosage (P, lin) vs antagonism (P, log)
if int(sys.argv[1]) == 2:
    p1_range = (0, 2)
    p2_range = (-3, 0)


    def spatial(pP, k):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.kAP = 10 ** k
        return m.initiate().run().polarised()


    def uniform(pP, k):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.kAP = 10 ** k
        return m.bistability_instability()

# 3. Antagonism (P, log) vs antagonism (A, log)
if int(sys.argv[1]) == 3:
    p1_range = (-3, 0)
    p2_range = (-3, 0)


    def spatial(kAP, kPA):
        m = copy.deepcopy(BaseS)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        return m.initiate().run().polarised()


    def uniform(kAP, kPA):
        m = copy.deepcopy(BaseU)
        m.kAP = 10 ** kAP
        m.kPA = 10 ** kPA
        return m.bistability_instability()

# 4. Dosage (P, lin) vs dosage (A, lin) - low antagonism
if int(sys.argv[1]) == 4:
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
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 5. Dosage (P, lin) vs dosage (A, lin) - medium antagonism
if int(sys.argv[1]) == 5:
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
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 6. Dosage (P, lin) vs dosage (A, lin) - high antagonism
if int(sys.argv[1]) == 6:
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
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

"""
PAR model (generic, nonlinear antag, positive feedback)

"""

# 7. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
if int(sys.argv[1]) == 7:
    BaseS.kAP = 0.001
    BaseS.kPA = 0.001
    BaseU.kAP = 0.001
    BaseU.kPA = 0.001
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 8. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
if int(sys.argv[1]) == 8:
    BaseS.kAP = 0.01
    BaseS.kPA = 0.01
    BaseU.kAP = 0.01
    BaseU.kPA = 0.01
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 9. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
if int(sys.argv[1]) == 9:
    BaseS.kAP = 0.1
    BaseS.kPA = 0.1
    BaseU.kAP = 0.1
    BaseU.kPA = 0.1
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 10. Dosage (P, lin) vs dosage A (lin) - low positive feedback
if int(sys.argv[1]) == 10:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.005
    BaseU.kposP = 0.005
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 11. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
if int(sys.argv[1]) == 11:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.01
    BaseU.kposP = 0.01
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 12. Dosage (P, lin) vs dosage A (lin) - high positive feedback
if int(sys.argv[1]) == 12:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.015
    BaseU.kposP = 0.015
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

"""
PAR model (generic, linear antag, positive feedback)

"""

if int(sys.argv[1]) in [13, 14, 15, 16, 17, 18]:
    BaseS.eAneg = 1
    BaseS.ePneg = 1
    BaseU.eAneg = 1
    BaseU.ePneg = 1

# 13. Dosage (P, lin) vs positive feedback (P, lin) - low antagonism
if int(sys.argv[1]) == 13:
    BaseS.kAP = 0.001
    BaseS.kPA = 0.001
    BaseU.kAP = 0.001
    BaseU.kPA = 0.001
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 14. Dosage (P, lin) vs positive feedback (P, lin) - medium antagonism
if int(sys.argv[1]) == 14:
    BaseS.kAP = 0.01
    BaseS.kPA = 0.01
    BaseU.kAP = 0.01
    BaseU.kPA = 0.01
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 15. Dosage (P, lin) vs positive feedback (P, lin) - high antagonism
if int(sys.argv[1]) == 15:
    BaseS.kAP = 0.1
    BaseS.kPA = 0.1
    BaseU.kAP = 0.1
    BaseU.kPA = 0.1
    p1_range = (0, 2)
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    p2_range = (0, BaseS.konP / y0)


    def spatial(pP, kposP):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.initiate().run().polarised()


    def uniform(pP, kposP):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        y0 = (m.konP * m.pP) / (m.psi * m.konP + m.koffP)
        m.konP = m.konP - (kposP * y0)
        m.kposP = kposP
        return m.bistability_instability()

# 16. Dosage (P, lin) vs dosage A (lin) - low positive feedback
if int(sys.argv[1]) == 16:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.005
    BaseU.kposP = 0.005
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 17. Dosage (P, lin) vs dosage A (lin) - medium positive feedback
if int(sys.argv[1]) == 17:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.01
    BaseU.kposP = 0.01
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

# 18. Dosage (P, lin) vs dosage A (lin) - high positive feedback
if int(sys.argv[1]) == 18:
    y0 = (BaseS.konP * BaseS.pP) / (BaseS.psi * BaseS.konP + BaseS.koffP)
    BaseS.kposP = 0.015
    BaseU.kposP = 0.015
    BaseS.konP = BaseS.konP - (BaseS.kposP * y0)
    BaseU.konP = BaseS.konP - (BaseS.kposP * y0)
    p1_range = (0, 2)
    p2_range = (0, 2)


    def spatial(pP, pA):
        m = copy.deepcopy(BaseS)
        m.pP = pP
        m.pA = pA
        return m.initiate().run().polarised()


    def uniform(pP, pA):
        m = copy.deepcopy(BaseU)
        m.pP = pP
        m.pA = pA
        return m.bistability_instability()

###############################################################################################


Bifurcation2D(uniform, p1_range=p1_range, p2_range=p2_range, cores=4, resolution0=50, resolution_step=2, n_iterations=1,
              direc=save_direc + (sys.argv[1]) + '/Uniform', parallel=True, crange=[1, 6]).run()

Bifurcation2D(spatial, p1_range=p1_range, p2_range=p2_range, cores=4, resolution0=50, resolution_step=2, n_iterations=1,
              direc=save_direc + (sys.argv[1]) + '/Spatial', parallel=True, crange=[1, 3]).run()
