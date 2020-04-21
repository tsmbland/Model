from Funcs import ParamSpaceQual2D
import copy
from Models.PDE.PAR import PAR as PARs
import matplotlib.pyplot as plt

BaseS = PARs(Da=0.1, Dp=0.1, konA=0.1, koffA=0.0101, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0.01, kPA=0.01,
             ePneg=2, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, L=50, psi=0.1, pA=1, pP=1)

p1_range = (0, 2)
p2_range = (-4, 0)


def spatial(pP, k):
    m = copy.deepcopy(BaseS)
    m.pP = pP
    m.kAP = 10 ** k
    m.initiate()
    m.run(kill_uni=True, kill_stab=False)
    return m.stateB()


# Bifurcation2D(spatial, p1_range=p1_range, p2_range=p2_range, cores=8, resolution0=10, resolution_step=2,
#               n_iterations=8, direc='2/Spatial', parallel=True, crange=[1, 5]).run()


m = copy.deepcopy(BaseS)
m.pP = 1.031250000000
m.kAP = 10 ** -1.604166666667
m.initiate()

# plt.plot(m.A)
# plt.plot(m.P)
# plt.show()

m.run(kill_uni=True, kill_stab=False)

plt.plot(m.A)
plt.plot(m.P)
plt.show()
