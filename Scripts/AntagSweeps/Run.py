import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/AntagSweeps/'

from Models.SimplePosFeedback import Model
from Funcs import ParamSweep
import numpy as np

print(sys.argv[1])

"""
1: Generic, linear

"""

if int(sys.argv[1]) == 1:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/1')
    a.sweep()

"""
2: Generic, non-linear feedback

"""

if int(sys.argv[1]) == 2:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/2')
    a.sweep()

"""
3: Generic, non-linear on (weak)

"""

if int(sys.argv[1]) == 3:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0.02, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/3')
    a.sweep()

"""
4: Generic, non-linear on (strong)

"""

if int(sys.argv[1]) == 4:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0.05, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/4')
    a.sweep()

"""
5: Generic, both (weak)

"""

if int(sys.argv[1]) == 5:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0.02, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/5')
    a.sweep()

"""
6: Generic, both (strong)

"""

if int(sys.argv[1]) == 6:
    base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                       koffP=0.01, kposP=0.05, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1.02)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/6')
    a.sweep()

"""
7: PAR, linear

"""

if int(sys.argv[1]) == 7:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                       koffP=koffP,
                       kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01, deltax=0.5,
                       psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/7')
    a.sweep()

"""
8: PAR, non-linear feedback

"""

if int(sys.argv[1]) == 8:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                       koffP=koffP,
                       kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01, deltax=0.5,
                       psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/8')
    a.sweep()

"""
9: PAR, non-linear on (ols)

"""

if int(sys.argv[1]) == 9:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP,
                       koffP=koffP,
                       kposP=8.6578099573 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/9')
    a.sweep()

"""
10: PAR, both (ols)

"""

if int(sys.argv[1]) == 10:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP,
                       koffP=koffP,
                       kposP=8.6578099573 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/10')
    a.sweep()

"""
11: PAR, non-linear on (odr), wide

"""

if int(sys.argv[1]) == 11:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/11')
    a.sweep()

"""
12: PAR, both (odr)

"""

if int(sys.argv[1]) == 12:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=1000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/12')
    a.sweep()

"""
13: PAR, both (odr), long, wide

"""

if int(sys.argv[1]) == 13:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/13')
    a.sweep()

"""
14: PAR, non-linear on (odr), long, wide

"""

if int(sys.argv[1]) == 14:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=10000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/14')
    a.sweep()

"""
15: PAR, non-linear on (odr), long

"""

if int(sys.argv[1]) == 15:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=10000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/15')
    a.sweep()

"""
16: PAR, linear, long, wide

"""

if int(sys.argv[1]) == 16:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                       koffP=koffP,
                       kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=10000, deltat=0.01, deltax=0.5,
                       psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/16')
    a.sweep()

"""
17: PAR (RING mutant), linear, long, wide

"""

if int(sys.argv[1]) == 17:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=10.436114976179985 * koffP,
                       koffP=koffP,
                       kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=1, xsteps=100, Tmax=10000, deltat=0.01, deltax=0.5,
                       psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/17')
    a.sweep()

"""
18: PAR, non-linear feedback, long

"""

if int(sys.argv[1]) == 18:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                       koffP=koffP,
                       kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01, deltax=0.5,
                       psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/18')
    a.sweep()

"""
19: PAR, both (odr), long

"""

if int(sys.argv[1]) == 19:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
                       koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, xsteps=100, Tmax=10000, deltat=0.01,
                       deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)

    a = ParamSweep(base_model, p1='kAP', p2='kPA', p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution=2,
                   n_iterations=10, direc=save_direc + '/19')
    a.sweep()
