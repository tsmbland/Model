import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Funcs import AntagSweep, Model
import numpy as np

print(sys.argv[1])

"""
1: Generic, linear

"""

if int(sys.argv[1]) == 1:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=1,
                       psi=0.1, pA=1, pP=1)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/1')
    a.sweep()

"""
2: Generic, non-linear feedback

"""

if int(sys.argv[1]) == 2:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2,
                       psi=0.1, pA=1, pP=1)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=8,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/2')
    a.sweep()

"""
3: Generic, non-linear on (weak)

"""

if int(sys.argv[1]) == 3:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.02, kAP=0, kPA=0, ePneg=1, eAneg=1,
                       psi=0.1, pA=1, pP=1)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/3')
    a.sweep()

"""
4: Generic, non-linear on (strong)

"""

if int(sys.argv[1]) == 4:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.05, kAP=0, kPA=0, ePneg=1, eAneg=1,
                       psi=0.1, pA=1, pP=1.02)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/4')
    a.sweep()

"""
5: Generic, both (weak)

"""

if int(sys.argv[1]) == 5:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.02, kAP=0, kPA=0, ePneg=1, eAneg=2,
                       psi=0.1, pA=1, pP=1)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/5')
    a.sweep()

"""
6: Generic, both (strong)

"""

if int(sys.argv[1]) == 6:
    base_model = Model(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.05, kAP=0, kPA=0, ePneg=1, eAneg=2,
                       psi=0.1, pA=1, pP=1)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/6')
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

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP, koffP=koffP, kposP=0, kAP=0,
                       kPA=0, ePneg=1, eAneg=1, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/7')
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

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP,
                       koffP=koffP, kposP=0, kAP=0, kPA=0, ePneg=1, eAneg=2, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/8')
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

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP, koffP=koffP,
                       kposP=8.6578099573 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/9')
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

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP, koffP=koffP,
                       kposP=8.6578099573 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-3, 0), p2_range=(-3, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/10')
    a.sweep()

"""
11: PAR, non-linear on (odr)

"""

if int(sys.argv[1]) == 11:
    koffA = 0.0092
    koffP = 0.0073
    psi = 0.10318684114244771
    dosP = 0.294005475663175
    dosA = 1.05143336288 * dosP

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP, koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=1, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-4, 0), p2_range=(-4, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/11')
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

    base_model = Model(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP, koffP=koffP,
                       kposP=12.7357711156 * koffP, kAP=0, kPA=0, ePneg=1, eAneg=2, psi=psi, pA=dosA, pP=dosP)

    a = AntagSweep(base_model, xrange=(0, 10), yrange=(0, 10), p1_range=(-4, 0), p2_range=(-4, 0), cores=32,
                   resolution=2,
                   n_iterations=10, direc=home_direc + '/Res/12')
    a.sweep()
