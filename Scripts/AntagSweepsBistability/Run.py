import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')
save_direc = home_direc + '/../../../../ModelData/AntagSweepsBistability/'

from Funcs import ParamSpaceQual2D

# print(sys.argv[1])
#
# """
# 1: Generic, linear
#
# """
#
# if int(sys.argv[1]) == 1:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=1, eAneg=1, psi=0.1, pA=1, pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100, resolution=2,
#                           n_iterations=9, direc=save_direc + '/1')
#
# """
# 2: Generic, non-linear feedback
#
# """
#
# if int(sys.argv[1]) == 2:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0, ePneg=1, eAneg=2, psi=0.1, pA=1, pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/2')
#
# """
# 3: Generic, non-linear on (weak)
#
# """
#
# if int(sys.argv[1]) == 3:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.02, ePneg=1, eAneg=1, psi=0.1, pA=1,
#                     pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/3')
#
# """
# 4: Generic, non-linear on (strong)
#
# """
#
# if int(sys.argv[1]) == 4:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.05, ePneg=1, eAneg=1, psi=0.1, pA=1,
#                     pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/4')
#
# """
# 5: Generic, both (weak)
#
# """
#
# if int(sys.argv[1]) == 5:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.02, ePneg=1, eAneg=2, psi=0.1, pA=1,
#                     pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/5')
#
# """
# 6: Generic, both (strong)
#
# """
#
# if int(sys.argv[1]) == 6:
#     b = Bistability(konA=0.1, koffA=0.01, kposA=0, konP=0.1, koffP=0.01, kposP=0.05, ePneg=1, eAneg=2, psi=0.1, pA=1,
#                     pP=1)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/6')
#
# """
# 7: PAR, linear
#
# """
#
# if int(sys.argv[1]) == 7:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP, koffP=koffP, kposP=0, ePneg=1,
#                     eAneg=1, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/7')
#
# """
# 8: PAR, non-linear feedback
#
# """
#
# if int(sys.argv[1]) == 8:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=70.623414063 * koffP, koffP=koffP, kposP=0, ePneg=1,
#                     eAneg=2, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/8')
#
# """
# 9: PAR, non-linear on (ols)
#
# """
#
# if int(sys.argv[1]) == 9:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP, koffP=koffP,
#                     kposP=8.6578099573 * koffP, ePneg=1, eAneg=1, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/9')
#
# """
# 10: PAR, both (ols)
#
# """
#
# if int(sys.argv[1]) == 10:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP, koffP=koffP,
#                     kposP=8.6578099573 * koffP, ePneg=1, eAneg=2, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-3, 0), p2_range=(-3, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/10')
#
# """
# 11: PAR, non-linear on (odr), wide
#
# """
#
# if int(sys.argv[1]) == 11:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP, koffP=koffP,
#                     kposP=12.7357711156 * koffP, ePneg=1, eAneg=1, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/11')
#
# """
# 12: PAR, both (odr)
#
# """
#
# if int(sys.argv[1]) == 12:
#     koffA = 0.0092
#     koffP = 0.0073
#     psi = 0.10318684114244771
#     dosP = 0.294005475663175
#     dosA = 1.05143336288 * dosP
#
#     b = Bistability(konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP, koffP=koffP,
#                     kposP=12.7357711156 * koffP, ePneg=1, eAneg=2, psi=psi, pA=dosA, pP=dosP)
#     a = ParamSweepGeneral(b.bistability, p1_range=(-4, 0), p2_range=(-4, 0), log=True, cores=32, resolution0=100,
#                           resolution=2, n_iterations=9, direc=save_direc + '/12')

# a.run()


