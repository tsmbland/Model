"""
Runs original model with original parameters, saves results for comparison to other models
180729
m0000

"""

import M as x
import Models.m0000 as m0000

x.datadirec = '../../ModelData'

p0 = m0000.Params(Da=0.28, Dp=0.15, konA=0.0085, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, ePneg=1,
                  eAneg=2,
                  pA=1.56, pP=1, L=67.3, xsteps=500, psi=0.174, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5,
                  Peqmax=1)

x.alg_singlesim(m=m0000.Model(p0), jobid=9999, subjobid=0, simid=0, compression=0)
