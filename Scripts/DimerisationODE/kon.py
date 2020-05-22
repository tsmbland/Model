import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Models.ODE.DM import Model
import numpy as np
import copy
from Funcs import ParamSpace2D
from scipy.optimize import curve_fit
import scipy.odr as odr

print(sys.argv[1])

save_direc = home_direc + '/../../../../ModelData/Dimerisation/kon/'

# Base parameter set

BaseModel = Model(kon=1, kon2=1, koff=1, kd_f=1, kd_b=1, psi=0.174, dosage=1)

"""
Funcs

"""


def pf_profile(x, kon0n, E):
    psi = 0.174
    a = kon0n * (1 - E)
    b = E * (psi * kon0n + 1)
    y = (a * x) / (1 - (b * x))
    return y


def fitting(cyt, mem):
    # Ols fit
    popt, pcov = curve_fit(pf_profile, cyt, mem)
    a_min_0 = popt[0]
    b_min_0 = popt[1]

    # Odr fit
    def perform_odr(x, y):
        quadr = odr.Model(lambda B, x: pf_profile(x, B[0], B[1]))
        mydata = odr.Data(x, y)
        myodr = odr.ODR(mydata, quadr, beta0=[a_min_0, b_min_0])

        output = myodr.run()
        return output

    regression = perform_odr(cyt, mem)
    a_min = regression.beta[0]
    b_min = regression.beta[1]
    return a_min, b_min


def rundown(model):
    # Create results containers
    pm1 = np.zeros([100])
    pm2s = np.zeros([100])
    pm2d = np.zeros([100])
    pc1 = np.zeros([100])
    pc2 = np.zeros([100])

    # Perform simulations
    for i, d in enumerate(np.linspace(0, 1, 100)):
        m = copy.deepcopy(model)
        m.dosage = d

        sol = m.solve()
        pm1[i] = sol[0]
        pm2s[i] = sol[1]
        pm2d[i] = sol[2]
        pc1[i] = sol[3]
        pc2[i] = sol[4]

    # Compile
    cyt = pc1 + 2 * pc2
    mem = pm1 * 2 * pm2s + 2 * pm2d

    # Fit
    a, b = fitting(cyt, mem)

    return np.log10(a)


#########################################################################

"""
kd_f vs kd_b

"""

if int(sys.argv[1]) == 0:
    def func(kd_f, kd_b):
        m = copy.deepcopy(BaseModel)
        m.kd_f = 10 ** kd_f
        m.kd_b = 10 ** kd_b
        res = rundown(m)
        return res

if int(sys.argv[1]) == 1:
    BaseModel.kon2 = 10


    def func(kd_f, kd_b):
        m = copy.deepcopy(BaseModel)
        m.kd_f = 10 ** kd_f
        m.kd_b = 10 ** kd_b
        res = rundown(m)
        return res

if int(sys.argv[1]) == 2:
    BaseModel.kon2 = 100


    def func(kd_f, kd_b):
        m = copy.deepcopy(BaseModel)
        m.kd_f = 10 ** kd_f
        m.kd_b = 10 ** kd_b
        res = rundown(m)
        return res

"""
kon vs kon2

"""

if int(sys.argv[1]) == 3:
    BaseModel.kd_f = 0.1
    BaseModel.kd_b = 0.1


    def func(kon, kon2):
        m = copy.deepcopy(BaseModel)
        m.kon = 10 ** kon
        m.kon2 = 10 ** kon2
        res = rundown(m)
        return res

if int(sys.argv[1]) == 4:
    BaseModel.kd_f = 1
    BaseModel.kd_b = 1


    def func(kon, kon2):
        m = copy.deepcopy(BaseModel)
        m.kon = 10 ** kon
        m.kon2 = 10 ** kon2
        res = rundown(m)
        return res

if int(sys.argv[1]) == 5:
    BaseModel.kd_f = 10
    BaseModel.kd_b = 10


    def func(kon, kon2):
        m = copy.deepcopy(BaseModel)
        m.kon = 10 ** kon
        m.kon2 = 10 ** kon2
        res = rundown(m)
        return res

###########################################################################


ParamSpace2D(func, p1_range=[-2, 2], p2_range=[-2, 2], direc=save_direc + str(sys.argv[1]), parallel=True,
             cores=32, resolution0=20, save_fig=True, n_iterations=3, resolution_step=2, crange=[-2, 2]).run()
