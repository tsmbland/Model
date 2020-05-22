import sys
import os

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/../..')

from Models.ODE.DM import Model
import numpy as np
import copy
from Funcs import ParamSpace2D
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.odr as odr

# Base parameter set

BaseModel = Model(kon=1, kon2=10, koff=1, kd_f=10, kd_b=10, psi=0.174, dosage=1)

"""
Fitting 1

"""


# def func_wt(x, a, b):
#     y = (a * x) / (1 - (b * x))
#     return y
#
#
# def fitting(cyt, mem):
#     # Ols fit
#     popt, pcov = curve_fit(lambda x, a, b: (a * x) / (1 - (b * x)), cyt, mem)
#     a_min_0 = popt[0]
#     b_min_0 = popt[1]
#
#     # Odr fit
#     def perform_odr(x, y):
#         quadr = odr.Model(lambda B, x: (B[0] * x) / (1 - (B[1] * x)))
#         mydata = odr.Data(x, y)
#         myodr = odr.ODR(mydata, quadr, beta0=[a_min_0, b_min_0])
#         output = myodr.run()
#         return output
#
#     regression = perform_odr(cyt, mem)
#     a_min = regression.beta[0]
#     b_min = regression.beta[1]
#     return a_min, b_min


"""
Fitting 2

"""

psi = 0.174


def func_wt(x, kon0n, E):
    a = kon0n * (1 - E)
    b = E * (psi * kon0n + 1)
    y = (a * x) / (1 - (b * x))
    return y


def fitting(cyt, mem):
    # Ols fit
    popt, pcov = curve_fit(func_wt, cyt, mem)
    a_min_0 = popt[0]
    b_min_0 = popt[1]

    # Odr fit
    def perform_odr(x, y):
        quadr = odr.Model(lambda B, x: func_wt(x, B[0], B[1]))
        mydata = odr.Data(x, y)
        myodr = odr.ODR(mydata, quadr, beta0=[a_min_0, b_min_0])

        output = myodr.run()
        return output

    regression = perform_odr(cyt, mem)
    a_min = regression.beta[0]
    b_min = regression.beta[1]
    return a_min, b_min


"""
Rundown

"""


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

    # Compile results
    cyt = pc1 + 2 * pc2
    mem = pm1 + 2 * pm2d + 2 * pm2s

    # Fit
    a, b = fitting(cyt, mem)
    print(b)

    # Plot results
    plt.plot(cyt, mem)
    plt.plot(cyt, func_wt(cyt, a, b))
    plt.show()


rundown(BaseModel)


# m = copy.deepcopy(BaseModel)
# sol = m.solve()
