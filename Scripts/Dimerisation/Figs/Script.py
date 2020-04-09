from Models.pPAR_Dimerisation_NonSpatial import Model
import numpy as np
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil

direc = 'Res'

"""
Run

"""


def run(base_model, kd_f, kd_b, dosages, name):
    print(name)

    # Set up model
    model = copy.deepcopy(base_model)
    model.pc1 = dosages
    model.pc2 = 0 * dosages
    model.pm1 = 0 * dosages
    model.pm2s = 0 * dosages
    model.pm2d = 0 * dosages

    # Set rates
    model.kd_f = kd_f
    model.kd_b = kd_b

    # Simulation
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat

    # Save
    model.pool()

    # Save results
    if os.path.exists(direc + '/' + name):
        shutil.rmtree(direc + '/' + name)
    os.mkdir(direc + '/' + name)
    model.save(direc + '/' + name + '/')


base_model = Model(kon_p=1, kon_p2=10, koff_p=1, kd_f=1, kd_b=1, psi=0.1, Tmax=1000, deltat=0.001,
                   pm1_0=0, pm2s_0=0, pm2d_0=0, pc1_0=0, pc2_0=0)

run(base_model, 10, 10, dosages=np.linspace(0, 1, 100), name='1')
run(base_model, 0.1, 0.1, dosages=np.linspace(0, 1, 100), name='2')
run(base_model, 10, 0.1, dosages=np.linspace(0, 1, 100), name='3')
run(base_model, 0.1, 10, dosages=np.linspace(0, 1, 100), name='4')

"""
Fig 1: mem vs cyt

"""


def func(direc, c):
    mems = np.loadtxt(direc + '/pco.txt')
    cyts = np.loadtxt(direc + '/pcy.txt')
    plt.plot(cyts, mems, c=c)


plt.close()
func('Res/1/', c='r')
func('Res/2/', c='b')
func('Res/3/', c='g')
func('Res/4/', c='y')
plt.xlabel('Cytoplasmic concentration')
plt.ylabel('Cortical concentration')
fig = plt.gcf()
fig.set_size_inches(4, 4)
plt.tight_layout()
sns.despine()
plt.savefig('Fig1.png', dpi=300)

"""
Fig 2: dimerisation vs dosage

"""


def func(direc, c):
    pm1 = np.loadtxt(direc + '/pm1.txt')
    pm2s = np.loadtxt(direc + '/pm2s.txt')
    pm2d = np.loadtxt(direc + '/pm2d.txt')
    pc1 = np.loadtxt(direc + '/pc1.txt')
    pc2 = np.loadtxt(direc + '/pc2.txt')
    pmem_m = pm1
    pmem_d = 2 * pm2s + 2 * pm2d
    pcyt_m = pc1
    pcyt_d = 2 * pc2
    ptot_m = pcyt_m + 0.1 * pmem_m
    ptot_d = pcyt_d + 0.1 * pmem_d

    plt.plot(ptot_m + ptot_d, ptot_d / (ptot_m + ptot_d), c=c)


plt.close()
func('Res/1/', c='r')
func('Res/2/', c='b')
func('Res/3/', c='g')
func('Res/4/', c='y')
plt.xlabel('Dosage')
plt.ylabel('% Dimeric protein')
fig = plt.gcf()
fig.set_size_inches(4, 4)
plt.tight_layout()
sns.despine()
plt.savefig('Fig2.png', dpi=300)
