import sys
import os

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../')
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '../../..')
import M as x
import Models.m0000 as m
import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import shutil
import copy
from matplotlib import animation
import pandas as pd
from random import sample

direc = x.ddirec + 's1904__m0000_optimise_B/'

tmax = 1000
psets = 10

db = pd.read_csv(direc + '/db.csv')
samples = sample(list((db['Directory'][db['MSE'] < 0.5])), psets)
print(samples)

"""
WT

"""

plt.close()

for pset in samples:
    split = pset.split('/')
    params = pickle.load(open('%s/%s/%s/_params.pkl' % (direc, split[-2], split[-1]), 'rb'))
    model = m.Model(**dict(copy.deepcopy(params)))
    model.Tmax = tmax
    asis = []
    times = []
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat
        asis.append(abs(x.asi(model.pm)))
        times.append(model.time)
    plt.plot(np.array(times), np.array(asis))

plt.xlabel('Time (seconds)')
plt.ylabel('pPAR ASI')
sns.despine()
plt.ylim(-0.05, 0.55)
plt.savefig('asi_wt.png')

"""
kAP = 0

"""

plt.close()

for pset in samples:
    split = pset.split('/')
    params = pickle.load(open('%s/%s/%s/_params.pkl' % (direc, split[-2], split[-1]), 'rb'))
    model = m.Model(**dict(copy.deepcopy(params)))
    model.Tmax = tmax
    model.kAP = 0
    asis = [abs(x.asi(model.pm))]
    times = [0]
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat
        asis.append(abs(x.asi(model.pm)))
        times.append(model.time)
    plt.plot(np.array(times), np.array(asis))

plt.xlabel('Time (seconds)')
plt.ylabel('pPAR ASI')
sns.despine()
plt.ylim(-0.05, 0.55)
plt.savefig('asi_kAP.png')

"""
kPA = 0

"""

plt.close()

for pset in samples:
    split = pset.split('/')
    params = pickle.load(open('%s/%s/%s/_params.pkl' % (direc, split[-2], split[-1]), 'rb'))
    model = m.Model(**dict(copy.deepcopy(params)))
    model.Tmax = tmax
    model.kPA = 0
    asis = [abs(x.asi(model.pm))]
    times = [0]
    for t in range(int(model.Tmax / model.deltat)):
        model.react()
        model.time = (t + 1) * model.deltat
        asis.append(abs(x.asi(model.pm)))
        times.append(model.time)
    plt.plot(np.array(times), np.array(asis))

plt.xlabel('Time (seconds)')
plt.ylabel('pPAR ASI')
sns.despine()
plt.ylim(-0.05, 0.55)
plt.savefig('asi_kPA.png')
