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
import glob

direc = x.ddirec + 's1904__m0000_robustness'

print('g')

direcs = glob.glob('%s/*/' % direc)
print(len(direcs))
params = ['eAneg', 'kAP', 'kPA']
db = []

print('a')

for i, d in enumerate(direcs):
    print(100 * i / len(direcs))
    score = np.loadtxt(d + '/score.txt')
    p = pickle.load(open('%s/_params.pkl' % d, 'rb'))

    dic = {'Directory': d, 'Score': score}
    for param in params:
        dic[param] = p[param]
    db.append(dic)

pd.DataFrame(db).to_csv('db.csv')
