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

direc = x.ddirec + 's1904__m0000_optimise_B'
db = pd.read_csv(direc + '/db.csv')
params = ['konA', 'konP', 'kAP', 'kPA']

# """
# MSE hist
#
# """
#
# x = np.array(db['MSE'])
# x = x[~np.isnan(x)]
# plt.hist(x, bins=np.linspace(0, 5, 500))
# plt.show()

"""
Param choice

"""

thresh = 0.2
array = db.loc[db['MSE'] < thresh]['kAP']
array2 = db.loc[db['MSE'] < thresh]['kPA']
plt.hist(array / array2)
plt.show()
