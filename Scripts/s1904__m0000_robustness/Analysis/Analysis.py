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

db = pd.read_csv('db.csv')

db2 = db[db['eAneg'] == 3]

print(len(db2['kPA']))

plt.scatter(db2['kPA'][db2['Score'] == 1], db2['kAP'][db2['Score'] == 1], c='r', marker='s', s=10)
plt.scatter(db2['kPA'][db2['Score'] == -1], db2['kAP'][db2['Score'] == -1], c='b', marker='s', s=10)
plt.scatter(db2['kPA'][db2['Score'] == 0], db2['kAP'][db2['Score'] == 0], c='0.8', marker='s', s=10)
# plt.scatter(db2['kPA'], db2['kAP'], marker='s', s=1)
plt.xlim(0, 0.005)
plt.ylim(0, 0.001)
sns.despine()
plt.ylabel('kAP')
plt.xlabel('kPA')
plt.tight_layout()
plt.show()
