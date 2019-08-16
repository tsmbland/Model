import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
Diffusion

"""

d_params = [0.1, 1, 10]
koff_base = 0.001

for i, d_param in enumerate(d_params):
    plt.plot(np.loadtxt('Res/%s_%s/m.txt' % (d_param, koff_base)), label=d_params[i])
sns.despine()
plt.legend(frameon=False, title='D')
plt.ylim(bottom=0)
plt.xlabel('Position')
plt.ylabel('Cortical concentration')
plt.show()

"""
Off rates

"""
koff_params = [0.0001, 0.001, 0.01]
d_base = 1

for i, koff_param in enumerate(koff_params):
    plt.plot(np.loadtxt('Res/%s_%s/m.txt' % (d_base, koff_param)), label=koff_params[i])
sns.despine()
plt.legend(frameon=False, title='koff')
plt.ylim(bottom=0)
plt.xlabel('Position')
plt.ylabel('Cortical concentration')
plt.show()
