import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
A

"""

# d_params = [0.1, 0.5, 10]
#
# for i, d_param in enumerate(d_params):
#     plt.plot(np.linspace(0, 10, 100), np.loadtxt('Res2/%s_%s/m.txt' % ('A', d_param)), label=d_params[i])
# sns.despine()
# plt.legend(frameon=False, title='D (μm²/s)')
# plt.ylim(bottom=0)
# plt.xlabel('Position (μm)')
# plt.ylabel('Cortical concentration (/μm)')
# plt.show()

"""
B

"""

k_params = [0.005, 0.006, 0.007, 0.0075, 0.008, 0.01, 0.02, 0.05, 0.1, 0.5, 1]

for i, k_param in enumerate(k_params):
    plt.plot(np.linspace(0, 10, 100), np.loadtxt('Res2/%s_%s/m.txt' % ('B', k_param)), label=k_params[i])
sns.despine()
plt.legend(frameon=False, title='k (/s)')
plt.ylim(bottom=0)
plt.xlabel('Position (μm)')
plt.ylabel('Cortical concentration (/μm)')
plt.show()

for i, k_param in enumerate(k_params):
    plt.plot(np.linspace(0, 10, 100), np.loadtxt('Res2/%s_%s/c.txt' % ('B', k_param)), label=k_params[i])
sns.despine()
plt.legend(frameon=False, title='k (/s)')
plt.ylim(bottom=0)
plt.xlabel('Position (μm)')
plt.ylabel('Cytoplasmic concentration (/μm)')
plt.show()

"""
C

"""

# k_params = [0.001, 0.01, 0.1]
#
# for i, k_param in enumerate(k_params):
#     plt.plot(np.linspace(0, 10, 100), np.loadtxt('Res2/%s_%s/m.txt' % ('C', k_param)), label=k_params[i])
# sns.despine()
# plt.legend(frameon=False, title='k (/s)')
# plt.ylim(bottom=0)
# plt.xlabel('Position (μm)')
# plt.ylabel('Cortical concentration (/μm)')
# plt.show()

# for i, k_param in enumerate(k_params):
#     plt.plot(np.linspace(0, 10, 100), np.loadtxt('Res2/%s_%s/c.txt' % ('C', k_param)), label=k_params[i])
# sns.despine()
# plt.legend(frameon=False, title='k (/s)')
# plt.ylim(bottom=0)
# plt.xlabel('Position (μm)')
# plt.ylabel('Cytoplasmic concentration (/μm)')
# plt.show()
