import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
A

"""

d_params = [0.01, 0.05, 1]
for i, d_param in enumerate(d_params):
    ax = plt.subplot2grid((len(d_params), 1), (i, 0))
    ax.imshow(np.tile(np.loadtxt('Res3/%s_%s/m.txt' % ('A', d_param)), (20, 1)), cmap='Greens', vmin=0, vmax=5.4,
              interpolation='bicubic')
    # ax.set_xlabel('Position (μm)')
    ax.set_ylabel('D = %s' % d_param)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.savefig('A_Final.svg')

"""
B

"""

k_params = [0.007, 0.01, 0.1]

for i, k_param in enumerate(k_params):
    ax = plt.subplot2grid((len(k_params), 1), (i, 0))
    ax.imshow(np.tile(np.loadtxt('Res2/%s_%s/m.txt' % ('B', k_param)), (20, 1)), cmap='Greens', vmin=0, vmax=3.2,
              interpolation='bicubic')
    # ax.set_xlabel('Position (μm)')
    ax.set_ylabel('k = %s' % k_param)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.savefig('B_Final.svg')

# for i, k_param in enumerate(k_params):
#     ax = plt.subplot2grid((len(k_params), 1), (i, 0))
#     ax.imshow(np.tile(np.loadtxt('Res2/%s_%s/c.txt' % ('B', k_param)), (20, 1)), cmap='Greens', vmin=0, vmax=1.5,
#               interpolation='bicubic')
#     # ax.set_xlabel('Position (μm)')
#     ax.set_ylabel('k = %s' % k_param)
#     ax.set_yticks([])
#     ax.set_xticks([])
# sns.despine()
# plt.show()

"""
C

"""

k_params = [0.001, 0.01, 0.1]

plt.rcParams.update({'mathtext.default': 'regular'})
for i, k_param in enumerate(k_params):
    ax = plt.subplot2grid((len(k_params), 1), (i, 0))
    ax.imshow(np.tile(np.loadtxt('Res2/%s_%s/m.txt' % ('C', k_param)), (20, 1)), cmap='Greens', vmin=0, vmax=5.4,
              interpolation='bicubic')
    # ax.set_xlabel('Position (μm)')
    ax.set_ylabel('$k_{recy}$ = %s' % k_param)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.savefig('C_final.svg')

# # for i, k_param in enumerate(k_params):
# #     ax = plt.subplot2grid((len(k_params), 1), (i, 0))
# #     ax.imshow(np.tile(np.loadtxt('Res2/%s_%s/c.txt' % ('C', k_param)), (20, 1)), cmap='Greens', vmin=0, vmax=1.5,
# #               interpolation='bicubic')
# #     # ax.set_xlabel('Position (μm)')
# #     ax.set_ylabel('k = %s' % k_param)
# #     ax.set_yticks([])
# #     ax.set_xticks([])
# # sns.despine()
# # plt.show()
