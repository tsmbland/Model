import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# print(np.append(np.array(range(1, 100)), [100 - 1]))
# print(np.append([0], np.array(range(100 - 1))))


# x = np.array(range(100)) * (100 / 100)
# v = (x / np.exp(0.08 * (x ** 1)))[::-1]
# v[0] = 0
# v[-1] = 0
#
# plt.plot(v)
# plt.show()
#
# x = np.array(range(100)) * (100 / 100)
# v = (x / np.exp(0.00075 * (x ** 2)))[::-1]
# v[0] = 0
# v[-1] = 0
#
# # plt.plot(1 / (1 + np.exp(0.2 * (50-np.arange(99)))))
# plt.plot(v)
# plt.show()
#
# # - 2 * concs + concs[
# #             np.append([1], np.array(range(len(concs) - 1)))]


# plt.plot(np.linspace(0, 10, 100), 0.1 * np.arange(100) / 100)
# sns.despine()
# plt.xlabel('Position (μm)')
# plt.ylabel('Cortical flow velocity (μm/s)')
# plt.show()


# plt.plot(np.linspace(0, 10, 100), np.r_[0. * np.ones([80]), 1 * np.ones([20])])
# sns.despine()
# plt.xlabel('Position (μm)')
# plt.ylabel('kon (/s)')
# plt.ylim(bottom=0)
# plt.show()
