import numpy as np
import M as x
import matplotlib.pyplot as plt
import tifffile as tiff

"""


"""

# direc = '/Volumes/lab-goehringn/working/Tom/ModelData/DosageSweeps/2/'
#
# dos = np.linspace(0, 2, 20)
# ant = np.linspace(0, 2, 20)
#
# res = np.zeros([20, 20])
# for i, d in enumerate(dos):
#     for j, a in enumerate(ant):
#
#         name = '%s_%s' % (float('%.4g' % d), float('%.4g' % a))
#         print(name)
#
#         pm = np.loadtxt(direc + '/' + name + '/pm.txt')
#         am = np.loadtxt(direc + '/' + name + '/am.txt')
#
#         if sum(am > pm) == len(am):
#             res[i, j] = 1
#         elif sum(am > pm) == 0:
#             res[i, j] = -1
#         else:
#             res[i, j] = 0
#
# fig, ax = plt.subplots()
# ax.imshow(res, cmap='bwr', origin='lower', alpha=0.5)
# plt.show()

"""


"""

direc = '/Volumes/lab-goehringn/working/Tom/ModelData/DosageSweeps2/8/'

vals = np.linspace(0, 2, 100)

resA = np.zeros([100, 500])
resP = np.zeros([100, 500])

for i, v in enumerate(vals):
    name = '%s' % (float('%.4g' % v))

    resA[i, :] = np.loadtxt(direc + '/' + name + '/am.txt')
    resP[i, :] = np.loadtxt(direc + '/' + name + '/pm.txt')

resA /= np.amax(resA)
resP /= np.amax(resP)

rgb = np.dstack((resA, resP, resP))
plt.imshow(rgb, aspect='auto', origin='lower')
plt.show()

