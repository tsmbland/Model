import numpy as np
import matplotlib.pyplot as plt

n_particles = 10000
max_size = 4
n_bleaches = 100

a = np.random.rand(n_particles * max_size)
bleachvals = np.linspace(0, 0.99, n_bleaches)
totalcounts = np.zeros([max_size, n_bleaches])

for i, bleach in enumerate(bleachvals):
    b = np.array([a > bleach]).astype(int)[0]
    counts = np.zeros([max_size])

    for j in range(0, n_particles * max_size, max_size):
        c = sum(b[j: j + max_size])
        if c != 0:
            counts[c - 1] += 1
    totalcounts[:, i] = counts / sum(counts)


plt.plot(totalcounts.T)
plt.show()

# plt.plot(bleachvals * 100, totalcounts[0, :] * 100, label='1')
# plt.plot(bleachvals * 100, totalcounts[1, :] * 100, label='2')
# plt.plot(bleachvals * 100, totalcounts[2, :] * 100, label='3')
# plt.plot(bleachvals * 100, totalcounts[3, :] * 100, label='4')
# plt.legend(frameon=False)
# plt.xlabel('Bleach %')
# plt.ylabel('% of visible particles')
# plt.show()


# plt.plot(bleachvals * 100, totalcounts[0, :] * 100, label='1')
# plt.plot(bleachvals * 100, totalcounts[1, :] * 100, label='2')
# plt.plot(bleachvals * 100, totalcounts[2, :] * 100, label='3')
# plt.plot(bleachvals * 100, totalcounts[3, :] * 100, label='4')
# plt.plot(bleachvals * 100, totalcounts[4, :] * 100, label='5')
# plt.plot(bleachvals * 100, totalcounts[5, :] * 100, label='6')
# plt.plot(bleachvals * 100, totalcounts[6, :] * 100, label='7')
# plt.plot(bleachvals * 100, totalcounts[7, :] * 100, label='8')
# plt.plot(bleachvals * 100, totalcounts[8, :] * 100, label='9')
# plt.plot(bleachvals * 100, totalcounts[9, :] * 100, label='10')
# plt.legend(frameon=False)
# plt.xlabel('Bleach %')
# plt.ylabel('% of visible particles')
# plt.show()
