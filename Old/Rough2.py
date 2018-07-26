import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# molecules = ['a1m', 'a2m', 'a3m', 'a1c', 'a2c', 'a3c', 'p1m', 'p2m', 'p3m', 'p1c', 'p2c', 'p3c']

# molecules = ['am', 'ac', 'p1m', 'p2m', 'p3m', 'p1c', 'p2c', 'p3c']

corts = ['am', 'p1m', 'p2m', 'p3m']
cyts = ['ac', 'p1c', 'p2c', 'p3c']

corts0 = [np.zeros([500]), np.zeros([500]), np.zeros([500]), np.zeros([500])]
cyts0 = [3, 1, 1, 1]


# df = pd.DataFrame(np.zeros([len(molecules), len(molecules)]), columns=molecules, index=molecules)

# starts = [np.zeros([500]), 1, np.zeros([500]), 1, np.zeros([500]), 1, np.zeros([500]), 1]

# starts = [np.zeros([500]), 3, np.zeros([500]), np.zeros([500]), np.zeros([500]), 1, 1, 1]


class MiscParams:
    """
    General parameters shared between all models

    """

    def __init__(self, L, xsteps, psi, Tmax, deltat):
        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s


def cortical_reactions(m, k):
    r = np.zeros([4, 4, 500])
    r[1, 2, :] = k.kphos * m[0] * m[1]
    r[2, 3, :] = k.kphos * m[0] * m[2]
    return r


def cytoplasmic_reactions(m, k):
    r = np.zeros([4, 4])
    r[3, 2] = k.kdephos * m[3]
    r[2, 1] = k.kdephos * m[2]
    return r


def on_reactions(m, k):
    r = np.zeros([4, 500])
    r[0, :] = k.aon * m[0]
    r[1, :] = k.pona * m[1]
    r[2, :] = k.ponb * m[2]
    r[3, :] = k.ponb * m[3]
    return r


def off_reactions(m, k):
    r = np.zeros([4, 500])
    r[0, :] = k.aoff * m[0] + k.kAP * (m[1] + m[2] + m[3])
    r[1, :] = k.poffa * m[1]
    r[2, :] = k.poffa * m[2]
    r[3, :] = k.poffb * m[3]
    return r


def diffusion(concs, coeff, p):
    diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
        np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)
    return diff


def diffusion_reactions(m, k, p):
    r = np.zeros([4, 500])
    r[0, :] = diffusion(m[0], k.Da, p)
    r[1, :] = diffusion(m[1], k.Dp, p)
    r[2, :] = diffusion(m[2], k.Dp, p)
    r[3, :] = diffusion(m[3], k.Dp, p)
    return r


class k:
    pona = 0.01
    ponb = 0.1
    poffa = 0.01
    poffb = 1
    kphos = 1
    kdephos = 1
    aon = 0.01
    aoff = 0.01
    kAP = 0
    Da = 1
    Dp = 1


def run_model(mcort, mcyt, k, p):
    rescy = np.zeros([4, 10000])
    rescort = np.zeros([4, 500, 10000])

    for t in range(int(p.Tmax / p.deltat)):
        d = diffusion_reactions(mcort, k, p)

        on = on_reactions(mcyt, k)
        off = off_reactions(mcort, k)
        cy = cytoplasmic_reactions(mcyt, k)
        cor = cortical_reactions(mcort, k)

        mcort += (d + on - off + np.sum(cor, axis=0) - np.sum(cor, axis=1)) * p.deltat
        mcyt += ((np.sum(off, axis=1) * p.psi) - (np.sum(on, axis=1) * p.psi) + np.sum(cy, axis=0)
                 - np.sum(cy, axis=1)) * p.deltat
        rescy[:, t] = mcyt
        rescort[:, :, t] = mcort

    return rescy, rescort


p = MiscParams(L=50, xsteps=500, psi=0.3, Tmax=100, deltat=0.01)

r1, r2 = run_model(corts0, cyts0, k, p)

print(r1)
print(r2)

plt.plot(r1[0, :])
plt.plot(r1[1, :])
plt.plot(r1[2, :])
plt.plot(r1[3, :])
plt.show()
