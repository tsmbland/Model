import numpy as np

molecules = ['am', 'p1m', 'p2m', 'p3m', 'ac', 'p1c', 'p2c', 'p3c']


t0 = [np.zeros([500]), np.zeros([500]), np.zeros([500]), np.zeros([500]), 3, 1, 1, 1]


def diffusion(concs, coeff, p):
    diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
        np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)
    return diff



def Reactions(m, k):

    r = [[None] * 8] * 8

    # On reactions
    r[4][0] = k.aon * m[4]
    r[5][1] = k.pona * m[5]
    r[6][2] = k.ponb * m[6]
    r[7][3] = k.ponb * m[7]

    # Off reactions
    r[0][4] = k.aon * m[0] + k.kAP * (m[1] + m[2] + m[3])
    r[1][5] = k.poffa * m[1]
    r[2][6] = k.poffa * m[2]
    r[3][7] = k.poffb * m[3]

    # Cytoplasmic reactions
    r[7][6] = k.kdephos * m[7]
    r[6][5] = k.kdephos * m[6]

    # Cortical reactions
    r[1][2] = k.kphos * m[0] * m[1]
    r[2][3] = k.kphos * m[0] * m[2]

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





