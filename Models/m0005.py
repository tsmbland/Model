import numpy as np
import pickle

"""
Model with static A species, detailed mass action description of P species (simple scheme)

"""


class Params:
    def __init__(self, Dp, konpn, koffpn, konpp, koffpp, kphosp, kdephosp, L, xsteps, psi, Tmax, deltat, starts):
        ########## pPAR ############

        # Diffusion
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konpn = konpn  # um s-1
        self.koffpn = koffpn  # s-1
        self.konpp = konpp  # um s-1
        self.koffpp = koffpp  # s-1

        # Antagonism
        self.kphosp = kphosp  # um2 s-1
        self.kdephosp = kdephosp  # s-1

        ########## Misc ############

        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        ######### Starts ###########

        self.a_0 = starts[0]
        self.pmn_0 = starts[1]
        self.pmp_0 = starts[2]
        self.pcn_0 = starts[3]
        self.pcp_0 = starts[4]


class Model:
    def __init__(self, p):
        self.params = p
        self.res = self.Res(p)
        self.a = self.params.a_0
        self.pmn = self.params.pmn_0 * np.zeros([p.xsteps])  # membrane bound, not phosphorylated
        self.pmp = self.params.pmp_0 * np.zeros([p.xsteps])  # membrane bound, phosphorylated
        self.pcn = self.params.pcn_0  # cytoplasmic, not phosphorylated
        self.pcp = self.params.pcp_0  # cytoplasmic, phosphorylated

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_pmn(self):
        diff = self.diffusion(self.pmn, self.params.Dp)
        onn = self.params.konpn * self.pcn
        offn = self.params.koffpn * self.pmn
        phos = self.params.kphosp * self.a * self.pmn
        self.pmn += ((diff + onn - offn - phos) * self.params.deltat)

    def update_pmp(self):
        diff = self.diffusion(self.pmp, self.params.Dp)
        onp = self.params.konpp * self.pcp
        offp = self.params.koffpp * self.pmp
        phos = self.params.kphosp * self.a * self.pmn
        self.pmp += ((diff + onp - offp + phos) * self.params.deltat)

    def update_pcn(self):
        onn = (1 / self.params.psi) * (self.params.konpn * self.pcn)
        offn = (1 / self.params.psi) * np.mean(self.params.koffpn * self.pmn)
        dephos = self.params.kdephosp * self.pcp
        self.pcn += ((offn + dephos - onn) * self.params.deltat)

    def update_pcp(self):
        onp = (1 / self.params.psi) * (self.params.konpp * self.pcp)
        offp = (1 / self.params.psi) * np.mean(self.params.koffpp * self.pmp)
        dephos = self.params.kdephosp * self.pcp
        self.pcp += ((offp - dephos - onp) * self.params.deltat)

    def update_all(self):
        self.update_pcn()
        self.update_pcp()
        self.update_pmn()
        self.update_pmp()

    def get_all(self):
        return [self.a, self.pmn, self.pmp, self.pcn, self.pcp]

    def run(self):
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.update_all()
            self.res.update(t, self.get_all())
        return self.res

    class Res:
        def __init__(self, p):
            self.params = p

            self.a = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pmn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pmp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pcn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pcp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.a[t + 1] = c[0]
            self.pmn[t + 1, :] = c[1]
            self.pmp[t + 1, :] = c[2]
            self.pcn[t + 1] = c[3]
            self.pcp[t + 1] = c[4]

            self.aco[t + 1] = c[0]
            self.pco[t + 1, :] = c[1] + c[2]

        def compress(self):
            self.a = self.a[-1]

            self.pmn = np.asarray([self.pmn[-1, :], ])
            self.pmp = np.asarray([self.pmp[-1, :], ])
            self.pcn = self.pcn[-1]
            self.pcp = self.pcp[-1]

            self.aco = self.a[-1]
            self.pco = np.asarray([self.pco[-1, :], ])
