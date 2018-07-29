import numpy as np

"""
Model with dynamic A and P species, detailed mass action description (simple scheme)

"""


class Params:
    def __init__(self, Da, konan, koffan, konap, koffap, kphosa, kdephosa, Dp, konpn, koffpn, konpp, koffpp, kphosp,
                 kdephosp, L, xsteps, psi, Tmax, deltat, starts):
        ########## aPAR ############

        # Diffusion
        self.Da = Da  # um2 s-1

        # Membrane exchange
        self.konan = konan  # um s-1
        self.koffan = koffan  # s-1
        self.konap = konap  # um s-1
        self.koffap = koffap  # s-1

        # Antagonism
        self.kphosa = kphosa  # um2 s-1
        self.kdephosa = kdephosa  # s-1

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

        ######### Misc ############

        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        ######### Starts ############

        self.amn_0 = starts[0]
        self.amp_0 = starts[1]
        self.acn_0 = starts[2]
        self.acp_0 = starts[3]

        self.pmn_0 = starts[4]
        self.pmp_0 = starts[5]
        self.pcn_0 = starts[6]
        self.pcp_0 = starts[7]


class Model:
    def __init__(self, p):
        """
        amn:
        apm:
        acn:
        acp:

        pmn:
        pmp:
        pcn:
        pcp:


        """

        self.params = p
        self.res = self.Res(self.params)

        self.amn = self.params.amn_0 * np.zeros([self.params.xsteps])  # membrane bound, not phosphorylated
        self.amp = self.params.amp_0 * np.zeros([self.params.xsteps])  # membrane bound, phosphorylated
        self.acn = self.params.acn_0  # cytoplasmic, not phosphorylated
        self.acp = self.params.acp_0  # cytoplasmic, phosphorylated

        self.pmn = self.params.pmn_0 * np.zeros([self.params.xsteps])  # membrane bound, not phosphorylated
        self.pmp = self.params.pmp_0 * np.zeros([self.params.xsteps])  # membrane bound, phosphorylated
        self.pcn = self.params.pcn_0  # cytoplasmic, not phosphorylated
        self.pcp = self.params.pcp_0  # cytoplasmic, phosphorylated

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_amn(self, lb, ub):
        diff = self.diffusion(self.amn[lb:ub], self.params.Da)
        onn = self.params.konan * self.acn
        offn = self.params.koffan * self.amn[lb:ub]
        phos = self.params.kphosa * (self.pmn[lb:ub] + self.pmp[lb:ub]) * self.amn[lb:ub]
        self.amn[lb:ub] += ((diff + onn - offn - phos) * self.params.deltat)

    def update_amp(self, lb, ub):
        diff = self.diffusion(self.amp[lb:ub], self.params.Da)
        onp = self.params.konap * self.acp
        offp = self.params.koffap * self.amp[lb:ub]
        phos = self.params.kphosa * (self.pmn[lb:ub] + self.pmp[lb:ub]) * self.amn[lb:ub]
        self.amp[lb:ub] += ((diff + onp - offp + phos) * self.params.deltat)

    def update_acn(self, lb, ub):
        onn = (1 / self.params.psi) * (self.params.konan * self.acn)
        offn = (1 / self.params.psi) * np.mean(self.params.koffan * self.amn[lb:ub])
        dephos = self.params.kdephosa * self.acp
        self.acn += ((offn + dephos - onn) * self.params.deltat)

    def update_acp(self, lb, ub):
        onp = (1 / self.params.psi) * (self.params.konap * self.acp)
        offp = (1 / self.params.psi) * np.mean(self.params.koffap * self.amp[lb:ub])
        dephos = self.params.kdephosa * self.acp
        self.acp += ((offp - dephos - onp) * self.params.deltat)

    def update_pmn(self, lb, ub):
        diff = self.diffusion(self.pmn[lb:ub], self.params.Dp)
        onn = self.params.konpn * self.pcn
        offn = self.params.koffpn * self.pmn[lb:ub]
        phos = self.params.kphosp * (self.amn[lb:ub] + self.amp[lb:ub]) * self.pmn[lb:ub]
        self.pmn[lb:ub] += ((diff + onn - offn - phos) * self.params.deltat)

    def update_pmp(self, lb, ub):
        diff = self.diffusion(self.pmp[lb:ub], self.params.Dp)
        onp = self.params.konpp * self.pcp
        offp = self.params.koffpp * self.pmp[lb:ub]
        phos = self.params.kphosp * (self.amn[lb:ub] + self.amp[lb:ub]) * self.pmn[lb:ub]
        self.pmp[lb:ub] += ((diff + onp - offp + phos) * self.params.deltat)

    def update_pcn(self, lb, ub):
        onn = (1 / self.params.psi) * (self.params.konpn * self.pcn)
        offn = (1 / self.params.psi) * np.mean(self.params.koffpn * self.pmn[lb:ub])
        dephos = self.params.kdephosp * self.pcp
        self.pcn += ((offn + dephos - onn) * self.params.deltat)

    def update_pcp(self, lb, ub):
        onp = (1 / self.params.psi) * (self.params.konpp * self.pcp)
        offp = (1 / self.params.psi) * np.mean(self.params.koffpp * self.pmp[lb:ub])
        dephos = self.params.kdephosp * self.pcp
        self.pcp += ((offp - dephos - onp) * self.params.deltat)

    def update_all(self):
        self.update_acn(0, self.params.xsteps)
        self.update_acp(0, self.params.xsteps)
        self.update_amn(0, self.params.xsteps)
        self.update_amp(0, self.params.xsteps)

        self.update_pcn(0, self.params.xsteps)
        self.update_pcp(0, self.params.xsteps)
        self.update_pmn(0, self.params.xsteps)
        self.update_pmp(0, self.params.xsteps)

    def equilibrate_all(self):
        self.update_acn(0, self.params.xsteps // 2)
        self.update_acp(0, self.params.xsteps // 2)
        self.update_amn(0, self.params.xsteps // 2)
        self.update_amp(0, self.params.xsteps // 2)

        self.update_pcn(self.params.xsteps // 2, self.params.xsteps)
        self.update_pcp(self.params.xsteps // 2, self.params.xsteps)
        self.update_pmn(self.params.xsteps // 2, self.params.xsteps)
        self.update_pmp(self.params.xsteps // 2, self.params.xsteps)

    def get_all(self):
        return [self.amn, self.amp, self.acn, self.acp, self.pmn, self.pmp, self.pcn, self.pcp]

    def run(self):

        # Equilibrate
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.equilibrate_all()
        self.res.update(-1, self.get_all())

        # Run model
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.update_all()
            self.res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, params):
            self.params = params
            self.scores = {}

            self.amn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.amp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.acn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.acp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.pmn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pmp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pcn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pcp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.amn[t + 1, :] = c[0]
            self.amp[t + 1, :] = c[1]
            self.acn[t + 1] = c[2]
            self.acp[t + 1] = c[3]

            self.pmn[t + 1, :] = c[4]
            self.pmp[t + 1, :] = c[5]
            self.pcn[t + 1] = c[6]
            self.pcp[t + 1] = c[7]

            self.aco[t + 1, :] = c[0] + c[1]
            self.pco[t + 1, :] = c[4] + c[5]

        def compress(self):
            self.amn = np.asarray([self.amn[-1, :], ])
            self.amp = np.asarray([self.amp[-1, :], ])
            self.acn = self.acn[-1]
            self.acp = self.acp[-1]

            self.pmn = np.asarray([self.pmn[-1, :], ])
            self.pmp = np.asarray([self.pmp[-1, :], ])
            self.pcn = self.pcn[-1]
            self.pcp = self.pcp[-1]

            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])
