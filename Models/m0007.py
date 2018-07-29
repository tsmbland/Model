import numpy as np

"""
Multi phosphorylation model

"""


class Params:
    def __init__(self, Da, konan, koffan, konam, koffam, konap, koffap, kphosan, kphosam, kdephosap, kdephosam, Dp,
                 konpn, koffpn, konpm, koffpm, konpp, koffpp, kphospn, kphospm, kdephospp, kdephospm, L, xsteps, psi,
                 Tmax, deltat, starts):
        ########### aPARS ###########

        # Diffusion
        self.Da = Da  # um2 s-1

        # Membrane exchange
        self.konan = konan  # um s-1
        self.koffan = koffan  # s-1
        self.konam = konam  # um s-1
        self.koffam = koffam  # s-1
        self.konap = konap  # um s-1
        self.koffap = koffap  # s-1

        # Phosphorylation
        self.kphosan = kphosan  # um2 s-1
        self.kphosam = kphosam  # um2 s-1

        # Dephosphorylation
        self.kdephosap = kdephosap  # s-1
        self.kdephosam = kdephosam  # s-1

        ########### pPARS ###########

        # Diffusion
        self.Dp = Dp  # um2 s-1

        # Membrane exchange
        self.konpn = konpn  # um s-1
        self.koffpn = koffpn  # s-1
        self.konpm = konpm  # um s-1
        self.koffpm = koffpm  # s-1
        self.konpp = konpp  # um s-1
        self.koffpp = koffpp  # s-1

        # Phosphorylation
        self.kphospn = kphospn  # um2 s-1
        self.kphospm = kphospm  # um2 s-1

        # Dephosphorylation
        self.kdephospp = kdephospp  # s-1
        self.kdephospm = kdephospm  # s-1

        ########### Misc ###########

        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        ########## Starts ##########

        self.amn_0 = starts[0]
        self.amm_0 = starts[1]
        self.amp_0 = starts[2]
        self.acn_0 = starts[3]
        self.acm_0 = starts[4]
        self.acp_0 = starts[5]

        self.pmn_0 = starts[6]
        self.pmm_0 = starts[7]
        self.pmp_0 = starts[8]
        self.pcn_0 = starts[9]
        self.pcm_0 = starts[10]
        self.pcp_0 = starts[11]


class Model:
    def __init__(self, params):
        """
        amn:
        amm:
        amp:
        acn:
        acm:
        acp:

        pmn:
        pmm:
        pmp:
        pcn:
        pcm:
        pcp:

        """

        self.params = params
        self.res = self.Res(params)

        self.amn = self.params.amn_0 * np.zeros([self.params.xsteps])  # membrane bound, not phosphorylated
        self.amm = self.params.amm_0 * np.zeros([self.params.xsteps])  # membrane bound, medium phosphorylated
        self.amp = self.params.amp_0 * np.zeros([self.params.xsteps])  # membrane bound, phosphorylated
        self.acn = self.params.acn_0  # cytoplasmic, not phosphorylated
        self.acm = self.params.acm_0
        self.acp = self.params.acp_0  # cytoplasmic, phosphorylated

        self.pmn = self.params.pmn_0 * np.zeros([self.params.xsteps])  # membrane bound, not phosphorylated
        self.pmm = self.params.pmm_0 * np.zeros([self.params.xsteps])
        self.pmp = self.params.pmp_0 * np.zeros([self.params.xsteps])  # membrane bound, phosphorylated
        self.pcn = self.params.pcn_0  # cytoplasmic, not phosphorylated
        self.pcm = self.params.pcm_0
        self.pcp = self.params.pcp_0  # cytoplasmic, phosphorylated

    def diffusion(self, concs, coeff):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_amn(self, lb, ub):
        diff = self.diffusion(self.amn[lb:ub], self.params.Da)
        onan = self.params.konan * self.acn
        offan = self.params.koffan * self.amn[lb:ub]
        phosan = self.params.kphosan * (self.pmn[lb:ub] + self.pmm[lb:ub] + self.pmp[lb:ub]) * self.amn[lb:ub]
        self.amn[lb:ub] += ((diff + onan - offan - phosan) * self.params.deltat)

    def update_amm(self, lb, ub):
        diff = self.diffusion(self.amm[lb:ub], self.params.Da)
        onam = self.params.konam * self.acm
        offam = self.params.koffam * self.amm[lb:ub]
        phosan = self.params.kphosan * (self.pmn[lb:ub] + self.pmm[lb:ub] + self.pmp[lb:ub]) * self.amn[lb:ub]
        phosam = self.params.kphosam * (self.pmn[lb:ub] + self.pmm[lb:ub] + self.pmp[lb:ub]) * self.amm[lb:ub]
        self.amm[lb:ub] += ((diff + onam - offam - phosam + phosan) * self.params.deltat)

    def update_amp(self, lb, ub):
        diff = self.diffusion(self.amp[lb:ub], self.params.Da)
        onap = self.params.konap * self.acp
        offap = self.params.koffap * self.amp[lb:ub]
        phosam = self.params.kphosam * (self.pmn[lb:ub] + self.pmm[lb:ub] + self.pmp[lb:ub]) * self.amm[lb:ub]
        self.amp[lb:ub] += ((diff + onap - offap + phosam) * self.params.deltat)

    def update_acn(self, lb, ub):
        onan = (1 / self.params.psi) * (self.params.konan * self.acn)
        offan = (1 / self.params.psi) * np.mean(self.params.koffan * self.amn[lb:ub])
        dephosam = self.params.kdephosam * self.acm
        self.acn += ((offan + dephosam - onan) * self.params.deltat)

    def update_acm(self, lb, ub):
        onam = (1 / self.params.psi) * (self.params.konam * self.acm)
        offam = (1 / self.params.psi) * np.mean(self.params.koffam * self.amm[lb:ub])
        dephosam = self.params.kdephosam * self.acm
        dephosap = self.params.kdephosap * self.acp
        self.acm += ((offam + dephosap - dephosam - onam) * self.params.deltat)

    def update_acp(self, lb, ub):
        onap = (1 / self.params.psi) * (self.params.konap * self.acp)
        offap = (1 / self.params.psi) * np.mean(self.params.koffap * self.amp[lb:ub])
        dephosap = self.params.kdephosap * self.acp
        self.acp += ((offap - dephosap - onap) * self.params.deltat)

    def update_pmn(self, lb, ub):
        diff = self.diffusion(self.pmn[lb:ub], self.params.Dp)
        onpn = self.params.konpn * self.pcn
        offpn = self.params.koffpn * self.pmn[lb:ub]
        phospn = self.params.kphospn * (self.amn[lb:ub] + self.amm[lb:ub] + self.amp[lb:ub]) * self.pmn[lb:ub]
        self.pmn[lb:ub] += ((diff + onpn - offpn - phospn) * self.params.deltat)

    def update_pmm(self, lb, ub):
        diff = self.diffusion(self.pmm[lb:ub], self.params.Dp)
        onpm = self.params.konpm * self.pcm
        offpm = self.params.koffpm * self.pmm[lb:ub]
        phospn = self.params.kphospn * (self.amn[lb:ub] + self.amm[lb:ub] + self.amp[lb:ub]) * self.pmn[lb:ub]
        phospm = self.params.kphospm * (self.amn[lb:ub] + self.amm[lb:ub] + self.amp[lb:ub]) * self.pmm[lb:ub]
        self.pmm[lb:ub] += ((diff + onpm - offpm - phospm + phospn) * self.params.deltat)

    def update_pmp(self, lb, ub):
        diff = self.diffusion(self.pmp[lb:ub], self.params.Dp)
        onpp = self.params.konpp * self.pcp
        offpp = self.params.koffpp * self.pmp[lb:ub]
        phospm = self.params.kphospm * (self.amn[lb:ub] + self.amm[lb:ub] + self.amp[lb:ub]) * self.pmm[lb:ub]
        self.pmp[lb:ub] += ((diff + onpp - offpp + phospm) * self.params.deltat)

    def update_pcn(self, lb, ub):
        onpn = (1 / self.params.psi) * (self.params.konpn * self.pcn)
        offpn = (1 / self.params.psi) * np.mean(self.params.koffpn * self.pmn[lb:ub])
        dephospm = self.params.kdephospm * self.pcm
        self.pcn += ((offpn + dephospm - onpn) * self.params.deltat)

    def update_pcm(self, lb, ub):
        onpm = (1 / self.params.psi) * (self.params.konpm * self.pcm)
        offpm = (1 / self.params.psi) * np.mean(self.params.koffpm * self.pmm[lb:ub])
        dephospm = self.params.kdephospm * self.pcm
        dephospp = self.params.kdephospp * self.pcp
        self.pcm += ((offpm + dephospp - dephospm - onpm) * self.params.deltat)

    def update_pcp(self, lb, ub):
        onpp = (1 / self.params.psi) * (self.params.konpp * self.pcp)
        offpp = (1 / self.params.psi) * np.mean(self.params.koffpp * self.pmp[lb:ub])
        dephospp = self.params.kdephospp * self.pcp
        self.pcp += ((offpp - dephospp - onpp) * self.params.deltat)

    def update_all(self):
        self.update_amn(0, self.params.xsteps)
        self.update_amm(0, self.params.xsteps)
        self.update_amp(0, self.params.xsteps)
        self.update_acn(0, self.params.xsteps)
        self.update_acm(0, self.params.xsteps)
        self.update_acp(0, self.params.xsteps)
        self.update_pmn(0, self.params.xsteps)
        self.update_pmm(0, self.params.xsteps)
        self.update_pmp(0, self.params.xsteps)
        self.update_pcn(0, self.params.xsteps)
        self.update_pcm(0, self.params.xsteps)
        self.update_pcp(0, self.params.xsteps)

    def update_all_e(self):
        self.update_amn(0, self.params.xsteps // 2)
        self.update_amm(0, self.params.xsteps // 2)
        self.update_amp(0, self.params.xsteps // 2)
        self.update_acn(0, self.params.xsteps // 2)
        self.update_acm(0, self.params.xsteps // 2)
        self.update_acp(0, self.params.xsteps // 2)
        self.update_pmn(self.params.xsteps // 2, self.params.xsteps)
        self.update_pmm(self.params.xsteps // 2, self.params.xsteps)
        self.update_pmp(self.params.xsteps // 2, self.params.xsteps)
        self.update_pcn(self.params.xsteps // 2, self.params.xsteps)
        self.update_pcm(self.params.xsteps // 2, self.params.xsteps)
        self.update_pcp(self.params.xsteps // 2, self.params.xsteps)

    def get_all(self):
        return [self.amn, self.amm, self.amp, self.acn, self.acm, self.acp, self.pmn, self.pmm, self.pmp, self.pcn,
                self.pcm, self.pcp]

    def run(self):

        # Equilibrate
        for t in range(int(self.params.Tmax / self.params.deltat)):
            self.update_all_e()
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
            self.amm = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.amp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.acn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.acm = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.acp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.pmn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pmm = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pmp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pcn = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pcm = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])
            self.pcp = np.zeros([int(self.params.Tmax / self.params.deltat) + 1])

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.amn[t + 1, :] = c[0]
            self.amm[t + 1, :] = c[1]
            self.amp[t + 1, :] = c[2]
            self.acn[t + 1] = c[3]
            self.acm[t + 1] = c[4]
            self.acp[t + 1] = c[5]

            self.pmn[t + 1, :] = c[6]
            self.pmm[t + 1, :] = c[7]
            self.pmp[t + 1, :] = c[8]
            self.pcn[t + 1] = c[9]
            self.pcm[t + 1] = c[10]
            self.pcp[t + 1] = c[11]

            self.aco[t + 1, :] = c[0] + c[1] + c[2]
            self.pco[t + 1, :] = c[6] + c[7] + c[8]

        def compress(self):
            self.amn = np.asarray([self.amn[-1, :], ])
            self.amm = np.asarray([self.amm[-1, :], ])
            self.amp = np.asarray([self.amp[-1, :], ])
            self.acn = self.acn[-1]
            self.acm = self.acm[-1]
            self.acp = self.acp[-1]

            self.pmn = np.asarray([self.pmn[-1, :], ])
            self.pmm = np.asarray([self.pmm[-1, :], ])
            self.pmp = np.asarray([self.pmp[-1, :], ])
            self.pcn = self.pcn[-1]
            self.pcm = self.pcm[-1]
            self.pcp = self.pcp[-1]

            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])
