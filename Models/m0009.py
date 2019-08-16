import numpy as np

"""
Model with pPAR membrane binding receptor
Not tested

"""


class Model:
    def __init__(self, ac_on, am_off, sc_on, sm_off, pc_on, pm_off, k_dim_f, k_dim_b, pc_on2, sc_on2, pm_off2, sm_off2,
                 r_am_p, r_sm_a, r_pm_a, r_spm10_a, r_spm11_a, r_spm01_a, d_am, d_sm, d_pm, d_spm10, d_spm11, d_spm01,
                 xsteps, psi, Tmax, deltat, deltal, radii, am_0, ac_0, pm_0, pc_0, sm_0, sc_0, spm10_0, spm11_0,
                 spm01_0, spc_0):
        """

        PARAMETERS

        p0: ac on rate                        ac_on
        p1: am off rate                       am_off
        p2: sc on rate                        sc_on
        p3: sm off rate                       sm_off
        p4: pc on rate                        pc_on
        p5: pm off rate                       pm_off
        p6: dimerisation                      k_dim_f
        p7: separation                        k_dim_b
        p8: p on 2                            pc_on2
        p9: s on 2                            sc_on2
        p10: p off 2                          pm_off2
        p11: s off 2                          sm_off2
        p12: antagonism from p to am          r_am_p
        p13: antagonism from a to sm          r_sm_a
        p14: antagonism from a to pm          r_pm_a
        p15: antagonism from a to spm10       r_spm10_a
        p16: antagonism from a to spm11       r_spm11_a
        p17: antagonism from a to spm01       r_spm01_a
        p18: diffusion of am                  d_am
        p19: diffusion of sm                  d_sm
        p20: diffusion of pm                  d_pm
        p21: diffusion of spm10               d_spm10
        p22: diffusion of spm11               d_spm11
        p23: diffusion of spm01               d_spm01


        SPECIES

        am: Membrane aPAR
        ac: Cytoplasmic aPAR
        pm: Membrane pPAR
        pc: Cytoplasmic pPAR
        sm: Membrane scaffold (unbound to p)
        sc: Cytoplasmic scaffold
        spm10: Membrane scaffold with non-membrane-bound p
        spm11: Membrane scaffold with membrane-bound p
        spm01: Non-membrane bound scaffold, membrane bound p
        spc: Cytoplasmic scaffold-p

        """

        # SPECIES
        self.am = am_0
        self.ac = ac_0
        self.pm = pm_0
        self.pc = pc_0
        self.sm = sm_0
        self.sc = sc_0
        self.spm10 = spm10_0
        self.spm11 = spm11_0
        self.spm01 = spm01_0
        self.spc = spc_0
        self.time = 0

        # PARAMETERS

        self.ac_on = ac_on
        self.am_off = am_off
        self.sc_on = sc_on
        self.sm_off = sm_off
        self.pc_on = pc_on
        self.pm_off = pm_off
        self.k_dim_f = k_dim_f
        self.k_dim_b = k_dim_b
        self.pc_on2 = pc_on2
        self.sc_on2 = sc_on2
        self.pm_off2 = pm_off2
        self.sm_off2 = sm_off2
        self.r_am_p = r_am_p
        self.r_sm_a = r_sm_a
        self.r_pm_a = r_pm_a
        self.r_spm10_a = r_spm10_a
        self.r_spm11_a = r_spm11_a
        self.r_spm01_a = r_spm01_a
        self.d_am = d_am
        self.d_sm = d_sm
        self.d_pm = d_pm
        self.d_spm10 = d_spm10
        self.d_spm11 = d_spm11
        self.d_spm01 = d_spm01

        ########### Misc ###########

        self.xsteps = int(xsteps)
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s
        self.deltal = deltal  # um
        self.radii = radii  # um

    def diffusion(self, concs):
        return concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]

    def reactions(self):
        """

        # Monomer binding
        r0: binding of ac to cortex
        r1: unbinding of am from cortex
        r2: binding of sc to cortex
        r3: unbinding of sm from cortex
        r4: binding of pc to cortex
        r5: unbinding of pm from cortex

        # Dimer binding
        r6: Single binding of spc to cortex -> spm10
        r7: Single  binding of spc to cortex -> spm01
        r8: unbinding of spm01 from cortex
        r9: unbinding of spm10 from cortex

        # Recruitment part 1
        r10: Dimerisation of sm with pc -> spm10
        r11: Separation of spm10
        r12: Dimerisation of pm with sc -> spm01
        r13: Separation of spm01

        # Recruitment part 2
        r14: binding of p in spm10 to cortex
        r15: binding of s in spm01 to cortex
        r16: unbinding of p in spm11 -> spm10
        r17: unbinding of s in spm11 -> spm01

        # Dimerisation
        r18: dimerisation of pm and sm
        r19: dimerisation of pc and sc
        r20: separation of spm11
        r21: separation of spc

        # Antagonism
        r22: antagonism from p to am
        r23: antagonism from am to sm
        r24: antagonism from am to pm
        r25: Antagonism from am to spm10
        r26: Antagonism from am to spm11
        r27: Antagonism from am to spm01

        # Diffusion
        r28: diffusion of am
        r29: diffusion of sm
        r30: diffusion of pm
        r31: diffusion of spm10
        r32: diffusion of spm11
        r33: diffusion of spm01


        """

        r = [None] * 34

        # Monomer binding
        r[0] = self.ac_on * self.ac
        r[1] = self.am_off * self.am
        r[2] = self.sc_on * self.sc
        r[3] = self.sm_off * self.sm
        r[4] = self.pc_on * self.pc
        r[5] = self.pm_off * self.pm

        # Dimer binding
        r[6] = self.sc_on * self.spc
        r[7] = self.pc_on * self.spc
        r[8] = self.sm_off * self.spm01
        r[9] = self.pm_off * self.spm10

        # Recruitment part 1
        r[10] = self.k_dim_f * self.sm * self.pc
        r[11] = self.k_dim_b * self.spm10
        r[12] = self.k_dim_f * self.pm * self.sc
        r[13] = self.k_dim_b * self.spm01

        # Recruitment part 2
        r[14] = self.pc_on2 * self.spm10
        r[15] = self.sc_on2 * self.spm01
        r[16] = self.pm_off2 * self.spm11
        r[17] = self.sm_off2 * self.spm11

        # Dimerisation
        r[18] = self.k_dim_f * self.pm * self.sm
        r[19] = self.k_dim_f * self.pc * self.sc
        r[20] = self.k_dim_b * self.spm11
        r[21] = self.k_dim_b * self.spc

        # Antagonism
        r[22] = self.r_am_p * (self.pm + self.spm11 + self.spm01) * self.am
        r[23] = self.r_sm_a * self.am * self.sm
        r[24] = self.r_pm_a * self.am * self.pm
        r[25] = self.r_spm10_a * self.am * self.spm10
        r[26] = self.r_spm11_a * self.am * self.spm11
        r[27] = self.r_spm01_a * self.am * self.spm01

        # Diffusion
        r[28] = self.d_am * self.diffusion(self.am)
        r[29] = self.d_sm * self.diffusion(self.sm)
        r[30] = self.d_pm * self.diffusion(self.pm)
        r[31] = self.d_spm10 * self.diffusion(self.spm10)
        r[32] = self.d_spm11 * self.diffusion(self.spm11)
        r[33] = self.d_spm01 * self.diffusion(self.spm01)

        return r

    def update_am(self, r):
        self.am += (+ r[0] - r[1] - r[22] + r[28]) * self.deltat

    def update_ac(self, r):
        self.ac += (-self.psi * r[0] + self.psi * np.mean(r[1]) + self.psi * np.mean(r[22])) * self.deltat

    def update_pm(self, r):
        self.pm += (+ r[4] - r[5] - r[12] + r[13] - r[18] + r[20] - r[24] + r[30]) * self.deltat

    def update_pc(self, r):
        self.pc += (
                       -self.psi * r[4] + self.psi * np.mean(r[5]) - self.psi * np.mean(r[10]) + self.psi * np.mean(
                           r[11]) -
                       r[19] + r[21] + self.psi * np.mean(r[24])) * self.deltat

    def update_sm(self, r):
        self.sm += (+ r[2] - r[3] - r[10] + r[11] - r[18] + r[20] - r[23] + r[29]) * self.deltat

    def update_sc(self, r):
        self.sc += (
                       -self.psi * r[2] + self.psi * np.mean(r[3]) - self.psi * np.mean(r[12]) + self.psi * np.mean(
                           r[13]) -
                       r[19] + r[21] + self.psi * np.mean(r[23])) * self.deltat

    def update_spm10(self, r):
        self.spm10 += (r[6] - r[8] + r[10] - r[11] - r[14] + r[16] - r[25] + r[31]) * self.deltat

    def update_spm11(self, r):
        self.spm11 += (r[14] + r[15] - r[16] - r[17] + r[18] - r[20] - r[26] + r[32]) * self.deltat

    def update_spm01(self, r):
        self.spm01 += (r[7] - r[9] + r[12] - r[13] + r[15] - r[17] - r[27] - r[33]) * self.deltat

    def update_spc(self, r):
        self.spc += (- self.psi * r[6] - self.psi * r[7] + self.psi * np.mean(r[8]) + self.psi * np.mean(r[9]) + r[19] -
                     r[21] + self.psi * np.mean(r[25]) + self.psi * np.mean(r[26]) + self.psi * np.mean(
            r[27])) * self.deltat

    def react(self):
        r = self.reactions()
        self.update_am(r)
        self.update_ac(r)
        self.update_pm(r)
        self.update_pc(r)
        self.update_sm(r)
        self.update_sc(r)
        self.update_spm10(r)
        self.update_spm11(r)
        self.update_spm01(r)
        self.update_spc(r)

    def run(self):
        for t in range(int(self.Tmax / self.deltat)):
            self.react()
            self.time = (t + 1) * self.deltat

    def save(self, direc):
        np.savetxt(direc + 'am.txt', self.am)
        np.savetxt(direc + 'ac.txt', [self.ac])
        np.savetxt(direc + 'pm.txt', self.pm)
        np.savetxt(direc + 'pc.txt', [self.pc])
        np.savetxt(direc + 'sm.txt', self.sm)
        np.savetxt(direc + 'sc.txt', [self.sc])
        np.savetxt(direc + 'spm10.txt', self.spm10)
        np.savetxt(direc + 'spm11.txt', self.spm11)
        np.savetxt(direc + 'spm01.txt', self.spm01)
        np.savetxt(direc + 'spc.txt', [self.spc])
        np.savetxt(direc + 'time.txt', [self.time])
