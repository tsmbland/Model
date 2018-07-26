import numpy as np
import seaborn as sns

sns.set()
sns.set_style("ticks")


class MiscParams:
    """
    General parameters shared between all models

    """

    def __init__(self, L, xsteps, psi, Tmax, deltat, Aeqmin, Aeqmax, Peqmin, Peqmax):
        # Misc
        self.L = L  # um
        self.xsteps = xsteps
        self.psi = psi  # um-1
        self.Tmax = Tmax  # s
        self.deltat = deltat  # s

        # Equilibration
        self.Aeqmin = Aeqmin
        self.Aeqmax = Aeqmax
        self.Peqmin = Peqmin
        self.Peqmax = Peqmax


def diffusion(concs, coeff, p):
    if hasattr(concs, '__len__'):
        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (p.L / p.xsteps)
    else:
        diff = 0
    return diff


class Molecule:
    def __init__(self, D, kon, koff, p):
        self.D = D
        self.kon = kon
        self.koff = koff
        self.cort = np.zeros([p.xsteps])
        self.cyt = 0


def update_mol(mol, p):
    diff = diffusion(mol.cort, mol.D, p)
    off = mol.koff * mol.cort
    on = mol.kon * mol.cyt
    mol.cort += ((diff + on - off) * p.deltat)
    mol.cyt += off


def react(molecules, ks):
    # Phosphorylation on cortex (A to P)
    # ant = param * (molecules['p'].cort ** param) * molecules['a'].cort

    # Phosphorylation on cortex (P to A)
    phos1 = ks.phos1k * (molecules['a'].cort ** ks.phos1e) * molecules['p1'].cort
    phos2 = ks.phos2k * (molecules['a'].cort ** ks.phos2e) * molecules['p2'].cort

    # Dephosphorylation in cytoplasm: 1
    dephos1 = ks.dephos1k * molecules['p3'].cyt
    dephos2 = ks.dephos2k * molecules['p2'].cyt

    # molecules['a'].cort -= ant
    # molecules['a'].cyt += ant
    molecules['p1'].cort -= phos1
    molecules['p2'].cort += phos1
    molecules['p2'].cort -= phos2
    molecules['p3'].cort += phos2
    molecules['p2'].cyt += dephos1
    molecules['p3'].cyt -= dephos1
    molecules['p1'].cyt += dephos2
    molecules['p2'].cyt -= dephos2


p = MiscParams(L=50, xsteps=500, psi=0.3, Tmax=1000, deltat=0.01, Aeqmin=0, Aeqmax=0.5, Peqmin=0.5, Peqmax=1)


class Ks:
    phos1k = 1
    phos2k = 1
    phos1e = 2
    phos2e = 2
    dephos1k = 1
    dephos2k = 1


Molecules = {'a': Molecule(D=1, kon=0.006, koff=0.005, p=p),
             'p1': Molecule(D=1, kon=0.006, koff=0.005, p=p),
             'p2': Molecule(D=1, kon=0.06, koff=0.005, p=p),
             'p3': Molecule(D=1, kon=0.06, koff=1, p=p)}


def run_model(Molecules, p):
    # Init
    Molecules['a'].cort[:] = 0
    Molecules['a'].cyt = 1
    Molecules['p1'].cort[:] = 0
    Molecules['p1'].cyt = 1
    Molecules['p2'].cort[:] = 0
    Molecules['p2'].cyt = 1
    Molecules['p3'].cort[:] = 0
    Molecules['p3'].cyt = 1

    # Run model
    for t in range(int(p.Tmax / p.deltat)):
        for key, value in Molecules.items():
            update_mol(value, p)

        react(Molecules, Ks)


run_model(Molecules, p)
