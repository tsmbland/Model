import numpy as np
import pickle

"""
Model description

"""


class Params:
    def __init__(self, param1):
        self.param1 = param1  # units


class Model:
    def __init__(self, params, starts):
        """

        :param params: object of Params class
        :param starts: array of starting values
        """

        """
        c0: component description
        """

        self.params = params
        self.starts = starts
        self.res = self.Res(params)
        self.c0 = starts[0] * np.ones([self.params.xsteps])

    def reactions(self, lb, ub):
        """

        :param lb: lower spatial bound
        :param ub: upper spatial bound
        :return: r: array of reaction rates
        """

        """
        r0: description of reaction
        """
        r = np.zeros([1])
        r[0] = ...

        return r

    def diffusion(self, concs, coeff):
        """

        :param concs: concentration array
        :param coeff: diffusion coefficient
        :return:
        """

        diff = coeff * (concs[np.append(np.array(range(1, len(concs))), [len(concs) - 2])] - 2 * concs + concs[
            np.append([1], np.array(range(len(concs) - 1)))]) / (self.params.L / self.params.xsteps)

        return diff

    def update_c0(self, r):
        """

        :param r: reaction rate array (returned from self.reactions(lb, ub)
        :return:
        """
        self.c0 += (r[0] + ...) * self.params.deltat

    ...

    def get_all(self):
        return [self.c0, ...]

    def run(self):

        # Equilibrate
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r1 = self.reactions(0, self.params.xsteps // 2)
            self.update_c0(r1)
            ...

            r2 = self.reactions(self.params.xsteps // 2, self.params.xsteps)
            ...

        # Run model
        for t in range(int(self.params.Tmax / self.params.deltat)):
            r = self.reactions(0, self.params.xsteps)
            self.update_c0(r)
            ...
            self.res.update(t, self.get_all())

        return self.res

    class Res:
        def __init__(self, params):
            self.params = params
            self.c0 = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            ...

            self.aco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])
            self.pco = np.zeros([int(self.params.Tmax / self.params.deltat) + 1, self.params.xsteps])

        def update(self, t, c):
            self.c0[t + 1, :] = c[0]
            ...

            self.aco[t + 1, :] = ...
            self.pco[t + 1, :] = ...

        def compress(self):
            self.c0 = np.asarray([self.c0[-1, :], ])
            ...
            self.aco = np.asarray([self.aco[-1, :], ])
            self.pco = np.asarray([self.pco[-1, :], ])
