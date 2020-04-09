import PyDSTool as dst
import numpy as np
from matplotlib import pyplot as plt

DSargs = dst.args(name='Wave pinning model')
DSargs.pars = {'k0': 0.02,
               'gama': 1,
               'K': 1,
               'delta': 1,
               'd': 0}
DSargs.pdomain = {'d': [0, 4],
                  'k0': [0, 0.1]}
DSargs.fnspecs = {'v': (['u'], 'd - u')}
DSargs.varspecs = {'u': 'v(u) * (k0 + (gama * (u ** 2)) / ((K**2) + (u ** 2))) - delta * u',
                   'w': 'u-w'}
DSargs.ics = {'u': 0, 'w': 0}
# DSargs.tdomain = [0, 30]
ode = dst.Generator.Vode_ODEsystem(DSargs)

# traj = ode.compute('traj')
# pts = traj.sample(dt=0.1)
#
# plt.plot(pts['t'], pts['u'])
# # plt.plot(pts['t'], pts['v'])
# plt.show()

PC = dst.ContClass(ode)  # Set up continuation class

PCargs = dst.args(name='EQ1',
                  type='EP-C')  # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
PCargs.freepars = ['d']  # control parameter(s) (it should be among those specified in DSargs.pars)
PCargs.MaxNumPoints = 500  # The following 3 parameters are set after trial-and-error
PCargs.MaxStepSize = 0.001
PCargs.MinStepSize = 1e-5
PCargs.StepSize = 2e-2
PCargs.LocBifPoints = 'LP'  # detect limit points / saddle-node bifurcations
PCargs.SaveEigen = True
PCargs.StopAtPoints = ['B']
PC.newCurve(PCargs)
PC['EQ1'].forward()
PC.display(['d', 'u'], stability=True)
plt.show()

PC['EQ1'].info()

# PCargs = dst.args(name='SN1', type='LP-C')
# PCargs.initpoint = 'EQ1:LP2'
# PCargs.freepars = ['d', 'k0']
# PCargs.MaxStepSize = 0.01
# PCargs.LocBifPoints = ['CP']
# PCargs.MaxNumPoints = 200
# PCargs.StopAtPoints = ['B']
# PC.newCurve(PCargs)
# PC['SN1'].forward()
# PC['SN1'].backward()
# PC['SN1'].display(['d', 'k0'])
# plt.show()
#
# PC['SN1'].info()
