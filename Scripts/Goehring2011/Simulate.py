from Goehring2011 import Model
import copy

"""
Set up model (default parameters)

"""

M = Model(Da=0.28, Dp=0.15, konA=0.00858, koffA=0.0054, konP=0.0474, koffP=0.0073, kAP=0.19, kPA=2, alpha=2, beta=1,
          xsteps=500, psi=0.127, Tmax=1000, deltat=0.01, L=67.3, pA=1.56, pP=1)

"""
Run, plot final result, save final time point

"""

m = copy.deepcopy(M)
m.initiate()
m.run()
m.save('Saved')
m.plot_res()

"""
Run, save all time points, animate

"""

# m = copy.deepcopy(M)
# m.initiate()
# m.run(save_direc='_temp', save_gap=10)
# m.animate(direc='_temp')
