from Models.SimplePosFeedback import Model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def animate(model, p1, framerate, name, ymax=None):
    """

    p1: model time between frames

    """

    # Initial equilibration (no antagonism)
    kAP, kPA = model.kAP, model.kPA
    model.kAP, model.kPA = 0, 0
    model.psi /= 2

    for t in range(10000):
        model.react()

    # Polarise
    model.am *= 1 * np.r_[np.ones([50]), np.zeros([50])]
    model.pm *= 1 * np.r_[np.zeros([50]), np.ones([50])]
    model.psi *= 2

    # Add antagonism
    model.kAP = kAP
    model.kPA = kPA

    # Seu up figure
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(4.5, 3.5)
    ax = fig.add_subplot(111)
    frames = range(0, int(model.Tmax / p1))

    def update_anim(i):
        ax.clear()

        # Plot current timepoint
        ax.plot(model.am, c='r')
        ax.plot(model.pm, c='b')

        # Configure plot
        ax.set_title('Time (s): {0:.0f}'.format(model.time))
        ax.set_ylim(bottom=0, top=ymax)
        ax.set_xticks([], [])
        ax.set_xlabel('Anterior                                                  Posterior')
        ax.set_ylabel('Membrane Concentration (a.u.)')

        # Run model until next save point
        if i != 0:
            for t in range(int(p1 / model.deltat)):
                model.react()
                model.time += model.deltat

    # Animate
    paranim = animation.FuncAnimation(fig, update_anim, frames=iter(frames), save_count=len(frames))
    writer = animation.writers['ffmpeg']
    writer = writer(fps=framerate, bitrate=2000)
    paranim.save(name, writer=writer, dpi=300)


# """
# Animation 1: Linear symmetric, A dominant
#
# """
#
# base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
#                    koffP=0.01, kposP=0, kAP=10 ** -2, kPA=10 ** -1.5, ePneg=1, eAneg=1, xsteps=100, Tmax=300,
#                    deltat=0.01, deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1)
# animate(base_model, p1=2, framerate=24, name='anim1.mp4', ymax=7.5)
#
# """
# Animation 2: Linear symmetric, P dominant
#
# """
#
# base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
#                    koffP=0.01, kposP=0, kAP=10 ** -1.5, kPA=10 ** -2, ePneg=1, eAneg=1, xsteps=100, Tmax=300,
#                    deltat=0.01, deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1)
# animate(base_model, p1=2, framerate=24, name='anim2.mp4', ymax=7.5)

# """
# Animation 3: Linear symmetric, equal
#
# """
#
# base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
#                    koffP=0.01, kposP=0, kAP=10 ** -1, kPA=10 ** -1, ePneg=1, eAneg=1, xsteps=100, Tmax=300,
#                    deltat=0.01, deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1)
# animate(base_model, p1=2, framerate=24, name='anim3.mp4', ymax=7.5)
#
# """
# Animation 4: Non-linear symmetric
#
# """
#
# base_model = Model(Da=0.1, Dp=0.1, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
#                    koffP=0.01, kposP=0, kAP=10 ** -1, kPA=10 ** -1, ePneg=2, eAneg=2, xsteps=100, Tmax=300,
#                    deltat=0.01, deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=1, pm_0=np.zeros([100]), pc_0=1)
# animate(base_model, p1=2, framerate=24, name='anim4.mp4', ymax=7.5)

# """
# Animation 5: Parameterised (ols)
#
# """
#
# koffA = 0.0092
# koffP = 0.0073
# psi = 0.10318684114244771
# dosP = 0.294005475663175
# dosA = 1.05143336288 * dosP
#
# base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=22.4738241751 * koffP,
#                    koffP=koffP, kposP=8.6578099573 * koffP, kAP=10 ** -0.638671875, kPA=10 ** -0.737890625, ePneg=1,
#                    eAneg=1,
#                    xsteps=100, Tmax=1000, deltat=0.01, deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA,
#                    pm_0=np.zeros([100]), pc_0=dosP)
# animate(base_model, p1=3, framerate=24, name='anim5.mp4', ymax=5)

# """
# Animation 6: Parameterised (odr)
#
# """
#
# koffA = 0.0092
# koffP = 0.0073
# psi = 0.10318684114244771
# dosP = 0.294005475663175
# dosA = 1.05143336288 * dosP
#
# base_model = Model(Da=0.095, Dp=0.15, konA=18.4 * koffA, koffA=koffA, kposA=0, konP=11.6824354189 * koffP,
#                    koffP=koffP,
#                    kposP=12.7357711156 * koffP, kAP=10 ** -1.03125, kPA=10 ** -1.3484375, ePneg=1, eAneg=1, xsteps=100,
#                    Tmax=1000,
#                    deltat=0.01,
#                    deltax=0.5, psi=psi, am_0=np.zeros([100]), ac_0=dosA, pm_0=np.zeros([100]), pc_0=dosP)
#
# animate(base_model, p1=3, framerate=24, name='anim6.mp4', ymax=5)
