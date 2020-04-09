from Models.SimplePosFeedbackCyt import Model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def animate(model, p1, framerate, name, ymax=None, ymax2=None):
    """

    p1: model time between frames

    """

    # Initial equilibration (no antagonism)
    kAP, kPA = model.kAP, model.kPA
    model.kAP, model.kPA = 0, 0
    model.psi /= 2

    for t in range(100000):
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
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    frames = range(0, int(model.Tmax / p1))

    def update_anim(i):
        ax1.clear()
        ax2.clear()

        # Plot current timepoint
        ax2.plot(model.ac, c=(1, 0.5, 0.5))
        ax2.plot(model.pc, c=(0.5, 0.5, 1))
        ax1.plot(model.am, c=(1, 0, 0))
        ax1.plot(model.pm, c=(0, 0, 1))


        # Configure plot
        ax1.set_title('Time (s): {0:.0f}'.format(model.time))
        ax1.set_xticks([], [])
        ax1.set_xlabel('Anterior                                                  Posterior')

        ax1.set_ylim(bottom=0, top=ymax)
        ax2.set_ylim(bottom=0, top=ymax2)
        ax2.tick_params(axis='y', colors=(0.5, 0.5, 0.5))

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


"""

"""

base_model = Model(Da=0.1, Dp=0.1, Da_cyt=2, Dp_cyt=2, konA=0.1, koffA=0.01, kposA=0, konP=0.1,
                   koffP=0.01, kposP=0, kAP=10 ** -1, kPA=10 ** -1, ePneg=2, eAneg=2, xsteps=100, Tmax=300,
                   deltat=0.01, deltax=0.5, psi=0.1, am_0=np.zeros([100]), ac_0=np.ones([100]), pm_0=np.zeros([100]),
                   pc_0=np.ones([100]))
animate(base_model, p1=2, framerate=24, name='animPAR.mp4', ymax=7.5, ymax2=1.5)
