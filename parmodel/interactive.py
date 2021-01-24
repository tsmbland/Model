import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

"""
To do: option to display or save movie

"""


def animatePAR(direc):
    """
    Display an animation for PAR model

    direc: directory to results

    """

    # Load data
    times = np.loadtxt(direc + '/times.txt')
    A = np.loadtxt(direc + '/A.txt')
    P = np.loadtxt(direc + '/P.txt')

    # Set up figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.25, wspace=0.5)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Time (s)', 0, times[-1], valinit=0, valfmt='%d')
    ymax = 1.1 * max([np.max(A), np.max(P)])

    def update(i):
        ax.clear()
        tpoint = np.argmin(abs(times - int(i)))
        a = A[tpoint, :]
        p = P[tpoint, :]
        ax.plot(a, c='tab:red')
        ax.plot(p, c='tab:blue')
        ax.set_ylim(bottom=0, top=ymax)
        ax.set_ylabel('Cortical concentration (a.u.)')
        ax.set_xlabel('Position')

    sframe.on_changed(update)
    plt.show()
