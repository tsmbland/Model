import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import animation
import seaborn as sns
import os
import pickle

sns.set()
sns.set_style("ticks")

os.chdir(os.path.expanduser('~/Desktop/ModelData'))


############################## ANALYSIS ##########################


def loaddata(jobid, subjobid, simid):
    data = open(
        '%s/%s/%s.pkl' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid), '{0:04}'.format(simid)),
        'rb')
    res = pickle.load(data)
    return res


def countsims(jobid, subjobid):
    count = 0
    for root, dirs, files in os.walk('%s/%s' % ('{0:04}'.format(jobid), '{0:04}'.format(subjobid))):
        for file in files:
            if file.endswith('.pkl'):
                count += 1
    return count


def countsubjobs(jobid):
    count = 0
    for root, dirs, files in os.walk('%s' % ('{0:04}'.format(jobid))):
        for d in dirs:
            count += 1
    return count


def stats_batch(job, subjobs):
    nsubjobs = len(subjobs)

    res = loaddata(jobid=job, subjobid=0)
    tsteps = int(res.Tmax / res.deltat) + 1

    class data:
        asi = np.zeros([nsubjobs, tsteps])
        a_cyt = np.zeros([nsubjobs, tsteps])
        a_mem = np.zeros([nsubjobs, tsteps])
        a_mem_cyt = np.zeros([nsubjobs, tsteps])
        a_size = np.zeros([nsubjobs, tsteps])
        p_cyt = np.zeros([nsubjobs, tsteps])
        p_mem = np.zeros([nsubjobs, tsteps])
        p_mem_cyt = np.zeros([nsubjobs, tsteps])
        p_size = np.zeros([nsubjobs, tsteps])
        subjob = np.zeros([nsubjobs, tsteps])

    class labels:
        asi = 'Asymmetry index'
        a_cyt = 'A cytoplasmic concentration [μm⁻³]'
        a_mem = 'A domain concentration [μm⁻²]'
        a_mem_cyt = 'A membrane:cytoplasmic ratio'
        a_size = 'A domain size [μm]'
        p_cyt = 'P cytoplasmic concentration [μm⁻³]'
        p_mem = 'P domain concentration [μm⁻²]'
        p_mem_cyt = 'P membrane;cytoplasmic ratio'
        p_size = 'P domain size [μm]'

    for subjob in range(nsubjobs):
        res = loaddata(jobid=job, subjobid=subjobs[subjob])
        data.asi[subjob, :] = (2 * np.sum((np.sign(res.aco - res.pco) + 1) / 2, axis=1) - res.p.xsteps) / res.p.xsteps
        data.a_mem[subjob, :] = np.amax(res.aco, axis=1)
        data.a_cyt[subjob, :] = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
        data.a_mem_cyt[subjob, :] = data.a_mem[subjob] / data.a_cyt[subjob]
        data.a_size[subjob, :] = np.sum(res.aco.transpose() > (0.5 * np.tile(data.a_mem[subjob], [res.p.xsteps, 1])),
                                        axis=0) * res.p.L / res.p.xsteps
        data.p_mem[subjob, :] = np.amax(res.pco, axis=1)
        data.p_cyt[subjob, :] = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
        data.p_mem_cyt[subjob, :] = data.p_mem[subjob] / data.p_cyt[subjob]
        data.p_size[subjob, :] = np.sum(res.pco.transpose() > (0.5 * np.tile(data.p_mem[subjob], [res.p.xsteps, 1])),
                                        axis=0) * res.p.L / res.p.xsteps
        data.subjob[subjob, :] = subjob

    return data, labels


def stats(res):
    tsteps = int(res.p.Tmax / res.p.deltat) + 1

    class data:
        asi = np.zeros([tsteps])
        a_cyt = np.zeros([tsteps])
        a_mem = np.zeros([tsteps])
        a_mem_cyt = np.zeros([tsteps])
        a_size = np.zeros([tsteps])
        p_cyt = np.zeros([tsteps])
        p_mem = np.zeros([tsteps])
        p_mem_cyt = np.zeros([tsteps])
        p_size = np.zeros([tsteps])
        subjob = np.zeros([tsteps])

    class labels:
        asi = 'Asymmetry index'
        a_cyt = 'A cytoplasmic concentration [μm⁻³]'
        a_mem = 'A domain concentration [μm⁻²]'
        a_mem_cyt = 'A membrane:cytoplasmic ratio'
        a_size = 'A domain size [μm]'
        p_cyt = 'P cytoplasmic concentration [μm⁻³]'
        p_mem = 'P domain concentration [μm⁻²]'
        p_mem_cyt = 'P membrane;cytoplasmic ratio'
        p_size = 'P domain size [μm]'

    data.asi = (2 * np.sum((np.sign(res.aco - res.pco) + 1) / 2, axis=1) - res.p.xsteps) / res.p.xsteps
    data.a_mem = np.amax(res.aco, axis=1)
    data.a_cyt = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
    data.a_mem_cyt = data.a_mem / data.a_cyt
    data.a_size = np.sum(res.aco.transpose() > (0.5 * np.tile(data.a_mem, [res.p.xsteps, 1])),
                         axis=0) * res.p.L / res.p.xsteps
    data.p_mem = np.amax(res.pco, axis=1)
    data.p_cyt = (res.p.pA - res.p.psi * np.mean(res.aco, axis=1))
    data.p_mem_cyt = data.p_mem / data.p_cyt
    data.p_size = np.sum(res.pco.transpose() > (0.5 * np.tile(data.p_mem, [res.p.xsteps, 1])),
                         axis=0) * res.p.L / res.p.xsteps

    return data, labels


############################## PLOTS ##############################


def parplot(i, ax, aco, pco, p):
    ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, aco[i, :], label='A', c='r')
    ax.plot(np.array(range(p.xsteps)) / (p.xsteps - 1) * p.L, pco[i, :], label='P', c='dodgerblue')
    ax.set_xlabel('x [μm]')
    ax.set_ylabel('Concentration [a.u.]')
    concmax = max(aco.max(), pco.max())
    ax.set_ylim(0, 1.1 * concmax)
    # ax.legend()


def plot_singlesim(jobid=0, subjobid=0, simid=0):
    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    parplot(-1, ax, res.aco, res.pco, res.p)
    plt.show()


def sliderplot(jobid=0, subjobid=0, simid=0):
    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Time (s)', 0, res.params.Tmax, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        i = int(sframe.val / res.params.deltat)
        parplot(i, ax, res.aco, res.pco, res.params)
        ax.set_title('Time (s): {0:.0f}'.format(i * res.params.deltat))

    parplot(0, ax, res.aco, res.pco, res.params)
    ax.set_title('Time (s): {0:.0f}'.format(0 * res.params.deltat))
    sframe.on_changed(update_slider)
    plt.show()


def anim(jobid=0, subjobid=0, simid=0, animrate=100, framerate=24):
    """
    :param animrate: seconds of model time per second of animation
    :param framerate: frames per second

    """

    res = loaddata(jobid, subjobid, simid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    frames = range(0, int(res.p.Tmax / res.p.deltat), int(animrate / (framerate * res.p.deltat)))

    def update_anim(i):
        ax.clear()
        parplot(i, ax, res.aco, res.pco, res.p)
        ax.set_title('Time (s): {0:.0f}'.format(i * res.p.deltat))

    paranim = animation.FuncAnimation(fig, update_anim, frames=iter(frames), save_count=len(frames))
    writer = animation.writers['ffmpeg']
    writer = writer(fps=framerate, bitrate=2000)
    paranim.save('animation.mp4', writer=writer)


def paramsliderplot(jobid, subjobid, param):
    sims = countsims(jobid, subjobid)
    xsteps = loaddata(jobid=jobid, subjobid=subjobid, simid=0).p.xsteps
    params = np.zeros(sims)
    aco = np.zeros([sims, xsteps])
    pco = np.zeros([sims, xsteps])

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        params[simid] = getattr(res.p, param)
        aco[simid, :] = res.aco[-1, :]
        pco[simid, :] = res.pco[-1, :]

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, param, 0, max(params), valinit=0)

    def update_slider(val):
        ax.clear()
        val = sframe.val
        subjob = (np.abs(params - val)).argmin()
        parplot(subjob, ax, aco, pco, res.p)

    parplot(0, ax, aco, pco, res.p)
    sframe.on_changed(update_slider)
    plt.show()


def paramsliderplot2(jobid, subjobid, param1, param2):
    sims = countsims(jobid, subjobid)
    xsteps = loaddata(jobid, subjobid, 0).p.xsteps
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    aco = np.zeros([sims, xsteps])
    pco = np.zeros([sims, xsteps])

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)
        aco[simid, :] = res.aco[-1, :]
        pco[simid, :] = res.pco[-1, :]

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    param1_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
    param1_slider = Slider(param1_slider_ax, param1, 0, max(param1vals), valinit=0)
    param2_slider_ax = plt.axes([0.25, 0.05, 0.65, 0.03])
    param2_slider = Slider(param2_slider_ax, param2, 0, max(param2vals), valinit=0)

    def update_slider(val):
        ax.clear()
        param1_val = param1_slider.val
        param2_val = param2_slider.val
        subjob = (np.abs(param1vals - param1_val) * np.abs(param2vals - param2_val)).argmin()
        parplot(subjob, ax, aco, pco, res.p)

    parplot(0, ax, aco, pco, res.p)
    param1_slider.on_changed(update_slider)
    param2_slider.on_changed(update_slider)
    plt.show()


def dataplot(jobid, subjobid, param, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    paramvals = np.zeros(sims)

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        paramvals[simid] = getattr(res.p, param)

    plt.plot(paramvals, getattr(d, var)[:, -1])
    plt.ylabel(getattr(labels, var))
    plt.xlabel(param)
    plt.show()


def dataplot2(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)

    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    vals = np.unique(param2vals)
    for param2val in range(len(vals)):
        print(param2vals[param2val])
        ax.plot(param1vals[param2vals == vals[param2val]],
                getattr(d, var)[:, -1][param2vals == vals[param2val]])
    # cbar = plt.colorbar(ax)
    # cbar.set_label(param2)
    plt.ylabel(getattr(labels, var))
    plt.xlabel(param1)
    plt.show()


def dataplot2_heatmap(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))
    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    sc = plt.scatter(param1vals, param2vals, c=getattr(d, var)[:, -1])
    cbar = plt.colorbar(sc)
    cbar.set_label(getattr(labels, var))
    plt.xlabel(param1)
    plt.ylabel(param2)

    plt.show()


def dataplot2_heatmap_180220(jobid, subjobid, param1, param2, var):
    sims = countsims(jobid, subjobid)
    (d, labels) = stats_batch(jobid, range(sims))

    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)

    asi = getattr(d, var)[:, -1]

    sc = plt.scatter(param1vals[abs(asi) < 0.25], param2vals[abs(asi) < 0.25])
    # cbar = plt.colorbar(sc)
    # cbar.set_label(getattr(labels, var))
    plt.xlabel(param1)
    plt.ylabel(param2)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.show()


def dataplot2_heatmap_180529(jobid, subjobid, param1, param2):
    sims = countsims(jobid, subjobid)

    param1vals = np.zeros(sims)
    param2vals = np.zeros(sims)
    pols = np.zeros(sims)
    for simid in range(sims):
        res = loaddata(jobid, subjobid, simid)
        param1vals[simid] = getattr(res.p, param1)
        param2vals[simid] = getattr(res.p, param2)
        pols[simid] = test_180529(res)

    plt.scatter(param1vals[abs(pols) == 0], param2vals[abs(pols) == 0], c='silver', marker='s')
    plt.scatter(param1vals[abs(pols) == 1], param2vals[abs(pols) == 1], c='dodgerblue', marker='s')
    plt.scatter(param1vals[abs(pols) == 2], param2vals[abs(pols) == 2], c='r', marker='s')

    plt.xlabel(param1)
    plt.ylabel(param2)
    # plt.xlim([0, 1])
    # plt.ylim([0, 1])
    plt.show()


def genalg_plot1(jobid, base=None):
    nsubjobs = countsubjobs(jobid)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Generation', 0, nsubjobs, valinit=0, valfmt='%d')

    def update_slider(val):
        ax.clear()
        subjobid = int(sframe.val)
        nsims = countsims(jobid, subjobid)

        # Plot base data
        if base is not None:
            res = loaddata(base[0], base[1], base[2])
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.aco[-1, :], c='g')
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.pco[-1, :], c='r')
            concmax = max(res.aco.max(), res.pco.max())
            ax.set_ylim(0, 1.5 * concmax)

        # Plot genentic algorithm simulations
        for simid in range(nsims):
            res = loaddata(jobid, subjobid, simid)
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.aco[-1, :], c='g', alpha=0.1)
            ax.plot(np.array(range(res.params.xsteps)) / (res.params.xsteps - 1) * res.params.L, res.pco[-1, :], c='r', alpha=0.1)

        ax.set_xlabel('Anterior - Posterior [μm]')
        ax.set_ylabel('Concentration [μm⁻²]')

    sframe.on_changed(update_slider)
    plt.show()


def genalg_spider(jobid, subjobid, params):
    # Set up plot
    angles = [n / float(len(params)) * 2 * np.pi for n in range(len(params))]
    angles += angles[:1]
    ax = plt.subplot(111, polar=True)
    plt.xticks(angles[:-1], params)
    ax.set_rlabel_position(0)

    # Plot data
    nsims = countsims(jobid, subjobid)
    for simid in range(nsims):
        values = np.zeros(len(params) + 1)
        res = loaddata(jobid, subjobid, simid)
        for param in range(len(params)):
            values[param] = getattr(res.p, params[param])
        values[-1] = values[0]
        ax.plot(angles, values)
    plt.show()


def genalg_spider_slider(jobid, params):
    nsubjobs = countsubjobs(jobid)

    plt.clf()
    fig = plt.figure()
    ax = plt.subplot(111, polar=True)
    axframe = plt.axes([0.25, 0.1, 0.65, 0.03])
    sframe = Slider(axframe, 'Generation', 0, nsubjobs, valinit=0, valfmt='%d')

    angles = [n / float(len(params)) * 2 * np.pi for n in range(len(params))]
    angles += angles[:1]

    def update_slider(val):
        ax.clear()
        subjobid = int(sframe.val)
        ax.set_xticks(angles[:-1], params)
        ax.set_rlabel_position(0)

        # Plot data
        nsims = countsims(jobid, subjobid)
        for simid in range(nsims):
            values = np.zeros(len(params) + 1)
            res = loaddata(jobid, subjobid, simid)
            for param in range(len(params)):
                values[param] = getattr(res.params, params[param])
            values[-1] = values[0]
            ax.plot(angles, values)

    sframe.on_changed(update_slider)
    plt.show()
