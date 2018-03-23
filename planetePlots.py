""" Plotting routines for planet population synthesis and analysis of results
from the planet formation Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set plot style
sns.set(context='notebook', style='whitegrid', font_scale=1., palette='colorblind',
        rc={
'text.usetex':True,
'text.latex.unicode':True,
'font.family' : 'sans-serif',
'font.style'         : 'normal',
'font.variant'        : 'normal',
'font.weight'         : 'normal',
'font.stretch'        : 'normal',
'savefig.dpi'         : 400,
'lines.linewidth'   : 1.0,
'lines.markersize'      : 3.,
'figure.subplot.left'    : 0.13,    # the left side of the subplots of the figure
'figure.subplot.right'   : 0.96,   # the right side of the subplots of the figure
'figure.subplot.bottom'  : 0.13,   # the bottom of the subplots of the figure
'figure.subplot.top'     : 0.96,    # the top of the subplots of the figure
'figure.subplot.hspace'  : 0.0,    # height reserved for space between subplots
'axes.xmargin' : 0.02,             # default margin for autoscale
'axes.ymargin' : 0.02,
'legend.handletextpad' : 0.5,
'legend.handlelength' : 0.75,
'xtick.minor.size'     : 2.,
'ytick.minor.size'     : 2.
})
sns.set_color_codes()


def normalize_rate(n_planet, n_star):
    """ normalize the occurrence rate to planets per 100 stars.

    Parameters
    ----------
    n_planet : int
        number of planets
    n_star : int
        number of stars

    Returns
    -------
    norm_rate : float
        normalized occurrence rate
    """
    norm_rate = 100*n_planet/n_star
    return norm_rate


def plot_occurrence(population, ax=None, xAxis='a', yAxis='r',*funcArgs, **funcKwargs):
    """Plot an occurrence map in two parameters.

    Parameters
    ----------
    population : pandas DataFrame
        planet population to plot
    keyword arguments

    Returns
    -------
    g : JointGrid
        seaborn JointGrid object with the plot on it
    """
    try:
        # if DataFrame has a column 'status', use only survived planets
        survivedPlanets = population[population['status'] == 0]
        print('using only planets with status "0"')
    except KeyError:
        survivedPlanets = population

    # # clip: do not allow negative values
    # g = sns.jointplot(xAxis, yAxis, data=survivedPlanets, kind="kde", color="m",
    #                   clip=((0.,1e12),(0.,1e12)), stat_func=None)
    if not ax:
        fig, ax = plt.subplots()

    # sort by systems
    Nsystems = survivedPlanets.systemNo.max()
    systems = survivedPlanets.set_index(['systemNo'])


    h, xedges, yedges, image = plt.hist2d(xAxis, yAxis,
        data=systems, bins=100)

    # H, xedges, yedges= np.histogram2d(x, y, bins=100)
    # implot = plt.imshow(H, interpolation='nearest', origin='low')



    # # overplot data points
    # g.plot_joint(plt.scatter, c="b", s=2, marker=".", alpha = .1)
    # g.ax_joint.collections[0].set_alpha(0)
    # g.set_axis_labels(xAxis,yAxis)
    # ax = g.ax_joint

    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # g.ax_marg_x.set_xscale('log')
    # g.ax_marg_y.set_yscale('log')
    return h, xedges, yedges


""" Plotting functions meant for single planet tracks.
"""
def plot_mass(tracks, ax):
    """plot a planet's total mass vs time"""
    ax.plot(tracks['t'],tracks['m'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_coreMass(tracks, ax):
    """plot core mass vs time"""
    ax.plot(tracks['t'],tracks['mCore'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('core mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_radius(tracks, ax):
    """plot radius vs time"""
    ax.plot(tracks['t'],tracks['r'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('radius [Jupiter radii]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_lum(tracks, ax):
    """ plot luminosity vs time"""
    ax.plot(tracks['t'], tracks['L'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('Luminosity [?]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax
