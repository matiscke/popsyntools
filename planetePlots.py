""" Plotting routines for planet population synthesis and analysis of results
from the planet formation Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import numpy as np
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

def compute_logbins(binWidth_dex, Range):
    """Compute the bin edges for a logarithmic grid.

    Parameters
    ----------
    binWidth_dex : float
        width of bins in log space (dex)
    Range : Tuple
        range for parameter

    Returns
    -------
    bins : array
        bins for one dimension
    """
    logRange = (np.log10(Range[0]), np.log10(Range[1]))
    return 10**np.arange(logRange[0], logRange[1], binWidth_dex)

def plot_occurrence(population, ax=None, xAxis='period', yAxis='r',
                    nBins=0, binWidth_dex=(0.25, 0.1), **funcKwargs):
    """Plot an occurrence map in two parameters.

    Parameters
    ----------
    population : pandas DataFrame
        planet population to plot
    ax : matplotlib axis
        axis to plot on
    xAxis : string
        parameter for the x axis
    yAxis : string
        parameter for the y axis
    nBins : integer
        number of bins for each axis. Only relevant if a positive integer is
        given, otherwise bins are defined via `binWidth_dex`.
    binWidth_dex : float or sequence of scalars
        width of each bin in dex for [xAxis, yAxis].
        If `binWidth_dex` is a scalar, it defines the bin width along both axes.
    **funcKwargs : keyword arguments
        kwargs to pass on to matplotlib


    Returns
    -------
    h : numpy array
        normalized 2D histogram
    xedges : numpy array
        bin edges along x axis
    yedges : numpy array
        bin edges along y axis
    ax : matplotlib axis
        axis with the plot
    """
    import scipy.ndimage as nd

    try:
        # if DataFrame has a column 'status', use only survived planets
        survivedPlanets = population[population['status'] == 0]
        print('using only planets with status "0"')
    except KeyError:
        survivedPlanets = population

    if not ax:
        fig, ax = plt.subplots()

    # define the bins
    xRange = (survivedPlanets[xAxis].min(), survivedPlanets[xAxis].max())
    yRange = (survivedPlanets[yAxis].min(), survivedPlanets[yAxis].max())
    if nBins:
        # logarithmic bins of equal width
        xBins = np.logspace(np.floor(np.log10(xRange[0])),
                            np.ceil(np.log10(xRange[1])), nBins)
        yBins = np.logspace(np.floor(np.log10(yRange[0])),
                            np.ceil(np.log10(yRange[1])), nBins)
    else:
        # define bins by their width
        if not np.iterable(binWidth_dex):
            # if only one number is given, use along both dimensions
            binWidth_dex = (binWidth_dex, binWidth_dex)
        xBins = compute_logbins(binWidth_dex[0], xRange)
        yBins = compute_logbins(binWidth_dex[1], yRange)


# DEBUGGING
    print('xRange {}'.format(xRange))
    print('yRange {}'.format(yRange))
    print('xBins {}'.format(xBins[-5:]))
    print(np.isnan(survivedPlanets[yAxis]).any())
    # print(xBins)
    # print(yBins)
    # return (xBins, yBins)

    # create 2D histogram
    h, xedges, yedges = np.histogram2d(survivedPlanets[xAxis],
                        survivedPlanets[yAxis], bins=(xBins, yBins))
    h = h.T

    # smooth out the contours and normalize to 1/100stars
    nd.gaussian_filter(h,(40,2))
    Nsystems = len(survivedPlanets)
    h = h*100/Nsystems

    X, Y = np.meshgrid(xedges, yedges)
    im = ax.pcolormesh(X, Y, h, **funcKwargs)

    # eyecandy
    plt.xscale('log')
    plt.yscale('log')
    cb = fig.colorbar(im)
    plt.xlabel(xAxis)
    plt.ylabel(yAxis)
    return h, xedges, yedges, ax


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
