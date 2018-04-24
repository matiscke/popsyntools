""" Plotting routines for planet population synthesis and analysis of results
from the planet formation Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plotstyle


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
    """ Compute the bin edges for a logarithmic grid.

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

    Example
    -------
    >>> binWidth_dex = 1.0
    >>> Range = (10., 1000.)
    >>> compute_logbins(binWidth_dex, Range)
    array([   10.,   100.,  1000.])
    """
    # add binWidth_dex to logrange to include last bin edge
    logRange = (np.log10(Range[0]), np.log10(Range[1]) + binWidth_dex)
    return 10**np.arange(logRange[0], logRange[1], binWidth_dex)


def r_Jup2r_Earth(r):
    """ Transform a radius given in Jupiter radii into Earth radii.
    """
    return r*10.973


def plot_occurrence(population, ax=None, xAxis='period', yAxis='r', nBins=0,
                    binWidth_dex=(0.25, 0.1), xRange=None, yRange=None,
                    smooth=False, normalize=True, discreteColors=False,
                    logColormap=False, annotated=False, **funcKwargs):
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
    xRange : sequence of scalars
        range of values to be considered in x direction
    yRange : sequence of scalars
        range of values to be considered in y direction
    smooth : Bool
        if True, apply Gaussian filter to the histogram
    normalize : Bool
        normalize occurrence to planets per 100 stars
    discreteColors : Bool
        use discrete color levels instead of a continuum colormap
    logColormap : Bool
        use logarithmic (log10) mapping of colors to data values. Has no effect
        if discreteColors = True.
    annotated : Bool
        show the value of each bin in the histogram
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

    # check existence of columns in the DataFrame
    if not (xAxis in population and yAxis in population):
        raise KeyError('population does not contain both columns for the histogram')

    # create new column with radius in R_Earth if not existent
    if (yAxis == 'r' or yAxis == 'r_rEarth'):
        if not 'r_rEarth' in population:
            population['r_rEarth'] = r_Jup2r_Earth(population.r)
        yAxis = 'r_rEarth'

    try:
        # if DataFrame has a column 'status', use only survived planets
        survivedPlanets = population[population['status'] == 0]
        print('using only planets with status "0"')
    except KeyError:
        survivedPlanets = population

    if not ax:
        fig, ax = plt.subplots()

    # define the bins
    if not xRange:
        xRange = (survivedPlanets[xAxis].min(), survivedPlanets[xAxis].max())
    if not yRange:
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

    # create 2D histogram
    h, xedges, yedges = np.histogram2d(survivedPlanets[xAxis],
                        survivedPlanets[yAxis], bins=(xBins, yBins))
    h = h.T

    if smooth:
        # smooth out the contours
        import scipy.ndimage as nd
        h = nd.gaussian_filter(h,(4,2))

    if normalize:
        # normalize to 1/100stars
        Nsystems = survivedPlanets.systemNo.nunique()
        print('Number of Systems: {}'.format(Nsystems))
        h = h*100/Nsystems
        cbarlabel = r"Planets per 100 Stars per $P-R_P$ interval"
    else :
        cbarlabel = r"Planets per $P-R_P$ interval"

    if logColormap:
        # logarithmic color mapping. Use linear scale around zero.
        import matplotlib.colors as colors
        threshold  = 0.001
        colorNorm = colors.SymLogNorm(vmin=h.min(), vmax=h.max(),
        linthresh=max(h.min(), threshold))
    else:
        colorNorm = None

    if annotated:
        import seaborn as sns
        ax = sns.heatmap(h, annot=True, fmt=".1f", xticklabels=xedges)
        ax.invert_yaxis()

        #get ticks right
        ax.set(xticks=range(len(xedges)), xticklabels=xedges,
               yticks=range(len(yedges)), yticklabels=yedges)
        return h, xedges, yedges, ax

    if discreteColors:
        """use discrete levels for occurrence. numbers are from
        Petigura et al. 2018
        """
        # levels = np.arange(-4, -1 + 1e-10, 0.25)
        cbarticklabels = [0.01, 0.03, 0.1, 0.3, 1, 3,10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        contourKwargs = dict(extend='min')
        im = plt.contourf(xedges[:-1], yedges[:-1], h, **contourKwargs,
                          **funcKwargs)
    else:
        cbarticks = None
        X, Y = np.meshgrid(xedges, yedges)
        im = ax.pcolormesh(X, Y, h, norm=colorNorm, **funcKwargs)

    # eyecandy
    plt.xscale('log')
    plt.yscale('log')
    cbar = fig.colorbar(im)
    cbar.set_label(cbarlabel, labelpad=15)
    if xAxis == 'period':
        ax.set_xlabel('Orbital Period [d]')
    else:
        plt.xlabel(xAxis)
    if yAxis == 'r_rEarth':
        ax.set_ylabel('Planet Size [$\mathrm{R_{Earth}}$]')
    else:
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
