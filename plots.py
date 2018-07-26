""" Plotting routines for planet population synthesis and analysis of results
from the planet formation Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

import plotstyle
import utils


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
                    kind='hist', smooth=False, normalize=True,
                    logColormap=False, **funcKwargs):
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
    kind : string
        the kind of plot to produce
        - 'hist' : 2D histogram
        - 'contour' : contour plot with discrete color levels
        - 'annotated' : 2D histogram with values written in the bins
    smooth : Bool
        if True, apply Gaussian filter to the histogram
    normalize : Bool
        normalize occurrence to planets per 100 stars
    logColormap : Bool
        use logarithmic (log10) mapping of colors to data values. Has no effect
        on a contour plot.
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

    # obtain number of systems
    if 'isystem' in survivedPlanets:
        Nsystems = survivedPlanets.isystem.nunique()
    else:
        print("""Column "isystem" missing - old ref_red file? Using column
              "systemNo" instead. Please check correct normalization.""")
        Nsystems = survivedPlanets.systemNo.nunique()
    print('Number of Systems: {}'.format(Nsystems))

    if normalize:
        # normalize to 1/100stars
        h = h*100/Nsystems
        cbarlabel = r"Planets per 100 Stars per $P-R_P$ interval"
    else:
        cbarlabel = r"Planets per $P-R_P$ interval"

    if logColormap:
        # logarithmic color mapping. Use linear scale around zero.
        import matplotlib.colors as colors
        threshold  = 0.01
        colorNorm = colors.SymLogNorm(vmin=h.min(), vmax=h.max(),
        linthresh=max(h.min(), threshold))
        cbar_kws = {'label' : cbarlabel, 'ticks' : [.01, .1, 1., 1e1, 1e2, 1e3]}
    else:
        colorNorm = None
        cbar_kws = {'label' : cbarlabel}

    if kind == 'annotated':
        # plot a histogram with values written in the bins.
        import seaborn as sns
        ax = sns.heatmap(h, annot=True, fmt=".1f", norm=colorNorm,
                         cbar_kws = cbar_kws, **funcKwargs)
        ax.invert_yaxis()

        # get axis ticks right
        try:
            ticklabels = [None, None]
            for dim, edges in enumerate([xedges, yedges]):
                ticks = [" " for i in range(len(edges))]
                keptTicks = edges[::int(len(edges)/10)]
                ticks[::int(len(edges)/10)] = ['{:.1f}'.format(t) for t in keptTicks]
                ticklabels[dim] = ticks
            ax.set(xticks=range(len(ticklabels[0])), xticklabels=ticklabels[0],
                   yticks=range(len(ticklabels[1])), yticklabels=ticklabels[1])
        except ValueError:
            pass
        plt.yticks(rotation=0)
    elif kind == 'contour':
        """use discrete levels for occurrence. numbers are from
        Petigura et al. 2018
        """
        # levels = np.arange(-4, -1 + 1e-10, 0.25)
        cbarticklabels = [0.01, 0.03, 0.1, 0.3, 1, 3,10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        contourKwargs = dict(extend='min')
        im = plt.contourf(xedges[:-1], yedges[:-1], h, **contourKwargs)
        if logColormap:
            print("logarithmic shading not possible with contours.")
    else:
        # plot a normal 2D histogram
        X, Y = np.meshgrid(xedges, yedges)
        im = ax.pcolormesh(X, Y, h, norm=colorNorm, **funcKwargs)

    # eyecandy
    if not kind == 'annotated':
        cbar = fig.colorbar(im)
        cbar.set_label(cbarlabel, labelpad=15)
        plt.xscale('log')
        plt.yscale('log')
        if xAxis == 'period':
            ax.set_xlabel('Orbital Period [d]')
        else:
            plt.xlabel(xAxis)
        if yAxis == 'r_rEarth':
            ax.set_ylabel('Planet Size [$\mathrm{R_{Earth}}$]')
        else:
            plt.ylabel(yAxis)

    return h, xedges, yedges, ax


def compare_surfaceDensity(disk1file, disk2file, sim1name="Type 2",
                           sim2name="Type 3"):
    """ compare the surface density as a function of time of two simulations.

    Parameters
    ----------
    disk1file : string
        path to the 'disk_structure' file of disk 1
    disk2file : string
        path to the 'disk_structure' file of disk 2
    sim1name : string
        Name of simulation 1 for the legend
    sim2name : string
        Name of simulation 2 for the legend

    Returns
    -------
    ax : array
        matplotlib axes with the plot
    """

    disk1 = pd.read_csv(disk1file, delim_whitespace=True,
                        header=None).rename(columns={9:'t', 1:'r', 2:'sigma'})
    disk2 = pd.read_csv(disk2file, delim_whitespace=True,
                        header=None).rename(columns={9:'t', 1:'r', 2:'sigma'})

    T1 = [t for t in disk1.t.unique()]
    T2 = [t for t in disk2.t.unique()]

    # compute difference of surface density
    disk1['sigDiff'] = disk2['sigma'] - disk1['sigma']

    fig, ax = plt.subplots(2, sharex=True)
    for t1 in T1:
        ax[0].plot(disk1[disk1['t'] == t1]['r'], disk1[disk1['t'] == t1]['sigma'],
                   c=cm.Blues(t1/.4/max(T1)), label='{:1.0E}'.format(t1))
        ax[0].plot(disk2[disk2['t'] == t1]['r'], disk2[disk2['t'] == t1]['sigma'],
                   c=cm.Greens(t1/.4/max(T1)))
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    # plot difference
    for t1 in T1:
        ax[1].plot(disk1[disk1['t'] == t1]['r'], disk1[disk1['t'] == t1]['sigDiff'],
                   c=cm.Reds(t1/.4/max(T1)))

    ax[1].set_xlabel('r [au]')
    ax[0].set_ylabel('$\Sigma$ [$\mathrm{g}/\mathrm{cm}^2$]')
    ax[1].set_ylabel('$\Delta \Sigma$ [$\mathrm{g}/\mathrm{cm}^2$]')

    # custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=cm.Blues(0.5), lw=2),
                    Line2D([0], [0], color=cm.Greens(0.5), lw=2),
                    Line2D([0], [0], color=cm.Reds(0.5), lw=2)]
    ax[0].legend(custom_lines[:2], ['$\Sigma$ ({})'.format(sim1name), '$\Sigma$ ({})'.format(sim2name)])
    ax[1].legend(custom_lines[2:], ['$\Sigma$ ({}) - $\Sigma$ ({})'.format(sim2name, sim1name)])

    return ax


def plot_planetTracks(simulation, truths=None, lwRange = (2., 40.)):
    """
    plot evolutionary tracks of a multiplanet system.

    The function plots semimajor axis as a function of time, displaying planet
    radii as varying line widths and planet mass as color code. The today-values
    of a real system can be added for comparison.

    Parameters
    ----------
    simulation : dictionary
        dictionary of pandas DataFrames, one for each planet, as returned by
        output.read_simFromFolder
    truths : pandas DataFrame
        DataFrame containing the planet parameters of a real system. Has to have
        columns for Semimajor axis 'a', planet radius 'r', and planet mass 'm'.
    lwRange : tuple
        the range of line widths for the plot. Widths vary along the tracks
        according to the current radius of a planet. The values given in this
        parameter correspond to the minimum and maximum planet radius in the
        whole simulation.

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    """
    import types
    import matplotlib.colors as colors
    from matplotlib.backend_bases import GraphicsContextBase, RendererBase

    # dirty hack to get round line caps
    class GC(GraphicsContextBase):
        def __init__(self):
            super().__init__()
            self._capstyle = 'round'

    def custom_new_gc(self):
        return GC()
    RendererBase.new_gc = types.MethodType(custom_new_gc, RendererBase)


    from matplotlib.collections import LineCollection
    fig, ax = plt.subplots()

    # find max planet mass for color code, radius range for line width
    maxM = max([max(simulation[df]['m']) for df in simulation])
    print(maxM)
    colorNorm = colors.LogNorm(0.01, maxM)
    rRange = (min([min(simulation[df]['r']) for df in simulation]),
              max([max(simulation[df]['r']) for df in simulation]))

    def radius2linewidth(r):
        return (r - rRange[0])/(rRange[1] - rRange[0])*(lwRange[1] - lwRange[0]) + lwRange[0]

    for df in simulation:
        planet = simulation[df]

        # exclude planets that never grew
        if max(planet['m']) < .2:
            print('exclude {}: too small'.format(planet.name))
            continue

        # restrict data to times with planet status = 0
        planet = planet[planet['status'] == 0]

        t = planet['t'][::5]
        a = planet['a'][::5]
        m = planet['m'][::5]
        r = planet['r'][::5]

        # linearly transform line widths into radius range
        lw = radius2linewidth(r)

        # create a set of line segments
        points = np.array([t, a]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, linewidths=lw, cmap='inferno_r', norm=colorNorm)
        lc.set_array(m)
        ax.add_collection(lc)

    if not truths is None:
        ax.scatter([max(t) + 1e6 for i in range(len(truths))], truths['a'], s=radius2linewidth(truths['r']),
                   c=truths['m'], cmap='inferno_r', norm=colorNorm)

    # add a colorbar
    cbar = fig.colorbar(lc)
    cbar.set_label('Planet Mass [$\mathrm{M_{Earth}}$]')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Semimajor Axis [au]')

    return ax


def plot_planetTypeBar(jointPop, ax=None, planetTypes = ['Earth', 'SuperEarth',
                   'Neptune', 'SubGiant', 'Giant', 'BrownDwarf']):
    """ Plot a categorical bar plot showing abundance of planet types.

    Parameters
    ----------
    jointPop : pandas multiindex DataFrame
        DataFrame with the host star mass as a multiindex.
    ax : matplotlib axis
        axis to draw the plot on
    planetTypes : list or array
        planet types to consider for the plot. The elements of this list must
        exist in config.massLimits.

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    """
    if ax == None:
        fig, ax = plt.subplots()

    M = [np.float(Mstar[:3]) for Mstar, df in jointPop.groupby(level=0)]
    xPos = [x for x in range(len(M))]
    barWidth = .95

    frequencies = []
    for n, planetType in enumerate(planetTypes):
        frequencies.append([len(subPop[subPop['planetType'] == planetType])
                            for i, subPop in jointPop.groupby(level=0)])
        if n == 0:
            ax.bar(xPos, frequencies[n], edgecolor='white', width=barWidth,
                   label=planetType)
        else:
            ax.bar(xPos, frequencies[n], bottom=np.sum(frequencies[:n], axis=0),
                   edgecolor='white', width=barWidth, label=planetType)

    plt.xticks(xPos, M)
    plt.legend()
    ax.set_xlabel('$M_{\star}$')
    ax.set_ylabel('Number of planets')

    return ax


def plot_diskFractions(diskFractions, ax=None, fractionLimits=(0.85, 0.01)):
    """Plot fractions of remaining systems with disks as a function of time.

    The fractions are fitted by an exponential function.

    Parameters
    ----------
    diskFractions : tuple or dict
        tuple containing arrays with times and corresponding disk fractions. If
        several populations are to be plotted, diskFractions has to be a
        dictionary containing a tuple for each population.
    ax : matplotlib axis
        axis to draw the plot on
    fractionLimits : tuple
        upper and lower bound for the range of disk fractions considered for the
        fit

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    diskFractions : tuple or dict
        updated array (fit results appended in the case of multiple populations)
    tau : float or array
        time constant(s) of disk lifetime
    std_tau : float array
        1 sigma error(s) of tau
    """

    def get_fitRangeMask(fractions, fractionLimits):
        """create mask to constrain range for the fit."""
        return np.ma.masked_inside(fractions, *fractionLimits)

    if ax == None:
        fig, ax = plt.subplots()

    tau = []
    std_tau = []
    if isinstance(diskFractions, dict):
        for Mstar, fractions in diskFractions.items():
            # mask times and fit exponential function
            fitRangeMask = get_fitRangeMask(fractions[1], fractionLimits)
            fitParams = utils.fit_diskFractions(fractions[0][fitRangeMask.mask],
                                                fractions[1][fitRangeMask.mask])
            diskFractions[Mstar] += fitParams

            tau.append(fitParams[0][1])
            std_tau.append(np.sqrt(np.diag(fitParams[1]))[1])

            # plot fractions
            ax.scatter(fractions[0], fractions[1], marker='x',
                       label=Mstar[:-4] + r'$\mathrm{M_\odot}$' +
                       r', $\tau = {:1.1f}$ Myr'.format(tau[-1]/1e6))

            # plot fits
            ax.plot(fractions[0][fitRangeMask.mask],
                    utils.exponential(fractions[0][fitRangeMask.mask],
                    fitParams[0][0], fitParams[0][1], fitParams[0][2]))
        plt.legend()

    else:
        fitRangeMask = get_fitRangeMask(diskFractions[1], fractionLimits)
        fitParams = utils.fit_diskFractions(diskFractions[0][fitRangeMask.mask],
                                            diskFractions[1][fitRangeMask.mask])
        tau = fitParams[0][1]
        std_tau = np.sqrt(np.diag(fitParams[1]))[1]

        ax.scatter(diskFractions[0], diskFractions[1], marker='x')
        ax.plot(diskFractions[0][fitRangeMask.mask],
                utils.exponential(diskFractions[0][fitRangeMask.mask],
                fitParams[0][0], fitParams[0][1], fitParams[0][2]))

    ax.axhline(.5, ls='--', c='gray')
    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Disk Fraction')
    return ax, diskFractions, tau, std_tau


def plot_diskLifetimeComparison(Mstar, tau, std_tau, ax=None, nSigma=1,
                                **funcKwargs):
    """Compare the time constants of several populations' disk lifetimes.

    Parameters
    ----------
    Mstar : array
        host star masses
    tau : array
        corresponding time constants
    std_tau : array
        1 sigma error of tau
    ax : matplotlib axis object, optional
        axis to plot on
    nSigma : int, optional
        number of standard deviations for confidence intervals
    **funcKwargs : keyword arguments
        kwargs to pass on to matplotlib

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    """
    from scipy.optimize import curve_fit

    if ax == None:
        fig, ax = plt.subplots()

    params, covariance = curve_fit(utils.linear, Mstar, tau, sigma=std_tau)
    pErr = np.sqrt(np.diag(covariance))*nSigma

    # prepare confidence intervals
    pHi = params + pErr
    pLo = params - pErr

    fit = utils.linear(Mstar, *params)
    fitHi = utils.linear(Mstar, *pHi)
    fitLo = utils.linear(Mstar, *pLo)

    (_, caps, _) = ax.errorbar(Mstar, tau, std_tau, lw=2, capsize=3, fmt='o',
                               ms=4, label=r'$\tau ( \mathrm{M_\star})$')

    # dirty hack to recover errorbar caps
    for cap in caps:
        cap.set_markeredgewidth(1)

    ax.plot(Mstar, fit, c='gray', label=r'$\tau = ({:1.1f}*\mathrm{{M_\star}} + {:1.1f}$) Myr'.format(*params/1e6))
    ax.fill_between(Mstar, fitHi, fitLo, alpha=.15,
                    label=r'$1\sigma$ confidence interval', **funcKwargs)


    ax.set_xlabel(r'Host Star Mass $[M_\odot]$')
    ax.set_ylabel(r'Time Constant $\tau$')
    ax.legend()

    return ax

################################################################################

""" Plotting functions meant for single planet tracks."""

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
