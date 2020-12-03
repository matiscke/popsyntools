""" Plotting routines for planet population synthesis and analysis of results
from the planet formation Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
import astropy
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

from popsyntools import plotstyle, config, utils


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


def fix_hist_step_vertical_line_at_end(ax):
    """remove vertical line at end of CDFs"""
    axpolygons = [poly for poly in ax.get_children() if isinstance(poly, mpl.patches.Polygon)]
    for poly in axpolygons:
        poly.set_xy(poly.get_xy()[:-1])

def log_col(pop, variables):
    """apply a log10 to columns and write result to new columns named <name>_log"""
    for v in variables:
        pop.loc[:,'{}_log'.format(v)] = np.log10(pop[v])
    return pop

def plot_occurrence(population, ax=None, xAxis='period', yAxis='r', nBins=0,
                    binWidth_dex=(0.25, 0.1), xRange=None, yRange=None,
                    zRange=None, kind='hist', smooth=False, normalize=True,
                    logColormap=False, showColorbar=True, **funcKwargs):
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
    zRange : sequence of scalars
        range of values for the color code. It's first value (lower bound)
        should be greater than zero for logarithmic color mapping.
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
    showColorbar : Bool
        option to hide the colorbar, e.g. if one bar for several subplots shall
        be shown.
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
    im : matplotlib AxesImage
        image for kind='contour', None otherwise
    """

    # check existence of columns in the DataFrame
    if not (xAxis in population and yAxis in population):
        raise KeyError('population does not contain both columns for the histogram')

    # create new column with radius in R_Earth if not existent
    if (yAxis == 'r' or yAxis == 'r_rEarth'):
        if not 'r_rEarth' in population:
            population['r_rEarth'] = utils.r_Jup2r_Earth(population.r)
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
        try:
            xBins = compute_logbins(binWidth_dex[0], xRange)
            yBins = compute_logbins(binWidth_dex[1], yRange)
        except ValueError:
            warnings.warn("""ValueError, probably a division by zero in log10.
                          Please check the input data if you want to plot on
                          logarithmic axes.""")
            return

    # create 2D histogram
    h, xedges, yedges = np.histogram2d(survivedPlanets[xAxis],
                        survivedPlanets[yAxis], bins=(xBins, yBins))
    h = h.T

    if smooth:
        # smooth out the contours
        import scipy.ndimage as nd
        h = nd.gaussian_filter(h,1)

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
        h = utils.normalize_rate(h, Nsystems)
        cbarlabel = r"Planets per 100 Stars"# per $P-R_\mathrm{p}$ interval"
    else:
        cbarlabel = r"Planets per $P-R_P$ interval"

    if logColormap:
        # logarithmic color mapping. Use linear scale around zero.
        threshold  = 0.01
        if zRange is not None:
            if zRange[0] > h.min() or zRange[1] < h.max():
                warnings.warn("""Values in the histogram exceed one or both bounds
                given in 'zRange'. The visual representation of your data might
                be corrupted!""")
            colorNorm = colors.SymLogNorm(vmin=zRange[0], vmax=zRange[1],
            linthresh=max(h.min(), threshold))
        else:
            colorNorm = colors.SymLogNorm(vmin=h.min(), vmax=h.max(),
            linthresh=max(h.min(), threshold))
        cbar_kws = {'label' : cbarlabel, 'ticks' : [.01, .1, 1., 1e1, 1e2, 1e3]}
    else:
        colorNorm = zRange
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
        im = None

    elif kind == 'contour':
        """use discrete levels for occurrence. numbers are from
        Petigura et al. 2018
        """
        # levels = np.arange(-4, -1 + 1e-10, 0.25)
        cbarticklabels = [0.01, 0.03, 0.1, 0.3, 1, 3,10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        contourKwargs = dict(extend='min')
        im = plt.contourf(xedges[:-1], yedges[:-1], h, **contourKwargs,
                          **funcKwargs)
        if logColormap:
            print("logarithmic shading not possible with contours.")
    else:
        # plot a normal 2D histogram
        X, Y = np.meshgrid(xedges, yedges)
        im = ax.pcolormesh(X, Y, h, norm=colorNorm, **funcKwargs)

    # eyecandy
    if not kind == 'annotated':
        if showColorbar:
            cbar = plt.colorbar(im, pad = 0.04)
            cbar.set_label(cbarlabel, labelpad=5)
        plt.xscale('log')
        plt.yscale('log')
        if xAxis == 'period':
            ax.set_xlabel('Orbital Period [d]')
        else:
            plt.xlabel(xAxis)
        if yAxis == 'r_rEarth':
            ax.set_ylabel('Planet Size [$\mathrm{R_\oplus}$]')
        else:
            plt.ylabel(yAxis)

    return h, xedges, yedges, ax, im


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


def plot_planetTracks(simulation, truths=None, lwRange = (2., 40.),
                      avgFactor=5):
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
    avgFactor : int
        window size of the running mean that is applied to the arrays.
        Represents number of timesteps included in each point to plot.

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    ax : matplotlib axis
        axis with the plot
    """
    import types
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

        # apply running mean on data
        def runningMean(array, avgFactor):
            return np.convolve(array, np.ones((avgFactor,))/avgFactor,
                               mode='valid')
        t = runningMean(planet['t'], avgFactor)
        a = runningMean(planet['a'], avgFactor)
        m = runningMean(planet['m'], avgFactor)
        r = runningMean(planet['r'], avgFactor)

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

    cbar = plt.colorbar(lc)
    cbar.set_label('Planet Mass [$\mathrm{M_{Earth}}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time [yr]')
    ax.set_ylabel('Semimajor Axis [au]')

    return fig, ax


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


def plot_multiHistogram(dataFrame, columnNames, ax=None, labels=None, **kwargs):
    """ Plot multiple histograms from a dataFrame into the same axis.

    The function produces histograms of either several columns of the same data
    frame or of the same column for each level-0-index of a multiindex data
    frame.

    Parameters
    ----------
    dataFrame : pandas DataFrame
        data frame containing column(s) stated in 'columnNames'
    columnNames : str or list
        name(s) of the column(s) to plot.
        If this is a string, the same column is used for each level-0-index
        (assuming a multiindex data frame).
        If this is a list, histograms for each column name in it are produced.
    ax : matplotlib axis object, optional
        axis to plot on
    labels : list
        list of labels for the legend
    **kwargs : dict
        additional keyword arguments to pass to matplotlib

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    """
    if ax == None:
        fig, ax = plt.subplots()
    mplKwargs = plotstyle.histKwargs(kwargs)

    if isinstance(columnNames, list):
        # plot different columns of the same dataFrame
        for i, name in enumerate(columnNames):
            if labels == None:
                mplKwargs['label'] = name
            else:
                mplKwargs['label'] = labels[i]

            if not 'bins' in mplKwargs.keys():
                # compute ideal bin width
                mplKwargs['bins'] = astropy.stats.knuth_bin_width(dataFrame[name],
                                        return_bins=True)[1]
            ax.hist(dataFrame[name], **mplKwargs)

    elif isinstance(columnNames, str):
        if not 'bins' in mplKwargs.keys():
            # compute ideal bin width
            mplKwargs['bins'] = astropy.stats.knuth_bin_width(dataFrame[columnNames],
                                        return_bins=True)[1]
        # assume that dataFrame is multiindexed, iterate over outermost level
        for label, subpopulation in dataFrame.groupby(level=0):
           ax.hist(subpopulation[columnNames], label=label, **mplKwargs)
    ax.legend()
    return ax


def plot_multiplicities(systemMultiplicities, ax=None, legend=True, **kwargs):
    """ Plot multiple histograms from a dataFrame into the same axis.

    The function produces histograms of either several columns of the same data
    frame or of the same column for each level-0-index of a multiindex data
    frame.

    Parameters
    ----------
    systemMultiplicities : dictionary
        dictionary containing for each subpopulation of planet types:
        dimension 0: list
            number of systems with multiplicity i+1, where i is the list's index
        dimension 1: int
            number of systems in the subpopulation
        dimension 2: str
            label for the legend
    ax : matplotlib axis object, optional
        axis to plot on
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    ax : matplotlib axis
        axis with the plot
    """
    if not ax:
        fig, ax = plt.subplots()
    print("""consistency check: sum of systems with all multiplicities (should be
            100 if sufficiently high multiplicity is considered)""")
    customLabel = False
    for key, arr in systemMultiplicities.items():
        nPlanets = [i+1 for i in range(len(arr[0]))]
        # normalize to "per 100 stars that contain this planet type"
        norm_rate = 100*np.array(arr[0])/arr[1]
        if 'label' in kwargs:
            legendLabel = kwargs.get('label')
            del kwargs['label']
            customLabel = True
        elif customLabel == False:
            legendLabel = arr[2]
        if (key == 'all') & ~('color' in kwargs):
            ax.plot(nPlanets, norm_rate, label=legendLabel, color='k', **kwargs)
        else:
            ax.plot(nPlanets, norm_rate, label=legendLabel, **kwargs)
        print('{}: sum = {}'.format(key, sum(norm_rate)))
    if legend:
        plt.legend()
    ticks = plt.xticks(np.arange(min(nPlanets), max(nPlanets)+1, 1.0))
    ax.set_xlabel('Number of Planets')
    ax.set_ylabel('Frequency per 100 Systems')
    return ax


def plot_histCDF(data, axes=None, axvlines=None, cdfbins=None, **kwargs):
    """ Plot a histogram and its empirical distribution function.

    Parameters
    ----------
    data : array, list, or series object
        the data for the statistics
    axes : tuple
        contains two Matplotlib axis objects
    axvlines : list
        list of x coordinates at which to plot vertical lines across both axes
    cdfbins : list or array
        bin edges to use for the cumulative histogram
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    ax : tuple or list
        contains two Matplotlib axis objects
    """
    if axes is None:
        fig, axes = plt.subplots(2, sharex=True, figsize=[6.4, 4.8])
    if 'histtype' in kwargs.keys():
            del kwargs['histtype']

    # plot histogram in first axis
    n, bins, patches = axes[0].hist(data, histtype='stepfilled',
                                    lw=1.2, **kwargs)
    if 'bins' in kwargs:
        del kwargs['bins']

    # plot vertical lines across all axes
    if isinstance(axvlines, list):
        for x in axvlines:
            [a.axvline(x, ls='--', color='gray') for a in axes]

    if cdfbins is None:
        # use bins of normal histogram
        cdfbins = bins
    axes[1].hist(data, cumulative=True, histtype='step',bins=cdfbins, lw=2, **kwargs)

    # fix collision of tick labels and write y labels
    axes[1].get_yticklabels()[-1].set_visible(False)
    axes[0].set_ylabel('Probability Density')
    axes[1].set_ylabel('Cumulative Fraction')
    axes[1].set_ylim(0,1)
    return axes


def plot_scatterColorCoded(x, y, z, fig=None, ax=None, diverging=True,
                           cbar=True, cbarlabel='', **kwargs):
    """ Produce a scatter plot with color code and color bar.

    Parameters
    ----------
    x : iterable
        x values
    y : iterable
        y values
    z : iterable
        values for color code
    fig : matplotlib figure object, optional
        figure to plot on
    ax : matplotlib axis object, optional
        axis to plot on
    diverging : bool
        use a diverging color palette
    cbar : bool
        if True, add a colorbar
    cbarlabel : str
        label for the color bar
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    ax : matplotlib axis
        axis with the plot
    sc : matplotlib PathCollection
        output of the scatter plot
    """
    if not ax:
        fig, ax = plt.subplots()
    if diverging:
        from matplotlib.colors import Normalize
        # norm = Normalize(vmin=min(z), vmax=max(z))
        cmap = sns.diverging_palette(h_neg=250, h_pos=9., s=99, l=50,
                                 n=19, center="dark", as_cmap=True)
    else:
        norm = None

    if kwargs is None:
        kwargs = {
            'edgecolors' : 'black',
            'linewidth' : .5
        }

    sc = ax.scatter(x, y, c=z, **kwargs)
    if cbar:
        cbar = fig.colorbar(sc, label=cbarlabel)
    return fig, ax, sc


def plot_correlationMap(pop, columns, fig=None, ax=None, **kwargs):
    """ Plot a map showing mutual correlations between population columns.

    Uses the parameter-less Spearman Rank Coefficient as a measure for
    correlation. The coefficient ranges from -1 to +1, where a positive
    (negative) coefficient denotes a positive (negative) rank correlation
    between two variables and $\rho = 0 $ corresponds to no correlation.

    Parameters
    ----------
    pop : Pandas dataframe
        planet population
    columns : list
        list of column names to use for the map
    fig : matplotlib figure object, optional
        figure to plot on
    ax : matplotlib axis object, optional
        axis to plot on
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    ax : matplotlib axis
        axis with the plot
    """
    if ax == None:
        fig, ax = plt.subplots()

    # Compute the correlation matrix
    corr = pop[columns].corr(method='spearman')

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Use Spearman, because it doesn't make assumptions on distributions
    # and works for both continuous and discrete variables
    ax.set_position([0.02, 0.08, 0.97, .99])
    cbar_ax = fig.add_axes([0.8, 0.084, 0.04, 0.86])
    mplKwargs = dict(cmap='seismic_r', annot_kws={'fontsize':13}, center=0,
                     cbar_ax=cbar_ax)
    hm = sns.heatmap(round(corr,2), mask=mask, square=True, annot=True, ax=ax,
                     fmt='.2f', linewidths=3, cbar_kws={"shrink": .8,
                     'fraction':0.15, 'pad': -.05, 'aspect':18,
                     'label' : 'Spearman Rank Coefficient'},
                     **dict(mplKwargs, **kwargs))

    # plot eyecandy
    colLabels = utils.columnLabels()
    hm.tick_params(axis=u'both', which=u'both',length=0)
    labels = hm.get_xticklabels()
    labels = [colLabels[l.get_text()] for l in labels]
    labels[-1] = ''
    hm.set_xticklabels(labels)
    labels = hm.get_yticklabels()
    labels = [colLabels[l.get_text()] for l in labels]
    labels[0] = ''
    l = hm.set_yticklabels(labels)

    # ax.figure.axes[-1].yaxis.label.set_size(18)
    # ax.figure.axes[0].tick_params(labelsize=16, labelrotation=0)
    # ax.figure.axes[-1].tick_params(labelsize=16)
    # ax.figure.axes[-1].yaxis.labelpad = 20
    return fig, ax


def plot_rollingMean(df, columns, labels, onCol, winSize=100,
                     errFun=utils.sqrtOfMean, fig=None, ax=None):
    """ plot a rolling mean with confidence intervals.

    Parameters
    ----------
    df : Pandas dataframe
        data frame containing column(s) with data for the rolling mean
    columns : list
        list of column names to use
    labels : list
        list with labels for the legend. Must have the same size as 'columns'.
    onCol : str
        name of the column along which the window moves
    winSize : int, optional
        size of the sliding window
    errFun : function handle
        function to apply on 'rolling' object for computation of uncertainties.
        Example: 'errFun=np.std'
    fig : matplotlib figure object, optional
        figure to plot on
    ax : matplotlib axis object, optional
        axis to plot on
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    ax : matplotlib axis
        axis with the plot
    """
    if ax == None:
        fig, ax = plt.subplots()

    # apply rolling mean along axis specified in 'onCol'
    roll = df.rolling(winSize, center=True, on=onCol).mean()
    for col, label in zip(columns, labels):
        ax.plot(roll[onCol], roll[col], label=utils.get_label(label))

    # plot confidence intervals (clip them below zero)
    roll_err = df.rolling(winSize, center=True, on=onCol).apply(errFun)
    for col in columns:
        conf_low = np.clip(roll[col] - roll_err[col], a_min=0, a_max=None)
        ax.fill_between(roll_err[onCol], conf_low, roll[col] + roll_err[col],
        alpha=.15)

    ax.legend()
    return fig, ax


def plot_initialConditionHist(simlist, columns, fig=None, axs=None, **kwargs):
    """
    Plot histograms of initial disk conditions.

    Parameters
    ----------
    simlist : pandas DataFrame
        simulation list
    columns : list
        list containing column names in simlist to plot
    fig : matplotlib figure object, optional
        figure to plot on
    axs : list (opional)
        list containing axis objects
    **kwargs : dict
        additional keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    axs : list
        list of subaxes
    """
    if axs is None:
        fig, axs = plt.subplots(1, len(columns), figsize=[12, 3], sharex=False)

    if 'cumulative' in kwargs:
        if kwargs['cumulative'] == True:
            mplKwargs = plotstyle.histKwargs(kwargs)
    else:
        mplKwargs = kwargs

    for ax, param in zip(axs, columns):
        ax.hist(simlist[param], density=False,
                **mplKwargs)
        ax.set_xlabel(utils.columnLabels()[param])
        ax.ticklabel_format(axis='both', style='sci', scilimits=(-2,3),
                            useMathText=True)
    return fig, axs


def plot_smaMassMetallicity(pop, fig=None, ax=None):
    """Plot mass-sma-metallicity scatter including orbital range."""
    if ax == None:
        fig, ax = plt.subplots()
    fig, ax = plot_scatterColorCoded(np.array(pop.a), np.array(pop.m),
                                     np.array(pop.metallicity), fig=fig, ax=ax,
                                     cbarlabel='Host Star Metallicity [Fe/H]',
                                     vmin=-.6, vmax=.6, zorder=100)
    ax.set_xlabel('Semi-major Axis [au]')
    ax.set_ylabel('Planet Mass [$\mathrm{M_{\oplus}}$]')

    # overlay orbital radius range
    pop.r_per, pop.r_apo = utils.get_ApoPeri(pop.a.values, pop.e.values)
    per2apoLines = [pop.a - pop.r_per, pop.r_apo - pop.a]
    ax.errorbar(np.array(pop.a), np.array(pop.m), xerr=per2apoLines,
                fmt='none', c='gray', lw=1., alpha=.5)
    return fig, ax


def plot_clusterScatter(pop, clusters, x='a', y='m', fig=None, ax=None, **kwargs):
    """ Make scatter plot of planets with colors defined by cluster affiliation.

    Parameters
    ----------
    pop : Pandas dataframe
        data frame containing a planet population. Must have boolean columns for
        each of the species in 'clusters'.
    clusters : list
        list of column names to use for cluster assignment.
    x : str
        column name of x variable
    y : str
        column name of y variable
    fig : matplotlib figure object, optional
        figure to plot on
    ax : matplotlib axis object, optional
        axis to plot on
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    ax : matplotlib axis
        axis with the plot
    """

    def rename_clusters(cluster):
        try:
            transl = {'HotRockies':'Hot Rockies','ColdStragglers':'Cold Stragglers',
                'WarmNeptunes':'Warm Neptunes', 'USPs':'USPs','Giants':'Giants'}
            return transl[cluster]
        except KeyError:
            return cluster

    if ax is None:
        fig, ax = plt.subplots(figsize=plotstyle.set_size('aa', scale=1.25))
    for cluster in clusters:
        ax.plot(pop[pop[cluster]][x], pop[pop[cluster]][y], 'o', ms=1,
                label=rename_clusters(cluster), **kwargs)
    plt.loglog()
    if x == 'a':
        ax.set_xlabel('Semi-major Axis [au]')
    if y == 'm':
        ax.set_ylabel('Mass [M$_\oplus$]')
    plt.legend()
    return fig, ax


def plot_singleSystemEvo(tpop, isystem, poptdisk=None, fig=None, axs=None,
                         nTime=5, times=None, survivorsOnly=True,
                         xAxis='a', yAxis='m', **kwargs):
    """ plot the temporal evolution of a system in a-m space.

    Based on ref_red times, it samples nTime geometrically spaced times.

    Parameters
    ----------
    tpop : Pandas multiindex dataframe
        data frame containing a planet population with simulation times as
        multi-indexes.
    isystem : int
        system number to plot
    poptdisk : Pandas dataframe, optional
        population of systems at their disk dispersal time
    fig : matplotlib figure object, optional
        figure to plot on
    axs : array, optional
        list of matplotlib axis objects containing axes to plot on
    nTime : int, optional
        number of times to sample
    times : iterable, optional
        list of times to plot
    survivorsOnly : bool
        should only surviving planets (status 0) be plotted?
    xAxis : string
        parameter for the x axis
    yAxis : string
        parameter for the y axis
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    axs : iterable
        list containing matplotlib axes with the plots
    """

    def plot_rvLimit(ax, K):
        periods = np.arange(.1,1e7, 100)
        M = utils.get_M_from_RVlim(K, periods, MstarRel=1.)
        sma = utils.get_sma(periods, MstarRel=1.)
        ax.plot(sma, M, '--', c='C6', alpha=.33)

        # annotate line, rotated according to the line's slope
        ax.annotate('K = {:.0f} m/s'.format(K), xy=(.97,.73), color='C6', xytext=(.97,.73), textcoords='axes fraction',
                    rotation=10, rotation_mode='anchor', horizontalalignment='right',
                    verticalalignment='bottom')
        return ax

    if not times:
        # sample some times
        times = tpop.index.levels[0]
        timeIdxs = np.rint(np.geomspace(1, 40, nTime)).astype(np.int32)
    else:
        nTime = len(times)
        timeIdxs = [i for i in range(nTime)]
    if axs is None:
        fig, axs = plt.subplots(1, nTime, figsize=plotstyle.set_size('aaDouble', subplot=[1, nTime], scale=2),
                                sharex=True, sharey=True)

    for i, time in enumerate([times[t] for t in timeIdxs]):
        if i == nTime:
            break

        if time == 'tdisk':
            # plot system at time of disk dispersal
            if survivorsOnly:
                poptdisk = poptdisk[poptdisk.status == 0]
            sys = poptdisk[poptdisk.isystem == isystem]
            axs[i].set_title('$t_\mathrm{disk}$')
            text = axs[i].annotate('t = {:1.1f} Myr'.format(sys.iloc[0]['t']/1e6), xy=(.04, .75),
                                   ha='left', textcoords='axes fraction', xytext=(.04, .75))

        elif time == 5e9:
            # 'observation' time. mark unobservable planets
            t = tpop.loc[time]
            axs[i].set_title('$t_\mathrm{obs} =$ 5 Gyr')
            if survivorsOnly:
                t = t[t.status == 0]
            sys = t[t.isystem == isystem]

            # mask unobservable planets
            sys = utils.get_orbitalPeriod(sys, MstarRel=1.)
            sys.loc[:,'K'] = utils.get_RVsemiamplitude(sys.m, sys.period)
            minRVamp = config.minRVamplitude()['SuperEarth']
            rvMask = sys.K >= minRVamp

            axs[i].scatter(sys[rvMask][xAxis], sys[rvMask][yAxis], s=8, **kwargs)
            axs[i].scatter(sys[~rvMask][xAxis], sys[~rvMask][yAxis], s=8, c='gray')

            # overlay orbital radius range
            r_per, r_apo = utils.get_ApoPeri(sys[xAxis], sys.e)
            per2apoLines = [sys[xAxis] - r_per, r_apo - sys[xAxis]]
            axs[i].errorbar(sys[xAxis], sys[yAxis], xerr=per2apoLines,
                            fmt='none', c='gray', lw=2., alpha=.5)

            axs[i] = plot_rvLimit(axs[i], minRVamp)
            axs[i].grid(which='major')
            break

        else:
            t = tpop.loc[time]
            axs[i].set_title('{:1.1f} Myr'.format(time/1e6))
            if survivorsOnly:
                t = t[t.status == 0]
            sys = t[t.isystem == isystem]

        axs[i].scatter(sys[xAxis], sys[yAxis], s=8, **kwargs)
        #         axs[i].set_xlabel('Semi-major Axis [au]')

        if xAxis == 'a':
            # overlay orbital radius range
            r_per, r_apo = utils.get_ApoPeri(sys[xAxis], sys.e)
            per2apoLines = [sys[xAxis] - r_per, r_apo - sys[xAxis]]
            axs[i].errorbar(sys[xAxis], sys[yAxis], xerr=per2apoLines,
                            fmt='none', c='gray', lw=2., alpha=.5)

    axs[i].grid(which='major')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlim([0.01, 1000])
    axs[0].set_ylim([0.1, 50000])
    try:
        axs[0].set_ylabel(utils.columnLabels()[yAxis])
    except:
        pass
    axs[0].yaxis.set_major_locator(plt.FixedLocator([1, 100, 10000]))
    text = axs[0].annotate('System {}'.format(isystem), xy=(.04, .75),
                           ha='left', textcoords='axes fraction', xytext=(.04, .75))
    plt.subplots_adjust(.04, .25, .99, .88)
    return fig, axs


def plot_randomSystemsEvo(pop, tpop, poptdisk=None, fig=None, axs=None, nTime=5, times=None,
                          nSystems=5, seed=None, xAxis='a', yAxis='m', **kwargs):
    """ sample random systems and plot their time evolution.

    draw some systems from the population 'pop', then plot their time evolution
    using the time-dependent full population 'tpop'.

    Parameters
    ----------
    pop : Pandas dataframe
        population from which to draw the systems
    tpop : Pandas multiindex dataframe
        data frame containing a planet population with simulation times as
        multi-indexes.
    poptdisk : Pandas dataframe
        population of systems at their disk dispersal time
    fig : matplotlib figure object, optional
        figure to plot on
    axs : array, optional
        list of matplotlib axis objects containing axes to plot on
    nTime : int, optional
        number of times to sample
    times : iterable, optional
        list of times to plot
    nSystems : int
        number of systems to sample
    seed : int
        seed for the random number generator
    xAxis : string
        parameter for the x axis
    yAxis : string
        parameter for the y axis
    **kwargs : keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    axs : iterable
        list containing matplotlib axes with the plots
    """
    if times:
        nTime = len(times)
    if axs is None:
        fig, axs = plt.subplots(nSystems, nTime,
                                figsize=plotstyle.set_size('aaDouble',
                                subplot=[nSystems, nTime], scale=1.5),
                                sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    np.random.seed(seed)
    for i, isys in enumerate(np.sort(np.random.choice(pop.isystem.unique(),
                                                      nSystems, replace=False))):
        fig, axs[i] = plot_singleSystemEvo(tpop, isystem=isys, poptdisk=poptdisk, nTime=nTime,
                                           times=times, fig=fig, axs=axs[i], xAxis=xAxis,
                                           yAxis=yAxis, **kwargs)
    # [ax.set_xlabel('Semi-major Axis [au]') for ax in axs[-1]]
    try:
        [ax.set_xlabel(utils.columnLabels()[xAxis]) for ax in axs[-1]]
    except:
        pass
    [ax.get_xticklabels()[1].set_visible(False) for ax in axs[-1][1:]]
    titles = [ax.title.set_visible(False) for ax in axs.flatten()[nTime:]]
    return fig, axs


def plot_ECDFofPops(pops, labels=None, columns=None, fig=None, axs=None,
                                bins=5000, **kwargs):
    """
    Plot empirical CDFs of initial conditions for different populations.

    Parameters
    ----------
    pops : list of pandas DataFrames
        data of different populations
    labels : list of strings
        labels for the legend
    columns : list
        list containing column names in pop to plot
    fig : matplotlib figure object, optional
        figure to plot on
    axs : list (opional)
        list containing axis objects
    bins : int
        number of bins
    **kwargs : dict
        additional keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    axs : list
        list of subaxes
    """
    if fig is None:
        fig, axs = plt.subplots(2,3, sharey=True,
                                          figsize=plotstyle.set_size('aaDouble', subplot=[2,3]))
        axs = axs.reshape(-1)

    if labels is None:
        labels= [None for i in range(len(pops))]
        showLegend = False
    else:
        showLegend = True

    if columns is None:
        columns = ['metallicity', 'Msolid', 'Mgas0','aCore','tdisk','aStart']

    for pop, label in zip(pops, labels):
        for ax, param in zip(axs, columns):
            ax.hist(pop[param], density=True, cumulative=True, histtype='step',
                       label=label, bins=bins, **kwargs)
            ax.set_xlabel(utils.columnLabels()[param])
            ax.ticklabel_format(axis='both', style='sci', scilimits=(-2,3),
                                useMathText=True)
            ax = fix_hist_step_vertical_line_at_end(ax)
    if showLegend:
        axs[0].legend(loc='lower left', ncol=99, bbox_to_anchor=(0., 1.),
                              frameon=False, columnspacing=1.6)

    axs[0].set_ylim([0,1])
    return fig, axs


def plot_iceMassFractionsPopulations(populations, labels, fig=None, axs=None, rugplot=True,
                                     seed=None):
    """ Make scatter plots showing ice mass fractions for different populations.

    This produces subplots for each population and an additional rug plot for the
    starting positions of super-Earths in the system, color-coded by final ice
    mass fraction. The number of planets is balanced to the size of the smallest
    population.
    """
    import matplotlib.cm as cmx

    N_rows = len(populations)
    N_samples = min([len(p[(p.status == 0) & (p.m > 0.5)]) for p in populations])
    # N_samples = min([len(p[(p.status==0) & (p.m > 1.)]) for p in populations])

    print('size of smallest population = N_samples = {}'.format(N_samples))

    fig = plt.figure(figsize=plotstyle.set_size(subplot=[round(1.5 * N_rows), 2]))
    gs = fig.add_gridspec(N_rows, 2, left=0.05, right=0.98, wspace=0.1,
                          hspace=0., width_ratios=[1, 0.07])
    axs = [fig.add_subplot(gs[i, 0]) for i in range(N_rows)]
    cax = fig.add_subplot(gs[:, -1])

    # normalize color range across populations
    vmin = min([min(p.fracIce) for p in populations])
    vmax = max([max(p.fracIce) for p in populations])

    for i, label, p in zip(range(N_rows), labels, populations):
        p = p[(p.status == 0) & (p.m > 0.5)].sample(N_samples,
                                                    random_state=seed, replace=False)

        fig, axs[i], sc = plot_scatterColorCoded(p.a, p.m, p.fracIce, fig=fig,
                                                       ax=axs[i], cmap='viridis_r',
                                                       diverging=False, cbar=False,
                                                       vmin=vmin, vmax=vmax, alpha=.75, s=3)

        if rugplot:
            '''add rugplot of initial starting positions. Consider only SE for rugplot'''
            p = p[p.planetType == 'SuperEarth']
            if len(p) > 0:
                cNorm = colors.Normalize(vmin=vmin, vmax=vmax)
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='viridis_r')
                axs[i] = sns.rugplot(p.aStart, ax=axs[i], c=scalarMap.to_rgba(p.fracIce),
                                     lw=.5, alpha=.99, height=0.075)
                text = axs[i].annotate('SE initial orbit', xy=(.65, .05), fontsize=8,
                                       xycoords='axes fraction',
                                       ha='right', va='center', textcoords='axes fraction',
                                       xytext=(.99, .05))
                text = axs[i].annotate('', xy=(.65, .05), xycoords='axes fraction',
                                       ha='right', va='center', textcoords='axes fraction',
                                       xytext=(.72, .05), arrowprops=dict(arrowstyle="->",
                                       connectionstyle="arc3"))

        axs[i].loglog()
        axs[i].set_xlim([1e-2, 1e3])
        axs[i].set_ylim([.15, 1e4])
        axs[i].set_ylabel('$\mathrm{M_P} [\mathrm{M}_{\oplus}$]')
        text = axs[i].annotate(label, xy=(.04, .83),
                               ha='left', textcoords='axes fraction', xytext=(.04, .83))
        # axs[i].grid(which='major')

    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label('Core Ice Mass Fraction')
    axs[-1].set_xlabel('Semi-major Axis [au]')

    [ax.xaxis.set_ticklabels([]) for ax in axs[:-1]]
    [ax.get_yticklabels()[-2].set_visible(False) for ax in axs[1:]]

    return fig, axs


def plot_massRadiusPopulations(populations, labels, fig=None, axs=None,
                               seed=None, **kwargs):
    """
    Make a mass-radius plot for different populations.

    Parameters
    ----------
    populations : list of pandas DataFrames
        data of different populations
    labels : list of strings
        labels for the legend
    fig : matplotlib figure object, optional
        figure to plot on
    axs : list (opional)
        list containing axis objects
    **kwargs : dict
        additional keyword arguments to pass to matplotlib

    Returns
    -------
    fig : matplotlib figure object
        figure with the plot
    axs : list
        list of subaxes
    sc : PathCollection
        as returned from plt.scatter
    """
    colors = ['C2', 'C4']

    def filterPop(p, massLim=(1., 47.), smaLim=(0., 0.3)):
        filtered = p[(p.status == 0) & (p.m.between(*massLim)) & (p.a.between(*smaLim))]
        return filtered

    N_samples = min([len(filterPop(p)) for p in populations])
    print('balanced sample: N={}'.format(N_samples))

    if axs is None:
        fig, axs = plt.subplots(3, 1, figsize=plotstyle.set_size(subplot=[3, 1], scale=1.),
                                          sharex=True)

    # normalize color range across populations
    vmin = min([min(p.fracIce) for p in populations])
    vmax = max([max(p.fracIce) for p in populations])

    # first two axes: the two populations separately, color-coded by ice mass fraction
    for ax, label, p, in zip(axs[:2], labels, populations):
        p = filterPop(p)
        sc = ax.scatter(p.m, p.r_rEarth, c=p.fracIce, cmap='viridis_r', alpha=.9,
                   s=3, label=label, vmin=vmin, vmax=vmax)
        text = ax.annotate(label, xy=(.04, .87),
                           ha='left', xycoords='axes fraction')

    try:
        # third axis: balanced samples of both populations in one scatter plot, colored by population.
        # reverse orders to have yellow on top.
        for label, p, col in zip(reversed(labels), reversed(populations), reversed(colors) ):
            p = filterPop(p)
            p = p.sample(N_samples, random_state=seed)
            axs[2].scatter(p.m, p.r_rEarth, edgecolors=col, alpha=.66,
                           s=9, label=label, facecolors='none')
            # reverse order of legend
            handles, labels = [hl[:-1] for hl in axs[2].get_legend_handles_labels()]
            axs[2].legend(handles[::-1], labels[::-1])
    except IndexError:
        # no third axis provided
        pass

    axs[-1].set_xlabel('Planet Mass [M$_\oplus$]')
    [ax.set_ylabel('Planet Radius [R$_\oplus$]') for ax in axs]
    [ax.loglog() for ax in axs]
    return fig, axs, sc


def plot_InitialsPairplot(data, variables, sample=True, samplesize=np.inf,
                          namedict=None):
    """ make a corner plot of disk features, color-coded by planet cluster affiliation.

    Parameters
    ----------
    data : pandas DataFrame
        DataFrame containing the data. Has to contain a numerical column 'labels'
        or 'labelname' for the color-coding.
    variables : list
        list of strings with column names to consider for the plot
    sample : bool
        sample from data?
    samplesize : int
        size of the sample

    Returns
    -------
    pp : seaborn PairGrid
        object with the pairplot
    """
    if sample:
        sample = data.sample(min(samplesize,len(data)), random_state=42)
    else:
        sample=data
    sample = log_col(sample, variables)

    pp = sns.pairplot(sample, vars=[v + '_log' for v in variables], hue='labels',
                      palette=sns.color_palette(),
                      plot_kws = {'alpha':.5,
                      's':2,
                      'linewidth':0}) # remove white marker edges)

    # remove upper triangle
    for i, j in zip(*np.triu_indices_from(pp.axes, 1)):
        pp.axes[i, j].set_visible(False)

    pp.fig.set_size_inches(8,8)

    # fix legend
    if namedict is None:
        namedict = {1: 'Neptunes', 2: 'icy cores', 3: 'giant planets', 4: '(super-)Earths'}
    handles = pp._legend_data.values()
    labels = [str(i) + ': ' + namedict[int(key)] for i, key in enumerate(pp._legend_data.keys())]
    pp._legend.remove()
    lgnd = pp.fig.legend(handles=handles, labels=labels, loc='upper center', ncol=1, title='cluster')

    #change the marker size manually for legend
    for i in range(len(lgnd.legendHandles)):
        lgnd.legendHandles[i]._sizes = [30]

    # fix labels
    for ax in pp.axes[-1,:]:
        xlabel = ax.xaxis.get_label_text()
        ax.xaxis.set_label_text(utils.get_plotlabels(xlabel))
    for ax in pp.axes[:,0]:
        ylabel = ax.yaxis.get_label_text()
        ax.yaxis.set_label_text(utils.get_plotlabels(ylabel))
    return pp