""" Config file with my personal preference of default matplotlibrc settings.
"""
from matplotlib import rcParams
from cycler import cycler
import seaborn as sns


def styleplots():
    # Set plot style
    sns.set(context='notebook', style='ticks', font_scale=1.,\
        rc={
    'text.usetex':False,
    # 'text.latex.unicode':True,
    'font.family' : 'sans-serif',
    #'font.serif':'Computer Modern',
    'font.style'         : 'normal',
    'font.variant'        : 'normal',
    'font.weight'         : 'normal',
    'font.stretch'        : 'normal',
    # Use 11pt font in plots, to match 11pt font in A&A documents
    "font.size": 11,
    "axes.labelsize": 11,
    'legend.fontsize': 8,
    'legend.handletextpad' : 0.5,
    'legend.handlelength' : 0.75,
    'xtick.labelsize' :8,
    'ytick.labelsize' :8,
    'xtick.minor.size'     : 2.,
    'ytick.minor.size'     : 2.,
    'savefig.dpi'         : 400,
    'lines.linewidth'     : 2.0,
    'figure.figsize'   : [6.4, 4.8],         # figure size in inches
    'figure.facecolor'      : 'white',
    'figure.subplot.left'    : 0.17,    # the left side of the subplots of the figure
    'figure.subplot.bottom'  : 0.18,   # the bottom of the subplots of the figure
    'figure.subplot.right'   : 0.96,   # the right side of the subplots of the figure
    'figure.subplot.top'     : 0.93,    # the top of the subplots of the figure
    'figure.subplot.hspace'  : 0.0,    # height reserved for space between subplots
    'axes.xmargin' : 0.02,             # default margin for autoscale
    'axes.ymargin' : 0.02,
    'image.cmap'   : 'inferno',
    # histograms
    'hist.bins' : 20,
    'patch.edgecolor' : 'black',
    })

    sns.set_color_codes()


def get_colorPalette():
    """ return a custom made color palette."""
    return ['#008fd5', '#fc4f30', '#e5ae38',
        '#810f7c', '#029e73', '#8b8b8b', '#00035b', '#fe828c','#005249']


def histKwargs(overrideDict=None):
    """get custom keyword arguments for matplotlib histograms.

    Values can be overwritten by specifying a dictionary 'overrideDict'.
    """
    histKwargs = {'alpha' : 0.6,
                     'histtype' : 'stepfilled',
                     'edgecolor' : 'black',
                     'linewidth' : 1.2,
                  }
    if overrideDict is not None:
        for key, val in zip(overrideDict.keys(), overrideDict.values()):
            histKwargs[key] = val
    return histKwargs


def set_size(width, fraction=1, subplot=[1,1]):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
        Width in pts
    fraction: float
        Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches

    Credits
    -------
    Adapted from https://jwalton.info/Embed-Publication-Matplotlib-Latex/
    """

    if width == 'aa':
        width_pt = 256.07
    elif width == 'aaDouble':
        width_pt = 523.53
    else:
        width_pt = width

    # Width of figure
    fig_width_pt = width_pt * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt

    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplot[0] / subplot[1])
    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


styleplots()
rcParams['axes.prop_cycle'] = cycler(color=get_colorPalette())
