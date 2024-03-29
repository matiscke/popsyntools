""" Config file with my personal preference of default matplotlibrc settings.
"""
import matplotlib as mpl
from matplotlib import rcParams
from cycler import cycler


def get_colorPalette():
    """ return a custom made color palette."""
    return ['#008fd5', '#fc4f30', '#e5ae38',
            '#810f7c', '#029e73', '#00035b',
            '#8b8b8b', '#fe828c','#005249']


def styleplots():
    # Set plot style
    spinewidth = 1.5

    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 9}                 # 9 pt corresponds to caption size of 1-col A&A
    legend = {'handlelength' : 0.75,
              'handletextpad' : 0.5}
    figure = {'dpi' : 100,
              'subplot.left'    : 0.16,   # the left side of the subplots of the figure
              'subplot.bottom'  : 0.21,   # the bottom of the subplots of the figure
              'subplot.right'   : 0.98,   # the right side of the subplots of the figure
              'subplot.top'     : 0.97,   # the top of the subplots of the figure
              'subplot.hspace'  : 0.0}    # height reserved for space between subplots

    mpl.rc('font', **font)
    mpl.rc('legend', **legend)
    mpl.rc('figure', **figure)
    mpl.rc('lines', linewidth = 2.5)

    mpl.rc('axes', linewidth = spinewidth)
    mpl.rc('xtick.major', width = spinewidth)
    mpl.rc('ytick.major', width = spinewidth)
    mpl.rc('xtick.minor', width = .5*spinewidth)
    mpl.rc('ytick.minor', width = .5*spinewidth)
    mpl.rc('image', cmap = 'inferno')
    mpl.rc('hist', bins = 20)
    mpl.rc('patch', edgecolor = 'black')
    mpl.rc('savefig', bbox = 'tight', dpi = 200)

    rcParams['axes.prop_cycle'] = cycler(color=get_colorPalette())


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


def set_size(width='aa', subplot=[1,1], scale=1):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float or string
        either width in pts or one of the following strings:
        'aa' : A&A column width
        'aaDouble' : A&A total text width
    subplot : list
        subplot dimensions in [rows, columns]
    scale: float
        fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches

    Credits
    -------
    Adapted from https://jwalton.info/Embed-Publication-Matplotlib-Latex/
    """
    aa_textwidth = 523.53       # textwidth of A&A template, as measured with "\the\textwidth" Latex command
    thesis_textwidth = 468.90   # measured in my PhD thesis
    if width == 'aaDouble':
        width_pt = aa_textwidth
    elif width == 'thesis':
        width_pt = thesis_textwidth
    else:
        width_pt = aa_textwidth/2

    # Width of figure
    fig_width_pt = width_pt * scale

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt

    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplot[0] / subplot[1])
    fig_dim = [fig_width_in, fig_height_in]

    return fig_dim


styleplots()
