""" Config file with my personal preference of default matplotlibrc settings.
"""
from matplotlib import rcParams
from cycler import cycler
import seaborn as sns


def styleplots():
    # Set plot style
    sns.set(context='notebook', style='ticks', font_scale=1.2, palette='Set2',\
        rc={
    # 'text.usetex':True,
    # 'text.latex.unicode':True,
    'font.family' : 'sans-serif',
    #'font.serif':'Computer Modern',
    'font.style'         : 'normal',
    'font.variant'        : 'normal',
    'font.weight'         : 'normal',
    'font.stretch'        : 'normal',
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
    'legend.handletextpad' : 0.5,
    'legend.handlelength' : 0.75,
    'xtick.minor.size'     : 2.,
    'ytick.minor.size'     : 2.,
    'image.cmap'   : 'inferno',
    # histograms
    'hist.bins' : 20,
    'patch.edgecolor' : 'black',
    })

    sns.set_color_codes()

    # 'fivethirtyeight' colors
    rcParams['axes.prop_cycle'] = cycler(color=['#008fd5', '#fc4f30', '#e5ae38',
        '#810f7c', '#029e73', '#8b8b8b', '#00035b', '#fe828c','#005249'])

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

styleplots()
