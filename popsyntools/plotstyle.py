""" Config file with my personal preference of default matplotlibrc settings.
"""

import seaborn as sns

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
'figure.subplot.left'    : 0.14,    # the left side of the subplots of the figure
'figure.subplot.bottom'  : 0.14,   # the bottom of the subplots of the figure
'figure.subplot.right'   : 0.96,   # the right side of the subplots of the figure
'figure.subplot.top'     : 0.95,    # the top of the subplots of the figure
'figure.subplot.hspace'  : 0.0,    # height reserved for space between subplots
'axes.xmargin' : 0.02,             # default margin for autoscale
'axes.ymargin' : 0.02,
'legend.handletextpad' : 0.5,
'legend.handlelength' : 0.75,
'xtick.minor.size'     : 2.,
'ytick.minor.size'     : 2.,

# histograms
'hist.bins' : 20,
'patch.edgecolor' : 'black',
})

sns.set_color_codes()

# 'fivethirtyeight' colors
from matplotlib import rcParams
from cycler import cycler
rcParams['axes.prop_cycle'] = cycler(color=['#008fd5', '#fc4f30', '#e5ae38',
    '#810f7c', '#029e73', '#8b8b8b', '#00035b', '#fe828c','#005249'])
