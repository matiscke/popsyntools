""" A configuration file with variables used for analysis of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""

import numpy as np

def massLimits(ZhuWu18=False):
    """ Provide mass limits for planet categories.

    Planet masses are given in Earth masses.

    Parameters
    ----------
    ZhuWu18 : Bool
        Flag to consider limits for comparison with Zhu & Wu 2018

    Returns
    -------
    massLimits : dict
        dictionary containing mass limits.
    """

    if ZhuWu18:
        massLimits = {
        'SuperEarth'  : (2., 20.),
        'ColdJupiter' : (95., np.inf),
        'HotJupiter'  : (95., np.inf),
        'WarmJupiter' : (95., np.inf)
        }
    else:
        massLimits = {
        'all' : (0., np.inf),
        'ltEarth' : (1., np.inf),
        'Earth' : (.5, 2.),
        'SuperEarth' : (2., 10.),
        'Neptune' : (10., 30.),
        'SubGiant' : (30., 300.),
        'Giant' : (300., 4322.),
        'BrownDwarf' : (4322., np.inf)
        }

    return massLimits

def periodLimits():
    """ Provide period limits (in days) for planet categories."""
    periodLim = {
    'SuperEarth'  : (0., 400.),
    'ColdJupiter' : (400., np.inf),
    'HotJupiter'  : (0., 10.),
    'WarmJupiter' : (10., 400.)
    }
    return periodLim

def minRVamplitude():
    """ Provide minimum RV amplitude 'K' (in m/s) for planet categories."""
    minRVamp = {
    'SuperEarth'  : 2.,
    'ColdJupiter' : 2.,
    'HotJupiter'  : 2.,
    'WarmJupiter' : 2.
    }
    return minRVamp
