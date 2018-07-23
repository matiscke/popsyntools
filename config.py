""" A configuration file with variables used for analysis of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""

import numpy as np

def massLimits():
    """ provide mass limits for planet categories"""
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
