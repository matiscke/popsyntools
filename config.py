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
    'SuperEarth' : (2., 10.)
    }

    return massLimits
