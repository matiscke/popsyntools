""" Contains routines for statistical analyses of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd

def categorize(population, Mgiant=300.):
    """ Sort planets into different categories.

    Parameters
    ----------
    population : pandas DataFrame
        planet population
    Mgiant : float
        minimum mass for a planet to be considered a giant
    """
    Nplanets = len(population[population['status'] == 0])
    Nplanets_ejected = len(population[population['status'] == 2])
    NltEarth = len(population[(population['status'] == 0) & (population['m'] > 1.)])
    NltEarth_ejected = len(population[(population['status'] == 2) & (population['m'] > 1.)])
    Ngiants = len(population[(population['status'] == 0) & (population['m'] > Mgiant)])
    Ngiants_ejected = len(population[(population['status'] == 2) & (population['m'] > Mgiant)])

    print('giant mass Mgiant = {}'.format(Mgiant))
    print('number of planets: {}'.format(Nplanets))
    print('number of ejected planets: {}'.format(Nplanets_ejected))
    print('Number of planets more massive than M_Earth: {}'.format(NltEarth))
    print('Number of planets more massive than M_Earth and ejected: {}'.format(NltEarth_ejected))
    print('Number of planets more massive than M_giant: {}'.format(Ngiants))
    print('Number of planets more massive than M_giant and ejected: {}'.format(Ngiants_ejected))
