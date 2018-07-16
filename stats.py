""" Contains routines for statistical analyses of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd


def filterPlanets(population, type):
    """ return a population with planets of a certain type.

    Parameters
    ----------
    population : pandas DataFrame
        planet population
    type : string
        type of planet

    Returns
    -------
    population_filtered : pandas DataFrame
        filtered population
    """
    if not type == 'giants_ejected':
        # first, keep only survived planets
        population = population[population['status'] == 0]

    if type == 'all':
        population_filtered = population[population['m'] > 0.]
    elif type == 'ltEarth':
        population_filtered = population[population['m'] > 1.]

    return population_filtered


def categorize(population, Mgiant=300.):
    """ Sort planets into different categories.

    Parameters
    ----------
    population : pandas DataFrame
        planet population
    Mgiant : float
        minimum mass for a planet to be considered a giant

    Returns
    -------
    categories : dictionary
        amount of planets in each category
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

    return {'Nplanets' : Nplanets, 'Nplanets_ejected' : Nplanets_ejected,
            'NltEarth' : NltEarth, 'NltEarth_ejected' : NltEarth_ejected,
            'Ngiants' : Ngiants, 'Ngiants_ejected' : Ngiants_ejected}
