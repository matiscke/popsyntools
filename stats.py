""" Contains routines for statistical analyses of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd

import utils


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

    Notes
    -----
    Supported planet types:
    all
    ltEarth
    Earth
    SuperEarth
    """
    if not type == 'giants_ejected':
        # first, keep only survived planets
        population = population[population['status'] == 0]

    if type == 'all':
        population_filtered = population[population['m'] > 0.]
    elif type == 'ltEarth':
        population_filtered = population[population['m'] > 1.]
    elif type == 'Earth':
        population_filtered = population[(population['m'] > .5)
            & (population['m'] <= 2.)]
    elif type == 'SuperEarth':
        population_filtered = population[(population['m'] > 2.)
            & (population['m'] <= 10.)]

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


def get_typeStats(population, population_filtered):
    """ Compute statistics concerning a certain planet type.

    Parameters
    ----------
    population : pandas DataFrame
        full planet population
    population_filtered : pandas DataFrame
        population of a certain planet type; could be a population that was
        filtered with the 'filterPlanets' function.

    Returns
    -------
    stats : dictionary
        the statistics for the planet type in question
    """
    stats = {}

    # Number of planets of this type
    stats['Nplanets'] = len(population_filtered)

    # Number of systems with min. 1 planet of this type
    stats['Nsystems'] = population_filtered.isystem.nunique()

    # fraction of systems with min. 1 planet of this type
    stats['fractionSystems'] = stats['Nsystems']/population.isystem.nunique()

    # occurrence rate per star: mean number of planets of this type per system
    stats['occurrence'] = stats['Nplanets']/len(population)

    # multiplicity: mean number of planets of this type per system that contains
    # this type
    stats['multiplicity'] = stats['Nplanets']/len(population_filtered)

    # metallicity of stars with min. 1 planet of this type: mean and std
    population_filtered = utils.convert_dgr2metallicity(population_filtered)
    stats['meanMetallicity'] = population_filtered.metallicity.mean()
    stats['stdMetallicity'] = population_filtered.metallicity.std()

    # eccentricity of planets of this type: mean and std
    stats['meanEccentricity'] = population_filtered.e.mean()
    stats['stdEccentricity'] = population_filtered.e.std()

    return stats
