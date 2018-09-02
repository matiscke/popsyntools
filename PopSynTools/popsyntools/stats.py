""" Contains routines for statistical analyses of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd

import utils, config


def print_categories(population, Mgiant=300.):
    """ Sort planets into different categories and print simple statistics.

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


def categorizePlanets(population, ZhuWu18=False):
    """ Label planets into different mass categories.

    Each planet is categorized according to its mass with limits specified in
    config.py. The planet category (e.g. "Earth"; "Earth_ejected"; ...) is
    written into a column "planetType", which is newly created if nonexistent.

    Parameters
    ----------
    population : pandas DataFrame
        planet population
    ZhuWu18 : Bool
        Flag to consider limits for comparison with Zhu & Wu 2018

    Returns
    -------
    population : pandas DataFrame
        categorized population
    """

    massLim = config.massLimits(ZhuWu18)

    if ZhuWu18:
        # consider other criteria than mass, too
        periodLim = config.periodLimits()
        minRVamp = config.minRVamplitude()

    # keep only survived and ejected planets
    mask_status = (population['status'] == 0) | (population['status'] == 2)

    for pType in massLim:
        if ZhuWu18:
            # consider m*sin(i) (with random isotropic i) instead of m
            while True:
                try:
                    masses = population['msini']
                except KeyError:
                    # apparently, there is no 'msini' column yet
                    print("computing m*sin(i) with random, isotropic \
                    inclination i. This will take a while.")
                    population = utils.add_msini(population)
                    masses = population['msini']
                    continue
                break
        else:
            masses = population['m']

        # assign planet type according to mass limits
        mask = (mask_status & (masses > massLim[pType][0])
                            & (masses <= massLim[pType][1]))

        if ZhuWu18:
            # consider other criteria than mass, too
            while True:
                try:
                    mask = (mask & (population['period'] > periodLim[pType][0])
                                 & (population['period'] <= periodLim[pType][1])
                                 & (population['K'] > minRVamp[pType]))
                except KeyError:
                    # has RV semi-amplitude not been calculated yet?
                    print("adding RV semi-amplitude using Mstar = 1.0 Msol.")
                    population['K'] = utils.get_RVsemiamplitude(population['m'],
                            population['period'], population['e'], Mstar=1.0)
                    continue
                break

        population.loc[mask, 'planetType'] = pType

    # label ejected planets
    population.loc[population['status'] == 2, 'planetType'] += '_ejected'

    return population


def filterPlanets(population, pType):
    """ return a population with planets of a certain type.

    Planet types and their mass limits are specified in config.py.

    Parameters
    ----------
    population : pandas DataFrame
        planet population
    pType : string
        type of planet

    Returns
    -------
    population_filtered : pandas DataFrame
        filtered population

    """
    import warnings

    massLim = config.massLimits()
    if not pType in massLim:
        warnings.warn("the given planet type '{}' is not known. Please choose "
        "one of the following types or specify a new one in 'config.py': {}".format(
            pType, [s for s in massLim.keys()]))
        return population

    population = categorizePlanets(population)
    population_filtered = population[population['planetType'] == pType]
    return population_filtered


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
    Nfiltered = len(population_filtered)
    if Nfiltered > 0:
        stats['multiplicity'] = stats['Nplanets']/Nfiltered
    else:
        stats['multiplicity'] = 0

    # metallicity of stars with min. 1 planet of this type: mean and std
    population_filtered = utils.convert_dgr2metallicity(population_filtered)
    stats['meanMetallicity'] = population_filtered.metallicity.mean()
    stats['stdMetallicity'] = population_filtered.metallicity.std()

    # eccentricity of planets of this type: mean and std
    stats['meanEccentricity'] = population_filtered.e.mean()
    stats['stdEccentricity'] = population_filtered.e.std()

    return stats
