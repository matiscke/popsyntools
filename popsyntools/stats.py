""" Contains routines for statistical analyses of planet populations.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd

from popsyntools import utils, config


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
    population['planetType'] = np.nan

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
            """Besides mass limits, consider period and RV semi-amplitude.
            We have to make sure that the ejected planets are not removed (they
            will be flagged in another step), thus we don't apply the period
            and RV limits for ejected planets.
            """
            while True:
                try:
                    # mask_per_RV = ~(population['status'] == 0 & ((population['period'] > periodLim[pType][1]) & (population['period'] <= periodLim[pType][0]) | (population['K'] < minRVamp[pType])))

                    mask_per_RV = (population['status'] == 2) | (
                                   (population['period'] > periodLim[pType][0])
                                 & (population['period'] <= periodLim[pType][1])
                                 & (population['K'] > minRVamp[pType]))
                    mask = mask & mask_per_RV

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


def multiplicityFreq(subPop, nMax):
    """ Obtain mean multiplicity and standard deviation of a (sub-)population
    as well as the frequency of each multiplicity. Important: 'subPop' has to be
    a subpopulation containing only one planet type (except you want to count
    all planets)

    Parameters
    ----------
    subPop : pandas DataFrame
        sub-population containing only the planet type to count
    nMax : int
        maximum multiplicity considered

    Returns
    -------
    mean : float
        mean multiplicity
    std : float
        standard deviation of multiplicity
    NsystemsPerMult : list
        number of systems for each multiplicity
    Nsystems : integer
        Number of unique systems in the sub-population
    """
    counts = subPop.groupby('isystem').planetType.count()
    Nsystems = subPop.isystem.nunique()
    # obtain number of systems with each multiplicity
    NsystemsPerMult = []
    for nPlanet in range(1, nMax+1):
        NsystemsPerMult.append(np.sum(subPop.groupby('isystem').planetType.count() == nPlanet))
    # print('median for {}: {}'.format(subPop.planetType.unique()[0], np.median(counts)))
    return np.mean(counts), np.std(counts), NsystemsPerMult, Nsystems, counts


def get_multiplicities(pop, pTypes=None, nMax=7, verbose=True):
    """Compute multiplicities for all distinct planet types in a population.

    Parameters
    ----------
    pop : pandas DataFrame
        planet population
    pTypes : list of strings
        the planet types in column "planetType" that should be considered. If
        this is 'all', survived planets of .5 Earth masses and higher that pass
        the detection limit are considered.
    nMax : int
        maximum multiplicity considered

    Returns
    -------
    systemMultiplicities : dict
        contains for each planet type:
            [0]: list
                No. of systems per multiplicity
            [1]: int
                No. of systems
            [2]: str
                label
            [3]: float
                mean multiplicity
            [4]: float
                standard deviation of multiplicity
    """
    systemMultiplicities = {}
    if pTypes == None:
        pTypes = pop.planetType.unique()
    for pType in pTypes:
        if (pd.isnull(pType)) | (pType == 'allTypes'):
            # select all systems including any considered planet type, but
            # exclude ejected planets
            subPop = pop[(pd.notnull(pop.planetType)) &
                         (~pop['planetType'].str.contains('ejected', na=False))]
        elif pType.lower() == 'all':
            # consider all survived planets more massive than 0.5 M_Earth and
            # with K > Kmin(SE)
            popAll = pop.copy()
            popAll['planetType'] = None     # should this be 'NaN' instead?
            minRVamp = config.minRVamplitude()
            massLim = config.massLimits()
            masses = popAll['msini']

            mask_status = popAll['status'] == 0
            mask = (mask_status & (masses > massLim['Earth'][0])
                                & (masses <= np.inf))
            mask = (mask & (popAll['K'] > minRVamp['SuperEarth']))
            popAll.loc[mask, 'planetType'] = 'all'
            subPop = popAll[popAll.planetType == 'all']
        else:
            subPop = pop[pop.planetType == pType]
        meanMul, std, NsystemsPerMult, Nsystems, counts = multiplicityFreq(subPop,
            nMax)
        if verbose:
            print("{}: mean multiplicity = {:1.2f}; std = {:1.2f}".format(pType,
                meanMul, std))
        systemMultiplicities[pType] = [NsystemsPerMult, Nsystems,
            utils.get_label(pType), meanMul, std]
    return systemMultiplicities


def add_multiplicity(pop):
    """ add a column with planet multiplicity per system.

    Considers 'all' planets based on criteria in stats.get_multiplicities().
    """
    for i in pop.isystem.unique():
        m = get_multiplicities(pop[pop.isystem == i], ['all'],
                               verbose=False)['all'][3]
        pop.loc[pop.isystem == i, 'multiplicity'] = m if not np.isnan(m) else 0
    return pop


def add_pTypeFrequency(pop):
    """ add a column with frequency per system of specific planet types.

    Realized planet types:
    'SuperEarth' => 'nSE'
    'ColdJupiter' => 'nCJ'
    """
    planetTypes = ['SuperEarth', 'ColdJupiter']
    for col, pt in zip(['nSE', 'nCJ'], planetTypes):
        for val in pop.isystem.unique():
            pop.loc[pop.isystem == val, col] = len(pop[(pop.planetType == pt)
                                                       & (pop.isystem == val)])
    return pop

def occurrenceTable(pop, onCol='Msolid'):
    """ create a table with occurrence rates depending on an initial parameter.

    One line per system. This version is specialized on Super Earths and Cold
    Jupiters.
    """
    isystem = []
    diskParam = []
    nSE = []
    nCJ = []
    for i, group in pop.groupby(onCol):
        isystem.append(int(group.isystem.unique()[0]))
        diskParam.append(group.median()[onCol])
        nSE.append(len(group[group.planetType == 'SuperEarth']))
        nCJ.append(len(group[group.planetType == 'ColdJupiter']))

    occ = pd.DataFrame({onCol : diskParam, 'nSE' : nSE, 'nCJ' : nCJ},
                       index=isystem)
    return occ


def finalFate(system, iplanet, iplanetOrig=None):
    """
    track the final fate of a planet that was accreted.

    goes iteratively through the accretion cascade until the planet that was
    *not* accreted anymore is found. The index of this final planet is returned.
    This function calls itself recursively.

    Parameters
    ----------
    system : pandas DataFrame
        frame containing only the planets of a single system
    iplanet : int
        index of the current planet
    iplanetOrig : int
        index of the planet in question, remains the same during recursive
        function calls.

    Returns
    -------
    system : pandas DataFrame
        frame containing only the planets of a single system
    """
    if iplanetOrig == None:
        # first call from outside, store the original planet index
        iplanetOrig = iplanet
    planet = system[system.iplanet == iplanet]
    if not (planet.status < 0).any():
        system.loc[system.iplanet == iplanetOrig, 'finalFate'] = iplanet
        return system
    return finalFate(system, -int(planet.status), iplanetOrig)


def get_finalFate(pop):
    """
    track down final planet fate of accreted planets for the whole population.

    Uses function 'finalFate()' to add the final fate of each planet (See
    docstring of that function).

    Parameters
    ----------
    pop : pandas DataFrame
        planet population

    Returns
    -------
    pop : pandas DataFrame
        planet population
    """

    # disable performance-killing calls of gc.collect()
    pd.set_option('mode.chained_assignment', None)

    for i, system in pop.groupby('isystem'):
        for iplanet in system.iplanet:
            system = finalFate(system, iplanet)
        pop.loc[pop.isystem == i, 'finalFate'] = system['finalFate']
    return pop
