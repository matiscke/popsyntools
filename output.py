""" Contains functions to analyze output of the planet formation Code 'Planete'
by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd
import tables

import stats


def rename_tracksColumns(planetTracks, ref_red=False):
    """Rename some of the columns of a planet tracks or ref_red table.

    Newer versions of ref_red files (from ~May 2018) contain additional
    columns that are renamed if the number of columns exceeds the size of the
    previous ref_reds.

    Parameters
    ----------
    planetTracks : Pandas DataFrame
        Table with the tracks of a single planet
    ref_red : bool
        rename additional columns for 'ref_red' files

    Returns
    -------
    dfOut : Pandas DataFrame
        DataFrame with changed column names
    """
    colnames = {
     0: 'n',
     1: 't',
     2: 'mCore',
     3: 'menv',
     4: 'm',
     5: 'L',
     6: 'Lcomp',
     7: 'mdotcore',
     8: 'mdotgasl',
     9: 'rcore',
     11: 'pcore',
     12: 'tcore',
     13: 'rhocen',
     14: 'r',
     18: 'a',
     20: 'mdiskg',
     21: 'mdiskp',
     22: 'rroche',
     23: 'racc',
     24: 'dt',
     25: 'sigmamean',
     26: 'rfeed',
     27: 'rcaptot',
     29: 'tneb',
     30: 'pneb',
     31: 'type_mig',
     32: 'mejetot',
     33: 'macctot',
     34: 'miso',
     35: 'sigmagas',
     39: 'dtpl',
     41: 'rhocore',
     42: 'mgazacc',
     43: 'mgazevap',
     44: 'pout',
     45: 'tout',
     51: 'lcont',
     52: 'enew',
     53: 'ediff',
     54: 'kenergdiff',
     55: 'lacccore',
     56: 'ep',
     57: 'ip',
     62: 'e',
     63: 'i',
     66: 'typemig',
     67: 'status',
     69: 'tmig',
     70: 'tmige',
     73: 'mdotgas',
     78: 'dtmode',
     100: 'lactual',
     101: 'corrlesti',
     102: 'correesti',
     104: 'mdotgasmax',
     107: 'lcontenv',
     108: 'lcontsum',
     109: 'dust2gas',
     116: 'systemNo',
     117: 'planetNo'}
     # 118: 'NoEmbryos'}
    dfOut = planetTracks.rename(columns=colnames)

    if ref_red:
        # check if it's a recent ref_red file and add additional columns
        if len(dfOut.columns) > 120:
            ref_redColumns = {
            118: 'isystem',
            119: 'iplanet',
            120: 'isystemorig',
            121: 'iplanetorig',
            122: 'nplanets',
            123: 'line'
            }
            dfOut = dfOut.rename(columns=ref_redColumns)
    return dfOut


def read_simFromFolder(foldername):
    """Read a simulation from a folder and return it as a pandas DataFrame.

    Parameters
    ----------
    foldername : string
        path of the folder containing the population data

    Returns
    -------
    simulation : dictionary
        contains a pandas DataFrame for each of the planets

    Example
    -------
    >>> simulation = read_simFromFolder(foldername)
    >>> planet005tracks = simulation['planet005']
    """
    import glob

    simulation = {}
    filenamepattern = foldername + "/tracks*.outputdat"
    print('Reading files: {}'.format(filenamepattern))
    for i, name in enumerate(sorted(glob.glob(filenamepattern))):
        if "tracks" in name and ".outputdat" in name:
            planetName = "planet{:03d}".format(i + 1)
            print('Reading {}...'.format(planetName))
            planetTracks = pd.read_csv(name, delim_whitespace=True, header=None)
            planetTracks = rename_tracksColumns(planetTracks)
            planetTracks.name = planetName
            simulation[planetName] = planetTracks
    return simulation


def read_popHdf5(filename, hierarchical=False, nSample=None):
    """Reads a population from a hdf5 file.

    The output is a dictionary of pandas panels that correspond to a simulation
    each. They contain the tracks of each planet in a DataFrame.

    Parameters
    ----------
    filename : string
        filename of the HDF5 file containing the population data
    hierarchical : bool
        if True, return a tree structure of systems and planets. Otherwise, all
        data is written into one data frame.
    nSample : integer
        number of samples randomly drawn from the simulations in the file.
        Recommended for large populations where it is not feasible to read
        everything into one table.
    Returns
    -------
    population : dict or pandas DataFrame
        depending on the 'hierarchical' flag: Dictionary of pandas panels or
        a single pandas DataFrame

    Example
    -------
    >>> population = read_popHdf5(filename)
    >>> SIM1planet005tracks = population['SIM1']['planet_005',:]
    """
    if hierarchical:
        # read into a hierarchical structure
        population = {}
        tab = tables.open_file(filename)
        for i, sim in enumerate(tab.walk_groups('/')):
            if i != 0:
                # ignore rootGroup
                print(sim)
                dfcontainer = {}
                for array in sim:
                    if "planet" in array.name:
                        # only planet tracks
                        df = pd.DataFrame(array.read())
                        df = rename_tracksColumns(df)
                        dfcontainer[array.name] = df
                population[sim._v_name] = pd.Panel.from_dict(dfcontainer)
    else:
        # read everything into one data frame
        import h5py
        h5file = h5py.File(filename, 'r')
        NoColumns = np.shape(h5file['SIM0001']['planet_001'])[1]
        population = pd.DataFrame(columns=range(NoColumns))

        if nSample:
            sample = np.sort(np.random.choice([sim for sim in h5file], nSample,
                          replace=False))
        else:
            sample = [sim for sim in h5file]

        for sim in sample:
            print(sim)
            sim = h5file.get(sim)
            for array in sim:
                if 'planet' in array:
                    population = pd.concat([population, pd.DataFrame(sim.get(array)[:])])
        population = rename_tracksColumns(population)

    return population


def read_ref_red(ref_redFilename):
    """Reads the content from a 'ref_redXeY' file into a pandas DataFrame.

    'ref_redXeY' files contain the data of all planets in a population at
    a specific time t = 'XeY' yr. This function transfers one such file
    into a DataFrame for further analysis.

    Parameters
    ----------
    ref_redFilename : string
        full path to the 'ref_red' file

    Returns
    -------
    tracks : pandas DataFrame
        DataFrame containing data of all planets
    """
    # Example
    # -------
    # >>> tracks = read_ref_red(ref_redFilename)
    # >>>
    # """
    tracks = pd.read_csv(ref_redFilename, delim_whitespace=True, header=None)
    tracks = rename_tracksColumns(tracks, ref_red=True)
    return tracks


def join_dataframes(simlist, ref_red):
    """ Join a simulation list with their resulting planet population.

    Parameters
    ----------
    simlist : Pandas DataFrame
        data frame containing the simulation list. Has to contain a column
        'simName'.
    ref_red : Pandas DataFrame
        data frame created from a 'ref_red' file. Contains the parameters of
        all systems of a population at a given time. Has to contain a column
        'isystemorig'.

    Returns
    -------
    jointDF : Pandas DataFrame
        data frame including data from both input tables. Sorted by a joint
        index from 'simName' and 'isystemorig' columns.
    """
    try:
        simlist.simName = simlist.simName.astype('int')
    except:
        pass

    return ref_red.set_index('isystemorig').join(simlist.set_index('simName'))


class Population():
    """ a planet population consisting of systems which in turn contain planets."""
    def __init__(self, data=None, name=None):
        if data is not None:
            self.data = self.readData(data)
        else:
            self.data = None
        self.name = name

    def __fileTypeWarning(self):
        print("""not able to determine file type. Please import manually by
                    using a suitable import function.""")

    def readData(self, populationFile):
        """ reads data into a pandas DataFrame.

        The method distinguishes between a single file and a list of files.
        In the latter case, data of all files are combined into one multiindex
        data frame, using parts of the filenames as level-0-index.

        Parameters
        ----------
        populationFile : string or list
            filename or list of filenames. The names should include 'hd5' or
            'ref_red' in order to distinguish between those file types.

        Returns
        -------
        data : pandas DataFrame
            data frame containing the population
        """
        if isinstance(populationFile, list):
            # create a multiindex data frame from multiple populations
            populationsData = [read_ref_red(f) if 'ref_red' in f else read_popHdf5(f) for f in populationFile]
            if not populationsData:
                self.__fileTypeWarning()
                return

            # guess names from filenames. suitable for e.g. '...ref_red0.5Msol.dat'
            populationsNames = [f[-11:-4] for f in populationFile]
            jointDF = pd.concat([p for p in populationsData], axis=0,
                                  keys=[name for name in populationsNames])
            self.data = jointDF
            return self.data

        if "ref_red" in populationFile:
            self.data = read_ref_red(populationFile)
            return self.data
        elif ("hd5" in populationFile) or ("hdf5" in populationFile):
            self.data = read_popHdf5(populationFile)
            return self.data
        else:
            # exception if file type could not be recognized
            self.__fileTypeWarning()

    def categorize(self):
        """ Sort planets into different categories."""
        self.categories = stats.categorize(self.data)

    def planetType(self, type):
        """ Restrict population to a certain planet type."""
        self.data = stats.filterPlanets(self.data, type)

    def get_typeStats(self):
        """ Compute statistics for specific planet types."""
        types = ['all', 'ltEarth', 'Earth', 'SuperEarth']
        typeStatsTable = {}

        for type in types:
            population_filtered = stats.filterPlanets(self.data, type)
            statistics = stats.get_typeStats(self.data, population_filtered)
            typeStatsTable[type] = statistics

        self.typeStats = pd.DataFrame(typeStatsTable)
        return self.typeStats
