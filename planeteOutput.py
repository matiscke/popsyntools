""" Contains functions to analyze output of the planet formation Code 'Planete'
by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import pandas as pd
import tables

# Define au, Msol, Gravitational Constant in cm, g, s
au = 1.496e13
Msol = 1.98855e33
G = 6.6740831e-8

def rename_tracksColumns(planetTracks):
    """Rename some of the columns of a planet tracks table.

    Parameters
    ----------
    planetTracks : Pandas DataFrame
        Table with the tracks of a single planet

    Returns
    -------
    planetTracks : Pandas DataFrame
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
     116: 'systemNo',
     117: 'planetNo',
     118: 'NoEmbryos'}
    return planetTracks.rename(columns=colnames)


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
            simulation[planetName] = planetTracks
    return simulation


def read_popHdf5(filename):
    """Reads a population from a hdf5 file.

    The output is a dictionary of pandas panels that correspond to a simulation
    each. They contain the tracks of each planet in a DataFrame.

    Parameters
    ----------
    filename : string
        filename of the HDF5 file containing the population data

    Returns
    -------
    population : dict
        Dictionary of pandas panels

    Example
    -------
    >>> population = read_popHdf5(filename)
    >>> SIM1planet005tracks = population['SIM1']['planet_005',:]
    """
    # read hdf5 file with pytables
    tab = tables.open_file(filename)
    population = {}
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
    tracks = rename_tracksColumns(tracks)
    return tracks
