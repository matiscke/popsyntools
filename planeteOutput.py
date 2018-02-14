""" Contains functions to analyze output of the planet formation Code 'Planete'
by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tables
import seaborn as sns

# Set plot style
sns.set(context='notebook', style='whitegrid', font_scale=1., palette='colorblind',
        rc={
'text.usetex':True,
'text.latex.unicode':True,
'font.family' : 'sans-serif',
'font.style'         : 'normal',
'font.variant'        : 'normal',
'font.weight'         : 'normal',
'font.stretch'        : 'normal',
'savefig.dpi'         : 400,
'lines.linewidth'   : 1.0,
'lines.markersize'      : 3.,
'figure.subplot.left'    : 0.13,    # the left side of the subplots of the figure
'figure.subplot.right'   : 0.96,   # the right side of the subplots of the figure
'figure.subplot.bottom'  : 0.13,   # the bottom of the subplots of the figure
'figure.subplot.top'     : 0.96,    # the top of the subplots of the figure
'figure.subplot.hspace'  : 0.0,    # height reserved for space between subplots
'axes.xmargin' : 0.02,             # default margin for autoscale
'axes.ymargin' : 0.02,
'legend.handletextpad' : 0.5,
'legend.handlelength' : 0.75,
'xtick.minor.size'     : 2.,
'ytick.minor.size'     : 2.
})
sns.set_color_codes()
#%%

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
    """Read a simulation from a folder and returns it as a pandas DataFrame.

    Parameters
    ----------
    foldername : string
        filename of the HDF5 file containing the population data

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


#%%
""" Some helper functions to interact with Planete and aid with calculations.
"""
def get_M0(rc, Sigma0, expo, r0=5.2):
    """Compute the total disk mass from initial condition parameters.

    Parameters
    ----------
    rc : float
        characteristic radius [au]
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    expo : float
        Power law slope
    r0 : float
        reference radius [au], in general 5.2 au

    Returns
    -------
    M0 : float
        total disk mass in solar masses
    """
    # Define au, Msol in cm and g
    au = 1.496e13
    Msol = 1.98855e33

    M0 = (2*np.pi)/(2-expo)*Sigma0*(r0*au)**expo*(rc*au)**(2-expo)
    return M0/Msol


def get_Sigma0(rc, M0, expo, r0=5.2):
    """Compute Sigma0, the gas surface density at the reference radius r0
    necessary for an initial total disk mass of M0.

    Parameters
    ----------
    rc : float
        characteristic radius [au]
    M0 : float
        total disk mass in solar masses
    expo : float
        Power law slope
    r0 : float
        reference radius [au], in general 5.2 au

    Returns
    -------
    Sigma0 : float
        gas surface density at 5.2 au [g/cm^2]
    """
    # Define au, Msol in cm and g
    au = 1.496e13
    Msol = 1.98855e33

    return (r0*au)**(-expo)*(rc*au)**(-2+expo)*(2-expo)*M0*Msol/(2*np.pi)


#%%
""" Plotting functions meant for single planet tracks.
"""
def plot_mass(tracks, ax):
    """plot a planet's total mass vs time"""
    ax.plot(tracks['t'],tracks['m'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_coreMass(tracks, ax):
    """plot core mass vs time"""
    ax.plot(tracks['t'],tracks['mCore'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('core mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_radius(tracks, ax):
    """plot radius vs time"""
    ax.plot(tracks['t'],tracks['r'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('radius [Jupiter radii]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_lum(tracks, ax):
    """ plot luminosity vs time"""
    ax.plot(tracks['t'], tracks['L'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('Luminosity [?]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax
