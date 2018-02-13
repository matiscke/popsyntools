""" Contains functions to handle output data of the Planet Population
Synthesis Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
# import matplotlib
# matplotlib.use('Qt5Agg')
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
#'font.serif':'Computer Modern',
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
### read simulation results from file
# out001 = np.genfromtxt('/media/martin/Daten/phd/planete/outputs/singleplanet/03dusttogasUnity/tracks_001.outputdat')
# out001 = np.genfromtxt('outputs/bernNov17/tracks_002.outputdat')
# tracks = pd.DataFrame(out001)
#tracks = pd.read_table('outputs/tracks_001.outputdat',header=None,sep='  ')
# tracks = tracks.rename(columns = {1:'t',2:'mCore',4:'m',5:'L',14:'r'})
#disk = pd.read_table('outputs/structure_disk2.outputdat',header=None,sep='  ')
#disk = disk.rename(columns = {9:'t',1:'r',2:'Sigma'})


### if planet parameters read by np.genfromtxt:
#t=out001[:,1]
#m=out001[:,4]
#mCore=out001[:,2]
#r=out001[:,14]
#L=out001[:,5]

#%%

"""
DEVELOPMENT HELPER FUNCTION:
The code below is used to generate a dictionary with ONLY column numbers and
column names from the quantities dict of synth.py
"""

# 'quantities' from synth.py
quantities = {'a': {'all': False,
  'column': 18,
  'label': 'Semi-major axis [AU]',
  'pmax': 1000.0,
  'pmin': 0.1,
  'scale': 'log'},
 'ae': {'all': False,
  'column': 18,
  'label': 'Semi-major axis [AU]',
  'scale': 'log',
  'test': True},
 'correesti': {'all': False,
  'column': 102,
  'label': 'corrEesti',
  'scale': 'log'},
 'corrlesti': {'all': False,
  'column': 101,
  'label': 'corrLesti',
  'scale': 'log'},
 'dt': {'all': False,
  'column': 24,
  'compute': '+dtall',
  'label': 'Timestep [1]',
  'max': 100000.0,
  'min': 0.1,
  'scale': 'log'},
 'dtall': {'all': True,
  'column': 1,
  'compute': 'diff',
  'label': 'timestep [1]',
  'scale': 'log'},
 'dtmode': {'all': False,
  'column': 78,
  'label': 'Timestep limitation factor',
  'scale': 'linear'},
 'dtpl': {'all': False,
  'column': 39,
  'label': 'Time step lim [1]',
  'scale': 'log'},
 'e': {'all': False,
  'column': 62,
  'label': 'Eccentricity',
  'pmax': 1.0,
  'pmin': 1e-08,
  'scale': 'log'},
 'ediff': {'all': False,
  'column': 53,
  'compute': 'ediff',
  'label': 'E_{diff}',
  'scale': 'log'},
 'enew': {'all': False, 'column': 52, 'label': 'E_{new}', 'scale': 'log'},
 'eold': {'all': False, 'column': 53, 'label': 'E_{old}', 'scale': 'log'},
 'ep': {'all': False, 'column': 56, 'label': 'e_p', 'scale': 'log'},
 'i': {'all': False,
  'column': 63,
  'label': 'Inclination',
  'pmax': 10.0,
  'pmin': 1e-15,
  'scale': 'log'},
 'ip': {'all': False, 'column': 57, 'label': 'i_p', 'scale': 'log'},
 'kenerg': {'all': False, 'column': 54, 'label': 'kenerg', 'scale': 'log'},
 'kenergdiff': {'all': False,
  'column': 54,
  'compute': 'diff',
  'label': 'kenerg',
  'max': 0.2,
  'min': -0.1,
  'scale': 'linear'},
 'lacccore': {'all': False,
  'column': 55,
  'label': 'L_{core,acc}',
  'scale': 'log'},
 'lactual': {'all': False,
  'column': 100,
  'label': 'L_{actual}',
  'scale': 'log'},
 'lcomp': {'all': False,
  'column': 6,
  'compute': 'lcomp',
  'label': 'L_{sum}',
  'scale': 'log'},
 'lcont': {'all': False,
  'column': 51,
  'label': 'L_{env,cont}',
  'scale': 'log'},
 'lcontcore': {'all': False,
  'column': 108,
  'label': 'L_{core,cont}',
  'scale': 'log'},
 'lcontenv': {'all': False,
  'column': 107,
  'label': 'L_{env,cont}',
  'scale': 'log'},
 'lcontsum': {'all': False,
  'column': 108,
  'compute': 'lcont',
  'label': 'L_{cont}',
  'scale': 'log'},
 'lcore': {'all': False,
  'column': 6,
  'label': 'L_{core} [L_1]',
  'scale': 'log'},
 'ltot': {'all': False, 'column': 5, 'label': 'L_{tot} [L_1]', 'scale': 'log'},
 'macctot': {'all': True,
  'column': 33,
  'label': 'Accreted solid mass [M_1]',
  'scale': 'linear'},
 'mcore': {'all': False,
  'column': 2,
  'label': 'Core mass [M_1]',
  'pmax': 500.0,
  'pmin': 0.01,
  'scale': 'log'},
 'mdiskg': {'all': True,
  'column': 20,
  'label': 'Gas disc mass [M_1]',
  'scale': 'linear'},
 'mdiskp': {'all': True,
  'column': 21,
  'label': 'Planetesimal disc mass [M_1]',
  'scale': 'linear'},
 'mdotcore': {'all': False,
  'column': 7,
  'label': 'Mdot core [M_1/1]',
  'max': 0.01,
  'min': 1e-08,
  'scale': 'log'},
 'mdotgas': {'all': False,
  'column': 73,
  'label': 'Gas acc rate [M_1/1]',
  'scale': 'log'},
 'mdotgasl': {'all': False,
  'column': 8,
  'label': 'Mdot gas [M_1/1]',
  'max': 0.01,
  'min': 1e-08,
  'scale': 'log'},
 'mdotgasmax': {'all': False,
  'column': 104,
  'label': 'Max gas acc rate [M_1/1]',
  'scale': 'log'},
 'mdotgasratio': {'all': False,
  'column': 8,
  'compute': 'mdotgasratio',
  'label': 'Mdotgas/Mdotgasmax',
  'max': 1.1,
  'min': 0.0001,
  'scale': 'log'},
 'mejetot': {'all': True,
  'column': 32,
  'label': 'Ejected solid mass [M_1]',
  'scale': 'linear'},
 'menv': {'all': False,
  'column': 3,
  'label': 'Envelope mass [M_1]',
  'pmax': 10000.0,
  'pmin': 0.01,
  'scale': 'log'},
 'mgazacc': {'all': False,
  'column': 42,
  'label': 'Accreted gas mass [M_1]',
  'scale': 'linear'},
 'mgazevap': {'all': False,
  'column': 43,
  'label': 'Evaporated gas mass [M_1]',
  'scale': 'linear'},
 'miso': {'all': False,
  'column': 34,
  'label': 'Isolation mass [M_1]',
  'scale': 'linear'},
 'mplan': {'all': False,
  'column': 4,
  'label': 'Planet mass [M_1]',
  'pmax': 10000.0,
  'pmin': 0.01,
  'scale': 'log'},
 'n': {'all': True,
  'column': 0,
  'label': 'n',
  'scale': 'linear',
  'scilimits': (0, 3)},
 'pcore': {'all': False,
  'column': 11,
  'label': 'Central pressure [bar]',
  'scale': 'log'},
 'pneb': {'all': False, 'column': 30, 'label': 'Pneb [bar]', 'scale': 'log'},
 'pout': {'all': False, 'column': 44, 'label': 'Pout [bar]', 'scale': 'log'},
 'racc': {'all': True,
  'column': 23,
  'label': 'Bondi radius [R_1]',
  'scale': 'log'},
 'rcap': {'all': False,
  'column': 27,
  'label': 'Capture radius [R_1]',
  'scale': 'log'},
 'rcapcore': {'all': False,
  'column': 27,
  'compute': 'rcapcore',
  'label': 'Rcap/Rcore',
  'scale': 'log'},
 'rcaptot': {'all': False,
  'column': 27,
  'compute': 'rcaptot',
  'label': 'Rcap/Rtot',
  'scale': 'log'},
 'rcore': {'all': False,
  'column': 9,
  'label': 'Core radius [R_1]',
  'scale': 'log'},
 'rfeed': {'all': False,
  'column': 26,
  'label': 'Feeding zone width [AU]',
  'scale': 'log'},
 'rhocen': {'all': False,
  'column': 13,
  'label': 'Central density [g/cm^3]',
  'scale': 'linear'},
 'rhocore': {'all': False,
  'column': 41,
  'label': 'Core density [g/cm^3]',
  'scale': 'linear'},
 'rratio': {'all': False,
  'column': 14,
  'compute': 'rratio',
  'label': 'Rtot/Rcore',
  'scale': 'log'},
 'rroche': {'all': True,
  'column': 22,
  'label': 'Hills radius [R_1]',
  'scale': 'log'},
 'rtot': {'all': False,
  'column': 14,
  'label': 'Outer radius [R_1]',
  'scale': 'log'},
 'sigmagas': {'all': False,
  'column': 35,
  'label': 'Gas surface density at planet location [g/cm^2]',
  'scale': 'log'},
 'sigmamean': {'all': False,
  'column': 25,
  'label': 'Mean solid surface density [g/cm^2]',
  'scale': 'log'},
 't': {'all': True,
  'column': 1,
  'label': 'Time [1]',
  'max': 10000000.0,
  'min': 1000.0,
  'pmax': 10000000.0,
  'pmin': 1000.0,
  'scale': 'log'},
 'tcore': {'all': False,
  'column': 12,
  'label': 'Central temperature [K]',
  'scale': 'log'},
 'tmig': {'all': False,
  'column': 69,
  'label': 'Migration timescale [1]',
  'linthresh': 10000.0,
  'scale': 'symlog'},
 'tmige': {'all': False,
  'column': 70,
  'label': 'Ecc. damp. timescale [1]',
  'pmax': 10000000000.0,
  'pmin': 0.1,
  'scale': 'log'},
 'tneb': {'all': False, 'column': 29, 'label': 'Tneb [K]', 'scale': 'log'},
 'tout': {'all': False, 'column': 45, 'label': 'Tout [K]', 'scale': 'log'},
 'type_mig': {'all': False,
  'column': 31,
  'label': 'Migration type',
  'scale': 'linear'},
 'typemig': {'all': False,
  'column': 66,
  'label': 'Type of migration',
  'max': 20.0,
  'min': 0.0,
  'scale': 'linear'}}

def get_colnames(quantities):
    """ function to extract column numbers and names from the "quantitites"
    dictionary in synth.py."""
    colnames = {}
    for varname, values in quantities.items():
        columnNo = values['column']
        columnName = varname
        colnames[columnNo] = columnName
    return colnames

colnames = get_colnames(quantities)

""" END DEVELOPMENT HELPER function
"""




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
    colnames = {1:'t',2:'mCore',4:'m',5:'L',14:'r'}
    return planetTracks.rename(columns=colnames)


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

def read_simFromFolder(foldername):
    """Read a simulation from a folder and returns it as a pandas DataFrame.

    Parameters
    ----------
    foldername : string
        filename of the HDF5 file containing the population data

    Returns
    -------
    simulation : pandas DataFrame
        DataFrame containing the tracks of the planets

    Example
    -------
    >>> simulation = read_simFromFolder(foldername)
    >>> planet005tracks = simulation['planet_005',:]
    """
    import glob

    simulation = {}
    filenamepattern = foldername + "/tracks*.outputdat"
    for i, name in enumerate(sorted(glob.glob(filenamepattern))):
        if "tracks" in name and ".outputdat" in name:
            simname = "SIM{:03d}".format(i + 1)
            planetTracks = pd.read_csv(name, delim_whitespace=True, header=None)
            planetTracks = rename_tracksColumns(planetTracks)
            simulation[simname] = planetTracks
    return simulation

# # read hdf5
# filename = '/home/schlecker/phd/planete/outputs/bernNov17/popu/popu.hdf5'
# population = read_popHdf5(filename)




#%%
''' The following plotting functions are meant for single planet tracks.
'''
def plot_mass(tracks, ax):
    # plot mass
    ax.plot(tracks['t'],tracks['m'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_coreMass(tracks, ax):
    ## plot core mass
    ax.plot(tracks['t'],tracks['mCore'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('core mass [$m_{Earth}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_radius(tracks, ax):
    ## plot radius
    ax.plot(tracks['t'],tracks['r'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('radius [Jupiter radii]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

def plot_lum(tracks, ax):
    ''' plot luminosity vs time'''
    ax.plot(tracks['t'], tracks['L'])
    ax.set_xlabel('time [yr]')
    ax.set_ylabel('Luminosity [?]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return ax

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
# fig,ax = plt.subplots(4)

### plot disk surface density
#ntimesteps=12
#step = int(len(disk['t'])/ntimesteps)
#print('tmin = {} yr'.format(min(disk['t'])))
#print('tmax = {} yr'.format(max(disk['t'])))
#for s in range(ntimesteps+1):
#    if s == ntimesteps:
#        # for last plot, use maximum time value
#        subset = disk.loc[disk['t'] == max(disk['t'])]
#    else:
#        subset = disk.loc[disk['t'] == disk.iloc[s*step]['t']]
#    ax.plot(subset['r'],subset['Sigma'],label='t={:.2E} yr'.format(Decimal(subset.iloc[0]['t'])))
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('radius [au]')
#ax.set_ylabel('gas surface density [g/cm2]')
#plt.legend()

### TEST for overlay plot of surface density (not working)
#ntimes=5
#step = int(len(disk['t'])/ntimes)
#sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#g = sns.FacetGrid(disk[::step], row="t", hue="t", aspect=15, size=.5, palette=pal)
#g.map(plt.plot, "r", "Sigma", alpha=1, lw=1.5)
#g.map(plt.plot, "r", "Sigma", color="w", lw=2)
##g.map(plt.axhline, y=0, lw=2, clip_on=False)

# ax[0] = plot_mass(tracks, ax[0])
# ax[1] = plot_coreMass(tracks,ax[1])
# ax[2] = plot_radius(tracks, ax[2])
# ax[3] = plot_lum(tracks,ax[3])
# plt.show()
