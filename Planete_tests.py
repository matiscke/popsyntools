""" Contains functions to handle output data of the Planet Population
Synthesis Code 'Planete' by the Bern planet formation group.

Written by: Martin Schlecker
schlecker@mpia.de
"""
#%%
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
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
out001 = np.genfromtxt('/media/martin/Daten/phd/planete/outputs/singleplanet/03dusttogasUnity/tracks_001.outputdat')
# out001 = np.genfromtxt('outputs/bernNov17/tracks_002.outputdat')
tracks = pd.DataFrame(out001)
#tracks = pd.read_table('outputs/tracks_001.outputdat',header=None,sep='  ')
tracks = tracks.rename(columns = {1:'t',2:'mCore',4:'m',5:'L',14:'r'})
#disk = pd.read_table('outputs/structure_disk2.outputdat',header=None,sep='  ')
#disk = disk.rename(columns = {9:'t',1:'r',2:'Sigma'})


### if planet parameters read by np.genfromtxt:
#t=out001[:,1]
#m=out001[:,4]
#mCore=out001[:,2]
#r=out001[:,14]
#L=out001[:,5]

#%%

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
    >>> SIM1planet001tracks = population['SIM1']['planet_005',:]
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
                    df = df.rename(columns = {1:'t',2:'mCore',4:'m',5:'L',14:'r'})
                    dfcontainer[array.name] = df
            population[sim._v_name] = pd.Panel.from_dict(dfcontainer)
    return population

# # read hdf5
# filename = '/media/martin/Daten/phd/planete/outputs/'
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
fig,ax = plt.subplots(4)
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

ax[0] = plot_mass(tracks, ax[0])
ax[1] = plot_coreMass(tracks,ax[1])
ax[2] = plot_radius(tracks, ax[2])
ax[3] = plot_lum(tracks,ax[3])
plt.show()
