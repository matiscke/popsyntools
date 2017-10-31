#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#from decimal import Decimal

# Set plot style
sns.set(context='notebook', style='whitegrid', font_scale=1., palette='colorblind',\
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
#out001 = np.genfromtxt('outputs/tracks_002.outputdat')
out001 = np.genfromtxt('outputs/bernNov17/tracks_002.outputdat')
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
''' let's plot something
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

