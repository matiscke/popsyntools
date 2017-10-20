import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Set plot style
sns.set(context='notebook', style='ticks', font_scale=1.3, palette='colorblind',\
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

### read simulation results from file
#out = pd.read_table('outputs/tracks_001.outputdat',header=None,sep='  ')
disk = pd.read_table('outputs/structure_disk.outputdat',header=None,sep='  ')
disk = disk.rename(columns = {9:'t',1:'r',2:'Sigma'})
## read planet parameters
#t=out001[:,1]
#m=out001[:,4]
#mCore=out001[:,2]
#r=out001[:,14]
#L=out001[:,5]


''' let's plot something
'''
#fig,ax = plt.subplots()
# plot mass
#ax.plot(t,m)
#ax.set_xlabel('time (arbitrary units)')
#ax.set_ylabel('mass (Earth masses)')

## plot core mass
#ax.plot(t,mCore)
#ax.set_xlabel('time (arbitrary units)')
#ax.set_ylabel('core mass (Earth masses)')

## plot radius
#ax.plot(t,r)
#ax.set_xlabel('time [yr]')
#ax.set_ylabel('radius [Jupiter radii]')



# plot disk surface density
ntimes=5
step = int(len(disk['t'])/ntimes)
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
g = sns.FacetGrid(disk[::step], row="t", hue="t", aspect=15, size=.5, palette=pal)
g.map(plt.plot, "r", "Sigma", alpha=1, lw=1.5)
g.map(plt.plot, "r", "Sigma", color="w", lw=2)
#g.map(plt.axhline, y=0, lw=2, clip_on=False)






#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('radius [au]')
#ax.set_ylabel('gas surface density [g/cm2]')


