
# coding: utf-8

# In[6]:


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

#sns.set_color_codes()

# read simulation results from file
out001 = np.genfromtxt('outputs/tracks_001.outputdat')
disk = np.genfromtxt('outputs/structure_disk.outputdat')


## read planet parameters
#t=out001[:,1]
#m=out001[:,4]
#mCore=out001[:,2]
#r=out001[:,14]
#L=out001[:,5]

# read disk parameters
t=disk[:,9]
r=disk[:,1]
Sigma=disk[:,2]


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
#for time in range(4):#t[:10:]:
#    ax.plot(r,np.sin(r))
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel('radius [au]')
#ax.set_ylabel('gas surface density [g/cm2]')
## plot it the fancy way ('joy plot')
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
step = int(len(t)/10)
g = np.tile(t[::step],step)
df = pd.DataFrame(dict(Sigma=Sigma, g=g))
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
g = sns.FacetGrid(df, row="g", hue="g", aspect=15, size=.5, palette=pal)
# Draw the densities in a few steps
g.map(sns.kdeplot, "Sigma", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g.map(sns.kdeplot, "Sigma", clip_on=False, color="w", lw=2, bw=.2)
g.map(plt.axhline, y=0, lw=2, clip_on=False)
plt.xscale('log')
plt.yscale('log')

# Define and use a simple function to label the plot in axes coordinates
def label(r, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color, 
            ha="left", va="center", transform=ax.transAxes)
g.map(label, "Sigma")

# Set the subplots to overlap
g.fig.subplots_adjust(hspace=-.25)

# Remove axes details that don't play will with overlap
g.set_titles("")
g.set(yticks=[])
g.despine(bottom=True, left=True)
#plt.show()