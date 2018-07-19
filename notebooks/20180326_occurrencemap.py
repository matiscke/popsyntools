# coding: utf-8
refred = '/home/schlecker/phd/planete/outputs/J31/ref_red5e9.dat_J31_11220'

population = output.read_ref_red(refred)
from output import *
refred = '/home/schlecker/phd/planete/outputs/J31/ref_red5e9.dat_J31_11220'

population = output.read_ref_red(refred)
refred = '/home/schlecker/phd/planete/outputs/J31/ref_red5e9.dat_J31_11220'

population = read_ref_red(refred)
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
popul
population
population.head
population.head()
population.set_index(['systemNo'])
population.set_index(['systemNo']).head()
systems = population.set_index(['systemNo'])
systems[1]
systems.count()
population.count()
grouped = population.groupby('systemNo').count()
grouped
grouped = population.groupby('systemNo')
systems
import matplotlib.pyplot as plt
systems.index()
systems.index
len([i for i in systems.index])
system.loc['systemNo']
systems.loc['systemNo']
systems.loc[3]
systems.loc[34]
systems.loc[34444]
systems.loc[-1]
systems[-1]
len(systems)
systems.systemNo()
systems.systemNo
systems.systemNo.describe()
population.systemNo
population.systemNo.max
population.systemNo.max()
get_ipython().run_line_magic('pwd', '')
get_ipython().run_line_magic('ll', '')
get_ipython().run_line_magic('save', 'notebooks/20180326_occurrencemap')
get_ipython().run_line_magic('save', 'notebooks/20180326_occurrencemap 1-9999')
