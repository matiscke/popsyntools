import pandas as pd
from utils import *
from output import *

def rmcf(a,system,aout=np.Infinity):
    if 'status' in system.columns:
        system=system.loc[system['status']==0]
    system=system.loc[system['a']<aout]
    if(system.empty):
        return 0
    else:
        sumM = sum(system['m'])
        sumDeno = sum(system['m']*(np.log10(a/system['a']))**2.)
        return sumM/sumDeno

def rmc(system,aout=np.Infinity):
    """
    Calculates the Radial Mass Concentration (RMC) of a given system as Dataframe
    RMC was introduced by Chambers, 2001 and used by Raymond, 2009 as a nondimensional measure for mass concentration in log of sma
    """    

    a = np.logspace(-3,1,100)
    #print([rmcf(asingle,system) for asingle in a])
    return max([rmcf(asingle,system,aout=aout) for asingle in a])


# Trappist for reference
t1=get_Trappist1data()
t1['rmc']=[rmc(t1) for pla in t1['a']]
t1['N'] = [t1['a'].count() for pla in t1['a']]

# Trappist is alreaedy sorted
t1['gradM'] = t1.shift(-1)['m'] - t1[:6]['m']
t1['gradIceCOMPLETO'] = t1.shift(-1)['ice'] - t1[:6]['ice']
t1['gradIceDorn'] = t1.shift(-1)['ice_dorn_2018_ucm'] - t1[:6]['ice_dorn_2018_ucm']
print(t1)

t1.to_csv('t1_with_system_properties.csv')

p = Population()
p.read_data(populationFile='/home0/mcwork/scratch/Julia/CD_J35_0p1Msun_50emb_10Myr_diffp2/ref_redevo_38850/ref_red5e9.dat')


p.data.rename(columns={95:'Rtransit'},inplace=True)
print(p.data['Rtransit'])
rmcvalues = []
rmcvaluesPlanets = []
NvaluesPlanets = []
NallvaluesPlanets = []
aout = 0.1
print(len(p.data))

data_red=p.data.loc[p.data['status']==0]

systems=data_red.groupby(['isystem'])

Nall = systems['m'].apply(lambda x : (x > 0.1).sum()) # sum counts "true"'s. here, x is a scalar
Na = systems[['m','a']].apply(lambda x : ( (x['m']>0.1) & (x['a']<aout) ).sum()) # here, x is a list
#print(data_red.loc[(data_red['a']<aout) & (data_red[95] < 2./11.21) & (data_red[95] >0.5/11.21)])
#NaRrange = data_red.loc[(data_red['a']<aout) & (data_red[95] < 2./11.21) & (data_red[95] >0.5/11.21)].groupby('isystem').size()
NaRrange = systems[['Rtransit','a']].apply(lambda x : ( (x['Rtransit'] < 2./11.21) & (x['Rtransit'] >0.5/11.21) & (x['a']<aout) ).sum())
#rmcs = systems.apply(rmc)
system_stats=systems[['a','m','Rtransit']].agg(np.mean,np.std)
print("\nNall:\n",Nall,"\n\nNa:\n",Na,"\n\nNaRrange:\n",NaRrange,"\n\nstats:\n",system_stats)
#print("\nNall:",Nall,"\nNa:",Na,"\nNaRrange:",NaRrange,"\nstats:",system_stats,"\n rmc:",rmcs)


import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import seaborn as sns
#sns.set(style="whitegrid",context="talk")
sns.set(style="ticks",context="talk")


ax=plt.subplot(313)
plt.hist(Na,np.arange(-0.5,30.5,1),label=r"$a < 0.1$ au and $m>0.1\,M_\oplus$",color='red')
#sns.distplot(Na,label=r"$a < 0.1$ au and $m>0.1\,M_\oplus$",color='yellow')
ax.set_ylim([0,156])
ax.grid(False)
ax.set_xlabel(r"Planets per system")
plt.legend()

ax2=ax.twinx()
ax2.set_ylim([0,156/777.])
ax2.set_yticks([0,0.05,0.1,0.15,0.2])
ax2.grid(True)

ax3=plt.subplot(312,sharex=ax)
plt.hist(NaRrange,np.arange(-0.5,25.5,1),label=r"$a < 0.1$ au and $0.5\,R_\oplus < r < 2\,R_\oplus$",color='green')
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.set_ylim([0,156])
#ax3.set_ylabel(r"Systems")
plt.legend()
ax3.grid(False)

ax4=ax3.twinx()
ax4.set_ylim([0,156/777.])
ax4.set_yticks([0,0.05,0.1,0.15,0.2])
ax4.grid(True)

ax5=plt.subplot(311,sharex=ax)
plt.hist(Nall,np.arange(-0.5,25.5,1), label=r"$m>0.1\,M_\oplus$")
plt.setp(ax5.get_xticklabels(), visible=False)
ax5.set_ylim([0,156])
plt.legend()
ax5.grid(False)

ax6=ax5.twinx()
ax6.set_ylim([0,156/777.])
ax6.set_yticks([0,0.05,0.1,0.15,0.2])
ax6.grid(True)

#plt.tight_layout(rect=[0, 0, 1., 1.],h_pad=0.3)

plt.gcf().set_size_inches(8,10)
plt.text(1.105,0.5,r"Fraction of systems",ha='left',va='center',rotation=90,transform=ax3.transAxes)
plt.text(-0.095,0.5,r"Systems",ha='right',va='center',rotation=90,transform=ax3.transAxes)
ax.set_xlim([-0.5,24.5])
#ax.tick_params(direction='out', length=6, width=2, colors='grey')

plt.tight_layout(pad=1.8,h_pad=0.8)
plt.savefig('../J35_Plots/number_of_planets.pdf')

exit()


# Count number of planets and assign back to same. Only works with .loc if the column does not exist
for i in range(1,data_red['isystem'].max()+1):
    data_red.loc[data_red['isystem']==i,'N'] = data_red.loc[(data_red['isystem']==i) & (data_red['a']>aout),'a'].count()
    data_red.loc[data_red['isystem']==i,'Ntot'] = data_red.loc[data_red['isystem']==i].count()
print(data_red[['a','m','N','Ntot','isystem','iplanet']])



for sys in [data_red[data_red['isystem']==i] for i in range(1,int(len(p.data)/50+1))]:
    N = sys.loc[(sys['a']<aout)]['a'].count()
    print("Number of planets in range:",N)
    for pla in sys['iplanet']: NvaluesPlanets.append(N)
    for pla in sys['iplanet']: NvaluesPlanets.append(N)
    if(len(sys.loc[(sys['a']<aout) & (sys['status']==0)]['a'])>1):
        r=rmc(sys,aout=aout)
        rmcvalues.append(r)
        for pla in sys['iplanet']:
            print(pla,r)
            rmcvaluesPlanets.append(r)
        print(sys['isystemorig'],r)
    else:
        for pla in sys['iplanet']:
            rmcvaluesPlanets.append(np.nan)
        print(sys['isystemorig'],r)
    

p.data['N'] = NvaluesPlanets

print(len(rmcvaluesPlanets))
p.data['rmc']=rmcvaluesPlanets
print(p.data)
print(p.data['rmc'])


data_red['gradM']=data_red.sort_values(['isystem','a']).shift(-1)['m'] - data_red.sort_values(['isystem','a']).loc[calc_gradient]['m'] # take the actual difference
data_red['gradIce']=data_red.sort_values(['isystem','a']).shift(-1)['115'] - data_red.sort_values(['isystem','a']).loc[calc_gradient]['115'] # 115 is ice mass fraction

p.data.to_csv('J35_5e9_with_rmc.csv')


'''
sys=p.data[p.data['isystemorig']==138]


import matplotlib.pyplot as plt

a = np.logspace(-4,4,1000)

rmcvalt1 =[rmcf(asingle,t1) for asingle in a]
rmcval = [rmcf(asingle,sys) for asingle in a]
rmcvalaout = [rmcf(asingle,sys,aout=0.1) for asingle in a]
aout = t1['a'][6]
aout = t1['a'][6]*1.5
print(aout)
print(rmc(sys))
print(rmc(sys,aout=aout))
print(rmc(t1))
#plt.scatter(sys[sys['status']==0]['a'],sys[sys['status']==0]['m'])
#plt.scatter(sys[sys['status']==0][sys['a']<aout]['a'],sys[sys['status']==0][sys['a']<aout]['m'], c = 'black',alpha=0.5)
#plt.scatter(t1['a'],t1['m'],c='red')
plt.scatter(sys[sys['status']==0]['a'],[-0.1 for i in sys[sys['status']==0]['m']], s = sys[sys['status']==0]['m']*10.,label='Whole System')
plt.scatter(sys[sys['status']==0][sys['a']<aout]['a'],[-0.7 for i in sys[sys['status']==0][sys['a']<aout]['m']], s=sys[sys['status']==0][sys['a']<aout]['m']*10., c = 'black',alpha=0.5,label='Cut System')
plt.scatter(t1['a'],[0.5 for i in t1['m']],s=t1['m']*10.,c='red',label='Trappist-1')
plt.semilogx()
plt.plot(a,rmcval)
plt.plot(a,rmcvalt1,c='red')
plt.plot(a,rmcvalaout,c='black',alpha=0.5)
plt.xlabel("a (au)")
plt.ylabel("RMC")
plt.legend()
plt.show()
'''
