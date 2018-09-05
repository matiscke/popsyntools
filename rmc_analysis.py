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
t1.to_csv('t1_with_rmc.csv')

p = Population()
p.read_data(populationFile='/home0/mcwork/scratch/Julia/CD_J35_0p1Msun_50emb_10Myr_diffp2/ref_redevo_38850/ref_red5e9.dat')

rmcvalues = []
rmcvaluesPlanets = []
NvaluesPlanets = []
aout = 0.1
print(len(p.data))
for sys in [p.data.loc[p.data['isystem']==i] for i in range(1,int(len(p.data)/50+1))]:
    N = sys.loc[(sys['a']<aout) & (sys['status']==0)]['a'].count()
    print("Number of planets in range:",N)
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
