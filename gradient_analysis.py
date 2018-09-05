import pandas as pd
from scipy import optimize

data = pd.read_csv('J35_5e9_with_rmc.csv')
import numpy as np

whereActive=data['status']==0
data_red=data.loc[whereActive]


# Thats where the pandas magic happens
calc_gradient=data_red.sort_values(['isystem','a'])['isystem'] == data_red.sort_values(['isystem','a'])['isystem'].shift(-1) # check where the systems are the same after shift
data_red['gradM']=data_red.sort_values(['isystem','a']).shift(-1)['m'] - data_red.sort_values(['isystem','a']).loc[calc_gradient]['m'] # take the actual difference
data_red['gradIce']=data_red.sort_values(['isystem','a']).shift(-1)['115'] - data_red.sort_values(['isystem','a']).loc[calc_gradient]['115'] # 115 is ice mass fraction

print(data_red[['a','iplanet','isystem','m','gradM','gradIce']].describe())

print(data_red.loc[data_red['isystemorig']==138].sort_values('a')[['a','m','gradM','gradIce']])
print(data_red.loc[data_red['a']<0.1]['gradM'].describe())

print(data_red.loc[data_red['a']<0.1]['gradIce'].describe())


print(data_red)


import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#plt.semilogx()
#plt.scatter(data_red.loc[data_red['a']<0.1]['a'],data_red.loc[data_red['a']<0.1]['gradM'],alpha=0.2)
xdatalt0p1 = data_red.loc[(data_red['a']<0.1) & (data_red['m']>0.1)]['a']
ydatalt0p1 = data_red.loc[(data_red['a']<0.1) & (data_red['m']>0.1)]['gradM']
xdatamt0p1 = data_red.loc[(data_red['a']>0.1) & (data_red['m']>0.1)]['a']
ydatamt0p1 = data_red.loc[(data_red['a']>0.1) & (data_red['m']>0.1)]['gradM']
xdata = data_red.loc[data_red['m']>0.1]['a']
ydata = data_red.loc[data_red['m']>0.1]['gradM']

idx = np.isfinite(xdata) & np.isfinite(ydata)
idxlt0p1 = np.isfinite(xdatalt0p1) & np.isfinite(ydatalt0p1)
idxmt0p1 = np.isfinite(xdatamt0p1) & np.isfinite(ydatamt0p1)

print(xdata[idx])
print(ydata[idx]) 
print(np.log10(xdata[idx]))
fit=np.polyfit(np.log10(xdata[idx]),ydata[idx],deg=1,full=True)
fitlt0p1=np.polyfit(np.log10(xdatalt0p1[idxlt0p1]),ydatalt0p1[idxlt0p1],deg=1,full=True)
fitmt0p1=np.polyfit(np.log10(xdatamt0p1[idxmt0p1]),ydatamt0p1[idxmt0p1],deg=1,full=True)

print(fit)
p = np.poly1d(fit[0])
print(fitlt0p1)
plt0p1 = np.poly1d(fitlt0p1[0])
print(fitmt0p1)
pmt0p1 = np.poly1d(fitmt0p1[0])

def piecewise_linear(x, x0, y0, k1, k2):
    #print("In piecewise_linear, x:",x[0:5],"x0:",x0,"y0:",y0,"k1:",k1,"k2:",k2)
    #print("for first element:")
    #if(x[0]<x0): print(k1*x[0]+y0-k1*x0)
    #else: print(k2*x[0] + y0-k2*x0)
    vals=[]
    for xv in x:
        if(xv<x0): r=k1*xv+y0-k1*x0
        else: r=k2*xv + y0-k2*x0
        vals.append(r)
    #print("values:",vals[0:5])
    return vals
    #return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])



print("xdata:",xdata[idx])
print("ydata:",ydata[idx])

#pbroken , e = optimize.curve_fit(piecewise_linear, xd, yd)
#pbroken , e = optimize.curve_fit(piecewise_linear, np.log10(xdata[idx]), ydata[idx])
pbroken , e = optimize.curve_fit(piecewise_linear, np.log10(xdata[idx]), ydata[idx], p0=[0.2,0,0.1,-0.1])
print(pbroken,e)

fig1=plt.figure(1)
plt.semilogx()
plt.scatter(xdata,ydata,alpha=0.2,label=r'$\nabla M$ data')
xp = np.linspace(-4, 2, 200)
plt.plot(10**(xp),p(xp), c='red', lw=3, zorder=10, label=r'fit for all data: $\propto 10^{{{:.2f}}}$'.format(fit[0][0]))
plt.plot(10**(xp),plt0p1(xp), c='green', lw=3, zorder=10,label=r'fit on $a < 0.1$ au: $\propto 10^{{{:.2f}}}$'.format(fitlt0p1[0][0]))
plt.plot(10**(xp),pmt0p1(xp), c='yellow', lw=3, zorder=10, label="fit on $a > 0.1$ au: $\propto 10^{{{:.2f}}}$".format(fitmt0p1[0][0]))
plt.plot(10**(xp), piecewise_linear(xp, *pbroken), c='black',lw=3, zorder=10, label="fit: $\propto 10^{{{:.3f}}}$ for $a < {{{:.2f}}}$ au, else: $\propto 10^{{{:.3f}}}$".format(pbroken[2],10**(pbroken[0]),pbroken[3]))
plt.legend()
plt.xlabel("a (au)")
plt.ylabel(r"Mass difference ($M_\oplus$)")
plt.xlim([5*10**(-4),2*10**2])
plt.ylim([-5,7])
plt.savefig("J35_dm_analysis.pdf")


# Ice Gradient Analysis (care, array-names reused)
xdata = data_red.loc[data_red['m']>0.1]['a']
ydata = data_red.loc[data_red['m']>0.1]['gradIce']
xdatalt0p1 = data_red.loc[(data_red['a']<0.3) & (data_red['m']>0.1)]['a']
ydatalt0p1 = data_red.loc[(data_red['a']<0.3) & (data_red['m']>0.1)]['gradIce']
xdatamt0p1 = data_red.loc[(data_red['a']>0.3) & (data_red['m']>0.1)]['a']
ydatamt0p1 = data_red.loc[(data_red['a']>0.3) & (data_red['m']>0.1)]['gradIce']

idx = np.isfinite(xdata) & np.isfinite(ydata)
idxlt0p1 = np.isfinite(xdatalt0p1) & np.isfinite(ydatalt0p1)
idxmt0p1 = np.isfinite(xdatamt0p1) & np.isfinite(ydatamt0p1)

fit=np.polyfit(np.log10(xdata[idx]),ydata[idx],deg=1,full=True)
fitlt0p1=np.polyfit(np.log10(xdatalt0p1[idxlt0p1]),ydatalt0p1[idxlt0p1],deg=1,full=True)
fitmt0p1=np.polyfit(np.log10(xdatamt0p1[idxmt0p1]),ydatamt0p1[idxmt0p1],deg=1,full=True)
p = np.poly1d(fit[0])
plt0p1 = np.poly1d(fitlt0p1[0])
pmt0p1 = np.poly1d(fitmt0p1[0])

print("ice linear a < 0.3 au:",fitlt0p1)
print("ice linear a > 0.3 au:",fitmt0p1)

print("ice linear fit:", fit)
pbroken , e = optimize.curve_fit(piecewise_linear, np.log10(xdata[idx]), ydata[idx], p0=[0.2,0,0.1,-0.1])
print("ice broken fit:",pbroken,e)


fig2 = plt.figure(2)
plt.semilogx()
plt.scatter(xdata,ydata,alpha=0.2,label=r'$\nabla$ ice data')
xp = np.linspace(-4, 2, 200)
#plt.plot(10**(xp),p(xp), c='red', lw=3, zorder=10, label=r'fit for all data: $\propto 10^{{{:.3f}}}$'.format(fit[0][0]))
#plt.plot(10**(xp),plt0p1(xp), c='green', lw=3, zorder=10,label=r'fit on $a < 0.3$ au: $\propto 10^{{{:.3f}}}$'.format(fitlt0p1[0][0]))
#plt.plot(10**(xp),pmt0p1(xp), c='yellow', lw=3, zorder=10, label="fit on $a > 0.3$ au: $\propto 10^{{{:.3f}}}$".format(fitmt0p1[0][0]))
plt.plot(10**(xp), piecewise_linear(xp, *pbroken), c='black',lw=3, zorder=10, label="fit: $\propto 10^{{{:.3f}}}$ for $a < {{{:.2f}}}$ au, else: $\propto 10^{{{:.3f}}}$".format(pbroken[2],10**(pbroken[0]),pbroken[3]))
plt.legend()
plt.xlabel("a (au)")
plt.ylabel("Ice mass fraction difference")
plt.xlim([5*10**(-4),2*10**2])
plt.ylim([-0.5,0.7])
plt.savefig("J35_dice_analysis_l.pdf")

plt.show()
