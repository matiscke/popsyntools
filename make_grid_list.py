import numpy as np

# constants (cgs)
au=1.495979e13
Msun=1.9891e33

# constant values
expo = 0.9
CDnumber= 0
a_start = 0
t_start = 2.e5

f= open("grid_list.dat","w")
simNumber = 1

#Grid parameters: mstar, windM, mdisk
for mstar in np.arange(0.1,1.2,0.15):
	for windM in np.logspace(-10,-5,20):
		for mdisk in np.logspace(-3,-0.5,10):
			mdisk=mdisk*mstar # Mdisk in solar masses, linear scaling
			a_out=10*(mdisk / 0.002)**(1./1.6) # Andrews 2010, Fig. 10
			a_in = 0.02222222*mstar+0.00777777778
			sigma=(5.2*au)**(-expo)*(a_out*au)**(-2+expo)*(2-expo)*mdisk*Msun/(2*np.pi)
			print("mstar:",mstar,"windM:",windM,"mdisk:",mdisk)
			print("sigma:",sigma, 'a_out',a_out,'a_in:',a_in)
			fgp = 10.**(-0.02)*0.02 # does this make a difference for lifetimes?
			print("CD_{:014d}".format(CDnumber)+"FP_{:.8E}".format(1./fgp)+"SI_{:.8E}".format(sigma)+"AI_{:.8E}".format(a_in)+"AO_{:.8E}".format(a_out)+"EX_{:.8E}".format(expo)+"MW_{:.8E}".format(windM)+"SIM{:014d}".format(simNumber)+"AS_{:.8E}".format(a_start)+"ST_{:.8E}".format(t_start)+"MS_{:.8E}".format(mstar),file=f)
			simNumber += 1
