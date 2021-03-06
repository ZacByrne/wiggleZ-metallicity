import numpy as np
#import os
import matplotlib, matplotlib.pyplot as plt
#import collections
#import math as math

#constant E(b-v)_gamma. Might later include calculations?
#EBV_gamma = 0.61113184

#adjust line flux for reddening as flux * 10^(-0.4*E(b-v)*k(lambda))  as per ly_2016
def adjustflux(flux, kred,EBVg):
	flux = flux * 10**(0.4*EBVg*kred)
	return flux



#load values
#import args from system? (Google)
# import line flux values from file into matrix.. NEED TO SEE OUTPUT FILE... CHANGE TO USECOL()
col, oiiiw, oiiiwe, oii, oiie, oiii, oiiie, hb, hbe, EW, hg, hge = np.loadtxt("fluxes.txt", dtype = {'names': ('column name', '4363', '4363 error ', '3728', '3278 error', '4959', '4959 error', 'hbeta', 'hbeta error', 'hgamma', 'hgamma error'), 'formats': ('S4', 'f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')}, skiprows=1, unpack=True)

col1, avm, avmg, avr = np.loadtxt("averages.txt", dtype = {'names':('column name', 'Av mass', 'Av Mg', 'Av Redshift'), 'formats': ('S4', 'f5', 'f5', 'f5')})

# calculate k values? Or simply import from matlab program?
k4363 = 4.1482
k3728 = 5.1632
k4959 = 3.5163
k4861 = 3.6092
kgamma = 4.17




# import mass, sfr, redshift? Or calc?np.random.normal(size = 4)

#Loop stats
count = 50

summet = 0
summet2 = 0
sumsfr = 0
sumsfr2 = 0 
#Need to show something else? av e temp? Av oII met?


#assume each variable is an array of variable named column
for x in range(0,count):
    
    #Uncertainity calculations:
    
    while not oiiiwu.all() > 0:
        oiiiwu = oiiiw + np.random.normal(size = (len(oiiiw))) * oiiiwe
        
    oiiiu = oiii + np.random.normal(size = (len(oiiiw))) * oiiie
    oiiu = oii + np.random.normal(size = (len(oiiiw))) * oiie
    hbu = hb + np.random.normal(size = (len(oiiiw))) * hbe
    hgu = hg + np.random.normal(size = (len(oiiiw))) * hge
    
    #calculate reddening value
    EBV = np.log10((hgu/hbu)/0.468)/(-0.4*(kgamma-k4861))
    
    print EBV    
    


    # Adjust line fluxes for reddening as per Ly 2016
    oiiiwa = adjustflux(oiiiwu, k4363, EBV)
    oiia = adjustflux(oiiu, k3728, EBV)
    oiiia = adjustflux(oiiiu, k4959, EBV)
    hba = adjustflux(hbu, k4861, EBV)

    # calc OIII temp as per Nichollis et al. 2014 
    oiiitemp = 13205 * (-1*np.log10((oiiiwa / (oiiia*4)))-0.92506)**(-0.98062)

    #Calc OII temp as per Andrews and Martini 2013
    oiitemp = 0.7 * oiiitemp + 0.17 * 10000


    #Calculate OII and OIII metallicity
    oiimet = np.log10(oiia/hba)+5.961+1.676/(oiitemp/10000) - 0.4*np.log10(oiitemp/10000)-0.034*oiitemp/10000
    oiiimet = np.log10(4*oiiia/hba)+6.2+1.251/(oiiitemp/10000) - 0.55*np.log10(oiitemp/10000)-0.014*oiiitemp/10000


    #Calculate O+ and O++ metallicity
    metal = np.log10(10**(oiimet-12)+10**(oiiimet-12))+12
    summet += metal
    summet2 += metal**2

    #CHECK CALC WITH MICHAEL FOR SFR

    #Continous Spectrum lum in AB system from av g band mag
    speclum = 10**(-0.4*(avmg-51.595))
    
    # change units
    speclum = (speclum * 3.0*10**(-10)) / (10**(-8) * 4861)**2
    
    
    # Lum Hydrogen beta
    LHB = EW * 10**(-8) * speclum
    
    #Lun Hydrogen Alpha from LHB
    LHA = LHB * 2.86
    
    #sfr from LHA and metallicty func from Ly 2016
    sfr = LHA * 10**(-41.34-0.39*(metal-12+3.31)+0.127*(metal-12+3.31)**2 )

    sumsfr +=sfr
    sumsfr2 += sfr**2

    
#end of loop thing   
summet = summet / count
summettest = np.sqrt(summet2/ count - (summet)**2)

print summet
print summettest

#umsfr = sumsfr/count
sumsfr2 = np.sqrt(sumsfr2/count - (sumsfr)**2)

print sumsfr
print sumsfr2

plt.errorbar(avm, summet, yerr = summettest, xerr = sumsfr2, linestyle="None") 
