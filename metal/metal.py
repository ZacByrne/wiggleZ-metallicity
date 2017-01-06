import numpy as np
#import os
#import matplotlib, matplotlib.pyplot as plt
#import collections
import math as math

#constant E(b-v)_gamma. Might later include calculations?
EBV_gamma = 0.61113184

#adjust line flux for reddening as flux * 10^(-0.4*E(b-v)*k(lambda))  as per ly_2016
def adjustflux(flux, kred):
	flux = flux * 10**(0.4*EBV_gamma*kred)
	return flux



#load values
#import args from system? (Google)
# import line flux values from file into matrix.. NEED TO SEE OUTPUT FILE... CHANGE TO USECOL()
col, oiiiw, oiiiwe, oii, oiie, oiii, oiiie, hb, hbe = np.loadtxt(".../data/fluxes.txt", dtype = {'names': ('column name', '4363', '4363 error ', '3728', '3278 error', '4959', '4959 error', 'hbeta', 'hbeta error'), 'formats': ('s4', 'f8','f8','f8','f8','f8','f8','f8','f8')})


# calculate k values? Or simply import from matlab program?
k4363 = 4.1482
k3728 = 5.1632
k4959 = 3.5163
k4861 = 3.6092

# import mass, sfr, redshift? Or calc?np.random.normal(size = 4)

#Loop stats
count = 10000

summet = 0
summet2 = 0

#Need to show something else? av e temp? Av oII met?


#assume each variable is an array of variable named column
for x in range(0,count):
    
    #Uncertainity calculations:
    oiiiw = oiiiw + np.random.normal(size = (len(oiiiw))) * oiiiwe
    oiii = oiii + np.random.normal(size = (len(oiiiw))) * oiiie
    oii = oii + np.random.normal(size = (len(oiiiw))) * oiie
    hb = hb + np.random.normal(size = (len(oiiiw))) * hbe
    
    # Adjust line fluxes for reddening as per Ly 2016
    oiiiw = adjustflux(oiiiw, k4363)
    oii = adjustflux(oii, k3728)
    oiii = adjustflux(oiii, k4959)
    hb = adjustflux(hb, k4861)

    # calc OIII temp as per Nichollis et al. 2014 
    oiiitemp = 13205 * (-1*math.log10((oiiiw / (oiii*4)))-0.92506)**(-0.98062)

    #Calc OII temp as per Andrews and Martini 2013
    oiitemp = 0.7 * oiiitemp + 0.17 * 10000

    #Calculate OII and OIII metallicity
    oiimet = math.log10(oii/hb)+5.961+1.676/oiitemp - 0.4*math.log10(oiitemp)-0.034*oiitemp
    oiiimet = math.log10(4*oiii/hb)+6.2+1.251/oiiitemp - 0.55*math.log10(oiitemp)-0.014*oiiitemp

    #Calculate O+ and O++ metallicity
    metal = math.log10(10**(oiimet-12)+10**(oiiimet-12))+12
    
    summet += metal
    summet += metal**2
    
#end of loop thing   
summet = summet / count
summet2 = math.sqrt(summet2/ count - (summet)**2)
