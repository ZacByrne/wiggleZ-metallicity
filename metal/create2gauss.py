# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:39:06 2017

@author: uqzbyrne

A python alternative to create2gauss. reads starlight files
"""
import numpy as np
#import os
import matplotlib, matplotlib.pyplot as plt
import linecache

import sys
ncol =  sys.argv[1]

#Get errors from base spectra
basefile = 'wig02bin' + str(ncol) + '.txt'

basedata = np.genfromtxt(basefile)
z = [basedata[:,0], basedata[:,2]];

#test1= data[1][np.where(data[0] == 5)[0][0]]

er4363 = z[1][(np.where(z[0] == 4300)[0][0]):(np.where(z[0] == 4500)[0][0])+1]

erhb = z[1][(np.where(z[0] == 4800)[0][0]):(np.where(z[0] == 5100)[0][0])+1]

er3727 = z[1][(np.where(z[0] == 3500)[0][0]):(np.where(z[0] == 4000)[0][0])+1]


#Get cont level for eqw for SFR 

#hbfile = 'wig02hb' + str(ncol) + '.BN'

hbdata = np.genfromtxt(basefile)
hbdata = [hbdata[:,0] , hbdata[:,1], hbdata[:,2]]

average1 = np.mean(hbdata[1][(np.where(hbdata[0] == 4800)[0][0]):(np.where(hbdata[0] == 4835)[0][0])+1])
average2 = np.mean(hbdata[1][(np.where(hbdata[0] == 4885)[0][0]):(np.where(hbdata[0] == 4925)[0][0])+1])

contlevel = (average1+average2)/2;

#Get norm constants
linenum = 16;



s = linecache.getline('wig02fit' + str(ncol) + '.BN', linenum-1)
normfull = float(s.split()[0])
s = linecache.getline('wig02hb' + str(ncol) + '.BN', linenum-1)
normhb = float(s.split()[0])
s = linecache.getline( str(ncol) + '4363.BN', linenum-1)
norm4363 = float(s.split()[0])


#read files, get fitted spectra, adjust with norm and cont level then write to file
outfile = 'wig02fitted' + str(ncol) + '.txt'

o2file = 'wig02fit' + str(ncol) + '.BN'

o2data = np.genfromtxt(o2file, skip_header = 978)
o2data = [o2data[:,0] , o2data[:,1], o2data[:,2]]

o2data = np.transpose(np.array((o2data[0],(o2data[1]-o2data[2])*normfull + contlevel, er3727)))

o3file = str(ncol) + '4363.BN'

o3data = np.genfromtxt(o3file, skip_header = 978)
o3data = [o3data[:,0] , o3data[:,1], o3data[:,2]]

o3data = np.transpose(np.array((o3data[0],(o3data[1]-o3data[2])*norm4363 + contlevel, er4363)))

hbfile = 'wig02hb' + str(ncol) + '.BN'

hbdata = np.genfromtxt(hbfile, skip_header = 978)
hbdata = [hbdata[:,0] , hbdata[:,1], hbdata[:,2]]
hbdata = np.transpose(np.array((hbdata[0],(hbdata[1]-hbdata[2])*normhb + contlevel, erhb)))

data = np.append(o2data,o3data, axis=0)
