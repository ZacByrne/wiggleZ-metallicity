# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:39:06 2017

@author: uqzbyrne

A python alternative to create2gauss. reads starlight files
"""
import numpy as np
#import os
import matplotlib, matplotlib.pyplot as plt

import sys
ncol =  sys.argv[0]

basefile = 'wig02bin' + str(ncol) + '.txt'

basedata = np.genfromtxt(basefile)
z = [basedata[:,0], basedata[:,2]];

test1= data[1][numpy.where(data[0] == 5)[0][0]]

er4363 = z((find(z == 4300)):(find(z==4500)),:);

erhb = z((find(z == 4800)):(find(z==5100)),:);

er3728 = z(1:(find(z==4000)),:);