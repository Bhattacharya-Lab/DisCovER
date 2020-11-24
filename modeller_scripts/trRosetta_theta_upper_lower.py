#!/usr/bin/python

#############################################################################################################
#     File Name           :     trRosetta_theta_upper_lower.py			                            #
#     Developed By        :     Rahmatullah Roche 					                    #
#     Creation Date       :     [2020-11-23 20:12]					                    #
#     Last Modified       :     [2020-11-23 20:12]                                                          #
#     Description         :     trRosetta predicted theta angles derived from the highest likelihood bins   #
#############################################################################################################

import operator
import os
import math
import sys
import optparse    # for option sorting
import random
from decimal import *
import numpy as np
getcontext().prec = 4
parser = optparse.OptionParser()
parser.add_option('-d', dest='dist',
        default = '',    # default empty
        help = 'name of trRosetta predicted file (.npz)')
parser.add_option('-a', dest='fasta',
        default = '',    # default empty
        help = 'name of fasta file')


parser.add_option('-r', dest='rr',
        default = 'out',    # default out
        help = 'name of output')

(options,args) = parser.parse_args()
dist = options.dist
fasta = options.fasta
rr = options.rr
try:
        ffasta = open(fasta, 'r')
        flines = ffasta.readlines()
except IOError:
        sys.exit()
        print ('skipping..')

frr = open(rr, 'w')
frr.write(flines[1].strip())
frr.write('\n')



#try:
#f = open(dist, 'r')
f = np.load(dist)
lines = f['theta']
rrList = []
for res_i in range(len(lines)):
	for res_j in range(len(lines[res_i])):
                #line = line.split()

                
                res1 = res_i + 1
                res2 = res_j + 1
                
                #if((res2-res1)<=5):
                #        continue                
                
		#if(res1 == 1 and res2 == 51):
			#print(lines[res_i][res_j])
                maxp = 0
                th = 24
                for i in range(1, 25): #2-25 bin for 15 degree angle interval
                        #th = i
                        count = float(lines[res_i][res_j][i]) #Decimal(float(line[i]))
                        if (maxp < count):

                                th = i
                                maxp = count
                if (th <= 13):
                        up = (13 - th) * (-15)
                        lo = up - 15
                        upr = up * 3.1416 / 180.0
                        lor = lo * 3.1416 / 180.0 
                else:
                        up = (th - 13) * 15
                        lo = up  - 15
                        upr = up * 3.1416 / 180.0
                        lor = lo * 3.1416 / 180.0
			
                if (maxp < 0.0001):
                        maxp = 0.0001
                rrList.append([str(res1), str(res2), str(lor), str(upr), str(maxp)])
                res_j += 1
	res_i += 1
dict1 = {}
for x in rrList:
	dict1[(x[0], x[1], x[2], x[3])] = x[4]
sorted_rr = sorted(dict1.items(), key=operator.itemgetter(1))
sorted_rr.reverse()
count = 0
for sr in sorted_rr:
	(i, j, lor, upr), p = sr[0], sr[1]
	
	frr.write(str(i))
	frr.write(' ')
	frr.write(str(j))
	frr.write(' ')
	frr.write(str(lor))
	frr.write(' ')
	frr.write(str(upr))
	frr.write(' ')
	frr.write(str(p))
	frr.write('\n')
	count += 1
frr.close()
