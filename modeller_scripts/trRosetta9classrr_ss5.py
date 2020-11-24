#!/usr/bin/python

##############################################################################################################
#     File Name           :     trRosetta9classrr_ss5.py 					             #
#     Developed By        :     Rahmatullah Roche 					                     #
#     Creation Date       :     [2020-11-23 19:30]					                     #
#     Last Modified       :     [2020-11-23 19:30]                                                           #
#     Description         :     trRosetta predicted distance thresholds corresponding to the 9 distance bins #                 
##############################################################################################################

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
        help = 'name of output contact map')

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
lines = f['dist']
rrList = []
for res_i in range(len(lines)):
	for res_j in range(len(lines[res_i])):
                #line = line.split()

                
                res1 = res_i + 1
                res2 = res_j + 1
                
                if((res2-res1)<=4):
                        continue                
                
		#if(res1 == 1 and res2 == 51):
			#print(lines[res_i][res_j])
                count = 0
                for i in range(1, 25): #2-25 bin for 0-14A threashold, 1st bin is for > 20A
                        th = i
                        count += (lines[res_i][res_j][i]) #Decimal(float(line[i]))
                        if (count >= 0.85):

                                break

                if (th <= 8): #2-9 bin for 6A
                        th = 6

                if (th > 8 and th <= 10): 
                        th = 7


                elif (th > 10 and th <= 12): #2-13 bin for 8A
                        th = 8

                elif (th > 12 and th <= 14): 
                        th = 9

                elif (th > 14 and th <= 16): #14-17 bin for 8-10A
                        th = 10
                elif (th > 16 and th <= 18): 
                        th = 11

                elif (th > 18 and th <= 20): #18-21 for 10-12A
                        th = 12
                elif (th > 20 and th <= 22): 
                        th = 13

                elif (th > 22 and th <= 24): #22-25 for 12-14A
                        th = 14

                if (count < 0.0001):
                        count = 0.0001
                rrList.append([str(res1), str(res2), str(th), str(count)])
                res_j += 1
	res_i += 1
dict1 = {}
for x in rrList:
	dict1[(x[0], x[1], x[2])] = x[3]
sorted_rr = sorted(dict1.items(), key=operator.itemgetter(1))
sorted_rr.reverse()
count = 0
for sr in sorted_rr:
	(i, j, th), p = sr[0], sr[1]
	
	frr.write(str(i))
	frr.write(' ')
	frr.write(str(j))
	frr.write(' ')
	frr.write('0')
	frr.write(' ')
	frr.write(th)
	frr.write(' ')
	frr.write(str(p))
	frr.write('\n')
	count += 1
frr.close()
