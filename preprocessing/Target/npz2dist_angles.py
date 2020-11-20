import os,sys
import numpy as np
import optparse 

#i/p: .npz file from trRosetta
#o/p: query.rr, omega.txt, theta.txt, phi.txt
#Credit:
#	1. https://yanglab.nankai.edu.cn/trRosetta/


parser = optparse.OptionParser()
parser.add_option('-i', dest='npz_file',
        default = '',    # default empty
        help = '.npz file from trRosetta (mandatory)')
parser.add_option('-o', dest='outDir',
        default = '',    # default empty
        help = 'existing output directory path (mandatory)')
(options,args) = parser.parse_args()
npz_file=options.npz_file
outDir=options.outDir

def get_angles_all_bin(npz_file, tarDir):
        NPZ=npz_file
        npz = np.load(NPZ)
        omega,theta,phi = npz['omega'],npz['theta'],npz['phi']
        nres=int(omega.shape[0])
        f=open(tarDir+"/omega.txt","w")
        for i in range(0, nres):
                for j in range(i+1, nres):
                        Praw = omega[i][j]
                        f.write(str(i+1)+" "+str(j+1))
                        cnt=1
                        for ii in range(1,25):
                                f.write(" "+str(Praw[ii]))

                        f.write("\n")
                        f.write(str(j+1)+" "+str(i+1))
                        for ii in range(1,25):
                                f.write(" "+str(Praw[ii]))
                        f.write("\n")
        f.close()
        nres=int(theta.shape[0])
        f=open(tarDir+"/theta.txt","w")
        for i in range(0, nres):
                for j in range(0, nres):
                        Praw = theta[i][j]
                        f.write(str(i+1)+" "+str(j+1))
                        cnt=1
                        for ii in range(1,25):
                                f.write(" "+str(Praw[ii]))
                        f.write("\n")
        f.close()
        nres=int(phi.shape[0])
        f=open(tarDir+"/phi.txt","w")
        for i in range(0, nres):
                for j in range(0, nres):
                        Praw = phi[i][j]
                        f.write(str(i+1)+" "+str(j+1))
                        cnt=1
                        for ii in range(1,13):
                                f.write(" "+str(Praw[ii]))
                        f.write("\n")
        f.close()

def run():
	get_angles_all_bin(npz_file, outDir)
	#dist
	out=open(outDir+"/query.rr","w")
	npz = np.load(npz_file)
        dat=npz['dist']
        nres=int(dat.shape[0])
        for i in range(0, nres):
                for j in range(i+5, nres):
			Praw = dat[i][j]
                        first_bin=5
                        pcont=0
                        out.write(str(i+1)+" "+str(j+1))
			for ii in range(1,37):
                                pcont += Praw[ii]
                                if (ii) > first_bin:
                                        out.write(" "+str(round(pcont,4)))
                        out.write("\n")
	out.close()

run()
