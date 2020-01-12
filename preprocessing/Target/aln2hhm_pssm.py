#deep alignment(.aln) to hhm and pssm to run spider3

import os,sys
import glob
from itertools import product
from subprocess import call
import pprint
from multiprocessing import Pool

pp = pprint.PrettyPrinter(indent=4)

inDir = "/home/project/conThreader/paper3/disThreader/data/test480_CATHER/allSeq/"

mtx_dir="/home/project/conThreader/paper3/disThreader/data/test480_CATHER/deepMSA/"
aln2a3m = "/home/szb0134/apps/cEthreader/CEthreader/lib/deepMSA/scripts/aln2a3m.sh"
buildProfile="/home/project/conThreader/scripts/c++_implementation/freeze_version/conTh_v3/lomets600/buildquery.pl"
hhmake="/home/project/conThreader/apps/hhsuite/hh-suite/bin/hhmake"

outDir="/home/project/conThreader/paper3/disThreader/data/test480_CATHER/prf_mtx_spd33/"

def run(file):
        tar=(file.split("/")[-1]).split(".")[0]
        if os.path.exists(outDir+tar):
                return
        os.system("mkdir "+outDir+tar)
        os.system("mkdir "+outDir+tar+"/temp")
	#os.system("cp "+inDir+tar+".fasta seq.fasta")
        os.chdir(outDir+tar+"/temp")
        os.system("cp "+inDir+tar+".seq seq.fasta")
        #aln --> a3m
        os.system("cp "+file+" seq.aln")
        os.system(aln2a3m+" "+file+" seq.a3m")
        os.system(buildProfile)
        #aln-->mtx
        os.system("perl /home/szb0134/apps/deepMSA/scripts/a3m2mtx.pl "+file+" seq.mtx -aln")
        os.system(hhmake+" -i seq.a3m -o seq.hhm -v 0 -M a3m")
        #os.system("cp seq.prf "+outDir+tar)
        #os.system("cp seq.mtx "+outDir+tar)
        os.system("cp protein.pssm "+outDir+tar+"/seq.pssm")
        os.system("cp seq.hhm "+outDir+tar)
	#os.system("cp seq.prf "+outDir+tar)
	#os.system("cp seq.mtx "+outDir+tar)
	os.chdir(outDir+tar)
        #os.system("rm -rf "+outDir+tar+"/temp")
	#run SPIDER3
	os.system("/home/project/conThreader/apps/SPIDER3-numpy-server/script/spider3_pred.py s*")
	os.system("cp temp/seq.mtx "+outDir+tar)
	os.system("cp temp/seq.prf "+outDir+tar)
	os.system("cp temp/seq.fasta "+outDir+tar)
	os.system("rm -rf "+outDir+tar+"/temp")


p = Pool(1)
p.map(run, glob.glob(mtx_dir+"*.aln"))
