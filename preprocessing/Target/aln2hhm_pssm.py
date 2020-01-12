# Input: FASTA Sequence and Alignment(.aln in PSICOV format) files of the target
# Output: .prf, .mtx
# Ref:
#	1.
#	2.
#	3.

import os,sys
##########################################
# Update paths 				
fasta=inDir = "" #Query sequence 
aln=mtx_dir= "" #Query alignment (.aln)
outDir= "" #Output directory
##########################################

aln2a3m = "aln2a3m.sh"
buildProfile="buildquery.pl"
a3m2mtx="a3m2mtx.pl"
hhmake="hhmake"



def run(file):
        tar=(file.split("/")[-1]).split(".")[0]
        os.system("mkdir "+outDir+tar)
        os.system("mkdir "+outDir+tar+"/temp")
        os.chdir(outDir+tar+"/temp")
        os.system("cp "+fasta+" seq.fasta")
        #aln --> a3m
        os.system("cp "+file+" seq.aln")
        os.system(aln2a3m+" "+file+" seq.a3m")
        os.system(buildProfile)
        #aln-->mtx
        os.system("perl "+a3m2mtx+file+" seq.mtx -aln")
        os.system(hhmake+" -i seq.a3m -o seq.hhm -v 0 -M a3m")
        #os.system("cp protein.pssm "+outDir+tar+"/seq.pssm")
        #os.system("cp seq.hhm "+outDir+tar)	
	os.chdir(outDir+tar)
	#run SPIDER3
	#os.system("/home/project/conThreader/apps/SPIDER3-numpy-server/script/spider3_pred.py s*")
	os.system("cp temp/seq.mtx "+outDir+tar)
	os.system("cp temp/seq.prf "+outDir+tar)
	os.system("cp temp/seq.fasta "+outDir+tar)
	os.system("rm -rf "+outDir+tar+"/temp")

run(aln)
