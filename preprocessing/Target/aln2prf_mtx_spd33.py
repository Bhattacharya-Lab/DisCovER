# Input: FASTA Sequence and Alignment(.aln in PSICOV format) files of the target
# Output: .prf, .mtx, .spd33 (Profiles and SPIDER3 predicted output)
# Credit:
#	1. Zheng,W., Wuyun,Q., et al. (2019) Detecting distant-homology protein structures by aligning deep neural-network based contact maps. PLOS Computational Biology 
#	2. Greener,J.G. et al. (2019) Deep learning extends de novo protein modelling cover- age of genomes using iteratively predicted structural constraints. Nat Commun.
#	3. https://github.com/soedinglab/hh-suite


import os,sys
import optparse 

parser = optparse.OptionParser()
parser.add_option('-f', dest='fasta',
        default = '',    # default empty
        help = 'FASTA file containing the amino acid sequence (mandatory)')
parser.add_option('-a', dest='aln',
        default = '',    # default empty
        help = 'Alignmnet (.ali) file containing the query (mandatory)')
parser.add_option('-o', dest='outDir',
        default = '',    # default empty
        help = 'existing output directory path (mandatory)')
(options,args) = parser.parse_args()
fasta=options.fasta
aln=options.aln
outDir=options.outDir

##########################################
# Update paths 	
DisCovER_path = "/home/project/conThreader/DisCovER/DisCovER-master/"			
spider3_script = "/home/project/conThreader/apps/SPIDER3-numpy-server/script/" #which contains "spider3_pred.py"
##########################################

aln2a3m = DisCovER_path+"preprocessing/Target/aln2a3m.sh"
buildProfile=DisCovER_path+"preprocessing/Target/buildquery.pl"
hhmake=DisCovER_path+"preprocessing/Target/hhmake"

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
        os.system(hhmake+" -i seq.a3m -o seq.hhm -v 0 -M a3m")
        os.system("cp protein.pssm "+outDir+tar+"/seq.pssm")
        os.system("cp seq.hhm "+outDir+tar)	
	os.chdir(outDir+tar)
	#run SPIDER3
	os.system(spider3_script+"spider3_pred.py s*")
	os.system("cp temp/protein.mtx "+outDir+tar+"/seq.mtx")
	os.system("cp temp/seq.prf "+outDir+tar)
	os.system("cp temp/seq.fasta "+outDir+tar)
	os.system("rm -rf "+outDir+tar+"/temp")
	os.remove("seq.hhm")
	os.remove("seq.pssm")
run(aln)
