import os,sys


pwd="/home/project/conThreader/paper3/disThreader/experiment_Sep5/scripts/casp13/adv_modeller/"
to_rr=pwd+"trRosetta9classrr_ss5.py"
to_omega=pwd+"trRosetta_omega_upper_lower.py"
to_phi=pwd+"trRosetta_phi_upper_lower.py"
to_theta=pwd+"trRosetta_theta_upper_lower.py"

modeller_script=pwd+"advanced_modeller_95_orientation_u_l_v14.py"
trRosetta="/home/project/conThreader/paper3/disThreader/experiment_Sep5/data/casp13/trRosetta_native_domain/"
outDir="/home/project/conThreader/paper3/disThreader/experiment_Sep5/scripts/casp13/trRosetta_rr_angles_for_adv_modeller/"
SEQ="/home/project/conThreader/paper3/disThreader/data/casp13/domain_based_exp/seq/"

for tar in os.listdir(trRosetta):
	if not tar.endswith(".npz"):
		continue
	target=tar[:-4]
	print target
	#print "python "+to_rr+" -d "+trRosetta+tar+" -a "+SEQ+target+".seq  -r "+outDir+target+".trRosetta9class.rr "
	#break
	os.system("python "+to_rr+" -d "+trRosetta+tar+" -a "+SEQ+target+".seq  -r "+outDir+target+".trRosetta9class.rr ")
	os.system("python "+to_omega+" -d "+trRosetta+tar+" -a "+SEQ+target+".seq  -r "+outDir+target+".omega")
	os.system("python "+to_theta+" -d "+trRosetta+tar+" -a "+SEQ+target+".seq  -r "+outDir+target+".theta")
	os.system("python "+to_phi+" -d "+trRosetta+tar+" -a "+SEQ+target+".seq  -r "+outDir+target+".phi")

	#break
