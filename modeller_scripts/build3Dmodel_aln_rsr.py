#############################################################################################################################################
#     File Name           :     build3Dmodel_aln_rsr.py 										    #
#     Developed By        :     Rahmatullah Roche 											    #
#     Creation Date       :     [2020-11-23 18:12]											    #
#     Last Modified       :     [2020-11-23 18:12] 											    #
#     Description         :     model building with additional restraints from predicted distances, orientations, and secondary structures  #
#############################################################################################################################################


from modeller import *
from modeller.automodel import *
import os,sys,optparse
import math
parser=optparse.OptionParser()
parser.add_option('-f', dest='alignment_file', 
        default= '',    #default empty!' 
        help= 'alignment file')
parser.add_option('--ss', dest='ss', 
        default= '',    #default empty!' 
        help= 'secondary structure file') 
parser.add_option('--rr', dest='rr', 
        default= '',    #default empty!' 
        help= 'contact map in standard rr format') 
parser.add_option('--target', dest='target', 
        default= '',    #default empty!' 
        help= 'name of target') 
parser.add_option('--template', dest='template', 
        default= '',    #default empty! 
        help= 'name of template') 
parser.add_option('--ctype', dest='ctype', 
        default= 'ca',    #default empty! 
        help= 'ca or cb. default ca')
parser.add_option('--phi', dest='phi',
        default= '',    #default empty!
        help= 'phi orientation')
parser.add_option('--theta', dest='theta',
        default= '',    #default empty!
        help= 'theta orientation')
parser.add_option('--omega', dest='omega',
        default= '',    #default empty!
        help= 'omega orientation')

 
(options,args) = parser.parse_args()
alignment_file = options.alignment_file
target = options.target
template = options.template
ss = options.ss
rr = options.rr
ctype = options.ctype
phi = options.phi
theta = options.theta
omega = options.omega

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(automodel):
       def special_restraints(self, aln): 
               rsr = self.restraints 
               at = self.atoms 
               cutoff = 0.85
               fss = open(ss, 'r') 
               self.secondary = fss.read() 
               self.secondary = self.secondary.split() #secondary structure is now in self.secondary[0] 
               fss.close() 
               print (self.secondary[0]) 
               frr = open(rr, 'r') 
               rrlines = frr.readlines() 
               faln = open(alignment_file, 'r') 
               alnlines = faln.readlines() 
               n = len(rrlines[0].strip())  
               pair_u = {} 
               pair_l = {} 
               self.CM = [[0 for x in range(n)] for y in range(n)] 
               max_upper = 0
               for line in rrlines[1:]: 
                       line = line.strip().split() 
                       self.CM[int(line[0]) - 1][int(line[1]) - 1] = self.CM[int(line[1]) - 1][int(line[0]) - 1] = 1   
                       if (float(line[4]) <= cutoff): 
                               self.CM[int(line[0]) - 1][int(line[1]) - 1] = self.CM[int(line[1]) - 1][int(line[0]) - 1] = -1 
                       pair_l[(int(line[0]) - 1, int(line[1]) - 1)] = pair_l[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[2]) 
                       pair_u[(int(line[0]) - 1, int(line[1]) - 1)] = pair_u[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[3]) 
	               if (max_upper < float(line[3])):
	                       max_upper = float(line[3]) 
               #mdl = template + '.B99990001.pdb' 
               #if os.path.exists(mdl): 
               #        print ('changing') 
               #        os.system('mv ' + template + '.B99990001.pdb ' + template + '.pdb') 
               fpdb = open(template + '.pdb', 'r') 
               fpdblines = fpdb.readlines() 
               pos = [] 
               if(ctype == 'cb' or ctype == 'CB'): 
                        ct = "CB" 
               else: 
                        ct = "CA" 
               for line in fpdblines: 
                        if(line[:4] == "ATOM" and line[12:16].strip() == ct): 
                                x = float(line[30:38].strip()) 
                                y = float(line[38:46].strip()) 
                                z = float(line[46:54].strip()) 
                                pos.append([x, y, z]) 
               for i in range(n): 
                       for j in range(i+6, n): 
                               if (abs(i - j) < 6): 
                                       continue 
                               #d = math.sqrt((pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2 + (pos[i][2] - pos[j][2]) ** 2) 
                               if (self.CM[i][j] == 1): # and d > pair_u[(i,j)]): 
                                       firstres = ct + ':' + str(i+1) 
                                       secondres = ct + ':' + str(j+1) 
                                       if(rrlines[0][i] == 'G'): 
                                               firstres = 'CA:' + str(i+1) 
                                       if(rrlines[0][j] == 'G'): 
                                               secondres = 'CA:' + str(i+1) 
                                       rsr.add(forms.lower_bound(group=physical.xy_distance, 
                                              feature=features.distance(at[firstres], 
                                                                at[secondres]), 
                                              mean=3.6, stdev=.1)) 
                                       rsr.add(forms.upper_bound(group=physical.xy_distance, 
                                                       feature=features.distance(at[firstres], 
                                                                at[secondres]), 
                                                      mean=pair_u[(i,j)], stdev=0.1)) 
                               elif (self.CM[i][j] == 0): # and d < 14): 
                                        firstres = ct + ':' + str(i+1) 
                                        secondres = ct + ':' + str(j+1) 
                                        if(rrlines[0][i] == 'G'): 
                                               firstres = 'CA:' + str(i+1) 
                                        if(rrlines[0][j] == 'G'): 
                                               secondres = 'CA:' + str(i+1) 
                                        rsr.add(forms.lower_bound(group=physical.xy_distance, 
                                               feature=features.distance(at[firstres], 
                                                                 at[secondres]), 
                                               mean=max_upper, stdev=0.1)) 
               fphi = open(phi, 'r')
               ftheta = open(theta, 'r')
               fomega = open(omega, 'r')
               philines = fphi.readlines()
               thetalines = ftheta.readlines()
               omegalines = fomega.readlines()
               
               for line in philines[1:]:
                       line = line.strip().split()
                       res1 = line[0]
                       res2 = line[1]
                       prob = float(line[4])
                       lb = float(line[2])
                       ub = float(line[3])
                       l_u_mean = (lb + ub) / 2
                       atm1 = 'CA:' + res1
                       atm2 = 'CB:' + res1
                       atm3 = 'CB:' + res2
                       if(rrlines[0][int(res1)-1] == 'G' or rrlines[0][int(res2)-1] == 'G'):
                               continue
                       if(abs(int(res1) - int(res2)) <= 5):
                               continue
		       if(prob <= 0.85):
                               break
                       #lb = max(0, lb-(30*3.1416/180))
                       #ub = min(180, ub + (30*3.1416/180))
                       rsr.add(forms.upper_bound(group=physical.angle,
                                               feature=features.angle(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3]),
                                               mean=ub, stdev=0.1))
                       rsr.add(forms.lower_bound(group=physical.angle,
                                               feature=features.angle(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3]),
                                               mean=lb, stdev=0.1))
               
               
               for line in thetalines[1:]:
                       line = line.strip().split()
                       res1 = line[0]
                       res2 = line[1]
                       prob = float(line[4])
                       lb = float(line[2])

                       ub = float(line[3])
                       l_u_mean = (lb + ub) / 2

                       atm1 = 'N:' + res1
                       atm2 = 'CA:' + res1
                       atm3 = 'CB:' + res1
                       atm4 = 'CB:' + res2
                       if(rrlines[0][int(res1)-1] == 'G' or rrlines[0][int(res2)-1] == 'G'):
                               continue
                       if(abs(int(res1) - int(res2)) <= 5):
                               continue

                       if(prob <= 0.85):
                               break
                       #lb = max(-180, lb-(15*3.1416/180))
                       #ub = min(180, ub + (15*3.1416/180))

                       rsr.add(forms.upper_bound(group=physical.dihedral,
                                               feature=features.dihedral(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3],
                                                                 at[atm4]),
                                               mean=ub, stdev=0.1))
                       rsr.add(forms.lower_bound(group=physical.dihedral,
                                               feature=features.dihedral(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3],
                                                                 at[atm4]),
                                               mean=lb, stdev=0.1))

               
               for line in omegalines[1:]:
                       line = line.strip().split()
                       res1 = line[0]
                       res2 = line[1]
                       prob = float(line[4])
                       lb = float(line[2])
                       ub = float(line[3])
                       l_u_mean = (lb + ub) / 2
                       atm1 = 'CA:' + res1
                       atm2 = 'CB:' + res1
                       atm3 = 'CB:' + res2
                       atm4 = 'CA:' + res2
                       if(rrlines[0][int(res1)-1] == 'G' or rrlines[0][int(res2)-1] == 'G'):
                               continue
                       if(abs(int(res1) - int(res2)) <= 5):
                               continue

                       if(prob <= 0.85):
                               break
                       #lb = max(-180, lb-(15*3.1416/180))
                       #ub = min(180, ub + (15*3.1416/180))

                       rsr.add(forms.gauss(group=physical.dihedral,
                                               feature=features.dihedral(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3],
                                                                 at[atm4]),
                                               mean=ub, stdev=0.1))
                       rsr.add(forms.gauss(group=physical.dihedral,
                                               feature=features.dihedral(at[atm1],
                                                                 at[atm2],
                                                                 at[atm3],
                                                                 at[atm4]),
                                               mean=lb, stdev=0.1))

              
               self.helices = [] 
               self.sheets = [] 
               i = 1 
               while(i < len(self.secondary[0])): 
                       if self.secondary[0][i] == 'H': 
                               j = i 
                               while (self.secondary[0][j] == 'H' and j < len(self.secondary[0]) - 1): 
                                       j = j + 1 
                               if (j - i >= 3): 
                                       self.helices.append((i, j)) 
                                       i = j 
                       i = i + 1 
               i = 1 
               while(i < len(self.secondary[0])): 
                       if self.secondary[0][i] == 'E': 
                               j = i 
                               while (self.secondary[0][j] == 'E' and j < len(self.secondary[0]) - 1): 
                                       j = j + 1 
                               if (j - i >= 3): 
                                       self.sheets.append((i, j-1)) 
                                       i = j 
                       i = i + 1 
               for (i,j) in self.helices: 
                       start = str(i) + ':' 
                       end = str(j) + ':' 
                       rsr.add(secondary_structure.alpha(self.residue_range(start, end))) 
                       print (start) 
                       print (end) 
               for (i,j) in self.sheets: 
                       start = str(i) + ':' 
                       end = str(j) + ':' 
                       rsr.add(secondary_structure.strand(self.residue_range(start, end))) 
                       print (start) 
                       print (end) 
log.verbose() 
env = environ() #consistant models each time 
a = MyModel(env,alnfile=alignment_file, 
                knowns=template,sequence=target, 
                assess_methods=(assess.DOPE)) 
a.starting_model= 1 
a.ending_model = 1 
a.make() 
#mdl = template + '.B99990001.pdb' 
#if not os.path.exists(mdl): 
#       os.system('mv ' + template + '.pdb ' + template + '.B99990001.pdb') 
