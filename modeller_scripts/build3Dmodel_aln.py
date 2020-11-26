#############################################################################################
#     File Name           :     build3Dmodel_aln.py 					    #
#     Developed By        :     Rahmatullah Roche 					    #
#     Creation Date       :     [2020-11-23 18:56]					    #
#     Last Modified       :     [2020-11-23 18:56]                                          #
#     Description         :     model building using standard automodel()                   #
#     Credit		  :     https://salilab.org/modeller/				    # 
#############################################################################################


from modeller import *
from modeller.automodel import *
import os,sys,optparse

parser=optparse.OptionParser()
parser.add_option('-f', dest='alignment_file',
        default= '',    #default empty!
        help= 'alignment file in MODELLER format')
parser.add_option('--target', dest='target',
        default= '',    #default empty!
        help= 'name of target')
parser.add_option('--template', dest='template',
        default= '',    #default empty!
        help= 'name of template')
(options,args) = parser.parse_args()

alignment_file = options.alignment_file
target = options.target
template = options.template

#log.verbose()
env = environ()
a = automodel(env,alnfile=alignment_file,
                knowns=template,sequence=target,
                assess_methods=(assess.DOPE))
		
#a.very_fast() # prepare for extremely fast optimization

a.starting_model= 1
a.ending_model = 1
a.make()
