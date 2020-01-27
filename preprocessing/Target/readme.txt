Prequisites:
	Install hh-suite (https://github.com/soedinglab/hh-suite)
	Install SPIDER3 (https://sparks-lab.org/downloads/)

1. Modify lines 29-30 of 'aln2prf_mtx_spd33.py'
2. Modify line 13 of 'buildquery.pl'
3. Set HHLIB path to the location where hh-suite is installed. 
	For ex: "$ export HHLIB=/home/apps/hh-suite"
3. Modify line 12 of 'a3m2psiblast_v2.pl'
4. Run 'aln2prf_mtx_spd33.py'. The script takes a FASTA sequence, an alignment file (.aln), and an output directory. Please provide the FULL path or absolute path for each argument as follows:
	"$ python aln2prf_mtx_spd33.py -f <full/path/of/FASTA/sequence>  -a <full/path/of/the/alignmnet/(.aln)/file>  -o <full/path/of/the/output/directory>"
	For example: 
	"$ python aln2prf_mtx_spd33.py -f /home/DisCovER-master/preprocessing/Target/ex/d4pv4a1.seq -a /home/DisCovER-master/preprocessing/Target/ex/d4pv4a1.aln -o /home/DisCovER-master/preprocessing/Target/ex/"
