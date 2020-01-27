# DisCovER

<h2>Distance-based Covariational threadER </h2>

### Download and Installation
```sh
$ git clone https://github.com/Bhattacharya-Lab/DisCovER
$ cd DisCovER
```
The DisCovER executable (named as `DisCovER`) is compiled and tested on x86_64 redhat linux system. 

## Usage

To see the usage instructions, run `$ ./DisCovER -h`
```
!!! HELP !!!

-----------------------------------------------------------------------------------
 DisCovER : distance-based covariational threadER (V1.0)
 For comments, please email to bhattacharyad@auburn.edu
-----------------------------------------------------------------------------------

Usage:   ./DisCovER [-option] [argument]

Options:  -T path to template library             		- input, required
          -L .txt file with template IDs 	 	    	- input, required
          -q query ID   		          	        	- input, required
          -o path to output directory              		- input, required
          -d path to query distance file 	  	    	- input, required 
          -m path to MODELLER       		     		(optional, but recommended)
          -n number of TOP templates to be selected  	(optional, default 50)
          -c sequence identity cutoff	             	(optional, default 0.30)


```
##### More details on some options:

* `-T` Path to template directory, containing `PDB, MTX, DEP, dssp, and seq` sub-directories.
          `PDB`: Non-redundant PDB structures,
          `MTX`: Sequence Profiles from PSI-BLAST search,
          `DEP`: Depth-dependent structure profiles,
         `dssp`: [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) generated output files,
          `seq`: FASTA sequence (Ref: https://zhanglab.ccmb.med.umich.edu/library/)         
* `-L` Path to the file containing template (pdb) IDs (sample `template_list.txt` file in [DisCovER_template_Lib.tar.gz](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example) provided for your reference.)
* `-d` Path to distance histogram (distogram) `rawdistpred.current` generated by [DMPfold](https://github.com/psipred/DMPfold)
* `-m` Path to [MODELLER's](https://salilab.org/modeller) `modpy.sh` script, e.g. `/home/aubszb/modpy.sh`


### Test DisCovER

To run DisCovER, go to the target directory, which contains target sequence file, standard [SPIDER3](https://sparks-lab.org/downloads/) predicted output file, and profiles. Particularly, DisCovER will assume `seq.fasta, seq.spd33, seq.prf, seq.mtx` files are prersent in the current directory.

```sh
$ chmod a+x DisCovER 						
$ ulimit -s 419430400 ./DisCovER             (To increase the stack limit, avoiding Segmentation Error)
$ cd Example
$ tar -zxvf DisCovER_template_Lib.tar.gz           (Uncompress a toy template library for the test run)
$ cd d1a9xb1/                     (Target directory contains respective input files as mentioned above)
#To run DisCovER without MODELLER
$ ../../DisCovER -T ../template_library/ -L ../template_library/template_list.txt -q d1a9xb1 -o ./ -d rawdistpred.current 
#To run DisCovER with MODELLER  
$ ../../DisCovER -T ../template_library/ -L ../template_library/template_list.txt -q d1a9xb1 -o ./ -d rawdistpred.current -m /home/aubszb/modpy.sh
```
The predicted 3D model is named as `seq_model1.pdb` and the corresponding alignment file is named as `top1.fasta`. If DisCovER is run without MODELLER (-m), then only `top1.fasta` file will be generated. Both the predicted 3D model and the alignmnet file for `d1a9xb1` are given [here](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example/Output/) for your reference. When you run DisCovER, the output screen looks like [it](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example/Output/d1a9xb1.log).

### Generating input files from Query sequence

* Generate an alignmnet file ([.aln](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/preprocessing/Target/ex/d4pv4a1.aln) format) of the query sequence. 
* Go through the steps discussed [here](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/preprocessing/Target).

### Generating template library

* You can either download the template library from https://zhanglab.ccmb.med.umich.edu/library/ (you need `PDB, MTX, DEP`) and generate [`dssp`](https://swift.cmbi.umcn.nl/gv/dssp/) and `seq` by yourself

* or build (`MTX and DEP`) by your own library as discussed in [I-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/download/).
 

## Data

- The template library can be found [here](http://sanger.cse.eng.auburn.edu/DisCovER/downloads/DisCovER_template_Lib.tar.gz) 
- The benhcmark datasets along with input files can be found [here](http://sanger.cse.eng.auburn.edu/DisCovER/downloads/DisCovER_dataset.tar.gz) 
- DisCovER predicted full-length 3D models can be found [here](http://sanger.cse.eng.auburn.edu/DisCovER/downloads/DisCovER_3D_models.tar.gz) 

## Cite
If you find DisCovER useful, please cite our paper:
