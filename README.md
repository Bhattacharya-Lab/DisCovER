# DisCovER

<h2>Distance- and orientation-based Covariational threadER </h2>

### Download and Installation
```sh
$ git clone https://github.com/Bhattacharya-Lab/DisCovER
$ cd DisCovER
$ g++ -O3 -std=c++0x -o DisCovER DisCovER.cpp 
```
The DisCovER executable (named as `DisCovER`) is compiled and tested on x86_64 Linux system. 

## Usage

To see the usage instructions, run `$ ./DisCovER -h`
```
!!! HELP !!!

-----------------------------------------------------------------------------------
 DisCovER : Distance- and orientation-based Covariational threadER (V1.0)
 For comments, please email to bhattacharyad@auburn.edu
-----------------------------------------------------------------------------------

Usage:   ./DisCovER [-option] [argument]

Options:  -T path to template library             	- input, required
          -L .txt file with template IDs 	 	- input, required
          -q query ID   		          	- input, required
          -o path to output directory              	- input, required
          -n number of TOP templates to be selected  	(optional, default 50)
          -c sequence identity cutoff	             	(optional, default 0.30)
```
##### More details on some options:

* `-T` Path to template directory, containing `PDB, MTX, DEP, dssp, and seq` sub-directories.
          `PDB`: Non-redundant PDB structures,
          `MTX`: Sequence Profiles from PSI-BLAST search,
          `DEP`: Depth-dependent structure profiles,
         `dssp`: [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) generated output files,
          `seq`: FASTA sequence (Ref: https://zhanglab.ccmb.med.umich.edu/library/). Please make sure your template directory path ends with "/" (e.g. /home/template_library/) to avoid Segmentation Error.        
* `-L` Path to the file containing template (pdb) IDs (sample `template_list.txt` file in [DisCovER_template_Lib.tar.gz](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example) provided for your reference.)


### Test DisCovER

To run DisCovER, go to the target directory, which contains target sequence file, standard [SPIDER3](https://sparks-lab.org/downloads/) predicted output file, profiles, and inter-residue distance and orientations. Particularly, DisCovER will assume `seq.fasta, seq.spd33, seq.prf, seq.mtx, query.rr, omega.txt, theta.txt, phi.txt` files are prersent in the current directory.

```sh
$ chmod a+x DisCovER 						
$ ulimit -s 419430400                        (To increase the stack limit, avoiding Segmentation Error)
$ cd Example
$ tar -zxvf DisCovER_template_Lib.tar.gz           (Uncompress a toy template library for the test run)
$ cd 6GW7_A/                     (Target directory contains respective input files as mentioned above)
# run DisCovER 
$ ../../DisCovER -T ../template_library/ -L ../template_library/template_list.txt -q 6GW7_A -o ./ -n 50 -c 1.0  

```
The predicted top-one alignment file is named as `top1.fasta`. The alignmnet file for `6GW7_A` is given [here](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example/Output/) for your reference. When you run DisCovER, the output screen looks like [it](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/Example/Output/6GW7_A.log).

### Generating input files from Query sequence

* Generate an alignmnet file ([.aln](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/preprocessing/Target/ex/d1a6qa1.aln) format) of the query sequence. 
* Get trRosetta output file (.npz).
* Go through the steps discussed [here](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/preprocessing/Target).

### Generating template library

* You can either download the template library from https://zhanglab.ccmb.med.umich.edu/library/ (you need `PDB, MTX, DEP`) and generate [`dssp`](https://swift.cmbi.umcn.nl/gv/dssp/) and `seq` by yourself

* or build (`MTX and DEP`) by your own library as discussed in [I-TASSER](https://zhanglab.ccmb.med.umich.edu/I-TASSER/download/).
 
### Building 3D models

* The standard automodel() class of MODELLER is available [here](https://github.com/Bhattacharya-Lab/DisCovER/blob/master/modeller_scripts/build3Dmodel_aln.py).
* The MODELLER with additional restraints from predicted distances, orientations, and secondary structures script is available [here](https://github.com/Bhattacharya-Lab/DisCovER/blob/master/modeller_scripts/build3Dmodel_aln_rsr.py). Go through the usage instruction [here](https://github.com/Bhattacharya-Lab/DisCovER/blob/master/modeller_scripts/readme). The helping script to convert trRosetta output file (.npz) to inter-residue geometries are available [here](https://github.com/Bhattacharya-Lab/DisCovER/tree/master/modeller_scripts).

## Data

- The template library can be found [here](http://sanger.cse.eng.auburn.edu/DisCovER/downloads/DisCovER_template_Lib.tar.gz) 
- The benchmark datasets along with input files can be found [here](/home/project/conThreader/paper3/disThreader/experiment_Sep5/public_data/DisCovER_dataset.tar.gz) 
- DisCovER predicted full-length 3D models can be found [here](/home/project/conThreader/paper3/disThreader/experiment_Sep5/public_data/DisCovER_3D_models.tar.gz) 

## Cite
If you find DisCovER useful, please cite our [paper in bioRxiv](https://www.biorxiv.org/content/10.1101/2020.01.31.923409v1) available as preprint.
