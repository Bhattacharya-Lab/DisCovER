# DisCovER

Distance-based Covariational Protein Threading for Low-homology Fold Recognition.

### Download and Installation
```sh
$ git clone https://github.com/Bhattacharya-Lab/DisCovER
$ cd DisCovER
```
The DisCovER executable is compiled and tested on x86_64 redhat linux system. 

## Usage

To see the usage instructions, run `$ ./discover -h`
```
!!! HELP !!!

-----------------------------------------------------------------------------------
                        DisCovER
 Distance-based Covariational Protein Threading  V2020.Jan.09
 For comments, please email to bhattacharyad@auburn.edu
-----------------------------------------------------------------------------------

Usage:   ./discover [-option] [argument]

Options:  -T path to template library                   - input, required
          -L list.txt with template IDs                 - input, required
          -q query ID                                   - input, required
          -o path to output director                    - input, required
          -d path to query distance file                - input, required
          -m path to MODELLER                           (optional, but recommended)
          -n number of TOP templates to be selected     (optional, default 50)
          -c sequence identity cutoff                   (optional, default 0.30)
```
##### More details on some options:

* `-T` Path to template directory, containing `PDB, MTX, DEP, dssp, and seq` sub-directories.
          `pdb`: Non-redundant PDB structures,
          `MTX`: Sequence Profiles from PSI-BLAST search,
          `DEP`: Depth-dependent structure profiles,
         `dssp`: [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) generated output file,
          `seq`: FASTA sequence (Ref: https://zhanglab.ccmb.med.umich.edu/library/)         
* `-L` Path to the file containing template (pdb) IDs
* `-d` Path to distance histogram (distogram) `rawdistpred.current` generated by [DMPfold](https://github.com/psipred/DMPfold)
* `-m` Path to [MODELLER](https://salilab.org/modeller), e.g. `/home/XXXX/bin/modeller9.20/bin/modpy.sh`


### Test DisCovER

To run DisCovER, go to the target directory, which contains target sequence file, standard [SPIDER3](https://sparks-lab.org/downloads/) predicted output file, and profiles. Particularly, DisCovER will assume `seq.fasta, seq.spd33, seq.prf, seq.mtx` files are prersent in the current directory.

```sh
$ ulimit -s 419430400 ./discover             (To increase the stack limit, avoiding Segmentation Error)
$ cd Example
$ tar -zxvf DisCovER_template_Lib.tar.gz           (Uncompress a toy template library for the test run)
$ cd d1a9xb1/                     (Target directory contains respective input files as mentioned above)
$ ../../discover -T ../template_library/ -L ../template_library/template_list.txt -q d1a9xb1 -o ./ -d rawdistpred.current -m /home/XXXX/bin/modeller9.20/bin/modpy.sh -n 50 -c 0.30 
```
The first-ranked predicted 3D model will be named as `seq_model1.pdb` and the corresponding alignment file is named as `top1.fasta`.

### Generating input files from Query sequence


### Generating template libraries


## Data

- The template library can be found [here](data/FRAGFOLD_150.txt) 
- The benhcmark datasets along with input files can be found [here](data/CASP12_13_FM.txt) 
- DisCovER predicted first-ranked 3D full-length models can be found [here](data/CASP12_13_FM.txt) 

## Cite
If you find DisCovER useful, please cite our paper at bioRxiv:
