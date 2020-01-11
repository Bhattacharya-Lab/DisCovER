# DisCovER

Distance-based Covariational Protein Threading for Low-homology Fold Recognition.

## Installation

Installing DisCovER is very straightforward. The following instructions should work for 64-bit Linux system:

- DisCovER has been tested on 64 bit Linux operating system.
- Install [MODELLER](https://salilab.org/modeller), a license key is required. This can be installed using command `conda install modeller -c salilab`. DisCovER has been tested on MODELLER version 9.20.

That's it! DisCovER is ready to be used.

## Usage

To see the usage instructions, run `./discover -h`
```
---------------------------------------------------------------------
                        DisCovER
 Distance-based Covariational Protein Threading  V2020.Jan.09
 For comments, please email to bhattacharyad@auburn.edu
---------------------------------------------------------------------

Usage:   ./discover [-option] [argument]

Options:  -T path to template library             - input, required
          -L list.txt with template IDs           - input, required
          -q query ID                             - input, required
          -o path to output director              - input, required
          -d path to query distance file          - input, required
          -m path to MODELLER                     - input, optional
          -n number of TOP templates to be selected  (optional, default 50)
          -c sequence identity cutoff                (optional, default 0.30)
```

### Test DisCovER

To run DisCovER, go to the target directory, which contains target sequence file, SPIDER3 predicted output file, profiles, and DMPfold predicted distance map. Particularly, DisCovER will assume `seq.fasta,seq.spd33,seq.prf,seq.mtx,rawdistpred.current` files are prersent in the current directory.

```sh
$ cd example
$ discover -a A.map -b B.map

Top predicted model will be named as `seq_model1.pdb`.
```

### Distance Map Format 





## Data

- The list of PDB chains of FRAGFOLD dataset can be found [here](data/FRAGFOLD_150.txt) 
- The list of Target IDs of CASP12 and CASP13 FM targets can be found [here](data/CASP12_13_FM.txt) 
- The list of PDB chains of 510 membrane protein dataset can be found [here](data/Membrane_510.txt) 
- The list of PDB chains of EVfold dataset can be found [here](data/EVfold_15.txt) 
