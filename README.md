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

To run DisCovER, go to the target directory, which contains target sequence file, standard SPIDER3 predicted output file, and profiles. Particularly, DisCovER will assume `seq.fasta, seq.spd33, seq.prf, seq.mtx` files are prersent in the current directory.

```sh
$ ulimit -s 419430400 ./discover          (To increase the stack limit, avoiding Segmentation Error)
$ cd example/d1a9xb1/                     (Target directory contains respective input files as mentioned above)
$ ./discover -T ../template_library/ -L ../template_library/template_list.txt -q d1a9xb1 -o ./ -d rawdistpred.current -m /home/XXXX/bin/modeller9.20/bin/modpy.sh -n 50 -c 0.30 
```
The first-ranked predicted 3D model will be named as `seq_model1.pdb` and the corresponding alignment file is named as `top1.fasta`.

## Data

- The template library can be found [here](data/FRAGFOLD_150.txt) 
- The benhcmark datastes along with input files can be found [here](data/CASP12_13_FM.txt) 
- DisCovER predicted first-ranked 3D full-length models can be found [here](data/CASP12_13_FM.txt) 
