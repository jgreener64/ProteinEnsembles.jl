# ProteinEnsembles

[![Build Status](https://travis-ci.org/jgreener64/ProteinEnsembles.jl.svg?branch=master)](https://travis-ci.org/jgreener64/ProteinEnsembles.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/flfqouj1otkuf1rk?svg=true)](https://ci.appveyor.com/project/jgreener64/proteinensembles-jl)
[![Coverage Status](https://coveralls.io/repos/github/jgreener64/ProteinEnsembles.jl/badge.svg?branch=master)](https://coveralls.io/github/jgreener64/ProteinEnsembles.jl?branch=master)

This [Julia](http://julialang.org/) package implements the ExProSE algorithm that takes two protein structures and generates an ensemble of protein structures. The ensembles span conformational space and can be used to predict allosteric sites. The method is described in:

Greener JG, Filippis I and Sternberg MJE, Predicting protein dynamics and allostery using multi-protein atomic distance constraints, *Structure*, 2017, 25, 546-558. [Link to paper](http://www.cell.com/structure/fulltext/S0969-2126(17)30008-4).


## Summary

Install using `add ProteinEnsembles` from within the Julia package mode. Run using

```bash
exprose --i1 input_1.pdb --d1 input_1.dssp \
    --i2 input_2.pdb --d2 input_2.dssp \
    -n 50 -o exprose_out
```

where `exprose` is in the `bin` directory of the package.


## Installation

Julia is required and can be downloaded [here](http://julialang.org/downloads). Install ProteinEnsembles.jl by running `add ProteinEnsembles` from the Julia package REPL, which is entered by pressing `]`. This will also automatically install a few other required Julia packages. If you want, the tests can be run using `test ProteinEnsembles`. If you wish to use the auto-parameterisation procedure (see below) you must also have [TM-score](https://zhanglab.ccmb.med.umich.edu/TM-score) installed.


## Requirements

To use ProteinEnsembles.jl you will need the following:
- PDB files of the protein of interest. Two is best, but one may be used (see the paper). They must have polar hydrogens only added; this can be done using tools such as [Chimera](https://www.cgl.ucsf.edu/chimera) or [pdbtools](https://github.com/harmslab/pdbtools). The chain labelling and residue numbering must be consistent between the files as this is used to find common atoms. Alternative atom locations are discarded. PDB files must also be a single model and not have any inserted residues. HETATM records are discarded by default.
- DSSP files corresponding to the PDB files above. These can be obtained using [dssp](http://swift.cmbi.ru.nl/gv/dssp).


## Usage

These instructions are tailored towards Mac/Unix. However they could be modified to work on Windows.

Although organised as a Julia package, ProteinEnsembles.jl is primarily designed for use from the command line. The `exprose` script in the `bin` directory implements this. For example, to see the command line options run

```bash
~/.julia/packages/ProteinEnsembles/xxxxx/bin/exprose -h
```

where `xxxxx` is the directory of the latest package version. For easy access to the `exprose` command you might like to add the following line to your profile:

```bash
export PATH=$PATH:~/.julia/packages/ProteinEnsembles/xxxxx/bin
```

Then, if all input files are in your current directory, run the program as follows:

```bash
# Generate an ensemble of 50 structures with an output directory exprose_out
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out

# Use a tolerance weighting of 0.5
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out -w 0.5

# Generate an ensemble from a single structure with a tolerance weighting of 1.0
exprose --i1 input_1.pdb --d1 input_1.dssp -n 50 -o exprose_out -w 1.0
```

The method may also be run from within Julia. The below Julia script does the same thing as the first example above:

```julia
using ProteinEnsembles
runpipeline(
    i1="input_1.pdb",
    d1="input_1.dssp",
    i2="input_2.pdb",
    d2="input_2.dssp",
    n_strucs=50,
    out_dir="exprose_out"
)
```

Or, to split it up a little into the constituent functions:

```julia
using ProteinEnsembles
constraints_com, constraints_one, constraints_two = interactions(
    "input_1.pdb",
    "input_1.dssp",
    "input_2.pdb",
    "input_2.dssp"
)
ensemble_com = generateensemble(constraints_com, 50)
runanalysis("exprose_out", ensemble_com, constraints_one, constraints_two)
```


### Selecting parameters

The auto-parameterisation procedure can select a more suitable tolerance weighting value (see the paper). [TM-score](https://zhanglab.ccmb.med.umich.edu/TM-score) must be installed to do this. For example:

```bash
exprose-param --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -o exprose_param -t TMscore
```

runs the auto-parameterisation procedure with the `-t` option specifying the command to run TM-score. The last line of the output gives a suggested tolerance weighting. This value is also written out to `suggested.tsv`. Use this value in a normal `exprose` run as above.


### Allosteric site prediction

To predict allosteric sites you should run [LIGSITEcs](http://projects.biotec.tu-dresden.de/pocket/download.html) on the *second* input structure (the one you give as `--i2`). You then need to run the `cluster-ligsite` script in `bin` to assign the points to pockets:

```bash
cluster-ligsite pocket_r.pdb pocket_all.pdb pocket_points.pdb
```

where `pocket_r.pdb` and `pocket_all.pdb` are in the LIGSITEcs output. Then carry out an `exprose` run with the `pocket_points.pdb` file (`-l`) and the number of pockets (e.g. top 4) to perturb at (`-m`) as parameters:

```bash
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out -l pocket_points.pdb -m 4
```

A tolerance weighting from an auto-parameterisation run can also be used here. View the `predictions.tsv` output file to get the order of allosteric pocket predictions. Note that other pocket prediction software can be used provided you can get the output into the same format as `pocket_points.pdb`, i.e. pocket cavity points with the pocket number in the residue number column.


### Output

The output directory contains the following:
- `input_1.pdb` and `input_2.pdb`: atoms used from the input structures are written back out and superimposed.
- `pdbs`: generated structures in PDB format. Superimposed to `input_1.pdb` and `input_2.pdb`.
- `pcs`: projections onto the principal components (PCs) from the principal component analysis of the generated structures. Contains files for generated (`pcs.tsv`) and input structures (`pcs_input_1.tsv` and `pcs_input_2.tsv`) - line n corresponds to structure n and column c corresponds to PC c. Has graphs of these for the first few PCs (`pc_x_y.png`). Also includes a list of PCs ordered by decreasing distance between the input structures (`pcs_input_dist.tsv`) and the percentage variation explained by each PC (`evals_spread.tsv`).
- `pymol`: [PyMol](https://www.pymol.org/) scripts to view PCs on `input_1.pdb`, e.g. run `pymol input_1.pdb pymol/view_pc_1.pml`.
- `rmsds_input_1.tsv` and `rmsds_input_2.tsv`: RMSDs of generated structures to the input structures. Line n corresponds to structure n.
- `rmsfs.tsv` and `rmsfs.png`: RMSFs of each residue over the ensemble of generated structures, and a plot of this. Line n corresponds to residue index n.
- `spe_scores.tsv`: SPE error scores of generated structures (see paper). Line n corresponds to structure n.

For allosteric site prediction there will be `pdbs_mod_n` and `mod_n` containing similar information for each perturbed ensemble, as well as the ratio of RMSF values to the unperturbed ensemble (`rmsfs_ratio.tsv`). There will also be the order of allosteric predictions (`predictions.tsv`) and the size of the perturbation on modulating each site (`perturbations.tsv`), which is the RMSD between the centroid structure of the perturbed and unperturbed ensembles.

The default plot colours are blue for generated structures, red for input structure 1, green for input structure 2 and orange for perturbed ensemble structures.


### Reproducing paper results

The results from the paper can be generated using the instructions in `paper_results`. Information on the ensemble and allosteric datasets is in `datasets`.


### Performance

ExProSE can generate 250 structures in ~20 minutes for T4-lysozyme (162 residues) on a 3.1 GHz Intel Core i7 processor.


## Reporting issues

If you find any bugs in the software or have a comment or feature request, please open an issue on GitHub or email Joe Greener (j.greener at ucl.ac.uk).


## Notes

- All default values for parameters used in the code can be found and modified in `src/defaults.jl`.
- Auto-parameterisation works fine on all OSs but the auto-parameterisation tests are disabled by default to make the CI build pass. If you want to run the parameterisation tests, set `run_param_test` in `test/runtests.jl` to `true`.
- Julia utilities to deal with protein structures and PDB files can be found in [Bio.jl](http://biojulia.github.io/Bio.jl/) and [MIToS.jl](http://diegozea.github.io/MIToS.jl/).
- ExProSE users might also like to try [tCONCOORD](http://wwwuser.gwdg.de/~dseelig/tconcoord.html) and [NMSim](http://cpclab.uni-duesseldorf.de/nmsim/).
