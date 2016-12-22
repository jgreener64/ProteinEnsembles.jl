# ProteinEnsembles

[![Build Status](https://travis-ci.org/jgreener64/ProteinEnsembles.jl.svg?branch=master)](https://travis-ci.org/jgreener64/ProteinEnsembles.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/flfqouj1otkuf1rk?svg=true)](https://ci.appveyor.com/project/jgreener64/proteinensembles-jl)
[![Coverage Status](https://coveralls.io/repos/github/jgreener64/ProteinEnsembles.jl/badge.svg?branch=master)](https://coveralls.io/github/jgreener64/ProteinEnsembles.jl?branch=master)

This [Julia](http://julialang.org/) package implements the ExProSE algorithm that takes two protein structures and generates an ensemble of protein structures. The ensembles span conformational space and can be used to predict allosteric sites. The method is described in:

Greener JG, Filippis I and Sternberg MJE, Predicting protein dynamics and allostery using multi-protein atomic distance constraints, *Structure* (Accepted)


## Summary

Install using `Pkg.clone("https://github.com/jgreener64/ProteinEnsembles.jl.git")` from within Julia v0.5. Run using

```bash
exprose --i1 input_1.pdb --d1 input_1.dssp \
    --i2 input_2.pdb --d2 input_2.dssp \
    -n 50 -o exprose_out
```

where `exprose` is in the `bin` directory.


## Installation

Julia v0.5 is required and can be downloaded [here](http://julialang.org/downloads). Install ProteinEnsembles.jl by running

```julia
Pkg.clone("https://github.com/jgreener64/ProteinEnsembles.jl.git")
```

from the Julia REPL. This will also automatically install a few other required Julia packages. If you want, the tests can be run using `Pkg.test("ProteinEnsembles")`.

If you wish to use the auto-parameterisation procedure you must also have [TM-score](http://zhanglab.ccmb.med.umich.edu/TM-score) installed (see below).


## Requirements

To use ProteinEnsembles.jl you will need the following:
- PDB files of the protein of interest. Two is best, but one may be used (see the paper). They must have polar hydrogens only added; this can be done using tools such as [Chimera](https://www.cgl.ucsf.edu/chimera) or [pdbtools](https://github.com/harmslab/pdbtools). The chain labelling and residue numbering must be consistent between the files as this is used to find common atoms. Alternative atom locations are discarded. PDB files must also be a single model and not have any inserted residues. HETATM records are discarded by default.
- DSSP files corresponding to the PDB files above. These can be obtained using [dssp](http://swift.cmbi.ru.nl/gv/dssp).


## Usage

Although organised as a Julia package, ProteinEnsembles.jl is primarily designed for use from the command line. The `exprose` script in the `bin` directory implements this. For example, to see the command line options run

```bash
~/.julia/v0.5/ProteinEnsembles/bin/exprose -h
```

For easy access to the `exprose` command you might like to add the following line to your profile:

```bash
export PATH=$PATH:~/.julia/v0.5/ProteinEnsembles/bin
```

Then, if all input files are in your current directory, run the program as follows:

```bash
# Generate an ensemble of 50 structures with an output directory ake_out
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out

# Use a tolerance weighting of 0.5
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out -w 0.5

# Generate an ensemble from a single structure with a tolerance weighting of 1.0
exprose --i1 input_1.pdb --d1 input_1.dssp -n 50 -o exprose_out -w 1.0

# Perturb the ensemble at 4 sites (see below)
exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -n 50 -o exprose_out -m 4 -l pocket_points.pdb
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

The auto-parameterisation procedure can select a more suitable tolerance weighting value (see the paper). [TM-score](http://zhanglab.ccmb.med.umich.edu/TM-score) must be installed to do this. For example:

```bash
# Run auto-parameterisation procedure; TM-score is run with the command TMscore
exprose-param --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
    --d2 input_2.dssp -o exprose_param -t TMscore
```
