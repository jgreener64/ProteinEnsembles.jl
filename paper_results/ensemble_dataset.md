Instructions to reproduce the ensemble generation results from the paper.

1. Download the apo and holo PDB files from `datasets/ensemble_dataset.tsv`. Add polar hydrogens.

2. Obtain DSSP files for the above PDB files.

3. For each protein run the auto-parameterisation procedure to obtain a suggested tolerance weighting (last line of the output):

    ```bash
    exprose-param --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
        --d2 input_2.dssp -o exprose_param
    ```

4. For each protein generate an ensemble:

    ```bash
    exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
        --d2 input_2.dssp -o exprose_param -w [tolerance_weight]
    ```

    where the tolerance weight is that suggested in the previous step.

5. Examine the RMSD files to see how close the generated structures came to the input structures.
