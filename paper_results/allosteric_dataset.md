Instructions to reproduce the allosteric prediction results from the paper.

1. Download the apo and holo PDB files from `datasets/allosteric_dataset.tsv` and retain the chains listed.

2. Obtain DSSP files for the above PDB files.

3. For each protein run LIGSITEcs with default arguments on the holo structure.

4. For each protein obtain pocket point PDB files from the LIGSITEcs output files:

    ```bash
    cluster-ligsite pocket_r.pdb pocket_all.pdb pocket_points.pdb
    ```

5. For each protein run the auto-parameterisation procedure to obtain a suggested tolerance weighting (last line of the output):

    ```bash
    exprose-param --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
        --d2 input_2.dssp -o exprose_param
    ```

6. For each protein generate an ensemble and perturbed ensembles:

    ```bash
    exprose --i1 input_1.pdb --d1 input_1.dssp --i2 input_2.pdb \
        --d2 input_2.dssp -o exprose_param -w [tolerance_weight] \
        -m 8 -l pocket_points.pdb
    ```

    where the tolerance weight is that suggested in the previous step and `pocket_points.pdb` is from step 4.

7. For each protein view the `predictions.txt` output file to get the order of allosteric pocket predictions.
