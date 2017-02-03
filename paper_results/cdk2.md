Instructions to reproduce the CDK2 results from the paper.

1. Obtain PDB files 1HCL (apo) and 3PXF (allosteric modulator bound). Add polar hydrogens. Save as `1HCL.pdb` and `3PXF.pdb`.

2. Obtain DSSP files for the above PDB files. Save as `1HCL.dssp` and `3PXF.dssp`.

3. Run

    ```bash
    exprose --i1 1HCL.pdb --d1 1HCL.dssp --i2 3PXF.pdb --d2 3PXF.dssp \
        -o exprose_out -n 250 -w 0.3 -m 8 -l 3PXF_pocket_points.pdb
    ```

    where `3PXF_pocket_points.pdb` is included in this directory, or can be generated manually using `cluster-ligsite` after running LIGSITEcs on 3PXF.

4. Compare the principal components analysis plots for each site.
