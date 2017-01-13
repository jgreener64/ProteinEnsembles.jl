Instructions to reproduce the T4-lysozyme results from the paper.

1. Obtain PDB files 169L (chain E only) and 2LZM. Add polar hydrogens. Save as `169L.pdb` and `2LZM.pdb`.

2. Obtain DSSP files for the above PDB files. Save as `169L.dssp` and `2LZM.dssp`.

3. Run

    ```bash
    exprose --i1 169L.pdb --d1 169L.dssp --i2 2LZM.pdb --d2 2LZM.dssp \
        -o exprose_out -n 250 -w 0.2
    ```
