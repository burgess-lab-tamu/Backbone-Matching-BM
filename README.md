# IBM

A novel strategy of loop mimicry to discover mimics by overlaying all-Ala cyclo-organopeptides on loop targets and find hits by assessing goodness of fit.

## Usage

dataset/
```python
scaffolds: SMILES strings of 86 rigid and cheap organic scaffolds mentioned in the paper.
cyclo-organopeptides: SMILES strings of 602 cyclo-organopeptides, comprised of Ala and organic fragments, catergorized by sequence length.
conformers: include a link to access generated conformers of 602 cyclo-organopeptides.
```

scripts/ 
```python
alignment:
align.py: thoroughly overlay conformers of cyclo-organopeptides with loop segments over 4-10 residue length.
precise_align.py: overlay conformers of cyclo-organopeptides with loop segments on a fixed length, i.e. 5 residue length, thus only cyclo-organopeptides with 5A will be used in the overlay.
reminder: 
1) conformers must be downloaded to run these scripts.
2) the loop segment is obtained directly from PDB crystals.  In align.py it is recommended to input a sequence with at least 10 residues; or researchers may change the script to make it adaptable to shorter loops.
3) in precise_align, it is assumed hot spots on the target loop is known, thus only shortest segment including all hot spots is required as input, but not the whole loop sequence.

postprocess:
1) copy_best_overlays.py: extract best-overlaid conformers with their mimic segments and save them in a new folder.
2) overlays_pymol.py: present results in pymol, overlaying best-overlaid conformers with their mimic segments.  the results are catergorized into different loop segments for convenient check.
```

illustrative output/
```python
two illustrative outputs were presented here:
3BT1_uPA_scan.pse:  a graphical presentation of thorough scan over uPA hot loop (20-31) by align.py
uPA_23_30_precise_scan.pse:  a graphical presentation of precise scan over uPA hot loop segment (23-30, including all 5 hot spots) by precise_align.py
```
