# Backbone Matching (BM)

streamlined protein-protein interface loop mimicry strategy<br>
a novel strategy of loop mimicry to discover mimics by overlaying all-Ala cyclo-organopeptides on loop targets and find hits by assessing goodness of fit.

## Usage

requirements
```python
language: Python
packages: Pymol, Pandas, Numpy
```


dataset
```python
scaffolds: SMILES strings of 86 rigid and cheap organic scaffolds mentioned in the paper.
cyclo-organopeptides: SMILES strings of 602 cyclo-organopeptides, comprised of Ala and organic fragments, categorized by sequence length.
conformers: include a link to access generated conformers of 602 cyclo-organopeptides.
```

scripts
```python
prerequisites: prepare following folders under the directory of scripts.
  1) loop_target: put loop structural files here
  2) conformer_for_matching: put unzipped conformer files of 602 cyclo-organopeptides here
  3) align_result: empty folder, to store the RMSD values from the alignment calculation
  4) mimic_result: empty folder, to store best hits from the alignment calculation

precise-alignment: assume a hot loop of length k, overlay this loop with cyclo-{-(Ala)k-organo-}.
prerequisite: cut the loop structure file (in PDB format) to the desired fragment
input: python precise-alignment.py 'loop-name.pdb' 'number of hits to store' 'length of input loop'
example input: python precise-alignment.py 'uPA.pdb' 50 8
output: RMSD values in folder 'align_result', overlays of best hits with the loop fragment in the same folder as the script in Pymol pse format

auto-slicing-alignment: slice a long loop to fragments with 4 - 10 amino acids, then precisely overlay each fragment with cyclo-{-(Ala)n-organo-}, n = length of the fragment.
prerequisite: extract the long loop from the original PDB file, then save in pdb format; if the length is shorter than 10, then modify the 'alignment' function in the script
input: python auto-slicing-alignment.py 'loop-name.pdb' 'number of hits to store'
example input: python auto-slicing-alignment.py 'uPA.pdb' 50
output: RMSD values in folder 'align_result', overlays of best hits with different loop fragments in the same folder as the script in Pymol pse format

virtual screening: screen a library of loops using precise alignment.
This is the script used to screen 1398 Kritzer hot loops.  The results are in a zip file shared in the published paper.
```

illustrative output
```python
two illustrative outputs are presented here:
3BT1_uPA_scan.pse:  a graphical presentation of a thorough scan over uPA hot loop (20-31) by auto-slicing alignment
uPA_23_30_precise_scan.pse:  a graphical presentation of precise scan over uPA segment (23-30, including all 5 hot spots) by precise alignment
```

## Citation
to be published
