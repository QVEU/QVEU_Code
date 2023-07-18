# QVEU Collaborative Markdown of ChimeraX Tips and Tricks

This **cheat sheet** provides a quick start and reference with tips specific to QVEU activities. Please feel free to contribute new sections to this guide, including links to any accessory scripts in github repositories. 

## ChimeraX
*ChimeraX is a macromolecular structure visualization software.*

![ChimeraX](https://www.cgl.ucsf.edu/chimerax/docs/quickstart/images/chimerax.png)

## Installation
ChimeraX can be installed from UCSF: [https://www.cgl.ucsf.edu/chimerax/download.html#release](https://www.cgl.ucsf.edu/chimerax/download.html#release)

(IT assistance will be needed to install on government device)

## Fetching structures
### Using the PDB ID
```
open 4WM7
```

## Appearance
```
show surface
hide surface

show atoms
hide atoms
```

## Colors
### Color by chain
```
color bychain
```
### Color by Sequence Conservation
Open alignment in ChimeraX containing sequences matching structure(s) of interest.
Structures should automatically associate with the relevant sequence. The header color in the alignment will match the structure color if this has happened. To manually associate a structure to a sequence, right click the alignment window, select "Structure" then "Associations" and match the sequence to the structure.
```
color byattr seq_conservation
```

Use sequence conservtion with conditions
```
color #1/A::seq_conservation=1 red #color fully conserved residues red in structure 1 chain A
```
## Biophysics
### Display hydrophobic surface
```
mlp
```
### Display electrostatic surface
```
coloumbic
```

## Compare Structures
```
match #2 to #1/a 
```
