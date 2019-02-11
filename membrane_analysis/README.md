# Membrane analysis scripts
Contains scripts to calculate the membrane thickness, Depletion-Enrichment index and the Depletion-Enrichment 2D map.
## Membrane thickess
A descritpion of the options availbe can be found using the following
```
python membrane_thickness.py -h
```
If a protein is present in the system the trajectory should be preprocessed to center the trajectory on the protein. This script will not work correctly if the trajectory has been rotationally fitted in any manner.

## Enrichment Map
A descritpion of the options availbe can be found using the following
```
python enrichment_map.py -h
```
This script's input is the output of either the `gmx density` command or the density script written by M.Castillo et.al (http://perso.ibcp.fr/luca.monticelli/tools/index.html).

## Enrichment Index
A descritpion of the options availbe can be found using the following
```
python enrichment_index.py -h
```
