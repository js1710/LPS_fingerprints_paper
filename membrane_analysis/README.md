# Membrane analysis scripts
Contains scripts to calculate the membrane thickness, Depletion-Enrichment index and the Depletion-Enrichment 2D map.
## Membrane thickess
A decritpion of the options availbe can be found using the follwing
```
python membrane_thickness.py -h
```
If a protein is present in the system the trajectory should be preprocessed to center the trajectory on the protein. This script will not work correctly if the trajectory has been rotationally fitted in any manner.
