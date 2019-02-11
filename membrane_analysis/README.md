# Membrane analysis scripts
Contains scripts to calculate the membrane thickness, Depletion-Enrichment index and the Depletion-Enrichment 2D map.
## Membrane thickess
A description of the options available can be found using the following
```
python membrane_thickness.py -h
```
If a protein is present in the system the trajectory should be preprocessed to center the trajectory on the protein. This script will not work correctly if the trajectory has been rotationally fitted in any manner.

## Enrichment Map
A description of the options available can be found using the following
```
python enrichment_map.py -h
```
This script's input is the output of either the `gmx densmap` command or the density script written by M.Castillo et.al (http://perso.ibcp.fr/luca.monticelli/tools/index.html).

## Enrichment Index
A description of the options available can be found using the following
```
python enrichment_index.py -h
```
The input for this script is quite complex:
1. We determine the number of each lipid type within `x` nm of the protein of interest. Below is is an example for coarse-grained POPE:
```
gmx select -s production.tpr -f production_whole.xtc -select "res_cog of resname POPE and within 0.7 of name BB SC1 SC2 SC3 SC4" -tu us -e 30 -os shell_7_POPE -pbc yes
```
Note that subsequent steps use all of the data outputted at this step.

2. An enrichment .inp file is required, an example of which is shown below
```
resname POPE POPG CDL2
POPE  shell_11_POPE.xvg  "resname POPE"  blue
POPG  shell_11_POPG.xvg  "resname POPG"  green
CDL2  shell_11_CDL2.xvg  "resname CDL2"  red
```
The first line contains the [MDanalysis style](https://www.mdanalysis.org/) selection of the membrane. The other lines have the format, `[lipid label] [gmx_select_ouput]  [selection for lipid] [colour of output bar in plot]`.Note that the `[lipid label]` is an arbitary name for the lipid that will be ouputted to the final plot.

3. Now we have all the files required to run this script. This script also outputs the enrichment index as a function of time or each lipid if the user wishes to do analysis further to what is availble here.

## Lipid survivability probability
A description of the options available can be found using the following
```
python lipid_lifetime.py -h
```
This script calculates two properties:
1. The average time contiguously spent at the surface of a protein by lipids of a given type.
2. The Probability of finding a lipid bound at t and t+dt (based on the implementation [here](https://github.com/GrossfieldLab/loos) 
