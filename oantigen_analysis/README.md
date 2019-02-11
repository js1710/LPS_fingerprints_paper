# O-Antigen analysis
This subdirectory contains all scripts related to the conformations of the O-Antigen chains.
## O-Antigen tilt
All options availble to the script can be found with the following command:
```
python oantigen_tilt.py -h
```
This script determines the tilt of O-Antigen chains wrt to the normal to the membrane (z-axis).

## O-Antigen pair correlation
All options availble to the script can be found with the following command:
```
python oantigen_pair_correlation.py -h
```
This script determines the realtive tilt of all pairs of O-Antigen chains within a given cutoff distance of each other. 
This script also outputs probability distributions of the pairwise tilt and relative pairwise tilt to 
`pair_prob_angles.dat` and `pair_prob_angles_dir.dat`, respectively.  

The relative pairwise tilt gives the angles between a pair of O-Antigen chains and if they tilt away from each other the angle is negative.
Thus the relative pairwise tilts takes values in the interval [-180, 180]. When working in units of cos(theta) then the relative pairwise tilt takes 
values in the interval [-2, 2].
