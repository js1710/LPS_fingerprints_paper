from __future__ import print_function
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import argparse

def largest_groups(atoms):
    '''
    From a list of sizes, find out the indices of the two largest groups. These should correspond to the two leaflets of the bilayer.
    Keyword arguments:
    atoms -- list of sizes of clusters indentified by LeafletFinder
    '''
    largest=0
    second_largest=0

    for i in atoms:
        if atoms[i]>largest:
            largest_index=i
            largest = atoms[i]

    for i in atoms:
        if atoms[i]>second_largest and i!=largest_index:
            second_largest_index=i
            second_largest = atoms[i]

    return (largest_index,second_largest_index)


def determine_leaflets(universe,phosphateSelection="name P*"):
    '''
    From a selection of phosphates, determine which belong to the upper and lower leaflets.
    Keyword arguments:
    universe -- an MDAnalysis Universe object
    phosphateSelection -- a string specifying how to select the phosphate atoms (e.g. "name P1")
    '''
    leaflets = {}
    #print(universe)
    # calculate the z value of the phosphates defining the bilayer (assumes bilayer is in x and y..)
    #po4=universe.atoms.select_atoms(phosphateSelection)
    #bilayerCentre = po4.center_of_geometry()[2]
    #print(bilayerCentre)

    # apply the MDAnalysis LeafletFinder graph-based method to determine the two largest groups which
    #  should correspond to the upper and lower leaflets
    phosphates = LeafletFinder(universe,phosphateSelection)


    # find the two largest groups - required as the first two returned by LeafletFinder, whilst usually are the largest, this is not always so
    (a,b) = largest_groups(phosphates.sizes())
    layer1_avz = np.mean(phosphates.group(a).center_of_geometry()[2])
    layer2_avz = np.mean(phosphates.group(b).center_of_geometry()[2])
    bilayerCentre = (layer1_avz + layer2_avz) /2.
    # check to see where the first leaflet lies
    if phosphates.group(a).centroid()[2] > bilayerCentre:
        #print(phosphates.group(b).residues.atoms)
        leaflets["upper"] = phosphates.group(a).residues.atoms
        leaflets["lower"] = phosphates.group(b).residues.atoms
    else:
        leaflets["lower"] = phosphates.group(a).residues.atoms
        leaflets["upper"] = phosphates.group(b).residues.atoms

    return leaflets

def user_interface():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("structure", help="coordinate file (*.gro or *.pdb)")
    parser.add_argument("-p", "--phosphates", default="name P*", help="MDAnalysis selection text for phosphates")
    parser.add_argument("-l", "--lipids", default="all", help="MDAnalysis selection for lipids to output to *.gro or"
                                                              "*.ndx files")
    parser.add_argument("-ps", "--protein", default=None, help="MDAnalysis selection for protein")
    parser.add_argument("-wl", "--whole_lipid", action="store_true", help="ouputs the whole lipid to *.ndx and *.gro files")
    parser.add_argument("-layer", "--layer", default="both", choices=["both", "upper", "lower"], help="select layer to output")
    parser.add_argument("-o", "--output", default="leaflets", help="prefix for all output files")
    parser.add_argument("-s", "--select", default="same", help="MDAnalysis selection for output. If 'same' then same as -p selection")
    return parser.parse_args()

if __name__ == '__main__':
    args = user_interface()
    u = mda.Universe(args.structure)
    leaflets = determine_leaflets(u, args.phosphates)
    prefix = args.output
    if args.protein != None:
        protein = u.select_atoms(args.protein)
    if args.layer == "both":
        for key in leaflets:

            if args.lipids != "all":
                leaflets[key] = leaflets[key].select_atoms("{0}".format(args.lipids))
            if not args.whole_lipid:
                leaflets[key] = leaflets[key].select_atoms("{0}".format(args.phosphates))
            elif args.select != "same":
                leaflets[key] = leaflets[key].select_atoms("{0}".format(args.select))


    else:
        key = args.layer
        if args.lipids != "all":
            leaflets[key] = leaflets[key].select_atoms("{0}".format(args.lipids))
        if not args.whole_lipid:
            leaflets[key] = leaflets[key].select_atoms("{0}".format(args.phosphates))
        elif args.select != "same":
            leaflets[key] = leaflets[key].select_atoms("{0}".format(args.select))


    with mda.selections.gromacs.SelectionWriter('{}.ndx'.format(prefix), mode='w') as ndx:
        if args.layer == "both":
            for key, leaflet in leaflets.iteritems():

                ndx.write(leaflet, name=key)
                leaflets[key] = leaflet

            if args.protein != None:
                ndx.write(protein, name="protein")
                leaflet += protein
                ndx.write(leaflet, name="all")
            leaflet.write(key + "_layer.gro")
        else:
            leaflet = leaflets[args.layer]
            ndx.write(leaflet, name=key)
            leaflets[key] = leaflet
            if args.protein != None:
                ndx.write(protein, name="Protein")
                leaflet += protein
                ndx.write(leaflet, name="all")
            leaflet.write(key + "_layer.gro")





