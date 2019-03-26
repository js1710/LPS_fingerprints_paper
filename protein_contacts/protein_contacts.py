from __future__ import print_function, division
import MDAnalysis as mda
from collections import OrderedDict
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
import pickle
import argparse

def get_hist(structure, trajectory, start=0, stop=-1, step=1, contact_sel="name PO1 PO2"):
    u = mda.Universe(structure, trajectory)
    protein = u.select_atoms("name BB SC1 SC2 SC3 SC4")
    contact_group = u.select_atoms(contact_sel)
    resname_hist = {}
    for frame in u.trajectory[start:stop:step]:
        tree = AtomNeighborSearch(protein, box=u.dimensions)
        for atom in contact_group:
            near = tree.search(atom, 6., level='R')
            if len(near) == 0:
                #print(near)
                continue
            else:
                resids = [x.resid for x in near]
                for resid in resids:
                    if resid in resname_hist:
                        resname_hist[resid] +=1
                    else:
                        resname_hist[resid] = 1
    resname_hist = OrderedDict(sorted(resname_hist.items()))
    return resname_hist

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", "--begin", type=int, default=0, help="starting step for analysis")
    parser.add_argument("-e", "--end", type=int, default=-1, help="final step for analysis") 
    parser.add_argument("-o", "--output", default="results", help="prefix for name of pickle file outputted") 
    parser.add_argument("-s", "--select", default="name P*", help="MDAnalysis selection string for contacts target group (origin group is the protein)") 
    parser.add_argument("-c", "--coord", required=True, help="Input structure file *.gro/*.pdb")
    parser.add_argument("-t", "--traj", required=True, help="Input trajectory *.xtc/*.trr" )
    args = parser.parse_args()
    hist_lpss = {}
    proteins = ["OMPA", "OMPX", "BTUB", "FHUA", "OMPF", "ESTA"]
    lps_all = ["RE", "RE_r"]
    for protein in proteins:
        print(protein)
        for lps in lps_all:
            print(lps)
            path = "./{0}/{1}/".format(protein, lps)
            hist = get_hist(path + args.coord, path + args.traj, start=args.begin, stop=args.end, contact_sel=args.select)
            if protein in hist_lpss:
                hist_lpss[protein][lps] = hist
            else:
                hist_lpss[protein] = {}
                hist_lpss[protein][lps] = hist
    pickle.dump( hist_lpss, open( args.output + ".p", "wb" ) )
