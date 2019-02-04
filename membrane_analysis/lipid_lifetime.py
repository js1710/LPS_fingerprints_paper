#!/home/jon/anaconda2/envs/py36/bin/python
import numpy as np
import tqdm
import MDAnalysis as mda
import MDAnalysis.lib.pkdtree as pkdtree
from MDAnalysis.core.AtomGroup import Atom
import argparse
from numba import jit


class NeighbourSearchFast(object):
    '''Class to rapility calculate the neighbours of 'atomgroup' with the 'searchgroup'. Note that all atomgroups
     are MDAanalysis objects'''
    def __init__(self, atomgroup, box=None):
        self.atomgroup = atomgroup
        self._u = atomgroup.universe
        self.box = box
        self.kdtree = pkdtree.PeriodicKDTree(box=self.box)


    def search(self, searchgroup, radius, level="A"):
        self.kdtree.set_coords(self.atomgroup.positions, cutoff=radius + 0.1)
        if isinstance(searchgroup, Atom):
            positions = searchgroup.position.reshape(1, 3)
        else:
            positions = searchgroup.positions

        unique_idx = self.kdtree.search(positions, radius)
        self.indices = unique_idx
        return self._index2level(unique_idx, level)

    def _index2level(self, indices, level):
        n_atom_list = self.atomgroup[indices]
        if level == 'A':
            if not n_atom_list:
                return []
            else:
                return n_atom_list
        elif level == 'R':
            return list({a.residue for a in n_atom_list})
        elif level == 'S':
            return list(set([a.segment for a in n_atom_list]))
        else:

            raise NotImplementedError('{0}: level not implemented'.format(level))

def distance_sq(pos1, pos2):
    diff = pos2 - pos1
    return np.sum(np.square(diff))

@jit(nopython=True)
def survival_prob(lipids, dt):
    bound = 0
    total = 0
    for i in range(len(lipids)):
        lipid = lipids[i]
        for i in range(len(lipid) - dt):
            if lipid[i] == 1:
                if lipid[i + dt] == 1:
                    bound += 1
                total += 1

    return bound, total

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Calculate on off rate of lipids with proteins",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("structure", help="*.gro/pdb file of system")
    parser.add_argument("trajectory", help="*.xtc/trr trajectory of system")
    parser.add_argument("--begin", type=int, default=0, help="starting frame number")
    parser.add_argument("--stop", type=int, default=-1, help="final frame number")
    parser.add_argument("-o", "--output", default="lifetime", help="Output files prefix")
    parser.add_argument("--proximity", type=float, default=6.0, help="contact distance for protein-lipid interactions")

    args = parser.parse_args()
    prefix = args.output
    proximity = args.proximity
    lipid_select = "resname REMP RAMP POPE POPG CDL2 DPPC"
    protein_select = "protein"
    u = mda.Universe(args.structure, args.trajectory)
    lipids = u.select_atoms(lipid_select)
    protein = u.select_atoms(protein_select)
    lipids_on_prev = []
    time_on = {}
    dt = (u.trajectory[1].time - u.trajectory[0].time) / 1000. # dt in ns
    u.trajectory[0] # reset iterator
    start = args.begin
    stop = args.stop
    step = 1
    start, stop, step = u.trajectory.check_slice_indices(start, stop, step)
    lifetimes_per_resid = {}

    lipid_frame_contacts = dict([(resid, np.zeros(len(u.trajectory[start:stop:step]))) for resid in lipids.residues.resids])
    for i, frame in enumerate(tqdm.tqdm(u.trajectory[start:stop:step])):
        frame_contacts = dict([(resid, 0) for resid in lipids.residues.resids])

        tree = NeighbourSearchFast(lipids, box=u.dimensions)
        lipids_near = tree.search(protein, proximity, level='R')
        resids_on = [lipid.resid for lipid in lipids_near]

        if not lipids_on_prev:
            lipids_on_prev = resids_on
            for resid in resids_on:
                time_on[resid] = dt
            continue

        for resid in resids_on:
            lipid_frame_contacts[resid][i] += 1
            if resid in lipids_on_prev:
                time_on[resid] += dt
            else:
                time_on[resid] = dt

        for resid in lipids_on_prev:
            if resid not in resids_on:
                try:
                    lifetimes_per_resid[resid].append(time_on[resid])

                except KeyError:
                    lifetimes_per_resid[resid] = [time_on[resid]]
                finally:
                    del time_on[resid]
        lipids_on_prev = resids_on




    for resid in time_on:
        try:
            lifetimes_per_resid[resid].append(time_on[resid])
        except KeyError:
            lifetimes_per_resid[resid] = [time_on[resid]]

    #print(time_on)

    #convert resids to resnames
    lifetimes = {}
    survival_rates = {}
    for resid in lifetimes_per_resid:
        resname = u.select_atoms("resid " + str(resid)).residues[0].resname
        try:
            survival_rates[resname].append(lipid_frame_contacts[resid])
            lifetimes[resname].extend(lifetimes_per_resid[resid])
        except KeyError:
            lifetimes[resname] = lifetimes_per_resid[resid]
            survival_rates[resname] = [lipid_frame_contacts[resid]]

    for resname in survival_rates:
        survival_rates[resname] = np.array(survival_rates[resname])

    #print(lifetimes["REMP"])
    for resname in survival_rates:
        with open("{0}_{1}_survival.dat".format(prefix, resname), "w") as out:
            template = "{0:^6d}{1:^8.3e}"

            for dt in tqdm.tqdm(range(0, 20001, 1)):
                bound, total = survival_prob(survival_rates[resname], dt)
                try:
                    prob = bound / total
                except ZeroDivisionError:
                    prob = 0.
                #print(template.format(dt, prob))
                print(template.format(dt, prob), file=out)

    with open("{0}_summary.dat".format(prefix), "w") as out:
        for resname, times in lifetimes.items():
            average = np.mean(times)
            error = np.std(times)
            std_error = error / np.sqrt(len(times))
            bootstraps = []
            for n in range(200):
                weights = np.random.random(len(times))
                weights = weights / np.sum(weights)
                resample = np.random.choice(times, len(times), p=weights)
                resample_av = np.mean(resample)
                bootstraps.append(resample_av)

            bootstraps = np.array(bootstraps)
            print(len(bootstraps), len(lifetimes))
            print("Lipid {0}: {1:8.3f} +\- {2:8.3f} ns".format(resname, average, std_error))
            print("Bootstrap: {0:8.3f} +\- {1:8.3f} ns".format(np.mean(bootstraps),
                                                               np.std(bootstraps)))
            print("Lipid {0}: {1:8.3f} +\- {2:8.3f} ns".format(resname, average, std_error), file=out)
            print("Bootstrap: {0:8.3f} +\- {1:8.3f} ns".format(np.mean(bootstraps),
                                                               np.std(bootstraps)), file=out)

    for resname, times in lifetimes.items():
        np.savetxt("{0}_{1}_all.dat".format(prefix, resname), np.array(times))
