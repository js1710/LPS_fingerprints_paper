from __future__ import print_function, division
import numpy as np
import MDAnalysis as mda
import multiprocessing as mp
import argparse
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch



def get_tilt(vec1, vec2):
    '''dot product method'''
    dot = np.dot(vec1, vec2)
    if dot < 1:
        tilt = np.arccos(dot)
    else:
        tilt = 0.
    tilt *= 0.5 * 360. / np.pi
    return tilt


def pair_angles(structure, trajectory, start=0, stop=-1, step=1, cutoff=6.0, costheta=False):
    '''calclates the angles between the end to end vectors of oantigen neighbours over a given trajectory'''
    u = mda.Universe(structure, trajectory)

    if stop == -1:
        stop = u.trajectory.n_frames
    start, stop, step = u.trajectory.check_slice_indices(
        start, stop, step)
    num_steps = len(range(start, stop, step))
    if num_steps == 0:
        raise ValueError("Selected time interval does not exist in the trajectory")

    oantigens = u.select_atoms("resname OANT").residues
    sugars_all = oantigens.atoms.select_atoms("name O*")
    sugars = sugars_all.split("residue")

    angles = []
    angles2 = []
    for frame in u.trajectory[start:stop:step]:
        vectors = {}
        head_endpoints = {}
        tail_cogs = {}
        for sugar in sugars:
            resid = sugar.resids[0]
            tail = sugar.select_atoms("name O51 to O54")
            head = sugar.select_atoms("name O97 to O99")
            head_cog = head.center_of_geometry()
            tail_cog = tail.center_of_geometry()
            vector = head_cog - tail_cog
            vector /= np.linalg.norm(vector)

            tail_cogs[resid] = tail_cog
            head_endpoints[resid] = tail_cog + vector
            vectors[resid] = vector
        tree = AtomNeighborSearch(sugars_all, box=u.dimensions)

        for sugar in sugars:
            near = tree.search(sugar, cutoff, level='R')
            for neighbour in near:
                resid1 = neighbour.resid
                resid2 = sugar.resids[0]
                if resid1 != resid2:
                    vec1 = vectors[resid1]
                    vec2 = vectors[resid2]

                    angle = get_tilt(vec1, vec2)

                    # determine if oantigen chains face each other
                    d1 = np.linalg.norm(tail_cogs[resid1] - tail_cogs[resid2])
                    d2 = np.linalg.norm(head_endpoints[resid1] - head_endpoints[resid2])
                    if d2 > d1:
                        angle2 = 0. - angle
                    else:
                        angle2 = angle

                    angles.append(angle)
                    angles2.append(angle2)
    return angles, angles2

def worker_unpack(args):
    '''allows multiple arguments to be passed to functions as a single argument'''
    return pair_angles(*args)

def mpi_wrapper(structure, trajectory, start=0, stop=-1, step=1, processes=-1, cutoff=6, costheta=False):
    '''preprocess trajectory into chunks which are threaded over with the multiprocessing module'''
    if processes == 1:
        angles, angles_dir = pair_angles(structure, trajectory, start=start, stop=stop, step=step)

    else:
        u = mda.Universe(structure, trajectory)
        if stop == -1:
            stop = u.trajectory.n_frames
        start, stop, step = u.trajectory.check_slice_indices(
            start, stop, step)
        # num_steps = np.ceil((stop - start)  / float(step))
        num_steps = len(range(start, stop, step))
        if num_steps == 0:
            raise ValueError("Selected time interval does not exist in the trajectory")
        if processes == -1:
            processes = mp.cpu_count() - 1
        times = np.linspace(start, stop, processes * 2 + 1).astype(int)
        process_bins =  [[] for x in times[1:]]
        #chunks = list(group(times))
        all_steps = np.arange(start, stop+step, step)
        bin_inds = np.digitize(all_steps, times[1:])
        bin_inds[-1] = bin_inds[-1] - 1
        for frame, bin_ind in zip(all_steps, bin_inds):
            #print(bin_ind, len(process_bins))
            process_bins[bin_ind].append(frame)
        chunks = [[np.min(x), np.max(x)+step] for x in process_bins]
        chunks[-1][1] = u.trajectory.n_frames
        p = mp.Pool(processes)
        angles_all = p.map(worker_unpack, [[structure, trajectory, x[0], x[1], step, cutoff] for x in chunks])
        angles = np.array([y for x in angles_all for y in x[0] ])
        angles_dir = np.array([y for x in angles_all for y in x[1] ])
    return angles, angles_dir

def make_histogram(data):

    bin_factor = 2.
    bin_num = int(360. * bin_factor)
    n, bins = np.histogram(data, bins=bin_num, range=(-180, 180), density=True)
    return bins[1:], n

def write_ouput(bins, probs, out_name="pair_prob_dist.dat"):
    with open(out_name, "w") as out:
        template = "{0:^12.56}{1:^12.6f}\n"
        for angle, prob in zip(bins, probs):
            out.write(template.format(angle, prob))

def write_angles(angles, out_name="pair_prob_angles.dat"):
    with open(out_name, "w") as out:
        template = "{0:^12d}{1:^12.6f}\n"
        for i, angle in enumerate(angles):
            out.write(template.format(i, angle))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--structure", required=True, help="input structure file *.gro or *.pdb")
    parser.add_argument("-f", "--trajectory", required=True, help="trajectory file *.xtc or *.trr")
    parser.add_argument("-b", "--begin", type=int, default=0, help="starting step for analysis")
    parser.add_argument("-e", "--end", type=int, default=-1, help="final step for analysis")
    parser.add_argument("-skip", type=int, default=1, help="analysis carried out every nth frame")
    parser.add_argument("-nt", "--num_threads", type=int, default=1, help="Number of threads to use.")
    parser.add_argument("-c", "--cutoff", type=float, default=0.6, help="cutoff for aontigen cahins to be defined "
                                                                        "as a pair in nm")
    parser.add_argument_group("--costheta", action='store_true', help="use to give all result in terms of cos(theta), "
                                                                      "where theta is an angle in radians")
    parser.add_argument("-o", "--out_prefix", default="pair_prob_dist", help="prefix of output *.dat file")
    args = parser.parse_args()

    angles, angles_dir = mpi_wrapper(args.structure, args.trajectory, start=args.begin, stop=args.end, step=args.skip,
                                     processes=args.num_threads, cutoff=args.cutoff * 10.)

    bins_dir, prob_dir = make_histogram(angles_dir)
    bins, prob = make_histogram(angles)

    write_ouput(bins, prob, out_name=args.out_prefix + ".dat")
    write_ouput(bins_dir, prob_dir, out_name=args.out_prefix + "_dir.dat")
    write_angles(angles, out_name="pair_prob_angles.dat")
    write_angles(angles_dir, out_name="pair_prob_angles_dir.dat")
