from __future__ import print_function, division
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.distances
import argparse


def get_tilte2e(atoms):
    '''measures the tilt of a vector between the two points supplied wrt the z axis'''
    z_axis = np.array([0., 0., 1.])
    z_axis_norm = np.linalg.norm(z_axis)

    protein_axis = atoms[1]- atoms[0]
    tilt = np.arccos( np.dot(z_axis, protein_axis) / (z_axis_norm * np.linalg.norm(protein_axis) ) )
    tilt *= 0.5 * 360. / np.pi
    return tilt


def calc_tilt(structure, trajectory, start=0, stop=-1, step=1, dr=25., rmax=50.):
    '''function that determines the tilt of the oantigen chain in CG lps relative to
    the z axis wrt to distance from a protein'''
    u = mda.Universe(structure, trajectory)
    if stop == -1:
        stop = u.trajectory.n_frames
    start, stop, step = u.trajectory.check_slice_indices(
        start, stop, step)
    num_steps = len(range(start, stop, step))
    if num_steps == 0:
        raise ValueError("Selected time interval does not exist in the trajectory")

    bin_width = dr  # angstromns
    rmin = 0.
    r_cut = rmax  # maximum binning distance
    bin_values = np.arange(rmin, r_cut, bin_width) + bin_width / 2.
    bins = [[] for x in bin_values]
    tilt_angles = []
    oantigens = u.select_atoms("resname OANT").residues
    protein = u.select_atoms("name BB SC1 SC2 SC3 SC4")
    lipids_base = u.select_atoms("resname OANT and not name O* S*").residues
    num_res = len(lipids_base)

    for frame in u.trajectory[start:stop:step]:
        frame_tilts = []
        cogs_lipids = np.zeros((num_res, 3), dtype=np.float32)
        for i in range(num_res):
            cogs_lipids[i, :] = lipids_base[i].atoms.center_of_geometry(pbc=False)
        distances = mda.analysis.distances.distance_array(protein.positions, cogs_lipids, box=u.dimensions,
                                                          backend="OpenMP")
        min_molecule_dist = np.min(distances, axis=0)
        bin_inds = (np.floor((min_molecule_dist - rmin) / bin_width)).astype(int)

        for i, oantigen in enumerate(oantigens):
            sugar = [oantigen.atoms.select_atoms("name O51 to O54").center_of_geometry(),
                     oantigen.atoms.select_atoms("name O97 to O99").center_of_geometry()]

            tilt = get_tilte2e(sugar)
            bin_ind = bin_inds[i]
            if bin_ind < len(bins):
                if bin_ind < 0:
                    bin_ind = 0
                bins[bin_ind].append(tilt)
            frame_tilts.append(tilt)
        tilt_angles.append(frame_tilts)
    tilt_angles = np.array(tilt_angles)
    r = []
    errors = []
    averages = []
    population = []
    # print(bins)
    new_bins = []
    for i in range(len(bins)):
        if len(bins[i]) != 0:
            new_bins.append(bins[i])
            population.append(len(bins[i]))
            averages.append(np.mean(bins[i]))
            r.append(bin_values[i])
            errors.append(np.std(bins[i]))

    return r, averages, errors, population, tilt_angles, new_bins

def make_histogram(data):

    bin_factor = 2
    bin_num = int(360 * bin_factor)
    n, bins = np.histogram(data, bins=bin_num, range=(-180, 180), density=True)
    return bins[1:], n

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--structure", required=True, help="input structure file *.gro or *.pdb")
    parser.add_argument("-f", "--trajectory", required=True, help="trajectory file *.xtc or *.trr")
    parser.add_argument("-b", "--begin", type=int, default=0, help="starting step for analysis")
    parser.add_argument("-e", "--end", type=int, default=-1, help="final step for analysis")
    parser.add_argument("-skip", type=int, default=1, help="analysis carried out every nth frame")
    parser.add_argument("-dr", "--bin_width", type=float, default=2.5, help="bin widths of histograms in nm")
    parser.add_argument("-rmax", "--rmax", type=float, default=5.0, help="bin cutoff in nm")
    parser.add_argument("-o", "--out_prefix", default="otilt", help="prefix of output *.dat file")
    args = parser.parse_args()

    r, averages, errors, population, tilt_angles, new_bins = calc_tilt(args.structure, args.trajectory, start=args.begin,
                                                                       stop=args.end, step=args.skip,
                                                                       dr=args.bin_width * 10., rmax=args.rmax * 10.)

    prefix = args.out_prefix
    for i, distance_range in enumerate(new_bins):
        bins, probs = make_histogram(distance_range)
        with open("{0}_b{1}_dr{2}_hist.dat".format(prefix, i, int(args.bin_width)), "w") as out:
            template = "{0:^12.56}{1:^12.6f}\n"
            for angle, prob in zip(bins, probs):
                out.write(template.format(angle, prob))
        with open("{0}_b{1}_dr{2}_angles.dat".format(prefix, i, int(args.bin_width)), "w") as out:
            template = "{0:^12d}{1:^12.6f}\n"
            for j, angle in enumerate(distance_range):
                out.write(template.format(j, angle))

    with open("{0}_dr{1}_dist.dat".format(prefix, int(args.bin_width)), "w") as out:
        layout = "{0:^10.3f}{1:^15.5f}{2:^15.5f}{3:^15d}\n"
        out.write("{0:^10s}{1:^15s}{2:^15s}{3:^15s}\n".format("r / Ang", "Thickness", "Error", "Population"))
        for i in range(len(r)):
            out.write(layout.format(r[i], averages[i], errors[i], population[i]))

