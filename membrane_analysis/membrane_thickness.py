from __future__ import print_function
import numpy as np
import MDAnalysis as mda
import .get_leaflets
from scipy import spatial
import argparse


def get_volume_shells(initial_dr, r_max):
    V_o = 4. * np.pi * (initial_dr ** 3) / 3.
    bin_values = []
    value = 0.
    count = 1.
    while value < r_max:
        if count > 1:
            bin_values.append(value)
        value = (3 * V_o * count / (4. * np.pi)) ** (1. / 3.)
        # print(value)
        count += 1
    return np.array(bin_values)

def protein_overlap_vectorised(xy_wrapped, z_wrapped, tri, protein_pos, tri_offset):
    insides = tri.find_simplex(protein_pos[:, :2])
    tri_pos_xy = xy_wrapped[tri.simplices[insides]]
    tri_pos_z = z_wrapped[tri.simplices[insides]]
    tri_pos = np.dstack((tri_pos_xy, tri_pos_z[:, :, np.newaxis]))
    tri2_pos = tri_pos + tri_offset
    tri3_pos = tri_pos - tri_offset

    prism = np.array([tri2_pos, tri3_pos])
    p1 = prism[:, :, 0]
    p2 = prism[:, :, 1]
    p3 = prism[:, :, 2]
    p1p2 = p2 - p1
    p1p3 = p3 - p1
    plane_normal = np.cross(p1p2, p1p3, axis=2)
    dot_products = np.einsum('ijk,ijk->ij', plane_normal, protein_pos - p1)

    signs = np.sign(dot_products)
    signs = signs.T
    b = signs[:, 0] != signs[:, 1]
    inds = np.where(b == True)[0]
    return inds

def projected_thickness(projected_pos, positions, z_wrap):
    zi = np.fabs(z_wrap - projected_pos[2])
    x = positions[:, 0]
    y = positions[:, 1]
    A = y[0] * (zi[1] - zi[2]) + y[1] * (zi[2] - zi[0]) + y[2] * (zi[0] - zi[1])
    B = zi[0] * (x[1] - x[2]) + zi[1] * (x[2] - x[0]) + zi[2] * (x[0] - x[1])
    C = x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1])
    D = -A * x[0] - B * y[0] - C * zi[0]

    z_projected = - (A / C) * projected_pos[0] - (B / C) * projected_pos[1] - D / C
    return z_projected




class MembraneThickness(object):
    def __init__(self, structure, trajectory, protein_sel=None):
        self.layer_choices = ["upper", "lower"]
        self.mulitple_universes = False
        self.protein_sel = protein_sel
        if len(structure) == 1:
            self.u = mda.Universe(structure[0], trajectory[0])  # universe object
            self.layers = get_leaflets.determine_leaflets(self.u, phosphateSelection="name P*")
            self.upper_heads = self.layers["upper"].select_atoms("name P*")
            self.lower_heads = self.layers["lower"].select_atoms("name P*")
            if self.protein_sel is not None:
                self.protein = self.u.select_atoms(protein_sel)
        else:
            self.mulitple_universes = True
            self.universes = [mda.Universe(struct, traj) for struct, traj in zip(structure, trajectory)]

    def __repr__(self):
        return "<Thickness Analysis object>"



    def get_thickness(self, dx=3.0, start=0, stop=-1, step=1, rmax=40., dr=1.0, layers="lower", volume_normlised=False):
        if layers not in self.layer_choices:
            raise IOError("Incorrect input for 'layers' - choose from: {}".format(", ".join(self.layer_choices)))


        if not self.mulitple_universes:
            bin_width = dr  # angstromns
            rmin = 0.
            r_cut = rmax  # maximum binning distance
            if volume_normlised:
                bin_values = get_volume_shells(bin_width, rmax)
            else:
                bin_values = np.arange(bin_width, r_cut, bin_width)
            bins = [[] for x in bin_values]
            if stop == -1:
                stop = self.u.trajectory.n_frames
            start, stop, step = self.u.trajectory.check_slice_indices(
                start, stop, step)

            num_steps = len(range(start, stop, step))
            if num_steps == 0:
                raise ValueError("Selected time interval does not exist in the trajectory")

            thicknesses = []

            grid_num = np.ceil(self.u.dimensions[:2] / float(dx)).astype(int)
            grid = np.zeros((grid_num), dtype=np.float64)
            grid_count = np.zeros((grid_num), dtype=np.int32)

            shape = np.array(grid.shape)
            if layers == "upper":
                atoms = self.upper_heads
                surface_org = self.lower_heads
            else:
                atoms = self.lower_heads
                surface_org = self.upper_heads

            surface = surface_org
            if self.protein_sel is not None:
                # surface_prot = surface + protein
                selection_keys = self.protein_sel.split()

            tri_offset = np.array([0., 0., 2.0])
            boxxy = []
            times = []
            for frame in self.u.trajectory[start:stop:step]:
                boxxy.append(self.u.dimensions[:2])
                time = self.u.trajectory.time / 1000.
                times.append(time)
                surface = surface_org

                voxel_size = self.u.dimensions[:2] / grid_num.astype(float)

                thickness_all = []  # thickness_delaunay(upper_heads, lower_heads)

                box_x, box_y, box_z = self.u.dimensions[:3]
                xy = surface.positions[:, :2]
                z = surface.positions[:, 2]
                # print(xy.shape, z.shape)
                # quit()
                xy_wrapped = np.vstack((xy, xy + [box_x, 0], xy - [box_x, 0], xy + [0, box_y], xy - [0, box_y],
                                        xy + [box_x, box_y], xy - [box_x, box_y], xy + [box_x, -box_y],
                                        xy - [box_x, -box_y]))
                z_wrapped = np.hstack((z, z, z, z, z, z, z, z, z))
                # print(xy_wrapped.shape)
                # tri = spatial.Delaunay(surface.positions[:, :2])

                tri = spatial.Delaunay(surface.positions[:, :2])

                # print("membrane", tri.simplices)
                if self.protein_sel is not None:
                    inds = protein_overlap_vectorised(xy_wrapped, z_wrapped, tri, self.protein.positions, tri_offset)

                    surface = surface_org + self.protein[inds]
                    xy = surface.positions[:, :2]
                    z = surface.positions[:, 2]
                    # print(xy_wrapped)
                    # break
                    xy_wrapped = np.vstack((xy, xy + [box_x, 0], xy - [box_x, 0], xy + [0, box_y], xy - [0, box_y],
                                            xy + [box_x, box_y], xy - [box_x, box_y], xy + [box_x, -box_y],
                                            xy - [box_x, -box_y]))
                    z_wrapped = np.hstack((z, z, z, z, z, z, z, z, z))
                    n_org = surface.names
                    names = np.concatenate((n_org, n_org, n_org, n_org, n_org, n_org, n_org, n_org, n_org,
                                            n_org, n_org, n_org), axis=0)
                    tri = spatial.Delaunay(xy_wrapped, incremental=False)

                if self.protein_sel is not None:
                    distances = mda.analysis.distances.distance_array(self.protein.atoms.positions, atoms.positions,
                                                                      box=self.u.dimensions,
                                                                      backend="OpenMP")
                    min_molecule_dist = np.min(distances, axis=0)

                    bin_inds = np.digitize(min_molecule_dist, bin_values)

                thickness_all = []
                simplex_IDS = tri.find_simplex(atoms.positions[:, :2])
                for i, pos in enumerate(atoms.positions):
                    simplex_ID = simplex_IDS[i]
                    simplex = tri.simplices[simplex_ID]
                    tri_atoms = xy_wrapped[simplex]
                    if self.protein_sel is not None:
                        tri_names = list(names[simplex])
                        if any(x in tri_names for x in selection_keys):
                            counts = np.sum([tri_names.count(x) for x in selection_keys])
                            if counts > 1:
                                continue

                    tri_atoms_z = z_wrapped[simplex]
                    thickness = projected_thickness(pos, tri_atoms, tri_atoms_z)
                    if self.protein_sel is not None:
                        bin_ind = bin_inds[i]
                        if bin_ind < len(bins):
                            if bin_ind < 0:
                                bin_ind = 0

                            bins[bin_ind].append(thickness)


                    i, j = np.floor(pos[:2] / voxel_size).astype(int)
                    if i >= shape[0]:
                        i = shape[0] - 1
                    elif i < 0:
                        i = 0

                    if j >= shape[1]:
                        j = shape[1] - 1
                    elif j < 0:
                        j = 0

                    grid[i, j] += thickness
                    grid_count[i, j] += 1
                    thickness_all.append(thickness)

                thicknesses.append(thickness_all)
            averages = []
            r = []
            errors = []
            population = []
            if self.protein_sel is not None:
                if not volume_normlised:
                    bin_values -= bin_width / 2.

                for i in range(len(bins)):
                    if len(bins[i]) != 0:
                        population.append(len(bins[i]))
                        averages.append(np.mean(bins[i]))
                        r.append(bin_values[i])
                        errors.append(np.std(bins[i]))
            boxxy = np.array(boxxy)
            boxxy = np.mean(boxxy, axis=0)
            boxxy /= 10.
            times = np.array(times)
            results = {"timeseries": None, "grid": None, "radial": None}
            thick_aver = np.array([np.mean(x) for x in thicknesses]) / 10.
            thick_errors = np.array([np.std(x) for x in thicknesses]) / 10.
            grid = np.divide(grid, grid_count) / 10.
            results["timeseries"] = (times, thick_aver, thick_errors)
            results["grid"] = (grid, grid_count, boxxy)
            if self.protein_sel is not None:
                results["radial"] = (r, averages, errors, population)
            self.results = results
        else:
            grid_num = np.ceil(self.universes[0].dimensions[:2] / float(dx)).astype(int)
            grid = np.zeros((grid_num), dtype=np.float64)
            grid_count = np.zeros((grid_num), dtype=np.int32)

            shape = np.array(grid.shape)

            bin_width = dr  # angstromns
            rmin = 0.
            r_cut = rmax  # maximum binning distance
            if volume_normlised:
                bin_values = get_volume_shells(bin_width, rmax)
            else:
                bin_values = np.arange(bin_width, r_cut, bin_width)
            bins = [[] for x in bin_values]
            boxxy = []
            thicknesses_muni = []
            for universe in self.universes:
                self.u = universe
                self.layers = get_leaflets.determine_leaflets(self.u, phosphateSelection="name P*")
                self.upper_heads = self.layers["upper"].select_atoms("name P*")
                self.lower_heads = self.layers["lower"].select_atoms("name P*")
                if self.protein_sel is not None:
                    self.protein = self.u.select_atoms(self.protein_sel)
                thicknesses = []

                if stop == -1:
                    stop = len(self.u.trajectory) - 1
                start, stop, step = self.u.trajectory.check_slice_indices(
                    start, stop, step)
                num_steps = len(range(start, stop, step))
                if num_steps == 0:
                    raise ValueError("Selected time interval does not exist in the trajectory")


                if layers == "upper":
                    atoms = self.upper_heads
                    surface_org = self.lower_heads
                else:
                    atoms = self.lower_heads
                    surface_org = self.upper_heads

                surface = surface_org
                if self.protein_sel is not None:
                    # surface_prot = surface + protein
                    selection_keys = self.protein_sel.split()

                tri_offset = np.array([0., 0., 2.0])

                times = []
                for frame in self.u.trajectory[start:stop:step]:
                    boxxy.append(self.u.dimensions[:2])
                    time = self.u.trajectory.time / 1000.
                    times.append(time)
                    surface = surface_org

                    voxel_size = self.u.dimensions[:2] / grid_num.astype(float)

                    thickness_all = []  # thickness_delaunay(upper_heads, lower_heads)

                    box_x, box_y, box_z = self.u.dimensions[:3]
                    xy = surface.positions[:, :2]
                    z = surface.positions[:, 2]
                    xy_wrapped = np.vstack((xy, xy + [box_x, 0], xy - [box_x, 0], xy + [0, box_y], xy - [0, box_y],
                                            xy + [box_x, box_y], xy - [box_x, box_y], xy + [box_x, -box_y],
                                            xy - [box_x, -box_y]))
                    z_wrapped = np.hstack((z, z, z, z, z, z, z, z, z))

                    tri = spatial.Delaunay(surface.positions[:, :2])

                    if self.protein_sel is not None:
                        inds = protein_overlap_vectorised(xy_wrapped, z_wrapped, tri, self.protein.positions, tri_offset)

                        surface = surface_org + self.protein[inds]
                        xy = surface.positions[:, :2]
                        z = surface.positions[:, 2]

                        xy_wrapped = np.vstack((xy, xy + [box_x, 0], xy - [box_x, 0], xy + [0, box_y], xy - [0, box_y],
                                                xy + [box_x, box_y], xy - [box_x, box_y], xy + [box_x, -box_y],
                                                xy - [box_x, -box_y]))
                        z_wrapped = np.hstack((z, z, z, z, z, z, z, z, z))
                        n_org = surface.names
                        names = np.concatenate((n_org, n_org, n_org, n_org, n_org, n_org, n_org, n_org, n_org,
                                                n_org, n_org, n_org), axis=0)
                        tri = spatial.Delaunay(xy_wrapped, incremental=False)

                    if self.protein_sel is not None:
                        distances = mda.analysis.distances.distance_array(self.protein.atoms.positions, atoms.positions,
                                                                          box=self.u.dimensions,
                                                                          backend="OpenMP")
                        min_molecule_dist = np.min(distances, axis=0)
                        bin_inds = np.digitize(min_molecule_dist, bin_values)

                    thickness_all = []
                    simplex_IDS = tri.find_simplex(atoms.positions[:, :2])
                    for i, pos in enumerate(atoms.positions):
                        simplex_ID = simplex_IDS[i]
                        simplex = tri.simplices[simplex_ID]
                        tri_atoms = xy_wrapped[simplex]
                        if self.protein_sel is not None:
                            tri_names = list(names[simplex])
                            if any(x in tri_names for x in selection_keys):
                                counts = np.sum([tri_names.count(x) for x in selection_keys])
                                if counts > 1:
                                    continue

                        tri_atoms_z = z_wrapped[simplex]
                        thickness = projected_thickness(pos, tri_atoms, tri_atoms_z)
                        if self.protein_sel is not None:
                            bin_ind = bin_inds[i]
                            if bin_ind < len(bins):
                                if bin_ind < 0:
                                    bin_ind = 0

                                bins[bin_ind].append(thickness)

                        i, j = np.floor(pos[:2] / voxel_size).astype(int)
                        if i >= shape[0]:
                            i = shape[0] - 1
                        elif i < 0:
                            i = 0

                        if j >= shape[1]:
                            j = shape[1] - 1
                        elif j < 0:
                            j = 0

                        grid[i, j] += thickness
                        grid_count[i, j] += 1
                        thickness_all.append(thickness)

                    thicknesses.append(thickness_all)
                thicknesses_muni.append(thicknesses)

            thicknesses = [ t1 + t2 for t1, t2, in zip(thicknesses_muni[0], thicknesses_muni[1])]
            averages = []
            r = []
            errors = []
            population = []
            if self.protein_sel is not None:
                if not volume_normlised:
                    bin_values -= bin_width / 2.
                # bin_values = [x / 2. if i ==0  else x  - ((x - bin_values[i-1]) / 2.) for i, x in enumerate(bin_values)]
                for i in range(len(bins)):
                    if len(bins[i]) != 0:
                        # print(len(bins[i]))
                        population.append(len(bins[i]))
                        averages.append(np.mean(bins[i]))
                        r.append(bin_values[i])
                        errors.append(np.std(bins[i]))
            boxxy = np.array(boxxy)
            boxxy = np.mean(boxxy, axis=0)
            boxxy /= 10.
            times = np.array(times)
            results = {"timeseries": None, "grid": None, "radial": None}
            thick_aver = np.array([np.mean(x) for x in thicknesses]) / 10.
            thick_errors = np.array([np.std(x) for x in thicknesses]) / 10.
            grid = np.divide(grid, grid_count) / 10.
            results["timeseries"] = (times, thick_aver, thick_errors)
            results["grid"] = (grid, grid_count, boxxy)
            if self.protein_sel is not None:
                results["radial"] = (r, averages, errors, population)
            self.results = results


if __name__ == '__main__':
    desc="Script that calculates the membrane thickness as a function of time and as a 2D grid in the xy plane. " \
         "Method of calculation utilises Delaunay Triangulation as Described in the AplVoro paper"
    parser = argparse.ArgumentParser(desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--structure",required=True,  nargs='+',help="input structure files *.gro or *.pdb")
    parser.add_argument("-f", "--trajectory", required=True, nargs='+', help="trajectory files *.xtc or *.trr")
    parser.add_argument("-b", "--begin", type=int, default=0, help="starting step for analysis" )
    parser.add_argument("-e", "--end", type=int, default=-1, help="final step for analysis")
    parser.add_argument("-skip", type=int, default=1, help="analysis carried out every nth frame")
    parser.add_argument("-p", "--protein", default=None, help="MDAnalysis selection string for proteins")
    parser.add_argument("-dx", type=float, default=3, help="Approximate size of each voxel in the grid (Angstroms). "
                                                           "Note that the number of voxels in the grid is "
                                                           "fixed and determined from the size of the box of the first "
                                                           "frame to be analysed ")
    parser.add_argument("-dr", type=float, default=0.2, help="Bin width for radial thickness plots in nm.")
    parser.add_argument("-rmax", type=float, default=5., help="Cutoff for radial thickness in nm.")
    parser.add_argument("-o", "--output_prefix", default="thickness", help="Prefix for all output files.")
    parser.add_argument("-vnorm", action='store_true', help="bins have equal volumes")
    args = parser.parse_args()

    layers = {"upper": None, "lower": None}
    for layer in layers:
        analysis = MembraneThickness(args.structure, args.trajectory,
                                     protein_sel=args.protein)
        analysis.get_thickness(layers=layer, dx=args.dx, start=args.begin, stop=args.end, step=args.skip,
                               dr=args.dr * 10., rmax=args.rmax * 10., volume_normlised=args.vnorm)
        layers[layer] = analysis


    for layer in layers:
        obj = layers[layer]
        prefix = args.output_prefix

        results = obj.results["timeseries"]
        #print(results)
        with open("{0}_{1}.dat".format(prefix, layer), "w") as out:
            template = "{0:^15.3f}{1:^15.3f}{2:^15.3f}\n"
            out.write("@xlabel Time (ns)\n")
            out.write("@ylabel Thickness (nm)\n")
            for i in range(results[0].shape[0]):
                out.write(template.format(results[0][i], results[1][i], results[2][i]))

        if args.protein is not None:
            results = obj.results["radial"]
            with open("{0}_{1}_dist.dat".format(prefix, layer), "w") as out:
                layout = "{0:^10.3f}{1:^15.5f}{2:^15.5f}{3:^15d}\n"
                out.write("{0:^10s}{1:^15s}{2:^15s}{3:^15s}\n".format("r / Ang", "Thickness", "Error", "Population"))
                for i in range(len(results[0])):
                    out.write(layout.format(results[0][i], results[1][i], results[2][i], results[3][i]))

        results = obj.results["grid"]
        values, grid_count, boxxy = results
        with open("{0}_{1}_grid.dat".format(prefix, layer), "w") as out:


            out.write("@xwidth {0:^10.3f}\n".format(boxxy[0]))
            out.write("@ywidth {0:^10.3f}\n".format(boxxy[1]))
            out.write("@xlabel X (nm)\n")
            out.write("@ylabel Y (nm)\n")
            out.write("@legend Thickness map (nm)\n")
            for row in values:
                print(" ".join([str(cell) for cell in row]), file=out)

        with open("{0}_{1}_sampling.dat".format(prefix, layer), "w") as out:
            out.write("@xwidth {0:^10.3f}\n".format(boxxy[0]))
            out.write("@ywidth {0:^10.3f}\n".format(boxxy[1]))
            out.write("@xlabel X (nm)\n")
            out.write("@ylabel Y (nm)\n")
            out.write("@legend Bin count\n")
            for row in grid_count:
                print(" ".join([str(cell) for cell in row]), file=out)
