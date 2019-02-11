from __future__ import print_function, division
import numpy as np
import MDAnalysis as mda
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pylab as pl


def get_enrichment(filename):
    template = {"xwidth": float, "ywidth": float,
                    "xlabel": str, "ylabel": str,
                    "legend": str}
    metadata = {}
    with open(filename) as inp:
        data = []
        for line in inp:
            splitted = line[:-1].split()
            if line[0] == "@":
                key = splitted[0][1:]
                value = template[key](" ".join(splitted[1:]))
                metadata[key] = value
            else:
                try:
                    values = map(float, line.split())
                    data.append(values)
                except ValueError:
                    pass

    data = np.array(data)
    try:
        metadata["xwidth"]
    except KeyError:
        data = np.loadtxt(args.dens)
        xgrid = data[0]
        ygrid = data[:, 0]
        data = np.delete(data, 0, 1)
        data = np.delete(data, 0, 0)

        metadata["xwidth"] = np.max(xgrid)
        metadata["ywidth"] = np.max(ygrid)

    zav = np.mean(data)
    print("Average density is: {0:<8.3f} nm^2".format(zav))
    print("Maximum density is: {0:<8.3f} nm^2".format(np.max(data)))
    #print(np.unravel_index(np.argmax(data), data.shape))
    znew = (data/zav - 1.) * 100
    return znew, metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Calculates and plots enrichment maps. Requires the 2D density map ouput from either gmx densmap or"
                                     "the density script by N. Castillo et.al. (see http://perso.ibcp.fr/luca.monticelli/tools/index.html)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("dens", help="2D density map .dat")
    #parser.add_argument("-od", "--output_data", default="enrich.dat", help="name of output")
    parser.add_argument("-og", "--output_graph", default="enrich_graph.pdf", help="name of output graph")
    parser.add_argument("--dpi", type=int, default=300, help="dpi of image")
    parser.add_argument("--zlim", type=float, nargs=2, default=(None, None), metavar=("zmin", "zmax"), help="Picture limits for colour scale.")
    parser.add_argument("--plot_protein", action='store_true', help="Plot protein position on graph")
    parser.add_argument("-s", "--structure", default=None, help="gro file containing protein coordinates")
    args = parser.parse_args()
    if args.plot_protein == True and args.structure is None:
        parser.error("if '--plot_protein' option is passed a structure must be provided to '-s'")

    mpl.rcParams.update({'font.size': 14})
    znew, metadata = get_enrichment(args.dens)
    zlim = args.zlim
    if zlim == (None, None):
        zlim = (np.min(znew), np.max(znew))

    extent = [0, metadata["xwidth"], 0, metadata["ywidth"]]

    print("Extent is x-dim {0:6.2f} and y-dim {1:6.2f}".format(metadata["xwidth"], metadata["ywidth"]))

    fill = plt.imshow(znew.transpose(), extent=extent,
                      interpolation='none',
                      vmin=zlim[0], vmax=zlim[1], origin='lower', cmap=cm.jet)
    plt.contour(znew.transpose(), colors='k', extent=extent,
                origin="lower")

    try:
        metadata["xlabel"]
    except KeyError:
        metadata["xlabel"] = "X (nm)"
        metadata["ylabel"] = "Y (nm)"

    plt.xlabel(metadata["xlabel"])
    plt.ylabel(metadata["ylabel"])

    plt.colorbar(fill, extend='both', label="Enrichment percentage")
    plt.xlim(0, metadata["xwidth"])
    plt.ylim(0, metadata["ywidth"])

    if args.plot_protein:
        u = mda.Universe(args.structure)
        protein = u.select_atoms("protein")
        coords = np.array([residue.atoms.center_of_geometry()[:2] for residue in protein.residues]) / 10.
        plt.plot(coords[:, 0], coords[:, 1], 'o', markerfacecolor="None",
         markeredgecolor='white', markeredgewidth=2.5, markersize=0.8,alpha=0.8)



    plt.savefig(args.output_graph, dpi=args.dpi)