from __future__ import print_function, division
import numpy as np
import shlex
import MDAnalysis as mda
import argparse
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import get_leaflets

def load_data(filename):
    lipids_data = OrderedDict({})
    lipids_select = {}
    colours = []
    with open(filename, "r") as inp:
        line = inp.readline()
        membrane_select = line.strip()
        for line in inp:
            name, lipid_file, lipid_select, colour = shlex.split(line)
            lipids_data[name] = read_tseries(lipid_file)
            lipids_select[name] = lipid_select
            colours.append(colour)
    return lipids_data, lipids_select, membrane_select, colours



def read_tseries(filename):
    with open(filename, "r") as inp:
        lines = inp.readlines()
        data = []
        for line in lines:
            try:
                values = map(float, line.split())
                data.append(values)
            except ValueError:
                continue
    data = np.array(data)
    r = data[:, 0]
    average = data[:, 1]
    if data.shape[1] != 2:
        print("incorrect input format")
        raise SyntaxError
    return r, average


if __name__ == '__main__':
    parser = argparse.ArgumentParser("calculate/plot depletion/enrichment indices", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="Input file where first column contains the name of the lipid and 2nd column the "
                                      "filename containing the output of 'gmx select' for said lipid.")
    parser.add_argument("structure", help="*.gro/pdb file of system which must contain all of the lipids of interest")
    parser.add_argument("--n-blocks", dest='n_blocks', type=int, default=None, help="number of blocks to average data over")
    parser.add_argument("--graph-name", dest="graph", default="DE_index_blocks.pdf", help="output filename of graph")
    parser.add_argument("--block-output", dest="out_block", default="DE_index_blocks.dat", help="Name of output file for block averaged results.")
    parser.add_argument("--output", dest="output", default="DE_index.dat",
                        help="Name of output file for timeseries of enrichment values for each lipid.")
    parser.add_argument("--dpi", type=int, default=300, help="dpi of image")
    parser.add_argument("--leaflet", choices=["all", "upper", "lower"], default="all", help="Membrane Leaflet to use for analysis."
                                                                                            "Note previous gmx select results must also use the leaflet selected here"
                                                                                            "or your results will be wrong.")
    args = parser.parse_args()

    lipids_data, lipids_select, membrane_select, colours = load_data(args.input)
    u = mda.Universe(args.structure)

    atoms = u.atoms
    if args.leaflet != "all":
        leaflets = get_leaflets.determine_leaflets(u)
        atoms = leaflets[args.leaflet]

    total_bulk_count = len(atoms.select_atoms(membrane_select).residues)

    enrichment = dict.fromkeys(lipids_data)

    total_local_count = None
    for key in lipids_data:
        if total_local_count is None:
            total_local_count = lipids_data[key][1]
        else:
            total_local_count = total_local_count + lipids_data[key][1]


    for key in lipids_data:
        try:
            residue_count = len(atoms.select_atoms(lipids_select[key]).residues)
            ratio_local = lipids_data[key][1] / total_local_count
            print(total_local_count)
            ratio_local = np.nan_to_num(ratio_local)
            print(ratio_local)
            ratio_bulk = residue_count / total_bulk_count
            enrichment[key] = ratio_local / ratio_bulk
        except ZeroDivisionError:
            raise Exception("No lipids are present in this leaflet")

    with open(args.output, "w") as out:
        template_header = "".join(["{" + str(x) + ":^15s}" for x in range(len(colours) + 1)])
        template_results = template_header.replace("s", ".4f")
        print("#" + template_header.format("time (ns)", *lipids_data.keys()), file=out)
        times = lipids_data.values()[0][0] # all time series should be the same length
        size = len(times)
        array = np.array(enrichment.values()).T
        for i in range(size):

            print(template_results.format(times[i] / 1000, *array[i]), file=out)

    mpl.rcParams.update({'font.size': 22})

    total_averages = []
    errors = []
    standard_errors = []
    if args.n_blocks is not None:
        with open(args.out_block, "w") as out:
            template_header = "{0:^15s}{1:^15s}{2:^15s}"
            template_results = template_header.replace("s", ".4f")
            print("#" + template_header.format("total_averages", "error", "standard_error"), file=out)
            for key in lipids_data:
                sections = np.array_split(enrichment[key], args.n_blocks)
                averages = np.array([np.mean(section) for section in sections])
                total_average = np.mean(averages)
                error = np.std(averages)
                se = error / np.sqrt(len(averages))
                print(template_results.format(total_average, error, se), file=out)

                total_averages.append(total_average)
                errors.append(error)
                standard_errors.append(se)

        fig, ax = plt.subplots()
        ax.bar(lipids_data.keys(), total_averages, yerr=standard_errors, color=colours)
        ax.set_ylabel("D-E index")
        fig.tight_layout()
        fig.savefig(args.graph, dpi=args.dpi)









