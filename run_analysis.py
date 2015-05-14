#!/usr/bin/env python
"""
Run the analysis on a given set of PDB files.
"""

import sys

from matplotlib import pyplot

from pdb_parser import find_pdb_files, parse_pdb

if __name__ == "__main__":
    # of_interest = ("TYR", "SER", "THR", "PTR", "SEP", "TPO")
    # counts = {x: 0 for x in of_interest}
    counts = []
    for pdb_file in find_pdb_files(sys.argv[1]):
        data = parse_pdb(pdb_file)
        phospo_sites = data.get_phoso_sites()
        # for x in of_interest:
        #     counts[x] += len(data.get_compounds([x]))
        print("{} has {} phosphorylation sites (types: {})".format(pdb_file, len(phospo_sites), list(set([x.name for x in phospo_sites]))))
        counts.append(len(phospo_sites))
    # print counts
    pyplot.hist(counts, 15)
    pyplot.xlabel("Number of phosphorylation sites")
    pyplot.ylabel("Number of occurances")
    pyplot.savefig("histogram.png")
