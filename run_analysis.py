#!/usr/bin/env python
"""
Run the analysis on a given set of PDB files.
"""

import sys

from matplotlib import pyplot

from pdb_parser import distance, find_pdb_files, parse_pdb
from wrappa import get_dehydrons

if __name__ == "__main__":
    # of_interest = ("TYR", "SER", "THR", "PTR", "SEP", "TPO")
    # counts = {x: 0 for x in of_interest}
    # counts = []
    # for pdb_file in find_pdb_files(sys.argv[1]):
    #     data = parse_pdb(pdb_file)
    #     phospo_sites = data.get_phoso_sites()
        # for x in of_interest:
        #     counts[x] += len(data.get_compounds([x]))
        # print("{} has {} phosphorylation sites (types: {})".format(pdb_file, len(phospo_sites), list(set([x.name for x in phospo_sites]))))
        # counts.append(len(phospo_sites))
    # print counts
    # pyplot.hist(counts, 15)
    # pyplot.xlabel("Number of phosphorylation sites")
    # pyplot.ylabel("Number of occurances")
    # pyplot.savefig("histogram.png")
    dehydrons = get_dehydrons("wrappers.txt")
    data = parse_pdb("./1a81H.pdb")
    phospo_sites = data.get_phoso_sites()
    print dehydrons
    for dehydron in dehydrons:
        residue1 = data.get_residue_by_id(dehydron[0][0], dehydron[0][1])
        residue2 = data.get_residue_by_id(dehydron[1][0], dehydron[1][1])
        position1 = residue1.get_atom("CA").position
        position2 = residue2.get_atom("CA").position
        for phospo_site in phospo_sites:
            phospo_position = phospo_site.get_atom("CA").position
            distance1 = distance(position1, phospo_position)
            distance2 = distance(position2, phospo_position)
            if distance1 < 6.5 and distance2 < 6.5:
                print(phospo_site, dehydron, distance1, distance2)
    # print data, dehydrons
