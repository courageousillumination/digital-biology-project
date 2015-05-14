#!/usr/bin/env python
"""
Run the analysis on a given set of PDB files.
"""

import sys

from pdb_parser import find_pdb_files, parse_pdb

if __name__ == "__main__":

    for pdb_file in find_pdb_files(sys.argv[1]):
        data = parse_pdb(pdb_file)
        phospo_sites = data.get_phoso_sites()
        print("{} has {} phosphorylation sites (types: {})".format(pdb_file, len(phospo_sites), list(set([x.name for x in phospo_sites]))))
