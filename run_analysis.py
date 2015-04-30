#!/usr/bin/env python
"""
Run the analysis on a given set of PDB files.
"""

from pdb_parser import find_pdb_files, parse_pdb

if __name__ == "__main__":
    data = parse_pdb("data/1a67.pdb")
