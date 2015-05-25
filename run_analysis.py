#!/usr/bin/env python
"""
Run the analysis on a given set of PDB files.
"""

import argparse
import logging
import os

from matplotlib import pyplot

from pdb_parser import distance, find_pdb_files, parse_pdb
from wrappa import get_dehydrons

RESIDUES_OF_INTEREST = ("TYR", "SER", "THR", "PTR", "SEP", "TPO")

def count_residues(pdb_data):
    """
    Count the number of different residues that occur in the PDB data.
    """

    return {x: len(pdb_data.get_compounds([x])) for x in RESIDUES_OF_INTEREST}

def phosphorylation_in_desolvation(pdb_data, sites, dehydrons):
    """
    This will identify all sites that fall within a desolvation domain of a dehydron.
    This will return a list of tuples where the first element represents the phosphorylation
    site and the second represents the dehydron which it is close to.
    """

    results = []
    for dehydron in dehydrons:
        residue1 = pdb_data.get_residue_by_id(dehydron[0][0], dehydron[0][1])
        residue2 = pdb_data.get_residue_by_id(dehydron[1][0], dehydron[1][1])

        position1 = residue1.get_atom("CA").position
        position2 = residue2.get_atom("CA").position
        for site in sites:
            phospo_position = site.get_atom("CA").position
            distance1 = distance(position1, phospo_position)
            distance2 = distance(position2, phospo_position)
            if distance1 <= 6.5 and distance2 <= 6.5:
                results.append((site, (residue1, residue2)))
    return results

def run_analysis(pdb_name, data_directory):
    """
    This will run the analysis on a single pdb file. The PDB name should be a raw
    PDB name (i.e. 1a81H) and the data directory should contain the following files:

        1) PDB_NAME.pdb
        2) PDB_NAME_wrappers.txt
        3) PDB_NAME_bonds.txt

    The analysis consists of the following substeps (all of which will be returned)

        1) Count each occurance of a residue of interest (see above). Returns a dictionary
           of the counts of each residue of interest.
        2) Counts the number of phosphorylation sites. Returns an integer containing
           the number of phosphorylation sites.
        3) Generate a list of phosprohylation sites that are within a desolvation sphere
           of a dehydron. Returns this list of tuples.
    """

    logging.info("Running analysis for PDB %s", pdb_name)

    # Load all of the data
    pdb_data = parse_pdb(os.path.join(data_directory, pdb_name + ".pdb"))
    dehydrons = get_dehydrons(os.path.join(data_directory, pdb_name + "_wrappers.txt"))

    # Get the sites of interest
    phospo_sites = pdb_data.get_phoso_sites()
    #non_phoso_sites = pdb_data.get_compounds(names=["TYR", "THR", "SER"])

    return count_residues(pdb_data), len(phospo_sites), phosphorylation_in_desolvation(pdb_data, phospo_sites, dehydrons)


def main():
    """
    Run the main functionality for this script.
    """

    parser = argparse.ArgumentParser(description="Analyze phosphorylation sites and dehydrons")
    parser.add_argument("data_directory", help="The directory where the data is stored")
    parser.add_argument("--limit", default=None, type=int,
                        help="Set a limit on the number of files to process")
    args = parser.parse_args()

    all_residue_counts = {x: 0 for x in RESIDUES_OF_INTEREST}
    num_phospo_sites = []
    all_pairs = {}

    for i, pdb_file in enumerate(find_pdb_files(args.data_directory)):
        if args.limit is not None and i >= args.limit:
            break

        residue_counts, num_sites, pairs = run_analysis(os.path.split(pdb_file[:-4])[1], args.data_directory)

        # Merge in the data
        for x in RESIDUES_OF_INTEREST:
            all_residue_counts[x] += residue_counts[x]
        num_phospo_sites.append(num_sites)
        all_pairs[os.path.split(pdb_file[:-4])] = pairs

    # Display the final data
    print(all_pairs)
    print(sum(len(x) for x in all_pairs.values()))
    #print(all_residue_counts)
    #pyplot.hist(num_phospo_sites)
    #pyplot.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
