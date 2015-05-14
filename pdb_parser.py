"""
This file provides an easy interface for reading from PDB files.
"""

import itertools
import logging
import os


def iCodeOrder(icname):
    """
    Shamelessly coppied from Professor Scott's code. This converts the icode to
    a number so we can do some actual manipulation with it.
    """

    ico=12
    if (icname==""):  ico=0
    if (icname=="A"): ico=1
    if (icname=="B"): ico=2
    if (icname=="C"): ico=3
    if (icname=="D"): ico=4
    if (icname=="E"): ico=5
    if (icname=="F"): ico=6
    if (icname=="G"): ico=7
    if (icname=="H"): ico=8
    if (icname=="I"): ico=9
    if (icname=="J"): ico=10
    if (icname=="X"): ico=11
    return ico

def compound_id_to_float(compound_id):
    """
    This will convert a sequence ID to a float (which allows us to actually compare
    them). Most of the time it just cast's it as an int. However, this may not
    work becasue of insertion codes. We simply give these codes a fractional id.
    """

    try:
        return float(int(compound_id))
    except ValueError:
        return float(int(compound_id[:-1])) + iCodeOrder(compound_id[-1]) / 12.0

class PDBData(object):
    """
    This is the container for all data pertaining to a PDB file. At it's core
    it is nothing other than a collection of atoms and secondary structure.
    However, this includes logic to do useful things with these atoms.
    """

    def __init__(self):
        self.atoms = []
        self.helixes = []
        self.sheets = []

    def add_atom(self, atom):
        """
        Add an atom to this set of data.
        """

        self.atoms.append(atom)

    def add_helix(self, chain_id, start_id, end_id):
        """
        Register an alpha helix with this data. Note that both start and end id
        should be in the form of a string and include insertion codes
        """

        self.helixes.append((chain_id, start_id, end_id))

    def add_sheet(self, chain_id, start_id, end_id):
        """
        Register a beta sheet with this data.
        """

        self.sheets.append((chain_id, start_id, end_id))

    def get_atoms(self, atom_names=None, in_residue=None):
        """
        Get all the atoms in this data that meet requirements. Possible requirements
        are membership in a list of atoms or membership in a specific residue.
        """

        result = self.atoms
        if atom_names is not None:
            result = [x for x in result if x.atom in atom_names]
        if in_residue is not None:
            result = [x for x in result if x.compound in in_residue]
        return result

    def get_phoso_sites(self):
        """
        Returns a list of phosphorolayted residues.
        """

        return [x for x in self.get_compounds() if x.name in ("PTR", "SEP", "TPO")]

    def get_compounds(self):
        """
        Breaks up the atoms into sets of compounds.
        """

        compounds = []
        for sequence_id, compound in itertools.groupby(self.atoms, lambda x: x.sequence_id):
            compounds.append(PDBCompound(sequence_id, list(compound)))
        return compounds

    def get_adjacent_compounds(self):
        """
        Returns a list of tuples of compounds which are adjacent to one another.
        This means that they are on the same chain and come right after one another
        in the sequence of atoms.

        There are a couple of gotchas here:

            1) Compounds on different chains are not adjacent. So we need to first
               group by chain.
            2) There are sometimes missing residues. In this case we just throw away
               the data.
            3) We have to deal with insertion codes. We assume that if there is an
               insertion code that the atoms are adjacent.
        """


        all_pairs = []
        for _, chain in itertools.groupby(self.get_compounds(), lambda x: x.chain_id):
            last_compound = None
            for compound in chain:
                if last_compound is not None:
                    try:
                        if int(last_compound.compound_id) + 1 == int(compound.compound_id):
                            all_pairs.append((last_compound, compound))
                    except ValueError: # We assume that insertion codes mean that they are adjacent
                        all_pairs.append((last_compound, compound))
                last_compound = compound
        return all_pairs

    def get_helixes(self):
        """
        Returns a list of lists each of which represents an alpah helix.
        """

        compounds = self.get_compounds()
        result = []
        for chain_id, start, end in self.helixes:
            helix = []
            for compound in compounds:
                if (compound.chain_id == chain_id and
                    compound_id_to_float(compound.compound_id) >= compound_id_to_float(start) and
                    compound_id_to_float(compound.compound_id) <= compound_id_to_float(end)):
                    helix.append(compound)
            result.append(helix)
        return result

class PDBCompound(object):
    """
    This represents a collection of atoms in the PDB. Will generally be a residue
    group but could possibly be something else.
    """

    def __init__(self, compound_id, atoms):
        self.compound_id = compound_id
        self.atoms = atoms
        # Copy some data off of the first atom (assumed to be true for the rest
        # of the atoms)
        self.chain_id = atoms[0].chain_id
        self.name = atoms[0].compound

    def get_atom(self, atom_name):
        """
        Get the first instance of atom_name in this compound. If no such atom
        is found this will return None
        """

        for atom in self.atoms:
            if atom.atom == atom_name:
                return atom
        return None

    def __str__(self):
        return "Compound {}".format(self.compound_id)

    def __repr__(self):
        return "Compound {}".format(self.compound_id)

class PDBAtom(object):
    """
    Represents a single atom in a PDB file.
    """

    def __init__(self):
        self.position = (0, 0, 0)
        self.compound = None
        self.atom = None
        self.sequence_id = None
        self.alt_id = None
        self.chain_id = 0

    def __str__(self):
        return "{} in compound {} at position {} (Sequence ID: {}, Chain ID: {})".format(self.atom, self.compound, self.position, self.sequence_id, self.chain_id)

    def __repr__(self):
        return "{} in compound {} at position {} (Sequence ID: {}, Chain ID: {})".format(self.atom, self.compound, self.position, self.sequence_id, self.chain_id)



def find_pdb_files(folder):
    """
    Finds all PDB files in a given folder.
    """

    return [os.path.join(folder, file_name) for file_name in os.listdir(folder) if file_name[-3:] == "pdb"]

def parse_pdb(file_name):
    """
    This is the function that should be called from external files.
    """

    return parse_pdb_text(file_name)

def parse_pdb_text(file_name):
    """
    Parse a text formatted PDB file.
    """

    data = PDBData()

    atom_count = 0
    processed_model = False
    with open(file_name) as pdb_file:
        for line in pdb_file:
            if "ATOM" == line[0:4] or "HETATM" == line[0:6]:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                compound = line[17:20].strip()
                atom = line[12:16].strip()
                sequence_id = line[22:26].strip()
                insertion_code = line[26:27].strip()
                alt_id = line[16:17].strip()
                chain_id = line[21:22].strip()

                pdb_atom = PDBAtom()
                pdb_atom.position = (x, y, z)
                pdb_atom.compound = compound
                pdb_atom.atom = atom
                pdb_atom.sequence_id = sequence_id + insertion_code
                pdb_atom.alt_id = alt_id
                pdb_atom.chain_id = chain_id
                data.add_atom(pdb_atom)
                atom_count += 1
            if "MODEL" == line[0:5]:
                if processed_model:
                    break
                processed_model = True
            if "HELIX" == line[0:5]:
                chain_id = line[19].strip()
                start_seq = (line[21:25].strip())
                start_insertion = line[25].strip()
                end_seq = (line[33:37].strip())
                end_insertion = line[37].strip()
                data.add_helix(chain_id, start_seq + start_insertion, end_seq + end_insertion)
            if "SHEET" == line[0:5]:
                chain_id = line[21].strip()
                start_seq = (line[22:26].strip())
                start_insertion = line[26].strip()
                end_seq = (line[33:37].strip())
                end_insertion = line[37].strip()
                data.add_sheet(chain_id, start_seq + start_insertion, end_seq +  end_insertion)

    logging.info("Loaded %s (%d atoms)", file_name, atom_count)
    return data
