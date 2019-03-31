#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import Bio.PDB as pdb
from Bio import pairwise2
import itertools
import copy
import __main__ as main
import pickle


# ---------------------------
# -- Resume realated methods
# ---------------------------
def resume(options):
    """
    Returns a list with data needed for running the program by reading the
    binary files created by pickle.
    """
    print("option_parsed input")
    print(main.options.input)
    print(options.input)
    prefix = 'resume/' + options.input.strip('/').split('/')[-1]
    chains = pickle.load(open(prefix + "_chains.p", "rb"))
    pairs = pickle.load(open(prefix + "_pairs.p", "rb"))
    similar_chains = pickle.load(open(prefix + "_similar_chains.p", "rb"))
    structures = pickle.load(open(prefix + "_structures.p", "rb"))

    if options.verbose:
        sys.stderr.write("The following structures have been recovered:\n")
        sys.stderr.write("\tChains: %s\n" % len(chains))
        sys.stderr.write("\tPaired chains: %s\n" % len(pairs))
        sys.stderr.write("\tSimilar Chains: %s\n" % len(similar_chains))
        sys.stderr.write("\tStructures: %s\n\n" % len(structures))

    return (chains, pairs, similar_chains, structures)


def get_information(input_files, options):
    """
    Gets the possible structures for the macrocomplex construction.

    It gets each pair of chains from each file, constructing a chain ID. Once
    completed, the following structures are created:

    - A dictionary with the ID as key and the chain as value.
    - A list with sublists that contain IDs by pairs.

    After that, similar chains are gotten from the dictionary, generating a
    dictionary with key: ID and values: similar chains' IDs.

    At the end, the unneeded chains are removed.
    """
    chains = {}
    pairs = []
    structures = {}
    chain_index = 1
    for file in input_files:
        paired_chains = []
        structure = get_structure_from_file(file, remove_het=True)
        for chain in get_chains_from_structure(structure, remove_het=True):
            chain_id = str(chain_index)+"_"+str(chain.id)
            chains[chain_id] = chain
            chain.id = chain_id
            paired_chains.append(chain_id)
        pairs.append(paired_chains)
        structures[tuple(paired_chains)] = structure
        chain_index += 1
    similar_chains = get_similar_chains(chains)
    (chains, similar_chains, pairs) = remove_useless_chains(chains,
                                                            similar_chains,
                                                            pairs)

    # Save everything to binary files to be able to resume.
    if not os.path.exists('resume'):
        os.makedirs('resume')
    prefix = 'resume/' + options.input.strip('/').split('/')[-1]
    chains_backup = open(prefix + "_chains.p", "wb")
    pairs_backup = open(prefix + "_pairs.p", "wb")
    similar_chains_backup = open(prefix + "_similar_chains.p", "wb")
    structures_backup = open(prefix + "_structures.p", "wb")
    pickle.dump(chains, chains_backup)
    pickle.dump(pairs, pairs_backup)
    pickle.dump(similar_chains, similar_chains_backup)
    pickle.dump(structures, structures_backup)

    if options.verbose:
        sys.stderr.write("The analysis of the input results in:\n")
        sys.stderr.write("\tChains: %s\n" % len(chains))
        sys.stderr.write("\tPaired chains: %s\n" % len(pairs))
        sys.stderr.write("\tSimilar Chains: %s\n" % len(similar_chains))
        sys.stderr.write("\tStructures: %s\n\n" % len(structures))
        sys.stderr.write("Resume files have been saved in %s/resume\n"
                         % os.getcwd())

    return (chains, pairs, similar_chains, structures)

# ---------------------------
#  Utils functions
# ---------------------------


def chain_to_fasta(chain):
    """
    Extracts the fasta sequence from a PDB file and returns a string
    containing the extracted sequence.
    """
    ppb = pdb.PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()


def get_structure_from_file(input_file, remove_het=False):
    """
    Extracts the structure from a pdb file.
    Returns a Structure object.
    If parameter remove_het is set to True, heteroatoms are removed.
    """
    parser = pdb.PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)
    if not remove_het:
        return structure
    else:
        for model in structure:
            for chain in model:
                remove_heteroatoms(chain)
        return structure


def get_chains_from_structure(structure, remove_het=False):
    """
    Creates a list of chains out of a given structure.
    If remove_het is set to true, the heteroatoms are removed from the chain.
    Returns the list of chains.
    """
    chains = []
    for model in structure:
        for chain in model:
            if not remove_het:
                chains.append(chain)
            else:
                remove_heteroatoms(chain)
                chains.append(chain)
    return chains


def remove_heteroatoms(chain):
    """Removes the heteroatoms of a given chain."""
    heteroatoms = list(filter(lambda x: x.id[0] != " ",
                              chain.get_residues()))
    for heteroatom in heteroatoms:
        chain.detach_child(heteroatom.id)

##???
def get_chain_permissive(structure, chain_id):
    """
    Extract a specific chain from a structure given a specified chain id.
    Returns a chain object.
    """
    for chain in get_chains_from_structure(structure):
        if chain.id == chain_id:
            return chain
        elif chain_id in chain.id:
            print("000000000000000000  COMING IN THIS ")
            return chain


def get_chain(structure, chain_id):
    """
    Extract a specific chain from a structure given a specified chain id.
    Returns a chain object.
    """
    for chain in get_chains_from_structure(structure):
        if chain.id == chain_id:
            return chain


def get_chain_ids_from_structure(structure):
    """
    Extracts all chain ids from a given sructure.
    Returns a list.
    """
    chains_ids = []
    for model in structure:
        for c in model:
            chains_ids.append(c.id)
    return chains_ids


def add_chain(structure, chain_real):
    """
    Add a given chain to a structure.
    Modifies the chain Id if necessary.
    """
    chain = copy.deepcopy(chain_real)
    chains_ids = get_chain_ids_from_structure(structure)
    while chain.id in chains_ids:
        chain.id = chain.id+"-"+chain.id
    structure[0].add(chain)


# ---------------------------------
# Macrocomplex building functions
# ---------------------------------
def get_similar_chains(chains, sequence_identity_threshold=0.95):
    """
    Compute which chains are similar to which ones given a dictionary having
    chain ID's as keys and chain objects as values.
    Returns a dictionary containing information about which chains are similar
    to which chains, according to the sequence_identity_threshold.
    Specifically the dictionary has as keys all chain ID's and as values all
    the chains that are similar to the key chain.
    """

    similar_chains = {}
    if main.options.verbose:
        sys.stderr.write("\nComputing similar chains ... \n")
    for chain_1, chain_2 in itertools.combinations(chains, 2):
        alignments = pairwise2.align.globalxx(chain_to_fasta(chains[chain_1]),
                                              chain_to_fasta(chains[chain_2]))
        pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
        if pairwise_sequence_identity >= sequence_identity_threshold:
            (superimposed, rmsd) = superimpose_chains(chains[chain_1],
                                                      chains[chain_2])
            if rmsd < 0.05:
                if chain_1 not in similar_chains:
                    similar_chains[chain_1] = []
                if chain_2 not in similar_chains:
                    similar_chains[chain_2] = []
                similar_chains[chain_1].append(chain_2)
                similar_chains[chain_2].append(chain_1)
    return similar_chains


def remove_useless_chains(chains, similar_chains, pairs):
    """
    Removes chains from the chains list that neither have any similar chain
    among the other ones nor they are paired with one chain that has
    similar chains among the others.
    """
    # Collect chains to be removed
    remove = []
    for chain in similar_chains:
        prefix = chain.split("_")[0]
        relevant_similars = list(filter(lambda chain:
                                        not chain.startswith(prefix),
                                        similar_chains[chain]))
        if len(relevant_similars) == 0:
            remove.append(chain)
    # Remove chains from similar chains and chains
    for k in remove:
        chain_id = similar_chains[k][0]
        for pair_list in pairs:
            if chain_id in pair_list:
                pair_list.remove(chain_id)
        pairs = list(filter(None, pairs))
        del similar_chains[k]
        del chains[k]

    return (chains, similar_chains, pairs)


def are_clashing(struct, chain, max_clashes=50, contact_distance=1.0):
    """
    Compares the CA atoms of two chains and checks for clashes according to the
    contact distance and a maximum number of clashes.
    Returns a boolean.
    """
    atoms_one = [atom for atom in struct.get_atoms() if atom.get_id() == 'CA']
    atoms_two = [atom for atom in chain.get_atoms() if atom.get_id() == 'CA']
    ns = pdb.NeighborSearch(atoms_one)
    clashes = 0
    for atom_two in atoms_two:
        for atom in ns.search(atom_two.get_coord(), contact_distance, 'A'):
            clashes += 1
            if clashes == max_clashes:
                if main.options.verbose:
                    sys.stderr.write("Clash Found!\n")
                return True
    return False


def superimpose_chains(chain_one_real, chain_two_real):
    """
    Superimposes two structuresor chains and returns the superimposer object.
    This will be applied to the chains that need to be rotated.
    """
    chain_one = copy.deepcopy(chain_one_real)
    chain_two = copy.deepcopy(chain_two_real)
    super_imposer = pdb.Superimposer()
    atoms_one = sorted(list(chain_one.get_atoms()))
    atoms_two = sorted(list(chain_two.get_atoms()))
    # Fix lengths so that they are the same
    min_len = min(len(atoms_one), len(atoms_two))
    atoms_one = atoms_one[:min_len]
    atoms_two = atoms_two[:min_len]
    super_imposer.set_atoms(atoms_one, atoms_two)
    return super_imposer


def get_possible_structures(chain_in_current_complex, similar_chains,
                            structures, used_pairs, clashing):
    """
    Selects the structures that have a possibility to be added to the current
    complex given informations about similarity of the chains and previous
    attempts.
    """
    possible_structures = {}
    if chain_in_current_complex in similar_chains:
        for sim_chain in similar_chains[chain_in_current_complex]:
            for tuple_key in structures:
                if sim_chain in tuple_key:
                    print(tuple_key)
                    print(used_pairs)
                    if tuple_key not in used_pairs:
                        if tuple_key not in clashing:
                            possible_structures[sim_chain] = (tuple_key)
    return possible_structures


def complex_is_ready(complex, stoichiometry):
    """
    Checks if the number of chains in the complex is the same as the maximum
    number of chains desired.
    """
    if len(get_chains_from_structure(complex)) == sum(stoichiometry.values()):
        return True
    else:
        return False


#?????
def complex_differ(structure_one, structure_two):
    chains_one_len = len(get_chain_ids_from_structure(structure_one))
    chains_two_len = len(get_chain_ids_from_structure(structure_two))
    if chains_one_len != chains_two_len:
        return True
    return False


# ----------------------------
# Writing and showing outputs
# ----------------------------
def write_structure_into_file(structure, name, format):
    """
    Writes the strcuture into a file. The file can be either a pdb or a mmcif.
    Format needs to be either "pdb" or "mmcif" depending on the desired output
    file.
    """
    if format == "pdb":
        io = pdb.PDBIO()
    elif format == "mmcif":
        io = pdb.MMCIFIO()
    io.set_structure(structure)
    io.save(name)


def open_in_chimera(directory, options):
    """
    Opens all the models in Chimera if the user specified it. Works with Linux
    and macOS.
    """
    for file in os.listdir(directory):
        if file.endswith('.pdb'):
            if options.verbose:
                sys.stderr.write('Opening model %s in Chimera' % file)
            if sys.platform == 'darwin':
                os.system('/Applications/Chimera.app/Contents/MacOS/chimera'
                          + file)
            else:
                os.system('chimera' + file)


def print_chain_in_structure(structure):
    """
    Prints all chains of a given structure in a human-readable way.
    Method used in development.
    """
    if structure is not None:
        for chain in get_chains_from_structure(structure):
            print("Chain id: " + str(chain.id) + " --> " + str(chain))
