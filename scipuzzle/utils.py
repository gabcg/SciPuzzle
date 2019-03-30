#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import Bio.PDB as pdb
from Bio import pairwise2
import itertools
import exceptions
import copy
import __main__ as main
import pickle
import os


def resume(options):
    """
    Returns a list with data needed for running the program by reading the
    binary files created by pickle.
    """
    prefix = 'resume/' + options.input
    chains = pickle.load(open(prefix + "_chains.p", "rb"))
    pairs = pickle.load(open(prefix + "_pairs.p", "rb"))
    similar_chains = pickle.load(open(prefix + "_similar_chains.p", "rb"))
    structures = pickle.load(open(prefix + "_structures.p", "rb"))

    if options.verbose:
        print("The following structures have been recovered:")
        print("Chains: \n" + str(chains))
        print("Paired chains: \n"+str(pairs))
        print("Similar Chains:\n"+str(similar_chains))
        print("structures: \n" + str(structures))

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
        structure = get_structure(file, remove_het=True)
        for chain in get_chain_from_structure(structure, remove_het=True):
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

    if options.verbose:
        print("Chains: \n" + str(chains))
        print("Paired chains: \n"+str(pairs))
        print("Similar Chains:\n"+str(similar_chains))
        print("structures: \n" + str(structures))

    # Save everything to binary files to be able to resume.
    if not os.path.exists('resume'):
        os.makedirs('resume')
    prefix = 'resume/' + options.input
    chains_backup = open(prefix + "_chains.p", "wb")
    pairs_backup = open(prefix + "_pairs.p", "wb")
    similar_chains_backup = open(prefix + "_similar_chains.p", "wb")
    structures_backup = open(prefix + "_structures.p", "wb")
    pickle.dump(chains, chains_backup)
    pickle.dump(pairs, pairs_backup)
    pickle.dump(similar_chains, similar_chains_backup)
    pickle.dump(structures, structures_backup)

    return (chains, pairs, similar_chains, structures)


def chain_to_fasta(chain):
    """
    Extracts the fasta sequence from a PDB file and returns a String
    containing the extracted sequence.
    """
    ppb = pdb.PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()


def get_structure(input_file, remove_het=False):
    """
    Extracts the structure from a pdb file.
    Returns the structure.
    """
    parser = pdb.PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)
    if not remove_het:
        return structure
    else:
        for model in structure:
            for chain in model:
                heteroatoms = list(filter(lambda x: x.id[0] != " ",
                                          chain.get_residues()))
                for heteroatom in heteroatoms:
                    chain.detach_child(heteroatom.id)
        return structure


def get_chains(input_file):
    """
    Extracts the chains from a pdb file and removes the heteroatoms from it.
    Returns a list containing n chain objects where n is the number of chains
    contained in the input file.
    """
    chains = []
    structure = get_structure(input_file)
    for model in structure:
        for chain in model:
            heteroatoms = list(filter(lambda x: x.id[0] != " ",
                                      chain.get_residues()))
            for heteroatom in heteroatoms:
                chain.detach_child(heteroatom.id)
            chains.append(chain)
    return chains




def get_similar_chains(chains, sequence_identity_threshold=0.95):
    """
    Compute which chains are similar to which ones given a dictionary having
    chain ids as keys and chain objects as values.
    Returns a dictionary containing informations about which chains are similar
    to which chains, according to the sequence_identity_threshold.
    Specifically the dictionary has as keys all chain ids and as values all the
    chains that are similar to the key chain.
    """
    print("Testing get similar chains")
    similar_chains = {}
    for chain_1, chain_2 in itertools.combinations(chains, 2):
        print("Chain one : " + str(chain_1))
        print("Chain two : " + str(chain_2))
        alignments = pairwise2.align.globalxx(chain_to_fasta(chains[chain_1]),
                                              chain_to_fasta(chains[chain_2]))
        pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
        if pairwise_sequence_identity >= sequence_identity_threshold:
            (superimposed, rmsd) = superimpose_chains(chains[chain_1], chains[chain_2])
            print(rmsd)
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


def are_clashing(chain_one, chain_two, max_clashes=50, contact_distance=1.0):
    """
    Compares the CA atoms of two chains and checks for clashes according to the
    contact distance.
    Returns a boolean.
    """
    atoms_one = [atom for atom in chain_one.get_atoms() if atom.get_id() == 'CA']
    atoms_two = [atom for atom in chain_two.get_atoms() if atom.get_id() == 'CA']
    ns = pdb.NeighborSearch(atoms_one)
    clashes = 0
    for atom_two in atoms_two:
        for atom in ns.search(atom_two.get_coord(), contact_distance, 'A'):
            clashes += 1
            if main.options.verbose:
                sys.stderr.write("Clash Found!\n")
            if clashes == max_clashes:
                return True
    return False


def get_chain_from_structure(structure, remove_het = False):
    chains = []
    for model in structure:
        for chain in model:
            if not remove_het:
                chains.append(chain)
            else:
                heteroatoms = list(filter(lambda x: x.id[0] != " ",
                                          chain.get_residues()))
                for heteroatom in heteroatoms:
                    chain.detach_child(heteroatom.id)
                chains.append(chain)
    return chains

def get_chain(structure, chain_id):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain

def get_possible_structures(chain_in_current_complex, similar_chains, structures):
    possible_structures = {}
    if chain_in_current_complex in similar_chains:
        for sim_chain in similar_chains[chain_in_current_complex]:
            for tuple_key in structures:
                if sim_chain in tuple_key:
                    possible_structures[sim_chain] = (tuple_key)
    return possible_structures


def get_chains_in_complex(used_pairs):
    chains = []
    for pair_in_current_complex in used_pairs:
        for chain_in_current_complex in pair_in_current_complex:
            chains.append(chain_in_current_complex)
    return chains


def remove_chain(structure, chain_id):
    for model in structure:
        model.detach_child(chain_id)
    return structure


def superimpose_chains_test(chain_one_real, chain_two_real):
    """
    Superimposes two structures and returns the superimposed structure and the rmsd
    of the superimposition.
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


def print_chain_in_structure(structure):
    if structure is not None:
        for chain in get_chain_from_structure(structure):
            print("Chain id: " + str(chain.id) + " --> " + str(chain))


def superimpose_chains(chain_one_real, chain_two_real):
    """
    Superimposes two structures and returns the superimposed structure and the rmsd
    of the superimposition.
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
    super_imposer.apply(list(chain_two.get_atoms()))
    return (chain_two, super_imposer.rms)




def superimpose(structure_one_real, structure_two_real):
    """
    Superimposes two structures and returns the superimposed structure and the rmsd
    of the superimposition.
    """
    structure_one = copy.deepcopy(structure_one_real)
    structure_two = copy.deepcopy(structure_two_real)
    super_imposer = pdb.Superimposer()
    atoms_one = list(structure_one.get_atoms())
    atoms_two = list(structure_two.get_atoms())
    # Fix lengths so that they are the same
    min_len = min(len(atoms_one), len(atoms_two))
    atoms_one = atoms_one[:min_len]
    atoms_two = atoms_two[:min_len]
    super_imposer.set_atoms(atoms_one, atoms_two)
    super_imposer.apply(list(structure_two[0].get_atoms()))
    return (structure_two, super_imposer.rms)


def get_all_similar_pairs(pair, similar_chains, structures):

    similar_pairs = []

    return similar_pairs


def superimpose_old(chain_one, chain_two):
    """
    Superimposes two chains and returns the superimposed chain and the rmsd
    of the superimposition.
    """
    super_imposer = pdb.Superimposer()
    atoms_one = list(chain_one.get_atoms())
    atoms_two = list(chain_two.get_atoms())
    # Fix lengths so that they are the same
    min_len = min(len(atoms_one), len(atoms_two))
    atoms_one = atoms_one[:min_len]
    atoms_two = atoms_two[:min_len]
    super_imposer.set_atoms(atoms_one, atoms_two)
    super_imposer.apply(chain_two.get_atoms())
    return (chain_two, super_imposer.rms)


def stoichiometry_is_not_ended(stoichiometry, current_chains_in_complex):
    """
    This function looks if it is possible to construct a complex with the
    files and the stoichiometry given by the user.
    """
    number_real_chains = sum(stoichiometry.values())
    if current_chains_in_complex < number_real_chains:
        return True
    else:
        return False


def write_structure_into_mmcif(structure, name):
    io = pdb.MMCIFIO()
    io.set_structure(structure)
    io.save(name)


def write_structure_into_pdb(structure, name):
    io = pdb.PDBIO()
    io.set_structure(structure)
    io.save(name)


def write_structure_into_pdb_old(chains, name):
    io = pdb.PDBIO()
    s = pdb.Structure.Structure('test')
    i = 1
    for chain in chains:
        s.add(pdb.Model.Model(i))
        s[i].add(chain)
        i += 1
    io.set_structure(s)
    io.save(name)


def complex_fits_stoich(complex, stoichiometry):
    if len(get_chain_from_structure(complex)) == sum(stoichiometry.values()):
        return True
    else:
        return False
