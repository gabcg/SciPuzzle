#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import Bio.PDB as pdb
from Bio import pairwise2
import itertools
import exceptions
import __main__ as main


def chain_to_fasta(chain):
    """
    Extracts the fasta sequence from a PDB file and returns a String
    containing the extracted sequence.
    """
    ppb = pdb.PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()


def get_chains(input_file):
    """
    Extracts the chains from a pdb file and removes the heteroatoms from it.
    Returns a list containing n chain objects where n is the number of chains
    contained in the input file.
    """
    chains = []
    parser = pdb.PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)
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
    similar_chains = {}
    for chain_one, chain_two in itertools.combinations(chains, 2):
        alignments = pairwise2.align.globalxx(chain_to_fasta(chains[chain_one])
                                              ,chain_to_fasta(chains[chain_two]))
        pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
        if pairwise_sequence_identity >= sequence_identity_threshold:
            if chain_one not in similar_chains:
                similar_chains[chain_one] = []
            if chain_two not in similar_chains:
                similar_chains[chain_two] = []
            similar_chains[chain_one].append(chain_two)
            similar_chains[chain_two].append(chain_one)
    return similar_chains


def remove_useless_chains(chains, similar_chains):
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
        del similar_chains[k]
        del chains[k]
    return (chains, similar_chains)


def are_clashing(chain_one, chain_two, contact_distance=1.4):
    """
    Compares the CA atoms of two chains and checks for clashes according to the
    contact distance.
    Returns a boolean.
    """
    atoms_one = [atom for atom in chain_one.get_atoms() if atom.name == 'CA']
    atoms_two = [atom for atom in chain_two.get_atoms() if atom.name == 'CA']
    ns = pdb.NeighborSearch(atoms_one)
    for atom_two in atoms_two:
        if len(ns.search(atom_two.get_coord(), contact_distance)) != 0:
            if main.options.verbose:
                sys.stderr.write("Clash Found!\n")
            return True
    return False


def superimpose(chain_one, chain_two):
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


def stoichiometry_is_possible(stoichiometry, chains):
    if chains == 1:
        return True
    else:
        raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)
