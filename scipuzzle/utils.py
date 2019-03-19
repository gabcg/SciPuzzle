#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import Bio.PDB as pdb
from Bio import pairwise2
import itertools
import os


def chain_to_fasta(chain):
    ppb = pdb.PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()

# TODO: do not include heteroatoms
def get_chains(input_file, no_hetatm=False):
    chains = []
    parser = pdb.PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)
    for model in structure:
        for chain in model:
            chains.append(chain)
    return chains


def filter_pairs(pairs):
    filtered_pairs = pairs
    for chain_one, chain_two in itertools.combinations(filtered_pairs, 2):
        sequence_identity_threshold = 0.95
        # aligner = pdb.Align.PairwiseAligner()
        alignments = pairwise2.align.globalxx(chain_to_fasta(chain_one),
                                              chain_to_fasta(chain_two))
        pairwise_sequence_identity = alignments[0][2] / len(alignments[0][0])
        if pairwise_sequence_identity >= sequence_identity_threshold:
            print(pairwise_sequence_identity)
    return filtered_pairs

def stoichiometry_is_possible(stoichiometry,chains):
    return True
