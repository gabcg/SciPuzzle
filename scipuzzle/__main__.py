#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.PDB import *
import os
from os import path
import itertools
from Bio import pairwise2





verbose = True

input_pdb = "../example_inputs/pair_his3_sc_XA.pdb"
input_folder = "/Users/luisasantus/Documents/GitHub/SciPuzzle/example_inputs/"



def chain_to_fasta(chain):
    ppb=PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()

def get_input_files():
    for i in os.listdir(input_folder):
        #TODO: remove the second part before submitting!
        if i.endswith('.pdb') and  not "chain" in i:
            pdb_files.append(os.path.join(input_folder, i))
    return pdb_files


def collect_chains(input_file):
    chains = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)
    for model in structure:
        for chain in model:
            chains.append(chain)
    return chains[]

# A-B B-C C-D

# A, B, B , C, C, D

pairs = []
for file in get_input_files():






for chain_one, chain_two  in itertools.combinations(chains, 2):
    sequence_identity_threshold = 0.95
    aligner = Align.PairwiseAligner()
    alignments = pairwise2.align.globalxx(chain_to_fasta(chain_one), chain_to_fasta(chain_two))
    pairwise_sequence_identity = alignments[0][2] /len(alignments[0][0])
    if pairwise_sequence_identity >= sequence_identity_threshold:
        print(pairwise_sequence_identity)
