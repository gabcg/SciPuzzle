#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Bio.PDB import *

parser = PDBParser()
structure = parser.get_structure('X', '../example_inputs/pair_his3_sc_XA.pdb')
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom)

ppb = CaPPBuilder()
for pp in ppb.build_peptides(structure):
    print(pp.get_sequence())
