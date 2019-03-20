#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import utils

# STEP 1: parse the arguments
options = arguments.read_args()
input_files = arguments.get_input_files(options.input)
stoichiometry = None
if options.stoichiometry is not None:
    stoichiometry = arguments.parse_stoichiometry(options.stoichiometry)
if options.verbose:
    utils.options = options
    sys.stderr.write("Input correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write(file+"\n")

# STEP2: Get possible structures for Macrocomplex construction and skip others
chains = {}
chain_index = 1
for file in input_files:
    for chain in utils.get_chains(file):
        chains[str(chain_index)+"_"+str(chain.id)] = chain
    chain_index += 1
similar_chains = utils.get_similar_chains(chains)
(chains, similar_chains) = utils.remove_useless_chains(chains, similar_chains)

print("Chains: \n"+ chains)
print("Similar Chains:\n"+similar_chains)
# STEP3: Check for stoichiometry requirements


# STEP4: Begin Macrocomplex reconstruction!
def construct_complex():
    test = utils.are_clashing(chains['1_C'], chains['1_D'])
    (test2, rmsd) = utils.superimpose(chains['1_C'], chains['1_D'])
    return complex
