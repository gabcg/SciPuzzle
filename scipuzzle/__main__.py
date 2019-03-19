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
    sys.err("Input correctly parsed.\nFiles used as input:")
    for file in input_files:
        sys.err(file)

# STEP2: Get possible structures for Macrocomplex construction and skip others
chains = {}
chain_index = 1
for file in input_files:
    for chain in utils.get_chains(file):
        chains[str(chain_index)+"_"+str(chain.id)] = chain
    chain_index += 1
similar_chains = utils.get_similar_chains(chains)
(chains, similar_chains) = utils.remove_useless_chains(chains, similar_chains)

# STEP3: Check for stoichiometry requirements
chains = []
if not utils.stoichiometry_is_possible(stoichiometry, chains):
    # TODO: change with trow exception!
    exit()

# STEP4: Begin Macrocomplex reconstruction!
