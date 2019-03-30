#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import interface
import utils
import exceptions
import copy
import random
import re

line = "-------------------------------------"


# STEP 1: parse the arguments
options = arguments.read_args()
if options.gui:
    options = interface.gui()

input_files = arguments.get_input_files(options.input)
stoichiometry = None
if options.stoichiometry is not None:
    stoichiometry = arguments.parse_stoichiometry(options.stoichiometry)
if options.verbose:
    utils.options = options
    sys.stderr.write("Input correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write("\t"+file+"\n")
    sys.stderr.write("\n")

# Step 2: get possible structures for macrocomplex construction and skip others
if options.resume:
    (chains, pairs, similar_chains, structures) = utils.resume(options)
else:
    (chains, pairs, similar_chains, structures) = utils.get_information(input_files, options)

# STEP3: Check for stoichiometry requirements
# if not utils.stoichiometry_is_possible(stoichiometry, chains, similar_chains):
#    raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)


complexes_found = []
# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex_real, chains,
                      similar_chains, stoichiometry, pairs,
                      structures, used_pairs_real):
    # bruteforce ending!
    current_complex = copy.deepcopy(current_complex_real)
    used_pairs = copy.deepcopy(used_pairs_real)

    # for the first round its going to be a random pair of chains.
    if current_complex is None:
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        print("First case: beginning with "+ str(random_choice_id))
        used_pairs.append(random_choice_id)
        construct_complex(random_choice, chains,
                          similar_chains, stoichiometry, pairs,
                          structures, used_pairs)
        return
    else:
        print("------------------------------------------")
        print("Beginning to try possible superimpositions")
        for chain_in_current_complex in utils.get_chains_in_complex(used_pairs):
            print("Comparing chain : " + str(chain_in_current_complex))
            ps = utils.get_possible_structures(chain_in_current_complex,
                                               similar_chains, structures, used_pairs)
            print("Possible structures : "+str(ps))
            if len(ps) == 0 and options.stoichiometry is None:
                complexes_found.append(current_complex)
                return
            for similar_chain_id in ps:
                structure_id = ps[similar_chain_id]
                structure_to_superimpose = structures[structure_id]
                other = [tuple_id for tuple_id in structure_id if tuple_id != similar_chain_id][0]
                print("similar_chain_id " + str(similar_chain_id))
                print("Other (to superimpose) : " + other)
                print("Superimposing")
                print("Structure to superimpose: "+ str(structure_to_superimpose))
                # Superimpose current complex with one of the possible structures.
                matrix = utils.superimpose_chains_test(utils.get_chain(current_complex,chain_in_current_complex)
                                                      ,utils.get_chain(structure_to_superimpose, similar_chain_id))
                matrix.apply(utils.get_chain(structure_to_superimpose,other))
                current_complex[0].add(utils.get_chain(structure_to_superimpose,other))
                used_pairs.append(structure_id)
                print("Used pairs : " + str(used_pairs))
                ## Check if finished!
                if current_complex is not None:
                    for chain in utils.get_chains_from_structure(current_complex):
                        print("Chain id: " + str(chain.id) + " --> " + str(chain))
                    if stoichiometry is not None and utils.complex_fits_stoich(current_complex,stoichiometry):
                        print("Recursion Finished!")
                        for chain in utils.get_chains_from_structure(current_complex):
                            print(chain)
                        complexes_found.append(current_complex)

                        return
                    else:

                        #elif --> add repeated!

                        construct_complex(current_complex, chains,
                                          similar_chains, stoichiometry, pairs,
                                          structures, used_pairs)
                        return

test_complex = construct_complex(None, chains, similar_chains,
                                 stoichiometry, pairs, structures, [])


# Step5: Filter the good ones

# Step6 : write output file
index_file = 0
for complex in complexes_found:
    index_file += 1
    utils.write_structure_into_file(complex, "results/output"+str(index_file)+".cif", "mmcif")

# Step 7

utils.write_structure_into_pdb(test_complex[0], 'test2.pdb')

# Step 7 (optional): open models in Chimera

# Need to store the outputs as a list
models = ['test.pdb']
if options.open_chimera:
    utils.open_in_chimera(models)

# if output is a directory:
# if options.open_chimera:
#     utils.open_in_chimera(options.output)
