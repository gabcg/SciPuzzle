#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import utils
import exceptions
import Bio.PDB as pdb
import copy
import random

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
pairs = []
structures = {}
chain_index = 1
for file in input_files:
    paired_chains = []
    structure = utils.get_structure(file)
    for chain in utils.get_chains(file):
        chain_id = str(chain_index)+"_"+str(chain.id)
        chains[chain_id] = chain
        paired_chains.append(chain_id)
    pairs.append(paired_chains)
    structures[tuple(paired_chains)] = structure
    chain_index += 1
similar_chains = utils.get_similar_chains(chains)
(chains, similar_chains, pairs) = utils.remove_useless_chains(chains, similar_chains, pairs)

print("Chains: \n" + str(chains))
print("Paired chains: \n"+str(pairs))
print("Similar Chains:\n"+str(similar_chains))
print("Sto: " + str(stoichiometry))
print("structures: \n" + str(structures))

# STEP3: Check for stoichiometry requirements
if not utils.stoichiometry_is_possible(stoichiometry, chains, similar_chains):
    raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)


# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex, chains,
                      similar_chains, stoichiometry, pairs, structures, used_pairs):
    # pairs_left =
    # test = utils.are_clashing(chains['1_C'], chains['1_D'])
    # (test2, rmsd) = utils.superimpose(chains['1_C'], chains['1_D'])

    # current_complex is a list of chains
    # for the first round its going to be a random pair of chains.
    if current_complex is None:
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        used_pairs.append(random_choice_id)
        construct_complex(random_choice, chains,
                          similar_chains, stoichiometry, pairs, structures, used_pairs)
        return

    print(current_complex)
    print(used_pairs)
    print("Wasdlakd")

    # Proceed recursion in all possible directions (based on similar_chains)
    possible_structures = set()
    for pair_in_current_complex in used_pairs:
        for chain_in_current_complex in pair_in_current_complex:
            print("Comparing chain : " + str(chain_in_current_complex))
            if chain_in_current_complex in similar_chains:
                for sim_chain in similar_chains[chain_in_current_complex]:
                    for tuple_key in structures:
                        if sim_chain in tuple_key:
                            print("Adding: " + str(sim_chain))
                            possible_structures.add(tuple_key)
    print(possible_structures)

    for possible_structure in possible_structures:
        print(possible_structure)
        # to be improved
        to_be_deleted = []
        for possible_chain in possible_structure:
            for pair_complex in used_pairs:
                for chain_complex in pair_complex:
                    for sim_to_used in similar_chains[chain_complex]:
                        if possible_chain in sim_to_used:
                            to_be_deleted.append(possible_chain)
        print("To be deleted")
        print(to_be_deleted)
        (superimposed, rmsd) = utils.superimpose(current_complex, structures[possible_structure])

        print(superimposed)
        print("Superimposing!\n")
        print("Possible structure : "+ str(possible_structure))
        print("used_pairs: "+ str(used_pairs))
        utils.write_structure_into_pdb(superimposed, "lala.pdb")
        exit(0)
        # check if are clashing
        #Call recursion


    # Each time we add a new chain - update current_complex

    # construct_complex(current_complex, chains, similar_chains, stoichiometry)

    # if current_complex fulfills stoichiometry - save in list and return
    if utils.complex_fits_stoich(current_complex, stoichiometry):
        # Output: List of lists of chain. Each list of chain corresponds to a
        # constructed complex
        return [current_complex]


test_complex = construct_complex(None, chains, similar_chains,
                                 stoichiometry, pairs, structures, [])
# Step5: Filter the good ones

# Step6 : write output file
#utils.write_structure_into_pdb(test_complex[0], 'test.pdb')
