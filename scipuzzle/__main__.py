#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import utils
import exceptions
import Bio.PDB as pdb
import copy
import random
import re


line = "-------------------------------------"

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
    for chain in utils.get_chain_from_structure(structure):
        chain_id = str(chain_index)+"_"+str(chain.id)
        chains[chain_id] = chain
        chain.id = chain_id
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
print(line)
print(line)

# STEP3: Check for stoichiometry requirements
# if not utils.stoichiometry_is_possible(stoichiometry, chains, similar_chains):
#    raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)



# def test(pair_one, pair_two):
#     (superimposed, rmsd) = utils.superimpose(pair_one,pair_two)
#     #pair_one.add()
#     utils.write_structure_into_mmcif(pair_one, "test.cif")
#
# print("Here i am ")

print("going further ")
# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex, chains,
                      similar_chains, stoichiometry, pairs,
                      structures, used_pairs):
    # bruteforce ending!
    if current_complex is not None:
        if len(utils.get_chain_from_structure(current_complex)) == 3:
            print(line)
            print("Recursion Finished!")
            for chain in utils.get_chain_from_structure(current_complex):
                print(chain)
            utils.write_structure_into_mmcif(current_complex, "lala.cif")
            return [current_complex]
    # current_complex is a list of chains
    # for the first round its going to be a random pair of chains.
    if current_complex is None:
        print(line)
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        print("First case: beginning with "+ str(random_choice_id))
        used_pairs.append(random_choice_id)
        construct_complex(random_choice, chains,
                          similar_chains, stoichiometry, pairs,
                          structures, used_pairs)
        return

    print("------------------------------------------")
    print("Beginning to try possible superimpositions")
    for chain_in_current_complex in utils.get_chains_in_complex(used_pairs):
        print("Comparing chain : " + str(chain_in_current_complex))
        print("similar_chains: "+ str(similar_chains))
        ps = utils.get_possible_structures(chain_in_current_complex,
                                           similar_chains, structures)
        print("Possible structures : "+ str(ps))
        for ps_id in ps:
            structure_id = ps[ps_id]
            structure_to_superimpose = structures[structure_id]
            print("Superimposing")
            print("Structure to superimpose: "+ str(structure_to_superimpose))
            # Superimpose current complex with one of the possible structures.
            (superimposed, rmsd) = utils.superimpose(current_complex,structure_to_superimpose)
            # print("ps id " + str(ps_id))
            # print("chain in cu co "+ str(chain_in_current_complex))
            for chain in utils.get_chain_from_structure(superimposed):
                if chain.id != ps_id:
                    print("Trying with: "+ str(chain.id))
                    # check if are clashing
                    if structure_id not in used_pairs:
                        #if not utils.are_clashing(current_complex, chains[chain.id]):

                            current_complex[0].add(utils.get_chain(superimposed, chain.id))
                            used_pairs.append(structure_id)
                            print("no clash found!")
                            print("Adding: "+ str(chain.id))
                            print("Used pairs "+str(used_pairs))
                            construct_complex(current_complex, chains,
                                              similar_chains, stoichiometry, pairs,
                                              structures, used_pairs)
                        #else:
                            #print("Clash was found :(")
                    else:
                        print("already in used pairs")


test_complex = construct_complex(None, chains, similar_chains,
                                 stoichiometry, pairs, structures, [])
# Step5: Filter the good ones

# Step6 : write output file
#utils.write_structure_into_pdb(test_complex[0], 'test.pdb')
