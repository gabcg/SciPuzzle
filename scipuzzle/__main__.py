#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import interface
import utils
import exceptions
import copy
import random
import pickle
from global_vars import *
import messages as msg
# STEP 1: parse the arguments
options = arguments.read_args()
if options.gui:
    options = interface.gui()

input_files = arguments.get_input_files(options.input)
stoichiometry = arguments.parse_stoichiometry(options.stoichiometry)
if options.verbose:
    utils.options = options
    msg.write_welcoming(input_files)

# Step 2: get possible structures for macrocomplex construction and skip others
if options.resume:
    (chains, pairs, similar_chains, structures) = utils.resume(options)
else:
    (chains, pairs, similar_chains, structures) = utils.get_information(input_files, options)


complexes_found = []
temp_complex = 0
if options.verbose:
    sys.stderr.write("\n# Beginning to construct the complex\n\n")


# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex_real,
                      similar_chains, stoichiometry, structures,
                      used_pairs_real=[], clashing_real=[],
                      old_complex_real=[]):
    # Making copies for avoiding modifying objects concurrently.
    current_complex = copy.deepcopy(current_complex_real)
    used_pairs = copy.deepcopy(used_pairs_real)
    clashing = copy.deepcopy(clashing_real)
    old_complex = copy.deepcopy(old_complex_real)

    # FIRST ROUND: for the first round its going to be a random pair of chains.
    if current_complex is None:
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        used_pairs.append(random_choice_id)
        if options.verbose:
            msg.beginning(random_choice_id)
        construct_complex(random_choice,
                          similar_chains, stoichiometry,
                          structures, used_pairs, clashing, old_complex)
        return
    # OTHER ROUNDS
    else:
        for chain_in_cc in utils.get_chain_ids_from_structure(current_complex):
            ps = utils.get_possible_structures(chain_in_cc,
                                               similar_chains, structures,
                                               used_pairs, clashing)
            print(" ==== Possible structures : "+str(ps))
            if len(ps) == 0 and options.stoichiometry is None:
                print("\t\t\t Stoichiometry not given - all possible structures consumed")
                if options.verbose:
                    msg.complex_built_no_stoich()
                complexes_found.append(current_complex)
                return
            elif len(ps) == 0:
                if len(utils.get_chain_ids_from_structure(current_complex)) == len(utils.get_chain_ids_from_structure(old_complex)):
                    #complexes_found.append(current_complex)
                    print(" -- ")
                    utils.print_chain_in_structure(current_complex)
                    print(" -- ")
                    utils.print_chain_in_structure(old_complex)
                    print(" -- ")
                    print("\t\t\t Found same complex twice! ")
                    print("Already found!")
                    if current_complex is not None:
                        temp_complex = current_complex
                    return
                else:
                    print("\t\t\t Trying OTHER process of recursion adding ALL possible chains again ")
                    old_complex = current_complex

                    if temp_complex != 0:
                        if not utils.complex_differ(current_complex, temp_complex):
                            return
                    clashing = []
                    print("calling myself again")
                    construct_complex(current_complex,
                                      similar_chains, stoichiometry,
                                      structures, [], clashing, old_complex)
                    return

            for similar_chain_id in ps:
                structure_id = ps[similar_chain_id]
                structure_to_superimpose = structures[structure_id]
                other = [tuple_id for tuple_id in structure_id if tuple_id != similar_chain_id][0]
                if options.verbose:
                    msg.trying_superimpose(other, structure_id)
                matrix = utils.superimpose_chains_test(utils.get_chain(current_complex,chain_in_cc)
                                                      ,utils.get_chain_permissive(structure_to_superimpose, similar_chain_id))
                utils.print_chain_in_structure(structure_to_superimpose)
                print(used_pairs)
                matrix.apply(utils.get_chain(structure_to_superimpose,other))
                chain_to_add = utils.get_chain_permissive(structure_to_superimpose,other)
                print(chain_to_add)
                if not utils.are_clashing(current_complex,chain_to_add):
                    if options.verbose:
                        sys.stderr.write("Not clashing -- adding chain! \n")
                    utils.add_chain(current_complex, chain_to_add)
                    used_pairs.append(structure_id)
                else:
                    clashing.append(structure_id)
                ## Check if finished!
                if current_complex is not None:
                    if stoichiometry is not None and utils.complex_fits_stoich(current_complex,stoichiometry):
                        if options.verbose:
                            print("Complex fits stoich!")
                            sys.stderr.write("Complex built!! :) \n")
                        complexes_found.append(current_complex)
                        return
                    else:
                        construct_complex(current_complex,
                                          similar_chains, stoichiometry,
                                          structures, used_pairs, clashing, old_complex)
                        if len(complexes_found) == 1:
                            print("One complex was already found! Break!")
                            utils.print_chain_in_structure(complexes_found)
                            return


test_complex = construct_complex(None, similar_chains,
                                 stoichiometry, structures, [], [], [])
# Step5: Filter the good ones

# Step 6: write output file
index_file = 0
if len(complexes_found) == 0:
    complexes_found.append(temp_complex)

for complex in complexes_found:
    index_file += 1
    outname = "results/"+str(options.output)+str(index_file)+".cif"
    utils.write_structure_into_file(complex, outname, "mmcif")


# Step 7 (optional): open models in Chimera
if options.open_chimera:
    utils.open_in_chimera(options)
