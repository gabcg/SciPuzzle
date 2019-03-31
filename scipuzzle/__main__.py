#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import interface
import utils
import copy
import random
import messages as msg

# STEP 1: parse the arguments
options = arguments.read_args()
if options.gui:
    options = interface.gui()

files = arguments.get_input_files(options.input)
nc = options.nc
if options.verbose:
    utils.options = options
    msg.write_welcoming(files)

# Step 2: get possible structures for macrocomplex construction and skip others
if options.resume:
    (chains, pairs, similar_chains, structures) = utils.resume(options)
else:
    (chains, pairs, similar_chains, structures) = utils.get_information(files,
                                                                        options)

# STEP3: Macrocomplex reconstruction!
complexes_found = []
temp_complex = 0
if options.verbose:
    sys.stderr.write("\n# Beginning to construct the complex\n\n")


def construct_complex(current_complex_real, similar_chains, nc,
                      structures, used_pairs_real=[], clashing_real=[],
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
                          similar_chains, nc,
                          structures, used_pairs, clashing, old_complex)
        return
    # OTHER ROUNDS
    else:
        for chain_in_cc in utils.get_chain_ids_from_structure(current_complex):
            ps = utils.get_possible_structures(chain_in_cc,
                                               similar_chains, structures,
                                               used_pairs, clashing)
            # Saves the current complex if it reached the max number of chains
            if len(ps) == 0 and options.nc is None:
                if options.verbose:
                    msg.complex_built_no_nc()
                complexes_found.append(copy.deepcopy(current_complex))
                return
            # CASE UNDER DEVELOPMENT
            # If code is executed as shown in tutorial, this will not be used.
            # Only for the special case described in under development section.
            # This code will not be executed if number of chains is not given
            # The purpose is to refresh all the input parameters of the
            # recursive function so that the same pair is allowed to be added
            # multiple time to the same complex
            elif len(ps) == 0:
                if not utils.complex_differ(current_complex, old_complex):
                    return
                clashing = []
                old_complex = copy.deepcopy(current_complex)
                construct_complex(current_complex, similar_chains,
                                  nc, structures, [],
                                  clashing, old_complex)
                return

            for similar_chain_id in ps:
                # get the new chain to be added to the structure
                structure_id = ps[similar_chain_id]
                structure_to_superimpose = structures[structure_id]
                other = [tuple_id for tuple_id in structure_id if tuple_id != similar_chain_id][0]
                if options.verbose:
                    msg.trying_superimpose(other, structure_id)
                # calculate the new coordinates accordin to the superimposition
                matrix = utils.superimpose_chains(utils.get_chain(current_complex,
                                                  chain_in_cc),
                                                  utils.get_chain(structure_to_superimpose,
                                                  similar_chain_id))
                matrix.apply(utils.get_chain(structure_to_superimpose, other))
                # adds the chain to the current complex
                chain_to_add = utils.get_chain(structure_to_superimpose, other)
                if not utils.are_clashing(current_complex, chain_to_add):
                    if options.verbose:
                        sys.stderr.write("Not clashing -- adding chain! \n")
                        sys.stderr.write("-----------------------------\n")
                    utils.add_chain(current_complex, chain_to_add)
                    used_pairs.append(structure_id)

                else:
                    clashing.append(structure_id)
                # Checks if finished!
                if current_complex is not None:
                    # if the current complex reached the max num of chains
                    # save the current complex
                    if nc is not None and utils.complex_is_ready(current_complex,nc):
                        if options.verbose:
                            sys.stderr.write("\n\nComplex built!! :) \n")
                        complexes_found.append(copy.deepcopy(current_complex))
                        return
                    else:
                        # call the recursion again
                        construct_complex(current_complex, similar_chains,
                                          nc, structures, used_pairs,
                                          clashing, old_complex)
                        # if one suitable model is found - stop!
                        if len(complexes_found) == 1:
                            return


macrocomplex = construct_complex(None, similar_chains, nc,
                                 structures)

# Step 6: write output file
index_file = 0
for complex in complexes_found:
    index_file += 1
    outname = "results/"+str(options.output)+str(index_file)+".cif"
    utils.write_structure_into_file(complex, outname, "mmcif")


# Step 7 (optional): open models in Chimera
if options.open_chimera:
    utils.open_in_chimera(options)
