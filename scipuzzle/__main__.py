#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import utils
import exceptions
import copy
import random
import pickle
import global_vars

# STEP 1: parse the arguments
options = arguments.read_args()
input_files = arguments.get_input_files(options.input)
stoichiometry = None
if options.stoichiometry is not None:
    stoichiometry = arguments.parse_stoichiometry(options.stoichiometry)
if options.verbose:
    utils.options = options
    sys.stderr.write("\n\n-------------------------------\n")
    sys.stderr.write("--  Welcome to SciPuzzle! ;) --\n")
    sys.stderr.write("-------------------------------\n")
    sys.stderr.write("\n\nInput correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write(file+"\n")

chains = {}
pairs = []
structures = {}
similar_chains = {}
directory = options.input.strip().split("/")[-2]
prefix = "resume/"+directory
if options.resume:
    try:
        chains = pickle.load(open(prefix + "_chains.p", "rb"))
        pairs = pickle.load(open(prefix + "_pairs.p", "rb"))
        similar_chains = pickle.load(open(prefix + "_similar_chains.p", "rb"))
        structures = pickle.load(open(prefix + "_structures.p", "rb"))
        # stoichiometry =pickle.load(open(options.input + "stoichiometry.p", "rb"))

        print("The following structures have been recovered:")
        print("Chains: \n" + str(chains))
        print("Paired chains: \n"+str(pairs))
        print("Similar Chains:\n"+str(similar_chains))
        print("Sto: " + str(stoichiometry))
        print("structures: \n" + str(structures))
    except FileNotFoundError:
        pass

else:
    # STEP2: Get possible structures for Macrocomplex construction and skip others
    chain_index = 1
    for file in input_files:
        paired_chains = []
        structure = utils.get_structure(file, remove_het = True)
        for chain in utils.get_chains_from_structure(structure, remove_het = True):
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
    chains_b = open(prefix + "_chains.p", "wb")
    pairs_b = open(prefix + "_pairs.p", "wb")
    similar_chains_b = open(prefix + "_similar_chains.p", "wb")
    structures_b = open(prefix + "_structures.p", "wb")
    # stioichiometry_b = open(options.input + "stoichiometry.p", "wb")
    pickle.dump(chains, chains_b)
    pickle.dump(pairs, pairs_b)
    pickle.dump(similar_chains, similar_chains_b)
    pickle.dump(structures, structures_b)
        # pickle.dump(stoichiometry, stioichiometry_b)

# STEP3: Check for stoichiometry requirements
# if not utils.stoichiometry_is_possible(stoichiometry, chains, similar_chains):
#    raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)


complexes_found = []
temp_complexes = []
if options.verbose:
    sys.stderr.write("\n# Beginning to construct the complex\n\n")
# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex_real,
                      similar_chains, stoichiometry,
                      structures, used_pairs_real, clashing_real,old_complex_real ):
    # bruteforce ending!
    current_complex = copy.deepcopy(current_complex_real)
    used_pairs = copy.deepcopy(used_pairs_real)
    clashing = copy.deepcopy(clashing_real)
    old_complex = copy.deepcopy(old_complex_real)

    print("HOLA")
    # for the first round its going to be a random pair of chains.
    if current_complex is None:
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        used_pairs.append(random_choice_id)
        if options.verbose:
            sys.stderr.write("Selecting starting pair: "+str(random_choice_id)
                             + "\n")
        construct_complex(random_choice,
                          similar_chains, stoichiometry,
                          structures, used_pairs,clashing, old_complex )
        return
    else:
        for chain_in_current_complex in utils.get_chain_ids_from_structure(current_complex):
            a = utils.get_chains_from_structure(current_complex)
            b = utils.get_chains_from_structure(old_complex)
            new_chain = set(a).difference(set(b))
            print("new chain "+ str(new_chain))
            ps = utils.get_possible_structures(chain_in_current_complex,
                                               similar_chains, structures, used_pairs, clashing, new_chain)
            print("Possible structures : "+ str(ps))
            if len(ps) == 0 and options.stoichiometry is None:
                if options.verbose:
                    sys.stderr.write("Complex built!! :)\n")
                    sys.stderr.write("If you are not satisfied with the number"+
                                     "chains used you can use the stoichiometry"+
                                     "parameter.\n")
                    sys.stderr.write("Do not forget to use --resume to avoid "+
                                     "unnecessary computation.\n")
                complexes_found.append(current_complex)
                return
            elif len(ps) == 0:
                if stoichiometry is None:
                    if len(utils.get_chain_ids_from_structure(current_complex)) == len(utils.get_chain_ids_from_structure(old_complex)):
                        complexes_found.append(current_complex)
                        utils.print_chain_in_structure(current_complex)
                        utils.print_chain_in_structure(old_complex)
                        print("Already found!")
                        return
                else:
                    old_complex = current_complex
                    clashing = []
                    print("calling myself again")
                    construct_complex(current_complex,
                                      similar_chains, stoichiometry,
                                      structures, [], clashing, old_complex)
                    return
                    # ps = utils.get_possible_structures(chain_in_current_complex,
                    #                                    similar_chains, structures, used_pairs, clashing)

            for similar_chain_id in ps:
                if stoichiometry is not None and utils.complex_fits_stoich(current_complex,stoichiometry):
                    "2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                    break
                    return

                structure_id = ps[similar_chain_id]
                print('chain in current complex '+ str(chain_in_current_complex))
                print("similar_chain_id "+ str(similar_chain_id))
                utils.print_chain_in_structure(current_complex)
                print(used_pairs)
                structure_to_superimpose = structures[structure_id]
                other = [tuple_id for tuple_id in structure_id if tuple_id != similar_chain_id][0]
                if options.verbose:
                    sys.stderr.write("\nTrying to add: "+str(other)+"\n")
                    sys.stderr.write("Superimposing structure: "+str(structure_id)+"\n")
                    print("-------------------")
                    utils.print_chain_in_structure(structure_to_superimpose)
                matrix = utils.superimpose_chains_test(utils.get_chain(current_complex,chain_in_current_complex)
                                                      ,utils.get_chain_permissive(structure_to_superimpose, similar_chain_id))
                print("other "+ str(other))
                utils.print_chain_in_structure(structure_to_superimpose)
                print(used_pairs)
                matrix.apply(utils.get_chain(structure_to_superimpose,other))
                chain_to_add = utils.get_chain_permissive(structure_to_superimpose,other)
                print(chain_to_add)
                if not utils.are_clashing(current_complex,chain_to_add):
                    if options.verbose:
                        sys.stderr.write("Not clashing -- adding chain! \n")
                    utils.add_chain(current_complex, chain_to_add)
                    #current_complex[0].add(utils.get_chain(structure_to_superimpose,other))
                    used_pairs.append(structure_id)
                else:
                    clashing.append(structure_id)
                ## Check if finished!
                if current_complex is not None:
                    if stoichiometry is not None and utils.complex_fits_stoich(current_complex,stoichiometry):
                        if options.verbose:
                            sys.stderr.write("Complex built!! :) \n")
                        for chain in utils.get_chains_from_structure(current_complex):
                            print(chain)
                        complexes_found.append(current_complex)
                        print("jkdfsjklfjhljhfdkslsdhjklsdfhjkl")
                        break
                        return
                    else:
                        print("getting here!")
                        construct_complex(current_complex, chains,
                                          similar_chains, stoichiometry, pairs,
                                          structures, used_pairs)
                        return

                        #elif --> add repeated!
                        if len(complexes_found) == 1:
                            print("skldjhgsalfjhgskdjhfgkajsdhfgkasdhjfgsdafkjhsad!!!!!222!")
                            return

                        construct_complex(current_complex,
                                          similar_chains, stoichiometry,
                                          structures, used_pairs, clashing, old_complex)
                        if len(complexes_found) == 1:
                            print("skldjhgsalfjhgskdjhfgkajsdhfgkasdhjfgsdafkjhsad!!!!!!")
                            return


test_complex = construct_complex(None, similar_chains,
                                 stoichiometry, structures, [], [], [])
# Step5: Filter the good ones

# Step6 : write output file
index_file = 0
for complex in complexes_found:
    index_file += 1
    utils.write_structure_into_file(complex, "results/output"+str(index_file)+".cif", "mmcif")
