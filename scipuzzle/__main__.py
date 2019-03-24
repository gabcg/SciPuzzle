#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments
import interface
import utils
import exceptions

# STEP 1: parse the arguments
if sys.argv[1] == '-gui':
    options = interface.gui()
else:
    options = arguments.read_args()

# Alternative using argparse
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
        sys.stderr.write(file+"\n")

# STEP2: Get possible structures for Macrocomplex construction and skip others
chains = {}
pairs = []
chain_index = 1
for file in input_files:
    paired_chains = []
    for chain in utils.get_chains(file):
        chain_id = str(chain_index)+"_"+str(chain.id)
        chains[chain_id] = chain
        paired_chains.append(chain_id)
    pairs.append(paired_chains)
    chain_index += 1
similar_chains = utils.get_similar_chains(chains)
(chains, similar_chains, pairs) = utils.remove_useless_chains(chains, similar_chains, pairs)

print("Chains: \n" + str(chains))
print("Paired chains: \n"+str(pairs))
print("Similar Chains:\n"+str(similar_chains))
print("Sto: " + str(stoichiometry))

if options.verbose:
    sys.stderr.write('Number of chains: %s' % (len(chains)))
    sys.stderr.write('Number of paired chains: %s' % (len(pairs)))
    sys.stderr.write('Number of similar chains: %s' % (len(similar_chains)))
    sys.stderr.write('Stoichiometry parsed as: %s' % (str(stoichiometry)))

# STEP3: Check for stoichiometry requirements
if not utils.stoichiometry_is_possible(stoichiometry, chains, similar_chains):
    raise exceptions.IncorrectStoichiometry(stoichiometry=stoichiometry)


# STEP4: Begin Macrocomplex reconstruction!
def construct_complex(current_complex, chains,
                      similar_chains, stoichiometry, pairs_left):
    test = utils.are_clashing(chains['1_C'], chains['1_D'])
    (test2, rmsd) = utils.superimpose(chains['1_C'], chains['1_D'])

    # current_complex is a list of chains
    # for the first round its going to be a random pair of chains.
    if len(current_complex) == 0:
        current_complex = chains[list(chains.keys())[0]]
        construct_complex(current_complex, chains,
                          similar_chains, stoichiometry, pairs_left)
    # remove pairs left!!
    print(current_complex)
    print("Wasdlakd")

    # Each time we add a new chain - update current_complex

    construct_complex(current_complex, chains, similar_chains, stoichiometry)

    # if current_complex fulfills stoichiometry - save in list and return
    if utils.complex_fits_stoich(current_complex, stoichiometry):
        # Output: List of lists of chain. Each list of chain corresponds to a
        # constructed complex
        return [current_complex]


test_complex = construct_complex([], chains, similar_chains, stoichiometry)

# Step5: Filter the good ones

# Step6 : write output file
utils.write_structure_into_pdb(test_complex[0], 'test2.pdb')

# Step 7 (optional): open models in Chimera
# Need to store the outputs as a list
models = ['test.pdb']
if options.open_chimera:
    utils.open_in_chimera(models)
