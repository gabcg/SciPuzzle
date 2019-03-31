#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re


def read_args():
    """Parse the commandline arguments and returns the namespace."""
    parser = argparse.ArgumentParser(description="Recreates a macrocomplex \
                                                  given different pdb files \
                                                  containing interacting \
                                                  protein pairs. \
                                                  \
                                                  If you need a graphic \
                                                  interface, use the -gui \
                                                  option by itself.")
    parser.add_argument('-i', '--input', dest="input",
                        action="store",
                        default=None,
                        help="Input directory containing PDB formatted files")
    parser.add_argument('-o', '--output', dest="output",
                        action="store",
                        default='model',
                        help="Ouput name")
    parser.add_argument('-s', '--stoichiometry', dest="stoichiometry",
                        action="store",
                        default=None,
                        help="Stoichiometry")
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action="store_true",
                        default=False,
                        help="Verbose mode")
    parser.add_argument('-r', '--resume', dest="resume",
                        action="store_true",
                        default=False,
                        help="Resume the program after a crash or when using \
                        a different stoichiometry")
    parser.add_argument('-c', '--chimera', dest="open_chimera",
                        action="store_true",
                        default=False,
                        help="Open models in Chimera when execution finishes")
    parser.add_argument('-gui', '--graphic_interface', dest="gui",
                        action="store_true",
                        default=False,
                        help="Graphic user interface mode")

    options = parser.parse_args()
    return options


def get_input_files(input):
    """
    Retrieves the input files from the input object and returns a list with the
    input files.

    The input object can either be a directory, in which case all the files
    terminating in ".pdb" are collected and returned in a list object, or it
    can be a list of files, in which case the input itself will be returned.
    """
    input_files = []
    if os.path.isdir(input):
        for filename in os.listdir(input):
            if filename.endswith(".pdb") and "chain" not in filename:
                input_files.append(os.path.join(input, filename))
    else:
        input_files = eval(input)
    return input_files


def parse_stoichiometry(stoichiometry):
    """
    Parses the String containing the Stioichiometry instructions and returns
    a dictionary.

    e.g. A3B12C1 --> {'A':3, 'B':12, 'C':1}
    """
    if stoichiometry is None:
        return None
    else:
        sto_dict = dict(re.findall("(\w)(\d*)", stoichiometry))
        for key in sto_dict:
            sto_dict[key] = int(sto_dict[key])
        return sto_dict
