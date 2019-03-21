#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from gooey import Gooey
from gooey import GooeyParser


@Gooey
def ui():
    """Parse the commandline arguments and returns the namespace."""
    parser = GooeyParser(description="Recreates a macrocomplex \
                                                  given different pdb files \
                                                  containing interacting \
                                                  protein pairs.")
    parser.add_argument('-i', '--input', dest="input",
                        action="store",
                        default=None,
                        help="Input Folder containing PDB formatted files",
                        required=True,
                        widget='DirChooser')
    parser.add_argument('-o', '--output', dest="output",
                        action="store",
                        default="reconstructed_macrocomplex.pdb",
                        help="PDB formatted outputfile",
                        widget='FileSaver')
    parser.add_argument('-s', '--stoichiometry', dest="stoichiometry",
                        action="store",
                        default=None,
                        help="Stoichiometry")
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action="store_true",
                        default=False,
                        help="Verbose Mode")
    options = parser.parse_args()
    return options
