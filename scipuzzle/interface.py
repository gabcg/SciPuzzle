#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from gooey import Gooey, GooeyParser

desc = """
Recreates a macrocomplex given different pdb
files containing interacting protein pairs."""


@Gooey(program_name='SciPuzzle')
def gui():
    """Parse the commandline arguments and returns the namespace."""
    parser = GooeyParser(description=desc)
    parser.add_argument('-i', '--input', dest="input",
                        action="store",
                        default=None,
                        help="Input Folder containing PDB formatted files",
                        required=True,
                        widget="DirChooser",
                        metavar="Input (-i)")
    parser.add_argument('-o', '--output', dest="output",
                        action="store",
                        default="reconstructed_macrocomplex.pdb",
                        help="PDB formatted outputfile",
                        widget="FileSaver",
                        metavar="Output (-o)")
    parser.add_argument('-s', '--stoichiometry', dest="stoichiometry",
                        action="store",
                        default=None,
                        help="Stoichiometry",
                        metavar="Stoichiometry (-s)")
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action="store_true",
                        default=False,
                        help="Enable verbose mode",
                        metavar="Verbose (-v)")
    parser.add_argument('-c', '--chimera', dest="open_chimera",
                        action="store_true",
                        default=False,
                        help="Open models in Chimera when execution finishes",
                        metavar="Open in Chimera (-c)")
    options = parser.parse_args()
    return options
