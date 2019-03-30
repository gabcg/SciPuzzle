#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from gooey import Gooey, GooeyParser

desc = """
Recreates a macrocomplex given different PDB
files containing interacting protein pairs."""


@Gooey(program_name='SciPuzzle')
def gui():
    """Parse the commandline arguments and returns the namespace."""
    parser = GooeyParser(description=desc)
    parser.add_argument('-i', '--input', dest="input",
                        action="store",
                        default=None,
                        help="Input directory containing PDB formatted files",
                        required=True,
                        widget="DirChooser",
                        metavar="Input (-i)")
    parser.add_argument('-o', '--output', dest="output",
                        action="store",
                        default="",
                        help="Output directory",
                        widget="DirChooser",
                        metavar="Output (-o)")
    # If we want the output to be a dir, then the model files need a prefix
    # parser.add_argument('-p', '--prefix', dest="prefix",
    #                     action="store",
    #                     default="model",
    #                     help="Prefix for your file models",
    #                     metavar="Prefix (-p)")
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
    parser.add_argument('-r', '--resume', dest="resume",
                        action="store_true",
                        default=False,
                        help="Resume the program after a crash or when using \
                        a different stoichiometry",
                        metavar="Resume (-r)")
    parser.add_argument('-c', '--chimera', dest="open_chimera",
                        action="store_true",
                        default=False,
                        help="Open models in Chimera when execution finishes",
                        metavar="Open in Chimera (-c)")
    options = parser.parse_args()
    return options
