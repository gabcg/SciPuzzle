import sys
def write_welcoming(input_files):
    sys.stderr.write("\n\n-------------------------------\n")
    sys.stderr.write("--  Welcome to SciPuzzle! ;) --\n")
    sys.stderr.write("-------------------------------\n")
    sys.stderr.write("\n\nInput correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write(file+"\n")


def complex_built_no_stoich():
    sys.stderr.write("Complex built!! :)\n")
    sys.stderr.write("If you are not satisfied with the number"+
                     "chains used you can use the stoichiometry"+
                     "parameter.\n")
    sys.stderr.write("Do not forget to use --resume to avoid "+
                     "unnecessary computation.\n")


def beginning(random_choice_id):
    sys.stderr.write("### Beginning")
    sys.stderr.write("Selecting starting pair: "+str(random_choice_id)
                     + "\n")

def trying_superimpose(other, structure_id):
    sys.stderr.write("\nTrying to add: "+str(other)+"\n")
    sys.stderr.write("Superimposing structure: "+str(structure_id)+"\n")
    print("-------------------")
