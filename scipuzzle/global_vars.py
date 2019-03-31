import sys
id_list = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
           'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
           'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
           '!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-',
           '.', '/', ':', ';', '<', '=', '>', '?', '@', '[', ']', '^', '_',
           '`', '{', '|', '}', '~']


def write_welcoming(input_files):
    sys.stderr.write("\n\n-------------------------------\n")
    sys.stderr.write("--  Welcome to SciPuzzle! ;) --\n")
    sys.stderr.write("-------------------------------\n")
    sys.stderr.write("\n\nInput correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write(file+"\n")
