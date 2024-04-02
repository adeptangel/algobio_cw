from NW_Part1 import Alignment
from Part2 import pretty_print_matrix

"""
Kimura model distance between pairs of sequences.
"""


def kimura_dist(a, b):
    """Kimura distance between sequences a and b."""
    align = Alignment(a, b)
    align.compute_alignment()

    # Calculate the number of columns (with both sequences justified
    # against one another) with no gaps
    n_gaps = 0
    for i in range(len(align.xa)):
        if align.xa[i] == "-" or align.ya[i] == "-":
            n_gaps += 1
    


def compute_distances_from_file(filename):
    print("Computing output for file " + filename)
    with open(filename, "r") as file:
        contents = file.readlines()
    # Strip newlines here as well
    reported_len = contents[0].strip()
    sequences = [seq.strip() for seq in contents[1:]]
    print("Found {0} sequences. Reports {1} sequences.".format(len(sequences), reported_len))

    matrix = [[kimura_dist(i, j) for j in sequences] for i in sequences]
    pretty_print_matrix(matrix)


if __name__ == "__main__":
    compute_distances_from_file("sequences/multiple3.txt")
