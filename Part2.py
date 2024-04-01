import itertools

from NW_Part1 import Alignment

"""
Running this script will compute a pairwise alignment between two given 
sequences X and Y. 

As you hopefully see when running this script, the score of an alignment is 
returned for each pair of inputs using the method score_alignment

After confirming that running this code indeed provides a visualisation of 
a sequence alignment along with a score, begin editing this script to return an 
array of score relating to each sequence found within one of the files in ./sequences

"""

# # Here is a sample output of working with the NW algorithm from Part1
# x = "GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL"
# y = "NNPELQAHAGKVFKLVYEAAIQLQVTGVVVTDATLKNLGSVHVSKG"
# 
# # Here we are calling the alignment method
# a = Alignment(x, y) # instantiates the alignment object, providing sequences X and Y
# a.compute_alignment() # run to do the actual computation
# print("Computed alignment:")
# a.display_alignment() # this can be displayed to help visualise what is going on
# 
# # a.score_alignment provides a float score for X and Y alignment
# print("Score of alignment: " + str(a.score_alignment()))

def score_alignment(seq_a, seq_b):
    agnmt = Alignment(seq_a, seq_b)
    agnmt.compute_alignment()
    return agnmt.score_alignment()


def pretty_print_matrix(mat):
    print("[")
    for row in matrix:
        print(row, end=",\n")
    print("]", end="")


if __name__ == "__main__":
    filename = "sequences/multiple3.txt"
    print("Computing output for file " + filename)
    with open(filename, "r") as file:
        contents = file.readlines()
    # Strip newlines here as well
    reported_len = contents[0].strip()
    sequences = [seq.strip() for seq in contents[1:]]
    print("Found {0} sequences. Reports {1} sequences.".format(len(sequences), reported_len))
    print("Sequences: {}".format(sequences))

    pairs = itertools.combinations(sequences, 2)
    copies = [(a, a) for a in sequences]
    print("Pairs:")
    print({pair: score_alignment(*pair) for pair in pairs})
    print("Copies:")
    print({copy: score_alignment(*copy) for copy in copies})

    matrix = [[score_alignment(i, j) for j in sequences] for i in sequences]
    pretty_print_matrix(matrix)
