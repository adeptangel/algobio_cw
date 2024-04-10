from math import log

from sklearn.cluster import AgglomerativeClustering

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
    n_nogap = 0
    for i in range(len(align.xa)): # len(align.xa) == len(align.ya)
        if align.xa[i] != "-" and align.ya[i] != "-":
            n_nogap += 1
    
    # Calculate the number of columns where both residues are identical
    n_identical = 0
    for i in range(len(align.xa)):
        if align.xa[i] == align.ya[i]:
            n_identical += 1
    
    s = n_identical / n_nogap
    d = 1 - s
    return -log(1 - d - 0.2 * (d ** 2)) # with one argument, math.log is the natural log
    

def compute_distances_from_file(filename):
    print("Computing output for file " + filename)
    with open(filename, "r") as file:
        contents = file.readlines()
    # Strip newlines here as well
    reported_len = contents[0].strip()
    sequences = [seq.strip() for seq in contents[1:]]
    print("Found {0} sequences. Reports {1} sequences.".format(len(sequences), reported_len))
    
    matrix = [[kimura_dist(i, j) for j in sequences] for i in sequences]
    return matrix


def make_cluster(distances):
    """
    Agglomerative clustering function. Expects a distance matrix as input.
    Returns the AgglomerativeClustering object.
    """
    print("Clustering...")

    cluster = AgglomerativeClustering(linkage="average", metric="precomputed", distance_threshold=None)
    cluster.fit(distances)
    # model is a reference to the same object as cluster
    # print("enumerate:", list(enumerate(cluster.children_)))
    # print("cluster.labels_:", cluster.labels_)
    # next_node = len(cluster.labels_)
    # for i, merge in enumerate(cluster.children_):
    #     print("Align {} with {} to give {}".format(merge[0], merge[1], next_node))
    #     next_node += 1
    print("Done.")
    return cluster


def align_from_clustering(cluster: AgglomerativeClustering, sequences: list):
    """Follow the guide tree produced by a clustering to align several sequences."""
    # Clustering isn't associated with list elements so we need to find the correct
    # mapping from list elements to their cluster nodes
    # named_sequences: dict = enumerate(sequences)
    
    for merger in cluster:
        # merger always has two elements
        


if __name__ == "__main__":
    distance_matrix = compute_distances_from_file("sequences/multiple3.txt")
    print("Python distance matrix: {}".format(pretty_print_matrix(distance_matrix)))
    cluster = make_cluster(distance_matrix)
    alignment_guide(cluster)
