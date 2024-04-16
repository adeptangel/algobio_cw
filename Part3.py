from math import log
from sys import argv
import statistics

from sklearn.cluster import AgglomerativeClustering

from NW_Part1 import Alignment
from Part2 import pretty_print_matrix
from multiple_sequence_alignment import MultipleAlignment

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


def read_seqs(filename):
    with open(filename, "r") as file:
        contents = file.readlines()
    # Strip newlines here as well
    reported_len = contents[0].strip()
    sequences = [seq.strip() for seq in contents[1:]]
    print("Found {0} sequences. Reports {1} sequences.".format(len(sequences), reported_len))
    return sequences
    

def compute_distances(sequences):
    print("Computing Kimura distances...")
    matrix = [[kimura_dist(i, j) for j in sequences] for i in sequences]
    print("Done.")
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
    print("Performing multiple alignment from guide tree...")
    align = MultipleAlignment()

    for merger in cluster.children_:
        # merger always has two elements
        profile = align.computeMultipleAlignment([sequences[merger[0]], sequences[merger[1]]])
        # Add the profile to sequences as though it is a new sequence
        # This enables sequence index access to the profiles
        sequences.append(profile)

    print("Done.")
    # Last elements of sequences is always the final alignment profile
    return sequences[-1]


def part3b_alternative():
    """Alternative multiple alignment technique (not using AgglomerativeClustering)"""
    if len(argv) < 2:
        print("Usage: Part3.py <filename>")
        return
    sequences: list = read_seqs(argv[1])
    align = MultipleAlignment()
    distance_matrix = compute_distances(sequences)
    
    # Find the two sequences in distance_matrix with the smallest mutual distance
    minimum_dist = distance_matrix[0][1]
    minimum_dist_index = (0, 1)
    for row in distance_matrix:
        diagonal_index = 0
        for item in row:
            if row.index(item) == diagonal_index:
                continue
            if item < minimum_dist:
                minimum_dist = item
                minimum_dist_index = (distance_matrix.index(row), row.index(item))
        diagonal_index += 1
    
    # Having found sequences with minimum mutual distance, align those sequences first
    profile = align.computeMultipleAlignment((sequences[minimum_dist_index[0]],
        sequences[minimum_dist_index[1]]))
    
    aligned_indices = [i for i in minimum_dist_index]
    # Repeat until all the sequences have been aligned
    while len(aligned_indices) < len(sequences):
        # Find the average distance of each sequence from the sequences already aligned
        # Make a virtual row based on the average of rows for sequences in the alignment
        average_distances = [statistics.mean([row[col] for row in distance_matrix])
                for col in range(len(distance_matrix))]
        # Choose the sequence from this row which has the lowest distance but isn't already in the alignment
        minimum_new_distance = min(filter(lambda el: average_distances.index(el) not in aligned_indices,
            average_distances))
        # We are only interested in which sequence this distance corresponds to
        minimum_new_distance_index = average_distances.index(minimum_new_distance)
        # Add the new nearest sequence to the alignment
        profile = align.computeMultipleAlignment((profile, sequences[minimum_new_distance_index]))
        # Update aligned_indices with the index of the newly-aligned sequence
        aligned_indices.append(minimum_new_distance_index)

    return profile

if __name__ == "__main__":
    # if len(argv) < 2:
    #     print("Usage: Part3.py <filename>")
    # else:
    #     sequences: list = read_seqs(argv[1])
    #     distance_matrix = compute_distances(sequences)
    #     # print("Python distance matrix: {}".format(pretty_print_matrix(distance_matrix)))
    #     cluster = make_cluster(distance_matrix)
    #     final_alignment: list = align_from_clustering(cluster, sequences)
    #     print("Final alignment: {}".format(final_alignment))
    #     # Score alignment
    #     m_align = MultipleAlignment()
    #     print("Scoring multiple alignment...")
    #     score: float = m_align.scoreMultipleAlignment(final_alignment)
    #     print("Done.")
    #     print("Score: {}".format(score))
    print("Profile: {}".format(part3b_alternative()))
