from sys import argv
import statistics

from Part3a import read_seqs, compute_distances, one_per_line
from multiple_sequence_alignment import MultipleAlignment


def align_by_nearest():
    """Alternative multiple alignment technique (not using AgglomerativeClustering)"""
    # Find the two sequences in distance_matrix with the smallest mutual distance
    # Ignore numbers along the diagonal of the matrix because they will always be zero
    minimum_dist = distance_matrix[0][1]
    minimum_dist_index = (0, 1)
    diagonal_index = 0
    for row in distance_matrix:
        for item in row:
            if row.index(item) == diagonal_index:
                continue
            if item < minimum_dist:
                minimum_dist = item
                minimum_dist_index = (distance_matrix.index(row), row.index(item))
        diagonal_index += 1
    
    # Having found sequences with minimum mutual distance, align those sequences first
    profile = ma.computeMultipleAlignment((sequences[minimum_dist_index[0]],
        sequences[minimum_dist_index[1]]))
    
    aligned_indices = [i for i in minimum_dist_index]
    # Repeat until all the sequences have been aligned
    while len(aligned_indices) < len(sequences):
        # Find the average distance of each sequence from the sequences already aligned
        # Extract from distance_matrix only the rows which correspond to already-aligned sequences
        aligned_rows = list(filter(lambda el: distance_matrix.index(el) in aligned_indices,
            distance_matrix))
        # Take the mean column value across the rows to produce a single "average distance" virtual row
        average_distances = [statistics.mean([row[col] for row in aligned_rows])
                for col in range(len(distance_matrix))]
        # Choose the sequence from this row which has the lowest distance but isn't already
        # in the alignment
        minimum_new_distance = min(filter(lambda el: average_distances.index(el) not in aligned_indices,
            average_distances))
        # We are only interested in which sequence this distance corresponds to
        minimum_new_distance_index = average_distances.index(minimum_new_distance)
        # Add the new nearest sequence to the alignment
        profile = ma.computeMultipleAlignment((profile, sequences[minimum_new_distance_index]))
        # Update aligned_indices with the index of the newly-aligned sequence
        aligned_indices.append(minimum_new_distance_index)

    return profile

if __name__ == "__main__":
    if len(argv) < 2:
        print("Usage: Part3b.py <filename>")
        exit()
    sequences: list = read_seqs(argv[1])
    distance_matrix = compute_distances(sequences)
    ma = MultipleAlignment()
    profile = align_by_nearest()
    score = ma.scoreMultipleAlignment(profile)
    print("Profile: {}\nScore: {}".format(one_per_line(profile), score))
