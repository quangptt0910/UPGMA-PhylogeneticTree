import numpy as np
from Bio.Phylo import BaseTree
from Bio import Phylo
from io import StringIO


class UPGMATreeBuilder:
    @staticmethod
    def upgma(dist_matrix, labels):
        n = len(labels)
        # Create initial clusters with names
        clusters = [BaseTree.Clade(name=label) for label in labels]
        heights = np.zeros(n)
        sizes = np.ones(n)

        D = dist_matrix.copy()
        np.fill_diagonal(D, 0)

        next_id = 0  # For internal node naming
        cluster_history = [] # track merge history

        while n > 1:
            # Find the closest clusters
            min_val = np.inf
            min_i, min_j = -1, -1
            for i in range(n):
                for j in range(i + 1, n):
                    if D[i, j] < min_val:
                        min_val = D[i, j]
                        min_i, min_j = i, j

            # Calculate branch lengths
            new_height = min_val / 2
            branch_i = new_height - heights[min_i]
            branch_j = new_height - heights[min_j]

            # Create new clade
            new_clade = BaseTree.Clade(name=f"N{next_id}")
            next_id += 1

            # Create new children with branch lengths
            left_child = BaseTree.Clade()
            left_child.branch_i = branch_i
            left_child.clades = [clusters[min_i]] if not clusters[min_i].clades else clusters[min_i].clades

            right_child = BaseTree.Clade()
            right_child.branch_i = branch_i
            right_child.clades = [clusters[min_i]] if not clusters[min_i].clades else clusters[min_i].clades

            new_clade.clades = [left_child, right_child]
            cluster_history.append(new_clade)

            # Update distance matrix
            new_size = sizes[min_i]
            new_dists = []
            for k in range(n):
                if k != min_i and k != min_j:
                    d_new = (D[min_i, k] * sizes[min_i] + D[min_j, k] * sizes[min_j]) / new_size
                    new_dists.append(d_new)

            # Update structures - preserve clusters not being merged
            remaining_indices = [i for i in range(n) if i not in (min_i, min_j)]

            # Create new distance matrix
            new_D = np.zeros((n - 1, n - 1))
            # Copy existing distances for remaining clusters
            for i, idx_i in enumerate(remaining_indices):
                for j, idx_j in enumerate(remaining_indices):
                    new_D[i, j] = D[idx_i, idx_j]

            # Add distances for the new merged cluster
            for i, idx in enumerate(remaining_indices):
                new_D[i, -1] = new_dists[i]
                new_D[-1, i] = new_dists[i]
            new_D[-1, -1] = 0  # Distance to itself

            D = new_D

            # Update other arrays
            clusters = [clusters[i] for i in remaining_indices] + [new_clade]
            sizes = np.concatenate((sizes[remaining_indices], [new_size]))
            heights = np.concatenate((heights[remaining_indices], [new_height]))
            n -= 1

        # Set root branch length to 0
        root = clusters[0]
        root.branch_length = 0
        return BaseTree.Tree(root=root, rooted=True)

    @staticmethod
    def tree_to_newick(tree):
        """Convert Bio.Phylo tree to Newick format string"""
        with StringIO() as handle:
            Phylo.write(tree, handle, "newick")
            newick_str = handle.getvalue().strip()
        return newick_str


# Test function with debugging
def test_upgma_detailed():
    # Test with your sequence data converted to distances
    # This simulates what your MSA alignment would produce
    dist_matrix = np.array([
        [0.0, 0.153846, 0.307692, 0.153846],  # Seq1 distances
        [0.153846, 0.0, 0.230769, 0.076923],  # Seq2 distances
        [0.307692, 0.230769, 0.0, 0.230769],  # Seq3 distances
        [0.153846, 0.076923, 0.230769, 0.0]  # Seq4 distances
    ])
    labels = ['Seq1', 'Seq2', 'Seq3', 'Seq4']

    print("Distance Matrix:")
    print(dist_matrix)
    print("\nRunning UPGMA...")

    tree = UPGMATreeBuilder.upgma(dist_matrix, labels)
    newick = UPGMATreeBuilder.tree_to_newick(tree)
    print("Newick:", newick)

    # Print branch lengths for debugging
    print("\nBranch Lengths:")
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            print(f"{clade.name or 'Root'}: {clade.branch_length:.6f}")

    return tree


# Simple test
def test_upgma():
    dist_matrix = np.array([
        [0, 2, 4],
        [2, 0, 4],
        [4, 4, 0]
    ])
    labels = ['A', 'B', 'C']
    tree = UPGMATreeBuilder.upgma(dist_matrix, labels)
    newick = UPGMATreeBuilder.tree_to_newick(tree)
    print("Simple test Newick:", newick)

    # Print branch lengths for debugging
    print("\nSimple test Branch Lengths:")
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            print(f"{clade.name or 'Root'}: {clade.branch_length:.6f}")

    return tree


def score_to_distance_matrix(score_matrix, method='max_score_normalize'):
    """
    Convert alignment score matrix to distance matrix for UPGMA

    Args:
        score_matrix: Matrix of pairwise alignment scores
        method: Method for conversion ('max_score_normalize', 'percent_identity', 'negative_log')

    Returns:
        Distance matrix suitable for UPGMA
    """
    n = score_matrix.shape[0]
    dist_matrix = np.zeros((n, n))

    if method == 'max_score_normalize':
        # Method 1: Normalize by maximum possible score
        # Distance = (max_score - actual_score) / max_score
        max_score = np.max(score_matrix)
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Convert score to distance: higher score = lower distance
                    dist_matrix[i][j] = (max_score - score_matrix[i][j]) / max_score
                    # Ensure non-negative distances
                    dist_matrix[i][j] = max(0, dist_matrix[i][j])

    elif method == 'percent_identity':
        # Method 2: Use percent identity approach
        # Assumes diagonal contains self-alignment scores
        diagonal_scores = np.diag(score_matrix)
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Calculate percent identity relative to self-scores
                    max_possible = (diagonal_scores[i] + diagonal_scores[j]) / 2
                    if max_possible > 0:
                        identity = score_matrix[i][j] / max_possible
                        dist_matrix[i][j] = max(0, 1 - identity)
                    else:
                        dist_matrix[i][j] = 1.0

    elif method == 'simple_inverse':
        # Method 3: Simple inverse transformation
        # Convert positive scores to small distances, negative scores to large distances
        max_score = np.max(score_matrix)
        min_score = np.min(score_matrix)
        score_range = max_score - min_score

        for i in range(n):
            for j in range(n):
                if i != j:
                    # Normalize score to 0-1, then invert
                    normalized_score = (score_matrix[i][j] - min_score) / score_range
                    dist_matrix[i][j] = 1 - normalized_score

    # Ensure diagonal is zero
    np.fill_diagonal(dist_matrix, 0)
    return dist_matrix


if __name__ == "__main__":
    print("=== Simple Test ===")
    test_upgma()
    print("\n=== Detailed Test (Similar to your sequences) ===")
    test_upgma_detailed()