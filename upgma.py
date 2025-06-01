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

            # Create new cluster
            new_size = sizes[min_i] + sizes[min_j]
            new_clade = BaseTree.Clade()
            new_clade.name = f"N{next_id}"  # Name internal node
            next_id += 1

            # Get the clusters to merge and set their branch lengths
            left_child = clusters[min_i]
            right_child = clusters[min_j]

            # Set branch lengths - this is the key fix
            left_child.branch_length = branch_i
            right_child.branch_length = branch_j

            # Set children for new clade
            new_clade.clades = [left_child, right_child]

            # Update distance matrix
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


if __name__ == "__main__":
    print("=== Simple Test ===")
    test_upgma()
    print("\n=== Detailed Test (Similar to your sequences) ===")
    test_upgma_detailed()