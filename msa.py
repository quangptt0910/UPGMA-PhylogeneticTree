
from itertools import combinations
from multiprocessing import Pool

import numpy as np
from typing import List, Tuple, Dict


from nw2 import scoring_path, traceback_alignment

def load_sequence_manual():
    """
    Load a sequence from manual input.

    Returns:
        tuple: (sequence, name)
    """
    name = input("Enter sequence name: ")
    seq = input(f"Enter sequence for {name}: ").strip().upper()
    return seq, name


def load_sequences_from_fasta(file_path):
    """
    Load multiple sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file

    Returns:
        tuple: (sequences, seq_names)
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        sequences = []
        seq_names = []
        current_seq = ""
        current_name = ""

        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Save the previous sequence if there was one
                if current_seq:
                    sequences.append(current_seq.upper())
                    seq_names.append(current_name)

                # Start a new sequence
                current_name = line[1:]  # Remove the '>' character
                current_seq = ""
            elif line:  # Only add non-empty lines
                current_seq += line

        # Add the last sequence
        if current_seq:
            sequences.append(current_seq.upper())
            seq_names.append(current_name)

        return sequences, seq_names
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return [], []
    except Exception as e:
        print(f"Error loading FASTA file: {str(e)}")
        return [], []

def worker(params):
    """
    Worker function for parallel computation of pairwise alignments

    Args:
        params: Tuple of ((i, seq_i), (j, seq_j), (match, mismatch, gap))

    Returns:
        Tuple: ((i, j), (aligned_i, aligned_j, score))
    """
    (i, seq_i), (j, seq_j), scoring_params = params
    match, mismatch, gap = scoring_params

    score_table, direction = scoring_path(seq_i, seq_j, match, mismatch, gap)
    aligned_i, aligned_j, score = traceback_alignment(seq_i, seq_j, direction, match, mismatch, gap)

    return ((i, j), (aligned_i, aligned_j, score))

def compute_all_pairwise_alignments(sequences, match=1, mismatch=-1, gap=-2):
    """
    Compute all pairwise alignments between sequences and their scores
    Args:
        sequences (list): List of sequences to align
        match (int): Score for matching
        mismatch (int): Score for mismatching
        gap (int): Score for gap penalty

    Returns:
        tuple: (score_matrix, alignment_pairs)
            - score_matrix: Matrix of alignment scores between sequences
            - alignment_pairs: Dictionary of pairwise alignments
    """
    n = len(sequences)
    score_matrix = np.zeros((n, n))
    alignment_pairs = {}

    # Create tasks for parallelization
    tasks = []
    for i, j in combinations(range(n), 2):
        tasks.append(((i, sequences[i]), (j, sequences[j]), (match, mismatch, gap)))

    # Parallel process
    with Pool() as pool:
        results = pool.map(worker, tasks)

    # Result process
    for (i, j), (aligned_i, aligned_j, score) in results:
        alignment_pairs[(i, j)] = (aligned_i, aligned_j)
        score_matrix[i, j] = score
        score_matrix[j, i] = score

    return score_matrix, alignment_pairs


def find_center_sequence(score_matrix):
    """
    Find the center sequence (most similar to all others)

    Args:
        score_matrix (np.array): Matrix of alignment scores
    Returns:
        int: Index of center sequence
    """
    row_sums = np.sum(score_matrix, axis=1)
    center_idx = np.argmax(row_sums)
    return int(center_idx)


def align_similar(s1, s2):
    """
    Find positions where gaps need to be inserted to make two sequences identical

    Args:
        s1 (str): First sequence
        s2 (str): Second sequence

    Returns:
        tuple: (changes_for_s1, changes_for_s2) - indices where gaps need insertion
    """
    change1, change2 = [], []
    i = 0

    # Copy sequences to avoid modifying originals
    s1_copy, s2_copy = s1, s2

    while s1_copy != s2_copy:
        if i >= len(s1_copy):
            # Need to extend s1 with remainder of s2
            s1_copy += s2_copy[i:]
            change1.extend(range(i, i + len(s2_copy[i:])))
            break

        if i >= len(s2_copy):
            # Need to extend s2 with remainder of s1
            s2_copy += s1_copy[i:]
            change2.extend(range(i, i + len(s1_copy[i:])))
            break

        if s1_copy[i] != s2_copy[i]:
            if s1_copy[i] == '-':
                # Insert gap in s2
                s2_copy = s2_copy[:i] + '-' + s2_copy[i:]
                change2.append(i)
            else:
                # Insert gap in s1
                s1_copy = s1_copy[:i] + '-' + s1_copy[i:]
                change1.append(i)

        i += 1

    return sorted(change1), sorted(change2)


def adjust(string_list, indices):
    """
    Insert gaps at specified positions in all strings in a list

    Args:
        string_list (list): List of strings to modify
        indices (list): Positions where gaps should be inserted
    """
    for idx in range(len(string_list)):
        s = string_list[idx]
        for pos in indices:
            if pos <= len(s):
                s = s[:pos] + '-' + s[pos:]
        string_list[idx] = s


def merge_alignments(center_idx, alignments, sequences, match=1, mismatch=-1, gap=-2):
    """
    Merge pairwise alignments into a multiple sequence alignment

    Uses center-star approach: progressively adds each sequence to the growing MSA by
    reconciling gaps between the current center profile and new pairwise alignment.

    Args:
        center_idx (int): Index of center sequence
        alignments (dict): Dictionary of pairwise alignments
        sequences (list): Original sequences
        match, mismatch, gap: Scoring parameters (unused here)

    Returns:
        list of str: Multiple sequence alignment ordered by original sequence indices
    """
    # Gather pairwise alignments with center
    pair_aligns = []  # tuples of (center_aln, seq_aln, idx)
    for idx in range(len(sequences)):
        if idx == center_idx:
            continue
        key = (center_idx, idx) if (center_idx, idx) in alignments else (idx, center_idx)
        a, b = alignments[key]
        # Ensure a is center alignment
        if key[0] == center_idx:
            center_aln, seq_aln = a, b
        else:
            center_aln, seq_aln = b, a
        pair_aligns.append((center_aln, seq_aln, idx))

    # If no other sequences, return solo center
    if not pair_aligns:
        return [sequences[center_idx]]

    # Initialize MSA with first pair
    first_center, first_seq, first_idx = pair_aligns[0]
    msa_rows = [first_center, first_seq]
    seq_order = [first_idx]

    # Add each remaining pair
    for center_aln, seq_aln, idx in pair_aligns[1:]:
        # Determine gap insertion positions to align current profile center and new center alignment
        ch1, ch2 = align_similar(msa_rows[0], center_aln)
        # Apply gaps to existing MSA rows
        adjust(msa_rows, ch1)
        # Apply gaps to the new sequence alignment
        new_seq_list = [seq_aln]
        adjust(new_seq_list, ch2)
        # Append adjusted new sequence to MSA
        msa_rows.append(new_seq_list[0])
        seq_order.append(idx)

    # Build final aligned list in original order
    aligned = [''] * len(sequences)
    aligned[center_idx] = msa_rows[0]
    for i, idx in enumerate(seq_order, start=1):
        aligned[idx] = msa_rows[i]
    return aligned



def center_star_alignment(sequences, match=1, mismatch=-1, gap=-2):
    """
    Perform multiple sequence alignment using the center star method.

    Args:
        sequences (list): List of sequences to align
        match (int): Score for matching characters
        mismatch (int): Score for mismatching characters
        gap (int): Score for gap penalty

    Returns:
        tuple: (msa, score_matrix, center_idx)
            - msa: List of aligned sequences in the MSA
            - score_matrix: Matrix of alignment scores
            - center_idx: Index of the center sequence
    """
    print("Computing pairwise alignments...")
    score_matrix, alignments = compute_all_pairwise_alignments(sequences, match, mismatch, gap)
    print(f"Score matrix:\n{score_matrix}")

    print("Finding center sequence...")
    center_idx = find_center_sequence(score_matrix)
    print(f"Center sequence is sequence #{center_idx + 1} ({sequences[center_idx]})")

    print("Merging alignments...")
    msa = merge_alignments(center_idx, alignments, sequences, match, mismatch, gap)

    return msa, score_matrix, center_idx

def compute_alignment_statistics(msa):
    """
    Compute statistics for the multiple sequence alignment (MSA).

    Args:
        msa (list of str): List of aligned sequences. Each sequence must have the same length.

    Returns:
        dict: A dictionary containing the statistics of the alignment.
    """
    if not msa or any(len(seq) != len(msa[0]) for seq in msa):
        raise ValueError("All sequences in MSA must have the same length and cannot be empty.")

    seq_length = len(msa[0])
    n = len(msa)

    # Initialize statistics
    stats = {
        "total_columns": seq_length,
        "identity_percentage": 0,
        "match_count": 0,
        "mismatch_count": 0,
        "gap_count": 0,
        "conserved_columns": 0,
        "column_stats": []
    }

    # Analyze each column in the alignment
    for col_idx in range(seq_length):
        column = [seq[col_idx] for seq in msa]
        gaps = column.count('-')
        residues = [res for res in column if res != '-']

        if residues:  # Ignore all-gap columns
            # Calculate matches and mismatches
            most_common_residue = max(set(residues), key=residues.count)
            matches = residues.count(most_common_residue)
            mismatches = len(residues) - matches

            # Update overall stats
            stats["gap_count"] += gaps
            stats["match_count"] += matches
            stats["mismatch_count"] += mismatches

            # Check for conserved column
            if matches == len(residues) and gaps == 0:
                stats["conserved_columns"] += 1

            # Record per-column statistics
            stats["column_stats"].append({
                "position": col_idx + 1,
                "gaps": gaps,
                "matches": matches,
                "mismatches": mismatches
            })

    # Calculate identity percentage
    total_pairs = stats["match_count"] + stats["mismatch_count"]
    stats["identity_percentage"] = (stats["match_count"] / total_pairs * 100) if total_pairs > 0 else 0

    return stats


def print_msa(msa, seq_names=None):
    """
    Print the multiple sequence alignment in a readable format

    Args:
        msa (list): List of aligned sequences
        seq_names (list): Optional list of sequence names
    """
    if not seq_names:
        seq_names = [f"Seq {i + 1}" for i in range(len(msa))]

    # Find max length of sequence names for proper formatting
    max_name_len = max(len(name) for name in seq_names)

    print()

    # Print the MSA
    for i, seq in enumerate(msa):
        print(f"{seq_names[i]:{max_name_len}} ", end="")
        for j, char in enumerate(seq):
            print(char, end="")
        print()


def compute_identity_matrix(msa, method='all_positions'):
    """
    Compute the percentage identity matrix from an MSA

    Args:
        msa (list): List of aligned sequences
        method (str): Method for computing identity
            - 'all_positions': Include all positions (gaps count as mismatches)
            - 'no_gaps': Only count positions where neither sequence has a gap
            - 'one_gap': Count positions where at most one sequence has a gap
            - 'normalized': Normalize by shorter sequence length

    Returns:
        np.array: Matrix of percent identities between sequences
    """
    n = len(msa)
    identity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            if i == j:
                identity_matrix[i, j] = 1.0  # 100% identity with itself
                continue

            if method == 'all_positions':
                # Count all positions, treat gaps as mismatches
                matches = sum(1 for k in range(len(msa[0])) if msa[i][k] == msa[j][k])
                total = len(msa[0])

            elif method == 'no_gaps':
                # Only count positions where neither sequence has a gap
                matches = 0
                total = 0
                for k in range(len(msa[0])):
                    if msa[i][k] != '-' and msa[j][k] != '-':
                        total += 1
                        if msa[i][k] == msa[j][k]:
                            matches += 1

            elif method == 'one_gap':
                # Count positions where at most one sequence has a gap
                matches = 0
                total = 0
                for k in range(len(msa[0])):
                    if not (msa[i][k] == '-' and msa[j][k] == '-'):  # Not both gaps
                        total += 1
                        if msa[i][k] == msa[j][k]:
                            matches += 1

            elif method == 'normalized':
                # Normalize by the length of the shorter original sequence
                seq_i_len = sum(1 for c in msa[i] if c != '-')
                seq_j_len = sum(1 for c in msa[j] if c != '-')
                min_len = min(seq_i_len, seq_j_len)

                matches = sum(1 for k in range(len(msa[0]))
                              if msa[i][k] == msa[j][k] and msa[i][k] != '-')
                total = min_len

            identity = matches / total if total > 0 else 0
            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity

    print(f"Identity matrix (method: {method}):\n{identity_matrix}")
    return identity_matrix


def save_alignment_to_file(msa: List[str], params: Dict, stats: Dict, filename: str, seq_names: List[str] = None):
    """
    Save alignment with parameters and statistics to file
    """
    with open(filename, 'w') as f:
        f.write("===== PROGRAM PARAMETERS =====\n\n")
        f.write(f"\tMatch score: {params['match']}\n")
        f.write(f"\tMismatch score: {params['mismatch']}\n")
        f.write(f"\tGap penalty: {params['gap']}\n\n")

        f.write("SEQUENCE NAMES: \n")

        if seq_names:
            for name in seq_names:
                f.write(f"\t{name}\n")
        else:
            for i in range(len(msa)):
                f.write(f"\t{msa[i]}\n")

        f.write("\nMULTIPLE SEQUENCE ALIGNMENT:\n")
        # Write each aligned sequence on a single line
        for seq in msa:
            f.write(f"{seq}\n")

        f.write("\nSTATISTICS:\n")
        f.write(f"\tIdentity percentage: {stats['identity_percentage']:.2f}%\n")
        f.write(f"\tTotal matches: {stats['match_count']}\n")
        f.write(f"\tTotal mismatches: {stats['mismatch_count']}\n")
        f.write(f"\tTotal gaps: {stats['gap_count']}\n")
        f.write(f"\tAlignment length: {stats['total_columns']}\n")