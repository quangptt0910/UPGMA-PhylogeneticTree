import numpy as np

def scoring_path(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Calculate the scoring matrix for sequence alignment (original single path logic)
    """
    len1, len2 = len(seq1), len(seq2)
    path = np.zeros((len1 + 1, len2 + 1), dtype=int) #score matrix
    direction = np.zeros_like(path, dtype=int)  # 0 diagonal, 1 Up, 2 left

    for i in range(1, len1 + 1):
        path[i, 0] = i * gap
        direction[i, 0] = 1
    for j in range(1, len2 + 1):
        path[0, j] = j * gap
        direction[0, j] = 2

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            diag_score = path[i - 1, j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up_score = path[i - 1, j] + gap
            left_score = path[i, j - 1] + gap

            max_score = max(diag_score, up_score, left_score)
            path[i, j] = max_score

            if max_score == left_score:
                direction[i, j] = 2
            if max_score == up_score:
                direction[i, j] = 1
            if max_score == diag_score:
                direction[i, j] = 0

    return path, direction


def traceback_alignment(seq1, seq2, direction, match=1, mismatch=-1, gap=-2):
    """
    Trace back through the direction matrix to get the optimal global alignment
    and its total score.

    Parameters:
    - seq1, seq2: input sequences to align
    - direction: matrix of directions from the scoring path
    - match, mismatch, gap: scoring parameters

    Returns:
    - aligned1, aligned2: aligned sequences
    - score: final alignment score
    """
    i, j = len(seq1), len(seq2)
    aligned1, aligned2 = [], []
    score = 0

    # Precompute length checks to avoid repeated evaluations
    seq1_idx = i - 1
    seq2_idx = j - 1

    # Main trace
    while i > 0 and j > 0:
        dir_val = direction[i, j]

        if dir_val == 0:  # diagonal: match or mismatch
            a, b = seq1[seq1_idx], seq2[seq2_idx]
            aligned1.append(a)
            aligned2.append(b)
            # Use boolean comparison directly for score adjustment
            score += match if a == b else mismatch
            i -= 1
            j -= 1
            seq1_idx -= 1
            seq2_idx -= 1

        elif dir_val == 1:  # up: gap in seq2
            aligned1.append(seq1[seq1_idx])
            aligned2.append("-")
            score += gap
            i -= 1
            seq1_idx -= 1

        else:  # left: gap in seq1
            aligned1.append("-")
            aligned2.append(seq2[seq2_idx])
            score += gap
            j -= 1
            seq2_idx -= 1

    # Handle remaining characters in seq1
    while i > 0:
        aligned1.append(seq1[seq1_idx])
        aligned2.append("-")
        score += gap
        i -= 1
        seq1_idx -= 1

    # Handle remaining characters in seq2
    while j > 0:
        aligned1.append("-")
        aligned2.append(seq2[seq2_idx])
        score += gap
        j -= 1
        seq2_idx -= 1

    # Reverse and join lists to get final strings
    return ''.join(aligned1[::-1]), ''.join(aligned2[::-1]), score



