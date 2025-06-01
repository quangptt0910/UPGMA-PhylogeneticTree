# UPGMA Phylogenetic Tree Builder with MSA

## Overview

This project implements the Unweighted Pair Group Method with Arithmetic Mean (UPGMA) algorithm for constructing phylogenetic trees, integrating multiple sequence alignment (MSA) to derive distance matrices when users input raw sequences. The implementation is split into:

* **`upgma.py`**: Contains the core UPGMA algorithm and utility functions to compute distance matrices from alignment scores and to convert the resulting tree into Newick format.
* **`upgma_ui.py`**: A graphical user interface (GUI) built using Tkinter, which allows users to provide either a set of sequences (in FASTA format) or a custom distance matrix, run the UPGMA algorithm, and visualize the resulting tree.
* **`nw2.py`** and **`msa.py`**: Needleman-wunsch for global alignments sequences and MSA algorithm in used for this UPGMA.

---

## Features

* **UPGMA clustering**: Hierarchical clustering that builds an ultra-metric tree based on pairwise distances.
* **Flexible input**: Accepts raw nucleotide/protein sequences, distance matrix or load from fasta (sequences), load from csv(distance matrix).
* **Multiple Sequence Alignment (MSA) integration**: When given raw sequences, performs a center-star alignment to compute pairwise percent-identity and converts that into a distance matrix.
* **Newick output**: Exports the tree in Newick format for compatibility with other phylogenetics tools.
* **Interactive GUI**:

  * Load sequences from a FASTA file or paste directly into a text box.
  * Specify or load a CSV distance matrix with custom labels.
  * Display the Newick string and render the phylogenetic tree using Matplotlib.
  * Save inputs and results to a text file.

---

## Prerequisites

* **Python (in use) 3.10+**
* **Biopython**: For tree data structures, reading/writing Newick format, and tree plotting.
* **NumPy**: For numerical operations on matrices.
* **Tkinter**: Standard Python GUI library (usually included with most Python distributions).
* **Matplotlib**: For rendering the phylogenetic tree in the GUI.
* **`msa` module**: A custom or third-party module providing `center_star_alignment` and `compute_identity_matrix`. Ensure this is installed or available in the project folder.

Install the required packages via pip:

```bash
pip install numpy biopython matplotlib
```

---

## Repository Structure

```
├── upgma.py           # Core UPGMA implementation
├── upgma_ui.py        # Tkinter-based GUI for input, execution, and visualization
├── msa.py             # (Provided separately) MSA functions (center_star_alignment, compute_identity_matrix)
└── README.md          # This documentation
├── test.fasta         # fasta test files with sequences
├── dist_matrix        # csv test files
```


## Usage

### 1. Command-Line / Programmatic (without GUI)


#### Testing UPGMA

A simple built-in test is provided at the bottom of `upgma.py`:

```bash
python upgma.py
```

This will run two tests:

* **`test_upgma()`**: A minimal 3-taxon example.
* **`test_upgma_detailed()`**: A 4-taxon example with a pre-specified distance matrix and debug printouts.

---

### 2. Graphical User Interface (GUI)

To launch the GUI:

```bash
python upgma_ui.py
```

#### GUI Workflow

1. **Choose Input Method**

   * **Sequences**: Enter/paste multiple sequences in FASTA format (or load from a `.fasta` file).
   * **Distance Matrix**: Enter a CSV of distances (with header labels) or load from a `.csv`.

2. **For Sequences**:

   * The tool performs a center-star multiple sequence alignment (`center_star_alignment` from `msa.py`).
   * Converts the MSA into a pairwise identity matrix (`compute_identity_matrix`).
   * Calculates a distance matrix as `1 - identity`.

3. **For Distance Matrix**:

   * Enter labels (comma-separated) and fill the distance matrix cells.
   * The diagonal is fixed to zero; only the upper triangle needs to be filled (the lower side auto-mirrors).

4. **Run UPGMA**

   * Click the **Run UPGMA** button.
   * The Newick representation appears in the text field.
   * A graphical phylogenetic tree is rendered. Branch lengths are labeled on internal nodes.

5. **Save Results**

   * Click the **Save Results** button to output a `.txt` file containing:

     * Input (FASTA sequences or distance matrix).
     * The Newick string.
   * Choose a filename and location via the file dialog.

---

## Example Usage

1. **Using Sequences**

   * Prepare a FASTA file (`example.fasta`):

     ```fasta
     >Seq1
     ATGCTAGCTAG
     >Seq2
     ATG-CAGCTAG
     >Seq3
     ATGCTAG-TAG
     ```
   * Run the GUI:

     ```bash
     python upgma_ui.py
     ```
   * Switch to **Sequences**, load `example.fasta`, then click **Run UPGMA**.
   * The Newick tree and plot window will appear.

2. **Using a Distance Matrix**

   * Create a CSV (`dist_matrix.csv`):

     ```csv
     ,A, B, C, D, E
     A, 0, 16, 6, 16, 6
     B, 16, 0, 16, 8, 16
     C, 6, 16, 0, 16 ,2
     D, 16, 8, 16, 0, 16
     E, 6, 16, 2, 16, 0
     ```
   * Launch the GUI.
   * Switch to **Distance Matrix**, load `dist_matrix.csv`, then click **Run UPGMA**.
   * Visualize and save as needed.

---

## Functionality Details

### `upgma(dist_matrix, labels)`

* **Inputs**:

  * `dist_matrix` (NumPy array): Symmetric matrix of pairwise distances.
  * `labels` (list of strings): Labels for each taxon/sequence.
* **Process**:

  1. Initializes each sequence as its own cluster (Biopython `Clade` objects).
  2. Iteratively finds the pair of clusters with the minimum distance.
  3. Computes the branch lengths for the two child node placements (half the distance minus their current height).
  4. Merges them into a new internal node (`Clade` named `N0`, `N1`, ...).
  5. Updates distances using average (arithmetic mean) weighted by cluster sizes.
  6. Repeats until a single root node remains.
  7. Sets the final root’s branch length to 0 and returns a rooted `Tree`.
* **Returns**: A Biopython `Tree` object.

### `upgma.tree_to_newick(tree)`

* Converts a Biopython `Tree` to a Newick-formatted string using `Phylo.write`.

### `compute_identity_matrix(msa, method)`

* Converts a pairwise alignment score matrix into a distance matrix.

