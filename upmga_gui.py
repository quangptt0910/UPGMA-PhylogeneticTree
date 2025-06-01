import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio import Phylo
from msa import *
from upgma import *

class UPGMA_GUI:
    def __init__(self, root):
        self.root = root
        self.root.title("UPGMA Phylogenetic Tree Builder")
        self.root.geometry("1000x800")

        # Input selection
        self.input_method = tk.IntVar(value=0)  # 0=Sequences, 1=Matrix
        input_frame = ttk.Frame(root)
        input_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(input_frame, text="Input Method:").grid(row=0, column=0, padx=5)
        ttk.Radiobutton(input_frame, text="Sequences", variable=self.input_method, value=0,
                        command=self.toggle_input).grid(row=0, column=1)
        ttk.Radiobutton(input_frame, text="Distance Matrix", variable=self.input_method, value=1,
                        command=self.toggle_input).grid(row=0, column=2)

        # Sequence input frame
        self.seq_frame = ttk.LabelFrame(root, text="Sequence Input (FASTA Format)")
        self.seq_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.seq_text = tk.Text(self.seq_frame, height=10, width=80)
        self.seq_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        seq_btn_frame = ttk.Frame(self.seq_frame)
        seq_btn_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(seq_btn_frame, text="Load FASTA File",
                   command=self.load_fasta).pack(side=tk.LEFT)
        ttk.Button(seq_btn_frame, text="Clear",
                   command=self.clear_sequences).pack(side=tk.LEFT, padx=5)

        # Matrix input frame
        self.matrix_frame = ttk.LabelFrame(root, text="Distance Matrix Input")
        self.matrix_frame.pack_forget()

        matrix_top = ttk.Frame(self.matrix_frame)
        matrix_top.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(matrix_top, text="Labels (comma separated):").pack(side=tk.LEFT)
        self.matrix_labels = ttk.Entry(matrix_top, width=50)
        self.matrix_labels.pack(side=tk.LEFT, padx=5)
        self.matrix_labels.insert(0, "Seq1,Seq2,Seq3")

        btn_frame = ttk.Frame(matrix_top)
        btn_frame.pack(side=tk.RIGHT)
        ttk.Button(btn_frame, text="Load CSV", command=self.load_csv).pack(side=tk.LEFT)
        ttk.Button(btn_frame, text="Add Row/Col", command=self.add_matrix_entry).pack(side=tk.LEFT, padx=5)

        # Matrix table frame
        self.table_frame = ttk.Frame(self.matrix_frame)
        self.table_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.matrix_entries = []
        self.init_matrix_table(3)

        # Control buttons
        ctrl_frame = ttk.Frame(root)
        ctrl_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(ctrl_frame, text="Run UPGMA", command=self.run_upgma).pack(side=tk.LEFT)
        ttk.Button(ctrl_frame, text="Save Results", command=self.save_results).pack(side=tk.LEFT, padx=5)

        # Output display
        self.output_frame = ttk.LabelFrame(root, text="Results")
        self.output_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Newick output
        newick_frame = ttk.Frame(self.output_frame)
        newick_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(newick_frame, text="Newick Format:").pack(side=tk.LEFT)
        self.newick_var = tk.StringVar()
        newick_entry = ttk.Entry(newick_frame, textvariable=self.newick_var, width=80)
        newick_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Tree visualization
        viz_frame = ttk.Frame(self.output_frame)
        viz_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.ax.set_title("Phylogenetic Tree Visualization")
        self.ax.axis('off')
        self.canvas.draw()

        # Initial layout
        self.toggle_input()

    def toggle_input(self):
        if self.input_method.get() == 0:
            self.matrix_frame.pack_forget()
            self.seq_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        else:
            self.seq_frame.pack_forget()
            self.matrix_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

    def clear_sequences(self):
        self.seq_text.delete("1.0", tk.END)

    def load_fasta(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        if not file_path:
            return

        try:
            with open(file_path, 'r') as f:
                fasta_content = f.read()
            self.seq_text.delete("1.0", tk.END)
            self.seq_text.insert("1.0", fasta_content)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")

    def parse_fasta_text(self):
        text = self.seq_text.get("1.0", tk.END).strip()
        if not text:
            return {}

        sequences = {}
        current_name = ""
        current_seq = []

        for line in text.split('\n'):
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            elif line.strip():
                current_seq.append(line.strip())

        if current_name:
            sequences[current_name] = ''.join(current_seq)

        return sequences

    def init_matrix_table(self, size):
        # Clear existing widgets
        for widget in self.table_frame.winfo_children():
            widget.destroy()

        self.matrix_entries = []

        # Create header row
        for j in range(size):
            ttk.Label(self.table_frame, text=f"Seq{j + 1}", width=8).grid(row=0, column=j + 1, padx=2, pady=2)

        # Create rows
        for i in range(size):
            row_entries = []
            ttk.Label(self.table_frame, text=f"Seq{i + 1}", width=8).grid(row=i + 1, column=0, padx=2, pady=2)

            for j in range(size):
                entry = ttk.Entry(self.table_frame, width=8)
                entry.grid(row=i + 1, column=j + 1, padx=2, pady=2)

                if i == j:
                    entry.insert(0, "0")
                    entry.config(state='disabled')
                elif j < i:
                    # Mirror upper triangle
                    entry.config(state='disabled')

                row_entries.append(entry)
            self.matrix_entries.append(row_entries)

    def add_matrix_entry(self):
        size = len(self.matrix_entries) + 1
        self.init_matrix_table(size)

        # Update labels
        current_labels = self.matrix_labels.get().split(',')
        if len(current_labels) < size:
            current_labels.extend([f"Seq{len(current_labels) + 1}" for _ in range(size - len(current_labels))])
            self.matrix_labels.delete(0, tk.END)
            self.matrix_labels.insert(0, ','.join(current_labels))

    def load_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not file_path:
            return

        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()

            # Parse labels from first line
            if lines:
                labels = [label.strip() for label in lines[0].split(',')]
                self.matrix_labels.delete(0, tk.END)
                self.matrix_labels.insert(0, ','.join(labels[1:]))

            # Resize matrix
            size = len(lines) - 1
            self.init_matrix_table(size)

            # Fill matrix values
            for i, line in enumerate(lines[1:], 1):
                values = line.split(',')[1:]
                for j, val in enumerate(values[:size]):
                    if j < i:
                        # Update lower triangle
                        self.matrix_entries[i - 1][j].config(state='normal')
                        self.matrix_entries[i - 1][j].delete(0, tk.END)
                        self.matrix_entries[i - 1][j].insert(0, val.strip())
                        # Mirror to upper triangle
                        if i - 1 != j:
                            self.matrix_entries[j][i - 1].config(state='normal')
                            self.matrix_entries[j][i - 1].delete(0, tk.END)
                            self.matrix_entries[j][i - 1].insert(0, val.strip())
                            self.matrix_entries[j][i - 1].config(state='disabled')
                        self.matrix_entries[i - 1][j].config(state='disabled' if j < i else 'normal')

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load CSV:\n{str(e)}")

    def parse_matrix(self):
        # Get labels
        labels = [label.strip() for label in self.matrix_labels.get().split(',')]
        size = len(self.matrix_entries)

        if len(labels) < size:
            labels.extend([f"Seq{len(labels) + 1}" for _ in range(size - len(labels))])

        # Create matrix
        dist_matrix = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                if i != j and j >= i:
                    try:
                        val = float(self.matrix_entries[i][j].get())
                        dist_matrix[i][j] = val
                        dist_matrix[j][i] = val
                    except ValueError:
                        messagebox.showerror("Error",
                                             f"Invalid number at position ({i + 1},{j + 1})")
                        return None, None

        return dist_matrix, labels

    def run_upgma(self):
        if self.input_method.get() == 0:  # Sequences
            sequences = self.parse_fasta_text()
            if not sequences:
                messagebox.showerror("Error", "No valid sequences found")
                return

            # Get sequences and labels
            labels = list(sequences.keys())
            seqs = list(sequences.values())

            try:
                msa, score_matrix, center_idx = center_star_alignment(
                    seqs,
                    match=1,
                    mismatch=-1,
                    gap=-2
                )
            except Exception as e:
                messagebox.showerror("Error", f"Failed to align sequences:\n{str(e)}")
                return

            # Compute identity matrix
            dist_matrix = score_to_distance_matrix(score_matrix, method='simple_inverse')

        else:  # Distance matrix+
            dist_matrix, labels = self.parse_matrix()
            if dist_matrix.size is None:
                return

        # Run UPGMA algorithm
        try:
            tree = UPGMATreeBuilder.upgma(dist_matrix, labels)
            # DEBUG: Print branch lengths
            print("\nBranch Lengths:")
            for clade in tree.find_clades():
                if clade.branch_length is not None:
                    print(f"{clade.name or 'Root'}: {clade.branch_length:.6f}")
            newick_str = UPGMATreeBuilder.tree_to_newick(tree)
            self.newick_var.set(newick_str)
            self.draw_tree(tree)
        except Exception as e:
            messagebox.showerror("Error", f"UPGMA failed: {str(e)}")

    def draw_tree(self, tree):
        self.ax.clear()
        self.ax.set_title("Phylogenetic Tree Visualization", fontsize=12)

        # First, layout the tree to calculate coordinates
        tree.ladderize()  # Improve tree layout

        # Draw tree using Bio.Phylo
        Phylo.draw(
            tree,
            axes=self.ax,
            branch_labels=lambda c: f"{c.branch_length:.4f}"
                if c.branch_length is not None
                   and c.branch_length > 0
                else "",
            label_func=lambda c: c.name if c.name else "",
            do_show=False
        )

        # Adjust layout to prevent cropping
        try:
            # Get all terminal nodes to calculate max label width
            terminal_nodes = tree.get_terminals()
            if terminal_nodes:
                max_label_length = max(len(str(clade.name)) for clade in terminal_nodes)

                # Calculate padding based on label length
                padding = max_label_length * 0.02  # Adjust based on font size

                # Adjust axis limits with padding
                x_min, x_max = self.ax.get_xlim()
                y_min, y_max = self.ax.get_ylim()

                # Add padding for labels
                self.ax.set_xlim(x_min, x_max + padding)
                self.ax.set_ylim(y_min - 0.5, y_max + 0.5)

        except Exception as e:
            print(f"Error adjusting tree layout: {e}")

        # Improve overall figure layout
        self.fig.tight_layout(pad=3.0)  # Add extra padding around the plot
        self.canvas.draw()

    def save_results(self):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if not file_path:
            return

        try:
            with open(file_path, 'w') as f:
                # Save input data
                if self.input_method.get() == 0:
                    f.write("=== INPUT SEQUENCES ===\n\n")
                    f.write(self.seq_text.get("1.0", tk.END))
                else:
                    f.write("=== DISTANCE MATRIX ===\n\n")
                    labels = self.matrix_labels.get()
                    f.write(f"Labels: {labels}\n\n")

                    size = len(self.matrix_entries)
                    for i in range(size):
                        row = []
                        for j in range(size):
                            row.append(self.matrix_entries[i][j].get())
                        f.write("\t".join(row) + "\n")

                # Save results
                f.write("\n\n=== UPGMA RESULTS ===\n\n")
                f.write(f"Newick Format: {self.newick_var.get()}\n")

            messagebox.showinfo("Success", f"Results saved to:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file:\n{str(e)}")

# Test UPGMATreeBuilder
def test_upgma():
    dist_matrix = np.array([
        [0, 2, 4],
        [2, 0, 4],
        [4, 4, 0]
    ])
    labels = ['A', 'B', 'C']
    tree = UPGMATreeBuilder.upgma(dist_matrix, labels)
    newick = UPGMATreeBuilder.tree_to_newick(tree)
    print("Newick:", newick)
    assert 'A' in newick and 'B' in newick and 'C' in newick

# if __name__ == "__main__":
#     test_upgma()

if __name__ == "__main__":
    root = tk.Tk()
    app = UPGMA_GUI(root)
    root.mainloop()