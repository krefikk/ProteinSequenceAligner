import os
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import matplotlib.pyplot as plt
import numpy as np
from alignment import slice_sequence, align_dp, format_alignment_text, parse_fasta

class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Protein Sequence Alignment Tool")
        self.root.geometry("950x800")
        
        self.line_thickness = tk.DoubleVar(value=0.5) 
        self.show_markers = tk.BooleanVar(value=False) 
        
        self.current_matrix = None
        self.current_path = None
        self.current_seq_a = None
        self.current_seq_b = None
        
        style = ttk.Style()
        style.theme_use('clam')
        
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Inputs
        input_frame = ttk.LabelFrame(main_frame, text="Alignment Parameters", padding="10")
        input_frame.pack(fill=tk.X, pady=(0, 10))
        
        # First Sequence
        seq_a_header_frame = ttk.Frame(input_frame)
        seq_a_header_frame.grid(row=0, column=0, columnspan=4, sticky=tk.W, pady=(0,2))
        ttk.Label(seq_a_header_frame, text="Sequence A:").pack(side=tk.LEFT)
        
        self.seq_a_file_label = ttk.Label(seq_a_header_frame, text="", foreground="blue", font=('Helvetica', 9, 'italic'))
        
        ttk.Button(seq_a_header_frame, text="Load from FASTA", command=lambda: self.load_fasta(self.seq_a_entry, self.seq_a_file_label)).pack(side=tk.LEFT, padx=10)
        self.seq_a_file_label.pack(side=tk.LEFT, padx=5)

        self.seq_a_entry = scrolledtext.ScrolledText(input_frame, height=3, width=80)
        self.seq_a_entry.grid(row=1, column=0, columnspan=4, pady=2, padx=5, sticky=tk.W)
        
        ttk.Label(input_frame, text="Bounds A (Start - End):").grid(row=2, column=0, sticky=tk.W, pady=2, padx=5)
        bound_a_frame = ttk.Frame(input_frame)
        bound_a_frame.grid(row=2, column=1, sticky=tk.W, pady=2)
        self.seq_a_start = ttk.Entry(bound_a_frame, width=6)
        self.seq_a_start.pack(side=tk.LEFT)
        ttk.Label(bound_a_frame, text=" to ").pack(side=tk.LEFT)
        self.seq_a_end = ttk.Entry(bound_a_frame, width=6)
        self.seq_a_end.pack(side=tk.LEFT)
        ttk.Label(bound_a_frame, text=" (Leave blank for full length)", foreground="gray").pack(side=tk.LEFT, padx=10)

        # Second Sequence
        seq_b_header_frame = ttk.Frame(input_frame)
        seq_b_header_frame.grid(row=3, column=0, columnspan=4, sticky=tk.W, pady=(10,2))
        ttk.Label(seq_b_header_frame, text="Sequence B:").pack(side=tk.LEFT)
        
        self.seq_b_file_label = ttk.Label(seq_b_header_frame, text="", foreground="blue", font=('Helvetica', 9, 'italic'))
        
        ttk.Button(seq_b_header_frame, text="Load from FASTA", command=lambda: self.load_fasta(self.seq_b_entry, self.seq_b_file_label)).pack(side=tk.LEFT, padx=10)
        self.seq_b_file_label.pack(side=tk.LEFT, padx=5)

        self.seq_b_entry = scrolledtext.ScrolledText(input_frame, height=3, width=80)
        self.seq_b_entry.grid(row=4, column=0, columnspan=4, pady=2, padx=5, sticky=tk.W)
        
        ttk.Label(input_frame, text="Bounds B (Start - End):").grid(row=5, column=0, sticky=tk.W, pady=2, padx=5)
        bound_b_frame = ttk.Frame(input_frame)
        bound_b_frame.grid(row=5, column=1, sticky=tk.W, pady=2)
        self.seq_b_start = ttk.Entry(bound_b_frame, width=6)
        self.seq_b_start.pack(side=tk.LEFT)
        ttk.Label(bound_b_frame, text=" to ").pack(side=tk.LEFT)
        self.seq_b_end = ttk.Entry(bound_b_frame, width=6)
        self.seq_b_end.pack(side=tk.LEFT)
        ttk.Label(bound_b_frame, text=" (Leave blank for full length)", foreground="gray").pack(side=tk.LEFT, padx=10)
        
        # Mode and gap penalty
        settings_frame = ttk.Frame(input_frame)
        settings_frame.grid(row=6, column=0, columnspan=4, sticky=tk.W, pady=(15, 0))

        ttk.Label(settings_frame, text="Mode:").pack(side=tk.LEFT, padx=(5,0))
        self.mode_var = tk.StringVar(value="global")
        ttk.Radiobutton(settings_frame, text="Global", variable=self.mode_var, value="global").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(settings_frame, text="Local", variable=self.mode_var, value="local").pack(side=tk.LEFT, padx=(5, 40))
        
        ttk.Label(settings_frame, text="Gap Penalty:").pack(side=tk.LEFT)
        self.gap_var = tk.StringVar(value="10")
        ttk.Entry(settings_frame, textvariable=self.gap_var, width=10).pack(side=tk.LEFT, padx=5)
        
        # Buttons
        btn_frame = ttk.Frame(input_frame)
        btn_frame.grid(row=7, column=0, columnspan=4, pady=15)
        
        self.run_btn = ttk.Button(btn_frame, text="Run Alignment", command=self.run_alignment)
        self.run_btn.pack(side=tk.LEFT, padx=5)
        
        self.plot_btn = ttk.Button(btn_frame, text="Show DP Matrix Plot", command=self.plot_matrix, state=tk.DISABLED)
        self.plot_btn.pack(side=tk.LEFT, padx=5)

        self.settings_btn = ttk.Button(btn_frame, text="⚙ Plot Settings", command=self.open_settings)
        self.settings_btn.pack(side=tk.LEFT, padx=20)
        
        # Results
        output_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        output_frame.pack(fill=tk.BOTH, expand=True)
        
        self.score_label = ttk.Label(output_frame, text="Alignment Score: -", font=('Helvetica', 12, 'bold'))
        self.score_label.pack(anchor=tk.W, pady=(0, 10))
        
        self.result_text = scrolledtext.ScrolledText(output_frame, font=("Courier", 10), bg="#f4f4f4", wrap=tk.NONE)
        self.result_text.pack(fill=tk.BOTH, expand=True)
        
        x_scroll = ttk.Scrollbar(output_frame, orient=tk.HORIZONTAL, command=self.result_text.xview)
        x_scroll.pack(fill=tk.X)
        self.result_text['xscrollcommand'] = x_scroll.set

    def load_fasta(self, target_entry_widget, target_label_widget):
        filepath = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=(("FASTA Files", "*.fasta *.fas *.fa *.txt"), ("All Files", "*.*"))
        )
        if not filepath:
            return # Canceled
            
        sequence = parse_fasta(filepath)
        
        if sequence:
            target_entry_widget.delete("1.0", tk.END)
            target_entry_widget.insert(tk.END, sequence)
            
            filename = os.path.basename(filepath)
            target_label_widget.config(text=f"Loaded: {filename}")
        else:
            messagebox.showerror("Error", "Could not read sequence from the selected file. Ensure it is a valid FASTA format.")

    def open_settings(self):
        settings_win = tk.Toplevel(self.root)
        settings_win.title("Plot Settings")
        settings_win.geometry("300x160")
        settings_win.resizable(False, False)
        settings_win.grab_set() 
        
        ttk.Label(settings_win, text="Traceback Line Thickness:").pack(pady=(15, 5))
        thickness_scale = tk.Scale(settings_win, from_=0.1, to_=5.0, resolution=0.1, orient=tk.HORIZONTAL, variable=self.line_thickness)
        thickness_scale.pack(fill=tk.X, padx=30)
        
        marker_check = ttk.Checkbutton(settings_win, text="Show Traceback Markers (Dots)", variable=self.show_markers)
        marker_check.pack(pady=15)

    def parse_bound(self, value):
        val = value.strip()
        if not val:
            return None
        return int(val)

    def run_alignment(self):
        seq_a = self.seq_a_entry.get("1.0", tk.END).strip().replace('\n', '').replace(' ', '').upper()
        seq_b = self.seq_b_entry.get("1.0", tk.END).strip().replace('\n', '').replace(' ', '').upper()
        mode = self.mode_var.get()
        
        if not seq_a or not seq_b:
            messagebox.showwarning("Missing Sequence", "Please provide both sequences. You can paste them or load from a FASTA file.")
            return

        try:
            gap_penalty = float(self.gap_var.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Gap penalty must be a number.")
            return

        try:
            start_a = self.parse_bound(self.seq_a_start.get())
            end_a = self.parse_bound(self.seq_a_end.get())
            start_b = self.parse_bound(self.seq_b_start.get())
            end_b = self.parse_bound(self.seq_b_end.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Bounds must be integer numbers (e.g., 5 to 30).")
            return

        def check_bounds(start, end, seq_len, seq_name):
            if start is not None:
                if start < 1:
                    return f"Start index for {seq_name} must be 1 or greater."
                if start > seq_len:
                    return f"Start index for {seq_name} ({start}) exceeds sequence length ({seq_len})."
            
            if end is not None:
                if end < 1:
                    return f"End index for {seq_name} must be 1 or greater."
                if end > seq_len:
                    return f"End index for {seq_name} ({end}) exceeds sequence length ({seq_len})."
            
            if start is not None and end is not None:
                if start > end:
                    return f"Start index ({start}) cannot be greater than End index ({end}) for {seq_name}."
            return None 

        error_a = check_bounds(start_a, end_a, len(seq_a), "Sequence A")
        if error_a:
            messagebox.showerror("Bound Error", error_a)
            return
            
        error_b = check_bounds(start_b, end_b, len(seq_b), "Sequence B")
        if error_b:
            messagebox.showerror("Bound Error", error_b)
            return

        self.run_btn.config(text="Calculating...", state=tk.DISABLED)
        self.plot_btn.config(state=tk.DISABLED)
        self.root.update()
        
        try:
            score, align1, lines, align2, matrix, path = align_dp(
                seq_a, seq_b, 
                mode=mode, 
                gap_penalty=gap_penalty,
                seq1_bounds=(start_a, end_a),
                seq2_bounds=(start_b, end_b)
            )
            
            self.current_matrix = matrix
            self.current_path = path
            self.current_seq_a = slice_sequence(seq_a, start_a, end_a)
            self.current_seq_b = slice_sequence(seq_b, start_b, end_b)
            
            bounds_text_a = f"[{start_a or 1}:{end_a or len(seq_a)}]"
            bounds_text_b = f"[{start_b or 1}:{end_b or len(seq_b)}]"
            
            self.score_label.config(text=f"Alignment Score: {score}  |  Seq A Sliced: {bounds_text_a}  |  Seq B Sliced: {bounds_text_b}")
            formatted_text = format_alignment_text(align1, lines, align2, chunk_size=160)
            
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, formatted_text)
            
            self.plot_btn.config(state=tk.NORMAL)
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during alignment:\n{str(e)}")
            
        finally:
            self.run_btn.config(text="Run Alignment", state=tk.NORMAL)

    def plot_matrix(self):
        if not self.current_matrix or not self.current_path:
            return
            
        mat = np.array(self.current_matrix)
        path = self.current_path
        seq1 = self.current_seq_a
        seq2 = self.current_seq_b
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        cax = ax.matshow(mat, cmap='coolwarm')
        fig.colorbar(cax, fraction=0.046, pad=0.04, label="Score")
        
        y_vals = [p[0] for p in path]
        x_vals = [p[1] for p in path]
        
        line_w = self.line_thickness.get()
        show_dots = self.show_markers.get()
        
        if show_dots:
            mark_style = 'o'
            mark_s = 4 if max(len(seq1), len(seq2)) <= 50 else 2 
        else:
            mark_style = ''
            mark_s = 0

        ax.plot(x_vals, y_vals, color='black', linewidth=line_w, marker=mark_style, markersize=mark_s, label='Traceback Path')
        
        if len(seq1) <= 50 and len(seq2) <= 50:
            ax.set_xticks(np.arange(len(seq2) + 1))
            ax.set_yticks(np.arange(len(seq1) + 1))
            ax.set_xticklabels(['-'] + list(seq2))
            ax.set_yticklabels(['-'] + list(seq1))
            
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    ax.text(j, i, str(int(mat[i, j])), va='center', ha='center', color='black', fontsize=8)
        else:
            ax.set_title(f"Dynamic Programming Matrix ({len(seq1)} x {len(seq2)})\nTraceback Path in Black", pad=20)
            ax.set_xlabel("Sequence B (Sliced)")
            ax.set_ylabel("Sequence A (Sliced)")
            
        plt.legend(loc="upper right")
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()