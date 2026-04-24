# Interactive Protein Sequence Alignment Tool

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![JavaScript](https://img.shields.io/badge/JavaScript-ES6+-yellow.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

ProAlign is an interactive educational and analytical bioinformatics tool designed to perform and visualize protein sequence alignments. It utilizes Dynamic Programming to implement both **Global (Needleman-Wunsch)** and **Local (Smith-Waterman)** alignment algorithms using the **BLOSUM62** substitution matrix.

This repository includes two fully functional versions of the tool:
1. **Desktop GUI Application:** Built with Python (Tkinter & Matplotlib).
2. **Serverless Web Application:** Built with pure HTML, CSS, and JavaScript (Plotly.js).

---

## Key Features

* **Dual Algorithms:** Choose between Global (Needleman-Wunsch) or Local (Smith-Waterman) alignment.
* **Interactive DP Matrix:** Visualizes the dynamic programming scoring matrix as a heat map and overlays the exact traceback path.
* **FASTA Support:** Instantly parse and load sequences from `.fasta` or `.txt` files directly into the UI.
* **Sequence Slicing (Bounds):** Option to align specific subsequences by defining Start and End boundaries (1-based biological indexing).
* **Dynamic Visualization Controls:** Adjust traceback line thickness and marker visibility for optimal viewing of large matrices.
* **Web Version:** All calculations, including the BLOSUM62 matrix processing, run directly in the browser via JavaScript.

---

## Web Version

Try the web version directly from your browser without installing anything!
**[Click here to open the Live Demo](https://krefikk.github.io/ProteinSequenceAligner/)**