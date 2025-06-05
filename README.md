# 🧬 EquilibraTor

**EquilibraTor** is a Python-based command-line tool that automates the setup and execution of molecular dynamics (MD) simulations for protein (and optionally ligand) systems using GROMACS. The pipeline runs from topology generation to energy minimization and equilibration, with customizable execution steps.

---

## 🚀 Overview

EquilibraTor streamlines MD preparation and execution in GROMACS by breaking the process into modular, automated steps. Users can run the entire pipeline or start/stop at specific steps — ideal for debugging, re-running specific stages, or customizing workflows.

---

## ⚙️ Usage

```text
python EquilibraTor.py -p example_protein.pdb -as

Available steps:
1: Generate topology for protein
2: Prepare to merge topology file(s) if ligand provided
3: Create the simulation box
4: Solvate the system
5: Add ions to neutralize the system
6: Run energy minimization
7: Plot potential energy
8: Obtain potential, backbone, and pressure xvgs
9: Plot panel of additional energy minimization results
10: Get final minimized pdb structure
11: Make refinement (Equilibrium)
12: Get refinement output
13: Run equilibration MD
14: Get equilibration MD output