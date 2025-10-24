# 🧬 EquilibraTor

EquilibraTor is a Python-based command-line tool that automates the setup and execution of molecular dynamics (MD) simulations for protein (and optionally ligand) systems using GROMACS. The pipeline runs from topology generation to energy minimization and equilibration, with customizable execution steps.

## 🛠️  Installation

---
### Option 1:

```bash
pip install EquilibraTor
```

### Option 2:

Clone this repository to your local machine:

```bash
# 1. clone the repository:
git clone https://github.com/Dias-Lab/EquilibraTor.git
cd EquilibraTor

# 2. create conda environment using yaml file and activate it. Use mamba instead of conda for faster installation:
   # with conda:
   conda env create -f equilibrator_env.yml
   conda activate EquilibraTor

# 3. install the python package
pip install -e .
```

---

## 🚀 Overview

<img src="paper/figures/equilibrator_workflow.png">

- **Protein PDB Preprocessing** — ⬛ 

- **Ligand PDB Preprocessing** *(Optional)* — ⬜ 

- **GROMACS Preprocessing** — 🟦 

- **Energy Minimization & Equilibration** — 🟧 

- **EquilibraTor Outputs** — 🟩 

---

## ⚙️ Usage

To show EquilibraTor arguments:

```text
EquilibraTor -h
```

```
usage: EquilibraTor [-h] [-l LIGANDS [LIGANDS ...]] -p PROTEIN [-aa {gaff,amber,gaff2,amber2}] [-an NET_CHARGE] [-cp]
                    [-gff {amber94,amber96,amber99,amber99sb,amber99sb-ildn,amber03}] [-gwm {spc,spce,tip3p,tip4p,tip5p}]
                    [-gbt {triclinic,cubic,dodecahedron,octahedron}] [-gd DISTANCE] [-gpi {NA,K,MG}] [-gni {CL,F,BR}] [-fs FIRST_STEP] [-ls LAST_STEP] [-as]

   ____          _ ___ __           ______        
  / __/__ ___ __(_) (_) /  _______ /_  __/__  ____
 / _// _ `/ // / / / / _ \/ __/ _ `// / / _ \/ __/
/___/\_, /\_,_/_/_/_/_.__/_/  \_,_//_/  \___/_/
      /_/
Equilibrator streamlines Molecular dynamics and equilibration simulations for proteins and protein-ligand complexes in a single execution
Developers: José D. D. Cediel-Becerra, Jose Cleydson F. Silva and Raquel Dias
Afiliation: Microbiology & Cell Science Deparment, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Version:v1.0.0

options:
  -h, --help            show this help message and exit

Input options:
  -l LIGANDS [LIGANDS ...], --ligands LIGANDS [LIGANDS ...]
                        Path(s) to the ligand file(s).
  -p PROTEIN, --protein PROTEIN
                        Path to the protein file.

Specify options for ligand topology generation with acpype:
  -aa {gaff,amber,gaff2,amber2}, --atom_type {gaff,amber,gaff2,amber2}
                        Specify the atom type supported by acpype: gaff, amber, gaff2 (default), amber2
  -an NET_CHARGE, --net_charge NET_CHARGE
                        net molecular charge (int), default is -an=0

Protein termini capping before protein topology generation:
  -cp, --cap_protein    Add ACE and NME terminal capping groups to the input protein PDB

Specify options for protein topology generation with gromacs:
  -gff {amber94,amber96,amber99,amber99sb,amber99sb-ildn,amber03}, --force_field {amber94,amber96,amber99,amber99sb,amber99sb-ildn,amber03}
                        Specify the Force Fields supported by GROMACS: amber94, amber96, amber99, amber99sb (default), amber99sb-ildn, amber03
  -gwm {spc,spce,tip3p,tip4p,tip5p}, --water_model {spc,spce,tip3p,tip4p,tip5p}
                        Specify the water model: spc, spce, tip3p (default), tip4p, tip5p

Specify options for box generation with gromacs:
  -gbt {triclinic,cubic,dodecahedron,octahedron}, --box_type {triclinic,cubic,dodecahedron,octahedron}
                        Specify the box type supported by GROMACS: triclinic, cubic (default), dodecahedron, octahedron
  -gd DISTANCE, --distance DISTANCE
                        Specify the distance between the solute and the box, default is -gd 1.2

Specify monoatomic cation/anion supported by the force field:
  -gpi {NA,K,MG}, --pos_ion {NA,K,MG}
                        Specify the monoatomic cation supported by the force field: NA (default), K, MG, etc
  -gni {CL,F,BR}, --neg_ion {CL,F,BR}
                        Specify the monoatomic anion supported by the force field: CL (default), F, BR, etc

Specify the steps for the execution:
  -fs FIRST_STEP, --first_step FIRST_STEP
                        Step number to start Equilibratior from (1-based)
  -ls LAST_STEP, --last_step LAST_STEP
                        Step number to end at (1-based)
  -as, --all_steps      List of Equilibrator steps and exit
```

## 📌 Note

- Only the protein file is required as input. 
- When the ligand is a polypeptide or another protein, it must be combined with the main protein into a single PDB file and provided using the -p option.
- Default .mdp files (ions.mdp, minim.mdp, equilibration.mdp, equilibration_2.mdp) are in equilibrator/flat. Modify them directly to change parameters.


To show the EquilibraTor steps to be performed for a protein file:

```Text
EquilibraTor -p example/example_protein.pdb -as

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
11: Run NVT equilibration
12: Get NVT equilibration output
13: Run NPT equilibration
14: Get NPT equilibration output
```
To run EquilibraTor for a protein file:

```Text
EquilibraTor -p example/example_protein.pdb
````

If you want to tweak certain parameters for your current protein file and avoid running the entire Equilibrator workflow, you can specify only the steps you wish to execute:

```Text
EquilibraTor -p example/example_protein.pdb -fs 10 -ls 14
```

To show the EquilibraTor steps to be performed when provided both protein and ligand files:

```Text
EquilibraTor -l example/example_ligand.pdb -p example/example_protein.pdb -as

Available steps:
1: Generate topology for protein
2: Convert ligand PDB to MOL2
3: Generate topology for ligand
4: Prepare to merge topology file(s) if ligand provided
5: Make a copy of protein if ligand provided
6: Merge topologies
7: Create the simulation box
8: Solvate the system
9: Add ions to neutralize the system
10: Run energy minimization
11: Plot potential energy
12: Obtain potential, backbone, and pressure xvgs
13: Plot panel of additional energy minimization results
14: Get final minimized pdb structure
15: Run NVT equilibration
16: Get NVT equilibration output
17: Run NPT equilibration
18: Get NPT equilibration output
```

To run EquilibraTor using this protein-ligand files:

```Text
EquilibraTor -l example/example_ligand.pdb -p example/example_protein.pdb
```


## 💾 Outputs

### Multiple panel figure
<img src="paper/figures/protein-ligand_equilibration.png">

### Equilibration PDB last frame
<img src="paper/figures/protein-ligand_eq_last_frame.png">
