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
- Multiple ligands can be provided, each into a separate PDB file
- When the ligand is a polypeptide or another protein, it must be combined with the main protein into a single PDB file and provided using the -p option.
- Default .mdp files (ions.mdp, minimization_stage.mdp, nvt_stage.mdp, npt_stage.mdp, production_stage.mdp) are in equilibrator/flat. Modify them directly to change parameters.


To show the EquilibraTor steps to be performed for a protein file:

```Text
EquilibraTor -p example/example_protein.pdb -as

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
Available steps:
1: Generating topology for the protein: example_protein
2: Checking wether merging topology file(s) is necessary
3: Combining and inserting unique atomtypes into main topology
4: Creating the simulation box
5: Solvating the system
6: Adding ions to neutralize the system
7: Running energy minimization
8: Plotting potential energy
9: Obtaining potential, backbone, and pressure xvgs
10: Plotting panel of additional energy minimization results
11: Getting final minimized pdb structure
12: Running NVT equilibration
13: Getting NVT equilibration output
14: Running NPT equilibration
15: Getting NPT equilibration output
16: Running Production stage
17: Getting Production output
```


To show the EquilibraTor steps to be performed for a protein file, when the protein capping feature is enabled:

```
EquilibraTor -p example/example_protein.pdb -cp -as
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
Available steps:
1: Adding ACE and NME terminal capping groups to the protein: example_protein
2: Generating topology for the protein: example_protein_capped.pdb
3: Checking wether merging topology file(s) is necessary
4: Combining and inserting unique atomtypes into main topology
5: Creating the simulation box
6: Solvating the system
7: Adding ions to neutralize the system
8: Running energy minimization
9: Plotting potential energy
10: Obtaining potential, backbone, and pressure xvgs
11: Plotting panel of additional energy minimization results
12: Getting final minimized pdb structure
13: Running NVT equilibration
14: Getting NVT equilibration output
15: Running NPT equilibration
16: Getting NPT equilibration output
17: Running Production stage
18: Getting Production output
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
Available steps:
1: Generating topology for the protein: example_protein
2: Converting example_ligand PDB to MOL2
3: Generating topology for the ligand: example_ligand
4: Checking wether merging topology file(s) is necessary
5: Making a copy of the protein: example/example_protein.pdb
6: Merging topologies
7: Combining and inserting unique atomtypes into main topology
8: Creating the simulation box
9: Solvating the system
10: Adding ions to neutralize the system
11: Running energy minimization
12: Plotting potential energy
13: Obtaining potential, backbone, and pressure xvgs
14: Plotting panel of additional energy minimization results
15: Getting final minimized pdb structure
16: Running NVT equilibration
17: Getting NVT equilibration output
18: Running NPT equilibration
19: Getting NPT equilibration output
20: Running Production stage
21: Getting Production output
```
To show the EquilibraTor steps to be performed when provided both protein and ligand files, enabling the protein capping feature:

```
EquilibraTor -l example/example_ligand.pdb -p example/example_protein.pdb -cp -as

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
Available steps:
1: Adding ACE and NME terminal capping groups to the protein: example_protein
2: Generating topology for the protein: example_protein_capped.pdb
3: Converting example_ligand PDB to MOL2
4: Generating topology for the ligand: example_ligand
5: Checking wether merging topology file(s) is necessary
6: Making a copy of the protein: example/example_protein.pdb
7: Merging topologies
8: Combining and inserting unique atomtypes into main topology
9: Creating the simulation box
10: Solvating the system
11: Adding ions to neutralize the system
12: Running energy minimization
13: Plotting potential energy
14: Obtaining potential, backbone, and pressure xvgs
15: Plotting panel of additional energy minimization results
16: Getting final minimized pdb structure
17: Running NVT equilibration
18: Getting NVT equilibration output
19: Running NPT equilibration
20: Getting NPT equilibration output
21: Running Production stage
22: Getting Production output

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
