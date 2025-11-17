# 🧬 EquilibraTor

EquilibraTor is a Python-based command-line tool that automates the setup and execution of molecular dynamics (MD) simulations for protein (and optionally ligand) systems using GROMACS. The pipeline runs from topology generation to energy minimization and equilibration, with customizable execution steps.

Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034

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

- **Energy Minimization, NVT/NPT stages, production stage, and convergence** — 🟧 

- **EquilibraTor Outputs** — 🟩 

---

## ⚙️ Usage

To show EquilibraTor arguments:

```text
EquilibraTor -h
```

```
usage: EquilibraTor [-h] [-l LIGANDS [LIGANDS ...]] -p PROTEIN [-aa {gaff,amber,gaff2,amber2}] [-an NET_CHARGE] [-cp] [-gff {amber94,amber96,amber99,amber99sb,amber99sb-ildn,amber03}] [-gwm {spc,spce,tip3p,tip4p,tip5p}] [-gbt {triclinic,cubic,dodecahedron,octahedron}] [-gd DISTANCE] [-gpi {NA,K,MG}]
                    [-gni {CL,F,BR}] [-oth OBS_THRESHOLDS [OBS_THRESHOLDS ...]] [-rao] [-bs BLOCK_SIZE] [-wsb WINDOW_SIZE_BLOCKS] [-ov OVERLAP] [-rc RMSD_CUTOFF] [-fs FIRST_STEP] [-ls LAST_STEP] [-as]

   ____          _ ___ __           ______        
  / __/__ ___ __(_) (_) /  _______ /_  __/__  ____
 / _// _ `/ // / / / / _ \/ __/ _ `// / / _ \/ __/
/___/\_, /\_,_/_/_/_/_.__/_/  \_,_//_/  \___/_/
      /_/
Equilibrator streamlines Molecular dynamics and equilibration simulations for proteins and protein-ligand complexes in a single execution
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
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
                        Specify the Force Fields supported by GROMACS
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

Automated checks to identify representative structures from converged trajectories:
  -oth OBS_THRESHOLDS [OBS_THRESHOLDS ...], --obs_thresholds OBS_THRESHOLDS [OBS_THRESHOLDS ...]
                        Relative drift thresholds for slope/drift equilibration detection (default: 10% drift, i.e. 0.10) Specify as observable=value pairs, e.g., rmsd=0.1 potential=0.2 rgyr=0.1 pressure=0.2. For relative thresholds, the value is multiplied by the standard deviation of the last half of the trajectory. If omitted, default values are used for all enabled observables. Example: rmsd=0.10 sets the RMSD drift tolerance to 10% of the reference noise.
  -rao, --req_all_obs   If set, fill in missing observables with defaults. Otherwise, only use user-specified observables.
  -bs BLOCK_SIZE, --block_size BLOCK_SIZE
                        The amount of time (in ps) for each block (default: 10)
  -wsb WINDOW_SIZE_BLOCKS, --window_size_blocks WINDOW_SIZE_BLOCKS
                        Number of consecutive blocks to analyze together for slope calculation in equilibration detection (default: 10)
  -ov OVERLAP, --overlap OVERLAP
                        The amount of time overlap (in ps) between consecutive blocks.(default: 0)
  -rc RMSD_CUTOFF, --rmsd_cutoff RMSD_CUTOFF
                        RMSD cutoff in nm for clustering structures using the GROMOS method (default: 0.2 nm)

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
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
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
10: Plotting additional energy minimization metrics
11: Getting final minimized pdb structure
12: Running NVT equilibration
13: Getting NVT equilibration output
14: Running NPT equilibration
15: Getting NPT equilibration output
16: Getting convergence and clustering
17: Running Production stage
18: Getting Production output
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
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
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
11: Plotting additional energy minimization metrics
12: Getting final minimized pdb structure
13: Running NVT equilibration
14: Getting NVT equilibration output
15: Running NPT equilibration
16: Getting NPT equilibration output
17: Running Production stage
18: Getting Production output
19: Getting convergence and clustering
```

To run EquilibraTor for a protein file:

```Text
EquilibraTor -p example/example_protein.pdb
````

If you want to tweak certain parameters for your EquilibraTor workflow (e.g., NPT equilibration), you can specify to run it from 14 step to the 15 step, without having to run again from topology to minimnization:

```Text
EquilibraTor -p example/example_protein.pdb -fs 14 -ls 15
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
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
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
14: Plotting additional energy minimization metrics
15: Getting final minimized pdb structure
16: Running NVT equilibration
17: Getting NVT equilibration output
18: Running NPT equilibration
19: Getting NPT equilibration output
20: Getting convergence and clustering
21: Running Production stage
22: Getting Production output
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
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
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
15: Plotting additional energy minimization metrics
16: Getting final minimized pdb structure
17: Running NVT equilibration
18: Getting NVT equilibration output
19: Running NPT equilibration
20: Getting NPT equilibration output
21: Getting convergence and clustering
22: Running Production stage
23: Getting Production output
```

To run EquilibraTor using these protein-ligand files:

```Text
EquilibraTor -l example/example_ligand.pdb -p example/example_protein.pdb
```

## 💾 Outputs

### Multiple panel figure
<img src="paper/figures/protein-ligand_equilibration.png">

### Equilibration PDB last frame
<img src="paper/figures/protein-ligand_eq_last_frame.png">

## Convergence analysis
EquilibraTor features an automated module for determining when a MD simulation reaches equilibrium. This analysis uses block-based statistical methods to detect stability in key observables, ensuring that only equilibrated segments are used for representative structure clustering and post-simulation analysis.

## Method Overview
EquilibraTor divides the trajectory into consecutive time blocks and computes the drift (slope) of several observables, such as RMSD, potential energy, pressure, and radius of gyration. The program checks whether these observables have stabilized by comparing their computed drift against user-defined relative thresholds.

Users define the acceptance thresholds via the `-oth` argument, for example:

-oth rmsd=0.1 potential=0.2 pressure=0.2

Each threshold specifies the maximum allowed drift as a fraction of the standard deviation of the last portion of the trajectory. When the measured drift of an observable falls below its threshold, the observable is considered stable, marking the onset of equilibration.

If multiple observables are specified, EquilibraTor can assess them either individually or jointly, depending on whether the `-rao` flag (`--req_all_obs`) is used.

## Output and Report Interpretation

Following analysis, EquilibraTor automatically generates an **Equilibration Analysis Report**, which summarizes the equilibration point for each observable and identifies the equilibrated region of the trajectory.  

Example report:

In this example, RMSD reaches its stability threshold (0.1) at 50 ps, marking the start of equilibration. The equilibrated segment spans from 50 ps to >

```
EQUILIBRATION ANALYSIS REPORT
============================================================

Analysis mode: Single observable (rmsd)
Block size: 10 ps
Block overlap: 0 ps
Total simulation time: 500.0 ps

EQUILIBRATION TIMES BY OBSERVABLE:
------------------------------------------------------------
  rmsd                                         50.0 ps (threshold=0.1)

------------------------------------------------------------
FINAL EQUILIBRATION START: 50.0 ps
EQUILIBRATED SEGMENT: 50.0 - 500.0 ps
EQUILIBRATED LENGTH: 450.0 ps
------------------------------------------------------------

```

After running EquilibraTor with the protein and ligand examples, the execution steps during the execution should be printed in the terminal:


```
   ____          _ ___ __           ______
  / __/__ ___ __(_) (_) /  _______ /_  __/__  ____
 / _// _ `/ // / / / / _ \/ __/ _ `// / / _ \/ __/
/___/\_, /\_,_/_/_/_/_.__/_/  \_,_//_/  \___/_/
      /_/
Equilibrator streamlines Molecular dynamics and equilibration simulations for proteins and protein-ligand complexes in a single execution
Developers: José D. D. Cediel-Becerra and Jose Cleydson F. Silva
Co-developer: Raquel Dias
Afiliation: Microbiology & Cell Science Department, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Find our paper here: https://doi.org/10.1016/j.csbj.2025.11.034
Version:v1.0.0

[1] - 2025-10-29 12:16:08,786 - INFO - Generating topology for the protein: example_protein
[2] - 2025-10-29 12:16:10,386 - INFO - Converting example_ligand PDB to MOL2
[3] - 2025-10-29 12:16:11,767 - INFO - Generating topology for the ligand: example_ligand
[4] - 2025-10-29 12:16:13,824 - INFO - Checking wether merging topology file(s) is necessary
[5] - 2025-10-29 12:16:13,854 - INFO - Making a copy of the protein: example/example_protein.pdb
[6] - 2025-10-29 12:16:14,396 - INFO - Merging topologies
[7] - 2025-10-29 12:16:14,400 - INFO - Combining and inserting unique atomtypes into main topology
[8] - 2025-10-29 12:16:14,407 - INFO - Creating the simulation box
[9] - 2025-10-29 12:16:15,008 - INFO - Solvating the system
[10] - 2025-10-29 12:16:16,075 - INFO - Adding ions to neutralize the system
[11] - 2025-10-29 12:16:17,516 - INFO - Running energy minimization
[12] - 2025-10-29 12:16:44,929 - INFO - Plotting potential energy
[13] - 2025-10-29 12:16:45,501 - INFO - Obtaining potential, backbone, and pressure xvgs
[14] - 2025-10-29 12:16:47,081 - INFO - Plotting additional energy minimization metrics
[15] - 2025-10-29 12:16:47,371 - INFO - Getting final minimized pdb structure
[16] - 2025-10-29 12:16:47,883 - INFO - Running NVT equilibration
[17] - 2025-10-29 12:21:46,878 - INFO - Getting NVT equilibration output
[18] - 2025-10-29 12:22:02,541 - INFO - Running NPT equilibration
[19] - 2025-10-29 12:27:27,176 - INFO - Getting NPT equilibration output
[20] - 2025-10-29 12:27:47,183 - INFO - Getting convergence and clustering
[21] - 2025-10-29 12:27:47,184 - INFO - Observables thresholds requested: {'rmsd': 0.1}
[22] - 2025-10-29 12:28:28,234 - INFO - Equilibration report written to: convergence_report.txt
[23] - 2025-10-29 12:28:28,237 - INFO - Extracting trajectory segment from 50 ps to the end of the simulation
[24] - 2025-10-29 12:28:32,255 - INFO - Clustering equilibrated trajectory with a RMSD cutoff value of 0.2
[25] - 2025-10-29 12:28:36,256 - INFO - Clustered representative structures saved at: /example_protein_example_ligand/clustering_anls/rep_structures.pdb
[26] - 2025-10-29 12:28:36,261 - INFO - Extracted MODEL 1 (cluster 1, 384.0 ps) to: /example_protein_example_ligand/clustering_anls/middle_repre_cluster1.pdb
[27] - 2025-10-29 12:28:36,261 - INFO - Extracting frame trajectory 384.0 ps to be used for production md
[28] - 2025-10-29 12:28:39,773 - INFO - Running Production stage
[29] - 2025-10-29 14:33:46,239 - INFO - Getting Production output
[30] - 2025-10-29 14:37:29,812 - INFO - Execution time: 2.36 hours

```
