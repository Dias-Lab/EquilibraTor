#! bin/env/ python

import os
import re
import sys
import argparse
import subprocess
import numpy as np
import equilibrator.flat
import matplotlib.pyplot as plt

## Constants
IONS_MDP = os.path.join(os.path.dirname(equilibrator.flat.__file__),'ions.mdp')
MINIM_MDP = os.path.join(os.path.dirname(equilibrator.flat.__file__),'minim.mdp')
EQUILIBRATION_MDP = os.path.join(os.path.dirname(equilibrator.flat.__file__),'equilibration.mdp')
EQUILIBRATION_2_MDP = os.path.join(os.path.dirname(equilibrator.flat.__file__),'equilibration_2.mdp')

VERSION = 'v0.1.2'

DESCRIPTION = """
   ____          _ ___ __           ______        
  / __/__ ___ __(_) (_) /  _______ /_  __/__  ____
 / _// _ `/ // / / / / _ \/ __/ _ `// / / _ \/ __/
/___/\_, /\_,_/_/_/_/_.__/_/  \_,_//_/  \___/_/
      /_/
Equilibrator streamlines Molecular dynamics and equilibration simulations for proteins and protein-ligand complexes in a single execution
Developers: José D. D. Cediel-Becerra, Jose Cleydson F. Silva and Raquel Dias
Afiliation: Microbiology & Cell Science Deparment, University of Florida
If you find any issues, please add a new issue in our GitHub repo (https://github.com/Dias-Lab/EquilibraTor)
Version:"""+VERSION

def run_equilibrator_steps(pipeline_steps, args):
    first_idx = args.first_step - 1
    last_idx = args.last_step

    if not (0 <= first_idx < last_idx <= len(pipeline_steps)):
        raise ValueError("Invalid step range: check --first_step and --last_step")

    for i, (name, func) in enumerate(pipeline_steps[first_idx:last_idx], start=first_idx + 1):
        print(f"[Step {i}] Running: {name}")
        func()

def list_equilibrator_steps(pipeline_steps):
    print("Available steps:")
    for i, (name, _) in enumerate(pipeline_steps, 1):
        print(f"{i}: {name}")

# Define utility functions
def run_command(command, cwd=None):
    """Run a shell command."""
    try:
        result = subprocess.run(command, shell=True, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr.decode()}")
        raise

def pdb_2_mol2(ligand_file, ligand_mol2):
    print("\n" + "="*100)
    print("[INFO]  Converting PDB to MOL2 format for the ligand.")
    print("="*100)
    print(f"obabel -ipdb -omol2 {ligand_file} -h > {ligand_mol2}")
    run_command(f"obabel -ipdb -omol2 {ligand_file} -h > {ligand_mol2}")

def generate_topology_ligand(ligand_file,ligand_name, output_dir):
    """Generate ligand topology using ACPYPE."""
    print("\n" + "="*100)
    print("[INFO]  Generating topology for the ligand.")
    print("="*100)
    print(f"acpype -i {ligand_file} -l -o gmx -b {ligand_name}")
    run_command(f"acpype -i {ligand_file} -l -o gmx -b {ligand_name}", cwd=output_dir)

def generate_topology_protein(protein_file,topology_file,protein_gro,output_dir):
    """Generate protein topology using GROMACS."""
    print("\n" + "="*100)
    print("[INFO]  Generating topology for the protein.")
    print("="*100)
    print(f"gmx pdb2gmx -f {protein_file} -o {protein_gro} -water tip3p -ff amber99sb -ignh -p {topology_file}")
    run_command(f"gmx pdb2gmx -f {protein_file} -o {protein_gro} -water tip3p -ff amber99sb -ignh -p {topology_file}", cwd=output_dir)


def prepare_to_merge_topologies(topology_file, ligand_itp, ligand_top, molecule_name, output_dir, ligand_provided):
    """
    Edits topology files to prepare for merging if ligand file provided.

    Parameters:
        topology_file (str): Path to the `topol.top` file.
        ligand_itp (str): Path to the `.itp` file.
        ligand_top (str): Path to the `.top` file.
        molecule_name (str): Name of the molecule (e.g., 'baricitinib').
    """
    print("\n" + "="*100)
    print("[INFO]  Preparing to merge topologies if ligand provided.")
    print("="*100)

    with open(topology_file, "r") as top_file:
        topology_lines = top_file.readlines()
    
    #; Include chain topologies
    if ligand_provided:
        include_lines = [
            f'; Include ligand topology\n',
            f'#include "{os.path.join(output_dir, ligand_itp)}"\n',
            f'#include "{os.path.join(output_dir, ligand_top)}"\n',   
        ]
        chain_includes_idx = next(
            (i for i, line in enumerate(topology_lines) 
             if line.strip() == '#include "amber99sb.ff/forcefield.itp"'),
            -1
        )
        if chain_includes_idx == -1:
            raise ValueError("Protein chain include lines not found in topol.top.")
        
        if not any(ligand_itp in line for line in topology_lines):
            # Insert the inclusion lines right after the strings
            topology_lines = (
                topology_lines[:chain_includes_idx + 1] +
                include_lines +
                ["\n"] +
                topology_lines[chain_includes_idx + 1:]
            )

        # Add the molecule information in the [ molecules ] section
        molecules_entry = f"{molecule_name}         1\n"
        molecule_section_idx = next(
            (i for i, line in enumerate(topology_lines) if line.strip().startswith("[ molecules ]")),-1)

        if molecule_section_idx != -1 and molecules_entry not in topology_lines[molecule_section_idx:]:
            topology_lines.append(molecules_entry)

    # Remove specific lines from the list
    topology_lines = [
        line.replace('#include "topol_Protein_chain_A.itp"', '').replace('#include "topol_Protein_chain_B.itp"', '')
        for line in topology_lines
    ]
    
    with open(topology_file, "w") as top_file:
        top_file.writelines(topology_lines)

    print(f"{topology_file} successfully updated")

    if ligand_provided:
        # Modify ligand top file
        with open(ligand_top, "r") as ligand_top_file:
            ligand_top_lines = ligand_top_file.readlines()

            modified_ligand_top = []
            in_defaults = False

            for line in ligand_top_lines:
                stripped_line = line.strip()

                # Ignore lines related to POSRES_LIG
                if stripped_line.startswith("#ifdef POSRES_LIG") or stripped_line.startswith("#endif") or 'posre_' in stripped_line:
                    modified_ligand_top.append(line)  
                    
                # Detect “[ defaults ]” section and comment out
                elif stripped_line.startswith("[ defaults ]"):
                    in_defaults = True
                    modified_ligand_top.append(f"; {line}")
                elif stripped_line.startswith("[ system ]"):
                    in_defaults = True
                    modified_ligand_top.append(f"; {line}")    
                elif in_defaults and stripped_line == "":
                    in_defaults = False
                elif in_defaults or stripped_line.startswith("#include") or stripped_line.startswith("[ molecules ]"):
                    modified_ligand_top.append(f"; {line}") 
                else:
                    modified_ligand_top.append(f"; {line}")

        # write modified file
        with open(ligand_top, "w") as ligand_top_file:
            ligand_top_file.writelines(modified_ligand_top)
        
        print(f"{ligand_top} modified succesfully")

def merge_topologies(protein_gro, ligand_gro, output_gro, ligand_file):
    """Merge protein and ligand topologies."""
    print("\n" + "="*100)
    print("[INFO] Merging topologies for the protein-ligand complex.")
    print("="*100)
    with open(protein_gro, 'r') as f1, open(ligand_gro, 'r') as f2, open(output_gro, 'w') as out:
        protein_lines = f1.readlines()
        ligand_lines = f2.readlines()
        total_atoms = int(protein_lines[1]) + int(ligand_lines[1])
        out.write(protein_lines[0])
        out.write(f"{total_atoms}\n")
        out.writelines(protein_lines[2:-1])
        out.writelines(ligand_lines[2:-1])
        out.write(protein_lines[-1])
        
def make_copy_of_protein(input_gro, output_gro, ligand_file):
    """make copy of protein"""
    print("\n" + "="*100)
    print("[INFO]  Making a copy of the protein structure.")
    print("="*100)
    run_command(f"gmx editconf -f {input_gro} -o {output_gro}")
    
def create_simulation_box(input_gro, output_gro):
    """Create a simulation box."""
    print("\n" + "="*100)
    print("[INFO] Creating the simulation box.")
    print("="*100)
    run_command(f"gmx editconf -f {input_gro} -o {output_gro} -c -d 1.2 -bt cubic")

def solvate_system(input_gro, output_gro, topology_file):
    """Add water to the system."""
    print("\n" + "="*100)
    print("[INFO] Solvating the system.")
    print("="*100)
    run_command(f"gmx solvate -cp {input_gro} -cs spc216.gro -o {output_gro} -p {topology_file}")

def modify_topology(atomtypes_file, topology_file):
    """
    Modifies the topology files to ensure that atom types are correctly defined.

    Parameters:
        atomtypes_file (str): Path to the file containing the [atomtypes] section (e.g., baricitinib_GMX.itp).
        topology_file (str): Path to the main topology file (e.g., topol.top).
    """
    # Read the atomtypes file
    with open(atomtypes_file, "r") as at_file:
        lines = at_file.readlines()

    # Extract the [atomtypes] section
    atomtypes_section = []
    in_atomtypes = False
    for line in lines:
        if line.strip().startswith("[ atomtypes ]"):
            in_atomtypes = True
        elif line.strip().startswith("[") and in_atomtypes:
            break  # Exit the section when another block definition is found
        if in_atomtypes:
            atomtypes_section.append(line)

    # Comment out the [atomtypes] section in the original file, preserving blank lines
    modified_lines = [
        f"; {line}" if line.strip() and line in atomtypes_section else line
        for line in lines
    ]
    with open(atomtypes_file, "w") as at_file:
        at_file.writelines(modified_lines)

    # Insert the [atomtypes] section into the beginning of topol.top, after the forcefield include
    with open(topology_file, "r") as top_file:
        topology_lines = top_file.readlines()

    forcefield_idx = next(
        (i for i, line in enumerate(topology_lines) if "forcefield.itp" in line), -1
    )
    if forcefield_idx == -1:
        raise ValueError("Could not find the forcefield.itp include in topol.top")

    # Update the topology file
    updated_topology = (
        topology_lines[:forcefield_idx + 1]
        + ["\n"] + atomtypes_section + ["\n"]
        + topology_lines[forcefield_idx + 1:]
    )
    with open(topology_file, "w") as top_file:
        top_file.writelines(updated_topology)

    print(f"Topology successfully modified: {topology_file}")

def add_ions_with_modifications(mdp_file, input_gro, output_gro, topology_file, atomtypes_file):
    """
    Modifica os arquivos de topologia e adiciona íons ao sistema para neutralizá-lo.

    Parameters:
        mdp_file (str): Caminho para o arquivo .mdp.
        input_gro (str): Arquivo .gro de entrada (e.g., sistema solventado).
        output_gro (str): Arquivo .gro de saída (e.g., sistema com íons adicionados).
        topology_file (str): Caminho para o arquivo de topologia principal.
        atomtypes_file (str): Caminho para o arquivo contendo a seção [atomtypes].
    """
    # Modify topology files
    modify_topology(atomtypes_file, topology_file)

    # Perform ion addition
    add_ions(mdp_file, input_gro, output_gro, topology_file, ions_tpr)
    print(f"Ions added successfully: {output_gro}")

def add_ions(mdp_file, input_gro, output_gro, topology_file, ions_tpr, output_dir):
    """Add ions to the system."""
    print("\n" + "="*100)
    print("[INFO] Adding ions to the system.")
    print("="*100)
    print (f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o {ions_tpr}")
    run_command(f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o {ions_tpr}", cwd=output_dir)
    print (f"echo SOL | gmx genion -s {ions_tpr} -o {output_gro} -p {topology_file} -pname NA -nname CL -neutral")
    run_command(f"echo SOL | gmx genion -s {ions_tpr} -o {output_gro} -p {topology_file} -pname NA -nname CL -neutral")

def minimize_energy(mdp_file, input_gro, output_gro, topology_file,em_tpr,em_edr,potential_xvg, output_dir):
    """Perform energy minimization."""
    print("\n" + "="*100)
    print("[INFO] Performing energy minimization.")
    print("="*100)
    run_command(f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o {em_tpr}", cwd=output_dir)
    run_command(f"gmx mdrun -v -deffnm {em_tpr.replace('.tpr','')}")
    run_command(f"echo '13 0' | gmx energy -f {em_edr} -o {potential_xvg}")

def plot_energy_results(xvg_file, output_pdf):
    """Generate plots from .xvg files."""
    print("\n" + "="*100)
    print("[INFO] Generating energy minimization plots.")
    print("="*100)
    data = np.loadtxt(xvg_file, comments=['@', '#'])
    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel('Time (ps)')
    plt.ylabel('Potential Energy (kJ/mol)')
    plt.title('Potential Energy vs Time')
    plt.savefig(output_pdf)

def get_potential_backbone_pressure_xvgs(em_edr, em_tpr, potential_xvg, rmsf_xvg, pressure_xvg, em_trr):
    print("\n" + "="*100)
    print("[INFO] Running GROMACS to calculate potential energy, RMSF for the backbone, and pressure")
    print("="*100)
    run_command(f"echo 'Potential' | gmx energy -f {em_edr} -o {potential_xvg}")
    run_command(f"echo 'Backbone' | gmx rmsf -s {em_tpr} -f {em_trr} -o {rmsf_xvg}")
    run_command(f"echo 'Pressure' | gmx energy -f {em_edr} -o {pressure_xvg}")

def plot_em_results(potential_xvg,pressure_xvg,rmsf_xvg,energy_minimization_results):
    # Load data from .xvg files
    print("\n" + "="*100)
    print("[INFO] Plotting energy minimization results for publication.")
    print("="*100)
    potential = np.loadtxt(potential_xvg, comments=['#', '@'])
    rmsf = np.loadtxt(rmsf_xvg, comments=['#', '@'])
    pressure = np.loadtxt(pressure_xvg, comments=['#', '@'])

    # Create plot panel
    fig, axs = plt.subplots(1, 3, figsize=(20, 6), sharex=False)

    # Potential Energy
    axs[0].plot(potential[:, 0], potential[:, 1])
    axs[0].set_ylabel('Potential Energy\n(kJ/mol)')
    axs[0].set_xlabel('Time (ps)')

    # Pressure
    axs[1].plot(pressure[:, 0], pressure[:, 1])
    axs[1].set_ylabel('Pressure (bar)')
    axs[1].set_xlabel('Time (ps)')

    # RMSF                                                                                                                                                                                             
    axs[2].plot(rmsf[:, 0], rmsf[:, 1])
    axs[2].set_ylabel('RMSF (nm)')
    axs[2].set_xlabel('Atom')

    plt.tight_layout()
    plt.savefig(energy_minimization_results, format='pdf', dpi=300)
    
def get_final_minimized_structure(em_tpr, em_trr, final_minimized):
    print("\n" + "="*100)
    print("[INFO]  Generating final minimized structure.")
    print("="*100)
    run_command(f"echo 'non-Water' | gmx trjconv -s {em_tpr} -f {em_trr} -o {final_minimized} -pbc nojump")

def load_xvg(filename):
    """Load data from an XVG file, ignoring comments."""
    return np.loadtxt(filename, comments=['#', '@'])

def plot_eq(eq_potential,eq_pressure_xvg,eq_temperature_xvg,eq_rmsd_xvg,eq_rmsf_xvg,eq_gyrate_xvg,equilibration_analysis):
    
    # Load data from XVG files
    potential = load_xvg(eq_potential)
    pressure = load_xvg(eq_pressure_xvg)
    temperature = load_xvg(eq_temperature_xvg)
    rmsd = load_xvg(eq_rmsd_xvg)
    rmsf = load_xvg(eq_rmsf_xvg)
    gyrate = load_xvg(eq_gyrate_xvg)

    # Create a 3x2 panel of plots
    fig, axs = plt.subplots(3, 2, figsize=(12, 15))
    fig.suptitle('Equilibration MD Analysis', fontsize=16)

    # Plot Potential Energy
    axs[0, 0].plot(potential[:, 0], potential[:, 1], label='Potential Energy', color='b')
    axs[0, 0].set_ylabel('Energy (kJ/mol)')
    axs[0, 0].set_xlabel('Time (ps)')
    axs[0, 0].legend()

    # Plot Pressure
    axs[0, 1].plot(pressure[:, 0], pressure[:, 1], label='Pressure', color='g')
    axs[0, 1].set_ylabel('Pressure (bar)')
    axs[0, 1].set_xlabel('Time (ps)')
    axs[0, 1].legend()

    # Plot Temperature
    axs[1, 0].plot(temperature[:, 0], temperature[:, 1], label='Temperature', color='r')
    axs[1, 0].set_ylabel('Temperature (K)')
    axs[1, 0].set_xlabel('Time (ps)')
    axs[1, 0].legend()

    # Plot RMSD
    axs[1, 1].plot(rmsd[:, 0], rmsd[:, 1], label='RMSD', color='c')
    axs[1, 1].set_ylabel('RMSD (nm)')
    axs[1, 1].set_xlabel('Time (ps)')
    axs[1, 1].legend()

    # Plot RMSF
    axs[2, 0].plot(rmsf[:, 0], rmsf[:, 1], label='RMSF', color='m')
    axs[2, 0].set_ylabel('RMSF (nm)')
    axs[2, 0].set_xlabel('Atom')
    axs[2, 0].legend()

    # Plot Radius of Gyration
    axs[2, 1].plot(gyrate[:, 0], gyrate[:, 1], label='Radius of Gyration', color='y')
    axs[2, 1].set_ylabel('Rg (nm)')
    axs[2, 1].set_xlabel('Time (ps)')
    axs[2, 1].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    plt.savefig(equilibration_analysis, format='pdf', dpi=300)

def get_last_frame_time(trr_file):
    """Extract the time of the last frame from the trajectory using `gmx check`."""
    result = subprocess.run(
        ["gmx", "check", "-f", trr_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        input="\n"
    )
    matches = re.findall(r"Reading frame\s+\d+\s+time\s+([\d.]+)", result.stdout)

    if matches:
        return float(matches[-1])  # Last one is the final frame
    else:
        raise RuntimeError("Could not determine last frame time from trajectory.")

def extract_last_frame_and_unwrap(equilibration_tpr, equilibration_trr,
                                   final_last_equilibrated_pdb, final_equilibrated_gro, final_equilibrated_pdb, group1="non-Water", group2="System"): 
    """Extract the last frame from a trajectory and also output a nojump .gro trajectory."""
    last_frame_time = get_last_frame_time(equilibration_trr)        
    get_last_frame = (
        f"echo '{group1}' | gmx trjconv -s {equilibration_tpr} -f {equilibration_trr} "
        f"-o {final_last_equilibrated_pdb} -dump {last_frame_time}"
    )
    run_command(get_last_frame)
    get_movie_pdb = (
        f"echo '{group1}' | gmx trjconv -s {equilibration_tpr} -f {equilibration_trr} "
        f"-o {final_equilibrated_pdb} -pbc nojump"
    )
    run_command(get_movie_pdb)
    get_movie_gro = (
        f"echo '{group2}' | gmx trjconv -s {equilibration_tpr} -f {equilibration_trr} "
	f"-o {final_equilibrated_gro} -pbc nojump"
    )
    run_command(get_movie_gro)
    
    
def make_refinement(topology_file, equilibration_tpr, em_gro, output_dir):
    print("\n" + "="*100)
    print("[INFO] Starting refinement (equilibrium) process.")
    print("="*100)

    print("\n" + "*"*100)
    print("[INFO] 1) Preparing input files for equilibration: Running `gmx grompp`...")
    print("\n" + "*"*100)
    run_command(f"gmx grompp -f {EQUILIBRATION_MDP} -c {em_gro} -p {topology_file} -o {equilibration_tpr}", cwd=output_dir)
    
    print("\n" + "*"*100)
    print("[RUNNING] 2) Equilibration simulation: Running `gmx mdrun`...")
    print("\n" + "*"*100)
    run_command(f"gmx mdrun -s {equilibration_tpr} -deffnm {equilibration_tpr.replace('.tpr','')}")

def get_refinement_output(equilibration_edr, eq_potential_xvg, eq_pressure_xvg, eq_temperature_xvg, equilibration_tpr, equilibration_trr, eq_rmsd_xvg, eq_rmsf_xvg, eq_gyrate_xvg, final_last_equilibrated_pdb, final_equilibrated_pdb, equilibration_analysis, final_equilibrated_gro):
    print("\n" + "*"*100)
    print("[INFO] 1) Extracting potential energy from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Potential' | gmx energy -f {equilibration_edr} -o {eq_potential_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 2) Extracting pressure data from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Pressure' | gmx energy -f {equilibration_edr} -o {eq_pressure_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 3) Extracting temperature data from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Temperature' | gmx energy -f {equilibration_edr} -o {eq_temperature_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 4) Calculating RMSD for the backbone from equilibration trajectory...")
    print("\n" + "*"*100)
    run_command(f"echo 'Backbone Backbone' | gmx rms -s {equilibration_tpr} -f {equilibration_trr} -o {eq_rmsd_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 5) Calculating RMSF for the backbone from equilibration trajectory...")
    print("\n" + "*"*100)
    run_command(f"echo 'Backbone' | gmx rmsf -s {equilibration_tpr} -f {equilibration_trr} -o {eq_rmsf_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 6) Calculating radius of gyration for the protein...")
    print("\n" + "*"*100)
    run_command(f"echo 'Protein' | gmx gyrate -s {equilibration_tpr} -f {equilibration_trr} -o {eq_gyrate_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 7)Extracting last frame and unwrapping trajectory...")
    print("\n" + "*"*100)
    extract_last_frame_and_unwrap(equilibration_tpr,equilibration_trr,final_last_equilibrated_pdb,final_equilibrated_gro, final_equilibrated_pdb) 
    
    print("\n" + "*"*100)
    print("[INFO] 8) Plotting equilibration analysis results...")
    print("\n" + "*"*100)
    plot_eq(eq_potential_xvg,eq_pressure_xvg,eq_temperature_xvg,eq_rmsd_xvg,eq_rmsf_xvg,eq_gyrate_xvg,equilibration_analysis)


def Run_NPT_Equilibration(topology_file, npt_tpr, nvt_gro, final_last_npt_pdb, output_dir):
    print("\n" + "="*100)
    print("[INFO] Starting additional refinement (equilibrium NPT) process.")
    print("="*100)

    print("\n" + "*"*100)
    print("[INFO] 1) Preparing input files for NPT equilibration: Running `gmx grompp`...")
    print("\n" + "*"*100)
    run_command(f"gmx grompp -f {EQUILIBRATION_2_MDP} -c {nvt_gro} -p {topology_file} -o {npt_tpr}", cwd=output_dir)

    print("\n" + "*"*100)
    print("[RUNNING] 2) NPT Equilibration simulation: Running `gmx mdrun`...")
    print("\n" + "*"*100)
    run_command(f"gmx mdrun -s {npt_tpr} -deffnm {npt_tpr.replace('.tpr','')}")

# Workflow execution
def main():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-l", "--ligand", 
        required=False,
        type=str, 
        help="Path to the ligand file."
    )
    parser.add_argument(
        "-p", "--protein", 
        required=True, 
        type=str, 
        help="Path to the protein file."
    )
    parser.add_argument(
        "-fs", "--first_step",
        required=False,
        type=int,
        default=1,
        help="Step number to start Equilibratior from (1-based)"
    )
    parser.add_argument(
        "-ls", "--last_step",
        required=False,
        type=int,
        help="Step number to end at (1-based)"
    )
    parser.add_argument(
        "-as", "--all_steps",
        required=False,
        action="store_true",
        help="List of Equilibrator steps and exit"
    )

    
    # Parse arguments
    args = parser.parse_args()

    # defining variable for input files 
    ligand_file = args.ligand
    protein_file = args.protein
    protein_name = os.path.splitext(os.path.basename(protein_file))[0]
    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0] if ligand_file else ''
    Project_dir = f"{protein_name}_{ligand_name}" if ligand_file else protein_name
            
    # Creating the directory to store outputs
    output_dir = os.path.join(os.getcwd(), Project_dir)
    os.makedirs(output_dir, exist_ok=True)
    ligand_mol2 = os.path.join(output_dir, ligand_name + '.mol2') if ligand_file else ''
    protein_gro = os.path.join(output_dir, f"{protein_name}_processed.gro")
    protein_gro_complex = protein_gro.replace('.gro','_complex.gro')
    merged_gro = os.path.join(output_dir, "merged.gro")
    protein_or_merged_gro = merged_gro if ligand_file else protein_gro        
    box_gro = os.path.join(output_dir, "box.gro")
    solvated_gro = os.path.join(output_dir, "solvated.gro")
    topology_file = os.path.join(output_dir, "topol.top")
    minimized_gro = os.path.join(output_dir, "minimized.gro")
    energy_plot = os.path.join(output_dir, f"{Project_dir}_potential.pdf")
    ions_tpr = os.path.join(output_dir, "ions.tpr")   
    
    # actype dir
    acpype_dir = os.path.join(output_dir, f"{ligand_name}.acpype") if ligand_file else ''
    ligand_itp = atomtypes_file = os.path.join(acpype_dir, f"{ligand_name}_GMX.itp") if ligand_file else ''
    ligand_top = os.path.join(acpype_dir, f"{ligand_name}_GMX.top") if ligand_file else ''
    ligand_acpype = os.path.join(acpype_dir, f"{ligand_name}_GMX.gro") if ligand_file else ''
    atomtypes_file = os.path.join(acpype_dir, "atomtypes.atp") if ligand_file else ''

    # Equilibration workflow
    equilibration_tpr = os.path.join(output_dir, "equilibration.tpr")
    equilibration_edr = os.path.join(output_dir, "equilibration.edr")
    equilibration_trr = os.path.join(output_dir, "equilibration.trr")
    eq_potential_xvg = os.path.join(output_dir, "eq_potential.xvg")
    eq_pressure_xvg = os.path.join(output_dir, "eq_pressure.xvg")
    eq_temperature_xvg = os.path.join(output_dir, "eq_temperature.xvg")
    eq_rmsd_xvg = os.path.join(output_dir, "eq_rmsd.xvg")
    eq_rmsf_xvg = os.path.join(output_dir, "eq_rmsf.xvg")
    eq_gyrate_xvg = os.path.join(output_dir, "eq_gyrate.xvg")
    final_equilibrated_pdb = os.path.join(output_dir, "final_minimized_equilibrated.pdb")
    final_last_equilibrated_pdb = os.path.join(output_dir, "final_minimized_equilibrated_last.pdb")
    energy_minimization_results = os.path.join(output_dir, f"{Project_dir}_energy_minimization_results.pdf")
    equilibration_analysis = os.path.join(output_dir, f"{Project_dir}_equilibration_analysis.pdf")
    final_equilibrated_gro = os.path.join(output_dir, "final_minimized_equilibrated.gro")

    # NPT Equilibration workflow
    npt_tpr = os.path.join(output_dir, "npt_equilibration.tpr")
    npt_edr = os.path.join(output_dir, "npt_equilibration.edr")
    npt_trr = os.path.join(output_dir, "npt_equilibration.trr")
    npt_potential_xvg = os.path.join(output_dir, "npt_potential.xvg")
    npt_pressure_xvg = os.path.join(output_dir, "npt_pressure.xvg")
    npt_temperature_xvg = os.path.join(output_dir, "npt_temperature.xvg")
    npt_rmsd_xvg = os.path.join(output_dir, "npt_rmsd.xvg")
    npt_rmsf_xvg = os.path.join(output_dir, "npt_rmsf.xvg")
    npt_gyrate_xvg = os.path.join(output_dir, "npt_gyrate.xvg")
    final_npt_pdb = os.path.join(output_dir, "final_npt_equilibrated.pdb")
    final_npt_gro = os.path.join(output_dir, "final_npt_equilibrated.gro")
    final_last_npt_pdb = os.path.join(output_dir, "final_npt_equilibrated_last.pdb")
    npt_analysis_pdf = os.path.join(output_dir, f"{Project_dir}_npt_equilibration_analysis.pdf")
    
    # Minimization workflow
    em_tpr = os.path.join(output_dir, "em.tpr")
    em_edr = os.path.join(output_dir, "em.edr")
    em_trr = os.path.join(output_dir, "em.trr")
    pressure_xvg = os.path.join(output_dir, "pressure.xvg")
    potential_xvg = os.path.join(output_dir, "potential.xvg")
    rmsf_xvg = os.path.join(output_dir, "rmsf.xvg")
    solv_ions = os.path.join(output_dir, "solv_ions.gro")
    em_gro = os.path.join(output_dir, "em.gro")
    final_minimized = os.path.join(output_dir, "final_minimized.pdb")
    
    equilibrator_steps = []
    # === Protein Topology ===
    equilibrator_steps.append(("Generate topology for protein", lambda: generate_topology_protein(protein_file, topology_file, protein_gro, output_dir)))

    # === Ligand Prep ===
    if ligand_file:
        equilibrator_steps.append(("Convert ligand PDB to MOL2", lambda: pdb_2_mol2(ligand_file, ligand_mol2)))
        equilibrator_steps.append(("Generate topology for ligand", lambda: generate_topology_ligand(ligand_mol2, ligand_name, output_dir)))

    # === Merge Prep ===
    equilibrator_steps.append(("Prepare to merge topology file(s) if ligand provided", lambda: prepare_to_merge_topologies(topology_file, ligand_itp, ligand_top, ligand_name, output_dir, ligand_file)))

    if ligand_file:
        equilibrator_steps.append(("Make a copy of protein if ligand provided", lambda: make_copy_of_protein(protein_gro, protein_gro_complex, ligand_file)))
        equilibrator_steps.append(("Merge topologies", lambda: merge_topologies(protein_gro_complex, ligand_acpype, merged_gro, ligand_file)))

    # === Simulation Setup ===
    equilibrator_steps.append(("Create the simulation box", lambda: create_simulation_box(protein_or_merged_gro, box_gro)))
    equilibrator_steps.append(("Solvate the system", lambda: solvate_system(box_gro, solvated_gro, topology_file)))
    equilibrator_steps.append(("Add ions to neutralize the system", lambda: add_ions(IONS_MDP, solvated_gro, solv_ions, topology_file, ions_tpr, output_dir)))

    # === Energy Minimization ===
    equilibrator_steps.append(("Run energy minimization", lambda: minimize_energy(MINIM_MDP, solv_ions, minimized_gro, topology_file, em_tpr, em_edr, potential_xvg, output_dir)))
    equilibrator_steps.append(("Plot potential energy", lambda: plot_energy_results(potential_xvg, energy_plot)))
    equilibrator_steps.append(("Obtain potential, backbone, and pressure xvgs", lambda: get_potential_backbone_pressure_xvgs(em_edr, em_tpr, potential_xvg, rmsf_xvg, pressure_xvg, em_trr)))
    equilibrator_steps.append(("Plot panel of additional energy minimization results", lambda: plot_em_results(potential_xvg, pressure_xvg, rmsf_xvg, energy_minimization_results)))
    equilibrator_steps.append(("Get final minimized pdb structure", lambda: get_final_minimized_structure(em_tpr, em_trr, final_minimized)))

    # === NVT Equilibration ===
    equilibrator_steps.append(("Run NVT equilibration", lambda: make_refinement(topology_file, equilibration_tpr, em_gro, output_dir
    )))
    equilibrator_steps.append(("Get NVT equilibration output", lambda: get_refinement_output(equilibration_edr, eq_potential_xvg, eq_pressure_xvg, eq_temperature_xvg, equilibration_tpr, equilibration_trr, eq_rmsd_xvg, eq_rmsf_xvg, eq_gyrate_xvg, final_last_equilibrated_pdb, final_equilibrated_pdb, equilibration_analysis, final_equilibrated_gro
    )))

    #===NPT Equilibration ===
    equilibrator_steps.append(("Run NPT equilibration", lambda: Run_NPT_Equilibration(topology_file, npt_tpr, final_equilibrated_gro, final_last_npt_pdb, output_dir
    )))
    equilibrator_steps.append(("Get NPT equilibration output", lambda: get_refinement_output(npt_edr, npt_potential_xvg, npt_pressure_xvg, npt_temperature_xvg, npt_tpr, npt_trr, npt_rmsd_xvg, npt_rmsf_xvg, npt_gyrate_xvg, final_last_npt_pdb, final_npt_pdb, npt_analysis_pdf, final_npt_gro
    )))

    if args.last_step is None:
        args.last_step = len(equilibrator_steps)
    if args.all_steps:
        list_equilibrator_steps(equilibrator_steps)
        return
    run_equilibrator_steps(equilibrator_steps, args)
    
if __name__ == "__main__":
    main()
