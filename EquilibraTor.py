import os
import sys
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Define utility functions
def run_command(command, cwd=None):
    """Run a shell command."""
    try:
        result = subprocess.run(command, shell=True, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr.decode()}")
        raise

def generate_topology_ligand(ligand_file,output_dir):
    """Generate ligand topology using ACPYPE."""
    print(f"acpype -i {ligand_file} -l -o gmx -b {output_dir}")
    run_command(f"acpype -i {ligand_file} -l -o gmx -b {output_dir}")

def generate_topology_protein(protein_file,topology_file,protein_gro):
    """Generate protein topology using GROMACS."""
    run_command(f"gmx pdb2gmx -f {protein_file} -o {protein_gro} -water tip3p -ff amber99sb -ignh -p {topology_file}")

def prepare_to_merge_topologies(topology_file, ligand_itp, ligand_top, molecule_name,output_dir):
    print (ligand_itp)
    print (ligand_top)
    """
    Edita os arquivos de topologia para preparar a fusão.

    Parameters:
        topology_file (str): Caminho para o arquivo `topol.top`.
        ligand_itp (str): Caminho para o arquivo `baricitinib_GMX.itp`.
        ligand_top (str): Caminho para o arquivo `baricitinib_GMX.top`.
        molecule_name (str): Nome da molécula (e.g., 'baricitinib').
    """
    #output_dir = output_dir.replace(f'{molecule_name}','')
    # Adicionar as linhas de inclusão ao `topol.top`
    with open(topology_file, "r") as top_file:
        topology_lines = top_file.readlines()
    
    #; Include chain topologies
    include_lines = [
        f'; Include ligand topology\n',
        f'#include "{os.path.join(output_dir, "topol_Protein_chain_A.itp")}"\n'
        f'#include "{os.path.join(output_dir, "topol_Protein_chain_B.itp")}"\n',
        f'#include "{os.path.join(output_dir, ligand_itp)}"\n',
        f'#include "{os.path.join(output_dir, ligand_top)}"\n',

    ]

    # Encontrar onde inserir as linhas de inclusão
    chain_includes_idx = next(
        (i for i, line in enumerate(topology_lines) 
         if line.strip() == '#include "amber99sb.ff/forcefield.itp"'),
        -1
    )

    if chain_includes_idx == -1:
        raise ValueError("Linhas de cadeia de proteínas não encontradas em topol.top.")

    if not any(ligand_itp in line for line in topology_lines):
        # Inserir as linhas de inclusão logo após as cadeias
        topology_lines = (
            topology_lines[:chain_includes_idx + 1] +
            include_lines +
            ["\n"] +
            topology_lines[chain_includes_idx + 1:]
        )

    # Adicionar a informação da molécula na seção [ molecules ]
    molecules_entry = f"{molecule_name}         1\n"
    molecule_section_idx = next(
        (i for i, line in enumerate(topology_lines) if line.strip().startswith("[ molecules ]")),-1)

    if molecule_section_idx != -1 and molecules_entry not in topology_lines[molecule_section_idx:]:
        topology_lines.append(molecules_entry)

    # Remover linhas específicas da lista
    topology_lines = [
        line.replace('#include "topol_Protein_chain_A.itp"', '').replace('#include "topol_Protein_chain_B.itp"', '')
        for line in topology_lines
    ]
    
    with open(topology_file, "w") as top_file:
        top_file.writelines(topology_lines)

    print(f"Arquivo {topology_file} atualizado com sucesso.")

    # Modificar `baricitinib_GMX.top`
    with open(ligand_top, "r") as ligand_top_file:
        ligand_top_lines = ligand_top_file.readlines()

    modified_ligand_top = []
    in_defaults = False

    for line in ligand_top_lines:
        stripped_line = line.strip()

        # Ignorar linhas relacionadas a POSRES_LIG
        if stripped_line.startswith("#ifdef POSRES_LIG") or stripped_line.startswith("#endif") or 'posre_' in stripped_line:
            modified_ligand_top.append(line)  # Não modificar essas linhas

        # Detectar seção "[ defaults ]" e comentar
        elif stripped_line.startswith("[ defaults ]"):
            in_defaults = True
            modified_ligand_top.append(f"; {line}")  # Comentar a linha
        # Detectar seção "[ system ]" e comentar
        elif stripped_line.startswith("[ system ]"):
            in_defaults = True
            modified_ligand_top.append(f"; {line}")  # Comentar a linha        
        # Detectar fim da seção de defaults
        elif in_defaults and stripped_line == "":
            in_defaults = False

        # Comentar outras linhas relevantes
        elif in_defaults or stripped_line.startswith("#include") or stripped_line.startswith("[ molecules ]"):
            modified_ligand_top.append(f"; {line}")  # Comentar a linha

        # Manter as demais linhas inalteradas
        else:
            modified_ligand_top.append(f"; {line}")

    # Escrever o arquivo modificado
    with open(ligand_top, "w") as ligand_top_file:
        ligand_top_file.writelines(modified_ligand_top)

    print(f"Arquivo {ligand_top} atualizado com sucesso.")


def merge_topologies(protein_gro, ligand_gro, output_gro):
    """Merge protein and ligand topologies."""
    with open(protein_gro, 'r') as f1, open(ligand_gro, 'r') as f2, open(output_gro, 'w') as out:
        protein_lines = f1.readlines()
        ligand_lines = f2.readlines()
        total_atoms = int(protein_lines[1]) + int(ligand_lines[1])
        out.write(protein_lines[0])
        out.write(f"{total_atoms}\n")
        out.writelines(protein_lines[2:-1])
        out.writelines(ligand_lines[2:-1])
        out.write(protein_lines[-1])

def make_copy_of_protein(input_gro, output_gro):
    """make copy of protein"""
    run_command(f"gmx editconf -f {input_gro} -o {output_gro}")
    
def create_simulation_box(input_gro, output_gro):
    """Create a simulation box."""
    run_command(f"gmx editconf -f {input_gro} -o {output_gro} -c -d 1.2 -bt cubic")

def solvate_system(input_gro, output_gro, topology_file):
    """Add water to the system."""
    run_command(f"gmx solvate -cp {input_gro} -cs spc216.gro -o {output_gro} -p {topology_file}")

def modify_topology(atomtypes_file, topology_file):
    """
    Modifica os arquivos de topologia para garantir que os atomtypes sejam definidos corretamente.

    Parameters:
        atomtypes_file (str): Caminho para o arquivo contendo a seção [atomtypes] (e.g., baricitinib_GMX.itp).
        topology_file (str): Caminho para o arquivo de topologia principal (e.g., topol.top).
    """
    # Lê o arquivo de atomtypes
    with open(atomtypes_file, "r") as at_file:
        lines = at_file.readlines()

    # Extrair a seção [atomtypes]
    atomtypes_section = []
    in_atomtypes = False
    for line in lines:
        if line.strip().startswith("[ atomtypes ]"):
            in_atomtypes = True
        elif line.strip().startswith("[") and in_atomtypes:
            break  # Sai da seção ao encontrar outra definição de bloco
        if in_atomtypes:
            atomtypes_section.append(line)

    # Remover ou comentar a seção [atomtypes] do arquivo original, preservando linhas em branco
    modified_lines = [
        f"; {line}" if line.strip() and line in atomtypes_section else line
        for line in lines
    ]
    with open(atomtypes_file, "w") as at_file:
        at_file.writelines(modified_lines)

    # Adicionar a seção [atomtypes] ao início de topol.top, após o include do forcefield
    with open(topology_file, "r") as top_file:
        topology_lines = top_file.readlines()

    forcefield_idx = next(
        (i for i, line in enumerate(topology_lines) if "forcefield.itp" in line), -1
    )
    if forcefield_idx == -1:
        raise ValueError("Não foi possível encontrar o include de forcefield.itp em topol.top")

    # Atualiza o arquivo de topologia
    updated_topology = (
        topology_lines[:forcefield_idx + 1]
        + ["\n"] + atomtypes_section + ["\n"]
        + topology_lines[forcefield_idx + 1:]
    )
    with open(topology_file, "w") as top_file:
        top_file.writelines(updated_topology)

    print(f"Topologia modificada com sucesso: {topology_file}")


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
    # Modificar os arquivos de topologia
    modify_topology(atomtypes_file, topology_file)

    # Executar a adição de íons (substitua pelo comando GROMACS real)
    add_ions(mdp_file, input_gro, output_gro, topology_file)
    print(f"Íons adicionados com sucesso: {output_gro}")

def add_ions(mdp_file, input_gro, output_gro, topology_file):
    """Add ions to the system."""
    print ("""Add ions to the system.""")
    print (f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o ions.tpr")
    run_command(f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o ions.tpr")
    print (f"echo SOL | gmx genion -s ions.tpr -o {output_gro} -p {topology_file} -pname NA -nname CL -neutral")
    run_command(f"echo SOL | gmx genion -s ions.tpr -o {output_gro} -p {topology_file} -pname NA -nname CL -neutral")

def minimize_energy(mdp_file, input_gro, output_gro, topology_file,em_tpr,em_edr,potential_xvg):
    """Perform energy minimization."""
    run_command(f"gmx grompp -f {mdp_file} -c {input_gro} -p {topology_file} -o {em_tpr}")
    run_command(f"gmx mdrun -v -deffnm {em_tpr.replace('.tpr','')}")
    run_command(f"echo '13 0' | gmx energy -f {em_edr} -o {potential_xvg}")

def plot_energy_results(xvg_file, output_pdf):
    """Generate plots from .xvg files."""
    data = np.loadtxt(xvg_file, comments=['@', '#'])
    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel('Time (ps)')
    plt.ylabel('Potential Energy (kJ/mol)')
    plt.title('Potential Energy vs Time')
    plt.savefig(output_pdf)

def plot_em_results(potential_xvg,pressure_xvg,rmsf_xvg,energy_minimization_results):

    # Load data from .xvg files
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
    axs[2, 0].set_xlabel('Residue')
    axs[2, 0].legend()

    # Plot Radius of Gyration
    axs[2, 1].plot(gyrate[:, 0], gyrate[:, 1], label='Radius of Gyration', color='y')
    axs[2, 1].set_ylabel('Rg (nm)')
    axs[2, 1].set_xlabel('Time (ps)')
    axs[2, 1].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to fit the title
    plt.savefig(equilibration_analysis, format='pdf', dpi=300)
    
    
def Make_Refinement(topology_file,equilibration_tpr,equilibration_edr,eq_potential_xvg,equilibration_trr,eq_pressure_xvg,eq_temperature_xvg,eq_rmsd_xvg,eq_rmsf_xvg,eq_gyrate_xvg,final_equilibrated_pdb,equilibration_analysis,em_gro,final_last_equilibrated_pdb):
    print("\n" + "*"*100)
    print("[INFO] 1) Preparing input files for equilibration: Running `gmx grompp`...")
    print("\n" + "*"*100)
    run_command(f"gmx grompp -f equilibration.mdp -c {em_gro} -p {topology_file} -o {equilibration_tpr}")
    
    print("\n" + "*"*100)
    print("[RUNNING] 2) Equilibration simulation: Running `gmx mdrun`...")
    print("\n" + "*"*100)
    run_command(f"gmx mdrun -s {equilibration_tpr} -deffnm {equilibration_tpr.replace('.tpr','')}")
    
    print("\n" + "*"*100)
    print("[INFO] 3) Extracting potential energy from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Potential' | gmx energy -f {equilibration_edr} -o {eq_potential_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 4) Extracting pressure data from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Pressure' | gmx energy -f {equilibration_edr} -o {eq_pressure_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 5) Extracting temperature data from equilibration results...")
    print("\n" + "*"*100)
    run_command(f"echo 'Temperature' | gmx energy -f {equilibration_edr} -o {eq_temperature_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 6) Calculating RMSD for the backbone from equilibration trajectory...")
    print("\n" + "*"*100)
    run_command(f"echo 'Backbone Backbone' | gmx rms -s {equilibration_tpr} -f {equilibration_trr} -o {eq_rmsd_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 7) Calculating RMSF for the backbone from equilibration trajectory...")
    print("\n" + "*"*100)
    run_command(f"echo 'Backbone' | gmx rmsf -s {equilibration_tpr} -f {equilibration_trr} -o {eq_rmsf_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 8) Calculating radius of gyration for the protein...")
    print("\n" + "*"*100)
    run_command(f"echo 'Protein' | gmx gyrate -s {equilibration_tpr} -f {equilibration_trr} -o {eq_gyrate_xvg}")
    
    print("\n" + "*"*100)
    print("[INFO] 9) Generating final equilibrated structure in PDB format (nojump)...")
    print("\n" + "*"*100)
    run_command(f"echo '17' | gmx trjconv -s {equilibration_tpr} -f {equilibration_trr} -o {final_equilibrated_pdb} -pbc nojump")
    
    print("\n" + "*"*100)
    print("[INFO] 10) Extracting specific frame (1000 ps) from equilibrated trajectory...")
    print("\n" + "*"*100)
    run_command(f"echo '17' | gmx trjconv -s {equilibration_tpr} -f {equilibration_trr} -o {final_last_equilibrated_pdb} -dump 1000")
    
    print("\n" + "*"*100)
    print("[INFO] 11) Plotting equilibration analysis results...")
    print("\n" + "*"*100)
    plot_eq(eq_potential_xvg,eq_pressure_xvg,eq_temperature_xvg,eq_rmsd_xvg,eq_rmsf_xvg,eq_gyrate_xvg,equilibration_analysis)

def PrepareMutantes(input_file):
    output_file = "individual_list.txt"  # Arquivo de saída

    # Abrir os arquivos
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Ignorar a primeira linha (equivalente a `tail -n +2`)
        lines = infile.readlines()[1:]

        for line in lines:
            # Substituir '>' por '\t' (equivalente a `tr '>' '\t'`)
            line = line.replace('>', '\t')

            # Dividir a linha em campos (equivalente a `cut -f 1,2`)
            fields = line.strip().split('\t')[:2]
            if len(fields) < 2:
                continue  # Pular linhas mal formatadas

            id_part, sequence = fields

            # Gerar os outputs A e B
            output_a = f"{sequence}A{id_part},"
            output_b = f"{sequence}B{id_part};"

            # Escrever no arquivo de saída
            outfile.write(output_a + output_b + "\n")

    return output_file
# Workflow execution
def main():

    # Define file paths
    # Configuração do argparse
    parser = argparse.ArgumentParser(
        description="Process files for ligand, protein, and mutants_uniprot."
    )
    parser.add_argument(
        "-l", "--ligand", 
        required=True, 
        type=str, 
        help="Path to the ligand file."
    )
    parser.add_argument(
        "-p", "--protein", 
        required=True, 
        type=str, 
        help="Path to the protein file."
    )

    
    # Parse os argumentos
    args = parser.parse_args()

    # Atribuição de arquivos
    ligand_file = args.ligand
    protein_file = args.protein
    
    molecule_name = ligand_file.replace('.pdb','')
    Project_dir = protein_file.replace('.pdb','') + '_' + ligand_file.replace('.pdb','')
    
    # Criação do diretório para armazenar os dados
    output_dir = os.path.join(os.getcwd(), Project_dir)
    os.makedirs(output_dir, exist_ok=True)   
    
    ligand_mol2 = os.path.join(output_dir, ligand_file.replace('.pdb','.mol2'))
    
    protein_gro = os.path.join(output_dir, f"{protein_file}_processed.gro")
    protein_gro_complex = protein_gro.replace('.gro','_complex.gro')
    merged_gro = os.path.join(output_dir, "merged.gro")
    box_gro = os.path.join(output_dir, "box.gro")
    solvated_gro = os.path.join(output_dir, "solvated.gro")
    topology_file = os.path.join(output_dir, "topol.top")
    minimized_gro = os.path.join(output_dir, "minimized.gro")
    energy_plot = os.path.join(output_dir, f"{Project_dir}_potential.pdf")
    
    # actype dir
    atomtypes_file = os.path.join(os.getcwd(), f"{molecule_name}.acpype/{molecule_name}_GMX.itp")
    ligand_itp = os.path.join(os.getcwd(), f"{molecule_name}.acpype/{molecule_name}_GMX.itp")
    ligand_top = os.path.join(os.getcwd(), f"{molecule_name}.acpype/{molecule_name}_GMX.top")
    ligand_acpype = os.path.join(os.getcwd(), f"{molecule_name}.acpype/{molecule_name}_GMX.gro")
        
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
    
    print("\n" + "="*100)
    print("[INFO] Setting paths for energy minimization output files.")
    print("="*100)
    
    # Convert PDB to MOL2
    print("\n" + "="*100)
    print("[INFO]  Converting PDB to MOL2 format for the ligand.")
    print("="*100)
    print(f"obabel -ipdb -omol2 {ligand_file} -h > {ligand_mol2}")
    run_command(f"obabel -ipdb -omol2 {ligand_file} -h > {ligand_mol2}")

    # Generate topologies
    print("\n" + "="*100)
    print("[INFO]  Generating topology for the ligand.")
    print("="*100)
    #generate_topology_ligand(ligand_mol2,molecule_name)
    generate_topology_ligand(ligand_file,molecule_name)
    
    print("\n" + "="*100)
    print("[INFO]  Generating topology for the protein.")
    print("="*100)
    generate_topology_protein(protein_file, topology_file,protein_gro)

    print("\n" + "="*100)
    print("[INFO]  Preparing to merge topologies.")
    print("="*100)
    prepare_to_merge_topologies(topology_file,ligand_itp,ligand_top,molecule_name,output_dir)
    
    print("\n" + "="*100)
    print("[INFO]  Making a copy of the protein structure.")
    print("="*100)
    # make_copy_of_protein
    make_copy_of_protein(protein_gro, protein_gro_complex)
    
    # Merge topologies
    print("\n" + "="*100)
    print("[INFO] Merging topologies for the protein-ligand complex.")
    print("="*100)
    merge_topologies(protein_gro_complex, ligand_acpype, merged_gro)
    
    # Create simulation box
    print("\n" + "="*100)
    print("[INFO] Creating the simulation box.")
    print("="*100)
    create_simulation_box(merged_gro, box_gro)

    # Solvate system
    print("\n" + "="*100)
    print("[INFO] Solvating the system.")
    print("="*100)
    solvate_system(box_gro, solvated_gro, topology_file)

    # Add ions (testar)
    print("\n" + "="*100)
    print("[INFO] ⚡ Adding ions to the system (with modifications).")
    print("="*100)
    add_ions_with_modifications("ions.mdp", solvated_gro, solv_ions, topology_file,atomtypes_file)

    # Minimize energy
    print("\n" + "="*100)
    print("[INFO] Performing energy minimization.")
    print("="*100)
    minimize_energy("minim.mdp", solv_ions, minimized_gro, topology_file,em_tpr,em_edr,potential_xvg)

    # Generate plot
    print("\n" + "="*100)
    print("[INFO] Generating energy minimization plots.")
    print("="*100)
    plot_energy_results(potential_xvg, energy_plot)
    
    # Generate additional plots for publication
    print("\n" + "="*100)
    print("[INFO] Running GROMACS to calculate potential energy.")
    print("="*100)
    run_command(f"echo 'Potential' | gmx energy -f {em_edr} -o {potential_xvg}")
    
    print("\n" + "="*100)
    print("[INFO] Running GROMACS to calculate RMSF for the backbone.")
    print("="*100)
    run_command(f"echo 'Backbone' | gmx rmsf -s {em_tpr} -f {em_trr} -o {rmsf_xvg}")
    
    print("\n" + "="*100)
    print("[INFO] Running GROMACS to calculate pressure.")
    print("="*100)
    run_command(f"echo 'Pressure' | gmx energy -f {em_edr} -o {pressure_xvg}")
    
    print("\n" + "="*100)
    print("[INFO] Plotting energy minimization results for publication.")
    print("="*100)
    plot_em_results(potential_xvg,pressure_xvg,rmsf_xvg,energy_minimization_results)
    
    print("\n" + "="*100)
    print("[INFO]  Generating final minimized structure.")
    print("="*100)
    run_command(f"echo 17 | gmx trjconv -s {em_tpr} -f {em_trr} -o {final_minimized} -pbc nojump")
    
    #Refinement (Equilibrium)
    print("\n" + "="*100)
    print("[INFO] Starting refinement (equilibrium) process.")
    print("="*100)
    Make_Refinement(topology_file,equilibration_tpr,equilibration_edr,eq_potential_xvg,equilibration_trr,eq_pressure_xvg,eq_temperature_xvg,eq_rmsd_xvg,eq_rmsf_xvg,eq_gyrate_xvg,final_equilibrated_pdb,equilibration_analysis,em_gro,final_last_equilibrated_pdb)
    

    
if __name__ == "__main__":
    main()

# https://jealous-jingle-3d3.notion.site/GROMACS-energy-minimization-131695fd261080caa8caf441a4031011
# https://jealous-jingle-3d3.notion.site/Using-GLM-score-to-predict-binding-affinity-and-interaction-between-protein-JAK2-and-inhibitors-14a695fd261080b08f7bc46fbe2244ee#
# https://jealous-jingle-3d3.notion.site/Mutate-protein-in-foldx-150695fd26108033bb48fa6158705b9c

