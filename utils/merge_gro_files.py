import sys

def parse_gro_line(line):
    """
    Parse a single line of a .gro file respecting its fixed-width format.
    """
    # Fixed-width fields based on .gro format
    residue_number = int(line[0:5].strip())  # Residue number (columns 1-5)
    residue_name = line[5:10].strip()        # Residue name (columns 6-10)
    atom_name = line[10:15].strip()          # Atom name (columns 11-15)
    atom_index = int(line[15:20].strip())    # Atom index (columns 16-20)
    x = float(line[20:28].strip())           # X coordinate (columns 21-28)
    y = float(line[28:36].strip())           # Y coordinate (columns 29-36)
    z = float(line[36:44].strip())           # Z coordinate (columns 37-44)
    return residue_number, residue_name, atom_name, atom_index, x, y, z

def format_gro_line(residue_number, residue_name, atom_name, atom_index, x, y, z):
    """
    Format a single line for a .gro file using fixed-width fields.
    """
    return f"{residue_number:>5}{residue_name:>5}{atom_name:>5}{atom_index:>5}{x:8.3f}{y:8.3f}{z:8.3f}\n"

def merge_gro_files(file1, file2, output_file):
    try:
        # Read the first file
        with open(file1, 'r') as f1:
            lines1 = f1.readlines()
        
        # Read the second file
        with open(file2, 'r') as f2:
            lines2 = f2.readlines()
        
        # Extract headers and atom data
        header1 = lines1[0].strip()
        atom_count1 = int(lines1[1].strip())
        atoms1 = lines1[2:-1]  # Exclude the last line (box dimensions)

        header2 = lines2[0].strip()
        atom_count2 = int(lines2[1].strip())
        atoms2 = lines2[2:-1]  # Exclude the last line (box dimensions)
        box_dimensions = lines2[-1].strip()  # Use box dimensions from the second file
        
        # Parse and reformat atoms in the first file
        parsed_atoms1 = [parse_gro_line(atom) for atom in atoms1]
        
        # Parse and reformat atoms in the second file with updated indices
        parsed_atoms2 = []
        for atom in atoms2:
            residue_number, residue_name, atom_name, atom_index, x, y, z = parse_gro_line(atom)
            new_atom_index = atom_index + atom_count1  # Update atom index
            parsed_atoms2.append((residue_number, residue_name, atom_name, new_atom_index, x, y, z))
        
        # Merge atoms and calculate new total atom count
        merged_atoms = parsed_atoms1 + parsed_atoms2
        total_atom_count = len(merged_atoms)
        
        # Write to output file
        with open(output_file, 'w') as out:
            out.write(f"{header1}\n")
            out.write(f"{total_atom_count}\n")
            for atom in merged_atoms:
                out.write(format_gro_line(*atom))
            out.write(f"{box_dimensions}\n")
        
        print(f"Files merged successfully into {output_file}")
    
    except Exception as e:
        print(f"Error: {e}")

# Main function to handle command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_gro_files.py <file1.gro> <file2.gro> <output.gro>")
    else:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        output_file = sys.argv[3]
        merge_gro_files(file1, file2, output_file)
