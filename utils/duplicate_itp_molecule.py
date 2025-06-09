import re
import sys


def process_file(input_file, output_file):
    # Define regex patterns for each section
    section_patterns = {
        "atoms": r"(\[ atoms \].*?)(?=\n\[|\Z)",
        "bonds": r"(\[ bonds \].*?)(?=\n\[|\Z)",
        "pairs": r"(\[ pairs \].*?)(?=\n\[|\Z)",
        "angles": r"(\[ angles \].*?)(?=\n\[|\Z)",
        "dihedrals": r"(\[ dihedrals \].*?)(?=\n\[|\Z)"
    }

    # Read the input file
    with open(input_file, 'r') as file:
        content = file.read()

    # Initialize variables to store the updated sections
    updated_content = content  # Start with the original content

    # Process each section
    for section, pattern in section_patterns.items():
        match = re.search(pattern, content, re.DOTALL)
        if match:
            original_section = match.group(1)
            lines = original_section.split("\n")
            header = lines[0]  # Keep the header line
            data_lines = lines[1:]  # Data lines to process

            # Prepare updated lines
            updated_lines = []
            for line in data_lines:
                if line.startswith(";"): #skip comments
                    continue
                if not line.strip():
                    updated_lines.append(line)  # Keep empty lines as is
                    continue

                if section == "atoms":
                    # Update the first column (nr) while preserving spaces
                    updated_line = re.sub(
                        r"^(\s*\d+)", 
                        lambda m: f"{int(m.group(1).strip()) + 37:>6}", 
                        line
                    )
                elif section in {"bonds", "pairs"}:
                    # Update ai and aj columns (first two numbers)
                    updated_line = re.sub(
                        r"^(\s*\d+)(\s+\d+)", 
                        lambda m: f"{int(m.group(1).strip()) + 37:>6}{int(m.group(2).strip()) + 37:>7}", 
                        line
                    )
                elif section == "angles":
                    # Update ai, aij, and ak columns (first three numbers)
                    updated_line = re.sub(
                        r"^(\s*\d+)(\s+\d+)(\s+\d+)", 
                        lambda m: f"{int(m.group(1).strip()) + 37:>6}{int(m.group(2).strip()) + 37:>7}{int(m.group(3).strip()) + 37:>7}", 
                        line
                    )
                elif section == "dihedrals":
                    # Update i, j, k, l columns (first four numbers)
                    updated_line = re.sub(
                        r"^(\s*\d+)(\s+\d+)(\s+\d+)(\s+\d+)", 
                        lambda m: f"{int(m.group(1).strip()) + 37:>6}{int(m.group(2).strip()) + 37:>7}{int(m.group(3).strip()) + 37:>7}{int(m.group(4).strip()) + 37:>7}", 
                        line
                    )
                else:
                    updated_line = line

                updated_lines.append(updated_line)

            # Concatenate the new data directly after the original data without extra empty lines
            updated_section = original_section.rstrip() + "\n" + "\n".join(updated_lines)

            # Replace the original section with the concatenated one
            updated_content = updated_content.replace(original_section, updated_section)

    # Write the modified content to the output file
    with open(output_file, 'w') as file:
        file.write(updated_content)


# Main execution
if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    process_file(input_filename, output_filename)

