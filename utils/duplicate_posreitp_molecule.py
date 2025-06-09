import re
import sys


def process_position_restraints(input_file, output_file):
    # Define regex pattern for the [ position_restraints ] section
    section_pattern = r"(\[ position_restraints \].*?)(?=\n\[|\Z)"

    # Read the input file
    with open(input_file, 'r') as file:
        content = file.read()

    # Process the [ position_restraints ] section
    match = re.search(section_pattern, content, re.DOTALL)
    if match:
        original_section = match.group(1)
        lines = original_section.split("\n")
        header = lines[0]  # Keep the header line
        data_lines = lines[1:]  # Data lines to process

        # Prepare updated lines
        updated_lines = []
        for line in data_lines:
            if line.startswith(";"):
                continue
            if not line.strip():
                updated_lines.append(line)  # Keep comments or empty lines as is
                continue

            # Update only the first column (atom number) while preserving spacing
            updated_line = re.sub(
                r"^(\s*\d+)", 
                lambda m: f"{int(m.group(1).strip()) + 37:>6}", 
                line
            )
            updated_lines.append(updated_line)

        # Concatenate new data directly after original data without extra blank lines
        updated_section = original_section.rstrip() + "\n" + "\n".join(updated_lines)

        # Replace the original section with the concatenated one
        content = content.replace(original_section, updated_section)

    # Write the modified content to the output file
    with open(output_file, 'w') as file:
        file.write(content)


# Main execution
if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    process_position_restraints(input_filename, output_filename)
