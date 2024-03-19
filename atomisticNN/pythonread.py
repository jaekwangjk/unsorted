# (input) Open the file
with open("diamond_Vacancy.data", "r") as file:
    # Read all lines
    lines = file.readlines()

# Extract lattice parameters
lattice_line = lines[1].split("Lattice=")[1].split("Properties")[0].strip()
lattice_params = list(map(float, lattice_line.strip('"').split()))

# Extract PBC information
pbc_line = lines[1].split("pbc=")[1].strip()
pbc_values = pbc_line.strip('"').split()
pbc_vector = tuple(1 if value == "T" else 0 for value in pbc_values)

print("pbc vector")
print(pbc_vector)

# Extract atom positions
atom_lines = lines[2:]
atom_data = []
for line in atom_lines:
    parts = line.split()
    atom_symbol = parts[0]
    atom_position = list(map(float, parts[1:]))
    atom_data.append((atom_symbol, atom_position))

# Extract atom numbers
num_atoms = int(lines[0])

# Display extracted data
print("Atom Numbers:", num_atoms)
print("Lattice Parameters:", lattice_params)
print("PBC Vector:", pbc_vector)

print("Atom Positions:")
for atom_symbol, atom_position in atom_data:
    print(atom_symbol, atom_position)

# (output) Write extracted data to the output file
with open("config.dat", "w") as outfile:
    outfile.write(f"{num_atoms}\n")
    
    #original configuration 
    outfile.write(f"{lattice_params[0]} {lattice_params[1]} {lattice_params[2]} \n")
    outfile.write(f"{lattice_params[3]} {lattice_params[4]} {lattice_params[5]} \n")
    outfile.write(f"{lattice_params[6]} {lattice_params[7]} {lattice_params[8]} \n")
    
    #deformed configuration 
    outfile.write(f"{lattice_params[0]} {lattice_params[1]} {lattice_params[2]} \n")
    outfile.write(f"{lattice_params[3]} {lattice_params[4]} {lattice_params[5]} \n")
    outfile.write(f"{lattice_params[6]} {lattice_params[7]} {lattice_params[8]} \n")
    
    
    outfile.write(f"{pbc_vector[0]} {pbc_vector[1]} {pbc_vector[2]} \n")

    for atom_symbol, atom_position in atom_data:
        outfile.write(f"{atom_symbol} {' '.join(str(coord) for coord in atom_position)} ")
        outfile.write(f"0.0 0.0 0.0 \n")

print("Data has been written to 'out_for_mdstress.dat'")
    