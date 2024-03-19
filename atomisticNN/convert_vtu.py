import numpy as np

# Define the file paths
input_file_path = "hardy1.stress"
output_file_path = "new.vtu"

def write_vtu_data(ngrid, data, file_name):
    nx=round(ngrid**(1/3)) # Grid dimensions
    ny=nx
    nz=nx

    with open(file_name, 'w') as file:
        file.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">\n')
        file.write('<UnstructuredGrid>\n')
        file.write(f'<Piece NumberOfPoints="{nx * ny * nz}" NumberOfCells="{(nx - 1) * (ny - 1) * (nz - 1)}">\n')
        
        # Point Data
        file.write('<PointData Scalars="scalars">\n')
        file.write('<DataArray type="Float32" Name="sigma_xx" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 3]
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')
        
        file.write('<DataArray type="Float32" Name="sigma_xy" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 4]
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')
        
        
        file.write('<DataArray type="Float32" Name="sigma_xz" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 5]
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')
        
  
        file.write('<DataArray type="Float32" Name="sigma_yy" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 7]
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')
        

        file.write('<DataArray type="Float32" Name="sigma_yz" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 8] 
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')

        
        file.write('<DataArray type="Float32" Name="sigma_zz" NumberOfComponents="1" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    value = data[i*ny*nx+j*nx+k, 11]
                    file.write(f'{value} ')
                file.write('\n')
        file.write('</DataArray>\n')
        
        file.write('</PointData>\n')
        
        
        # Coordinate Points
        file.write('<Points>\n')
        file.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in range(nz):
            for j in range(ny):
                for k in range(nx):
                    file.write(f'{data[i*ny*nx+j*nx+k, 0]} {data[i*ny*nx+j*nx+k, 1]} {data[i*ny*nx+j*nx+k, 2]}\n')
        file.write('</DataArray>\n')
        file.write('</Points>\n')
        
        # Cell Connectivity
        file.write('<Cells>\n')
        file.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for i in range(nz - 1):
            for j in range(ny - 1):
                for k in range(nx - 1):
                    org = i * nx * ny + j * nx + k
                    file.write(f'{org} {org + 1} {org + 1 + nx} {org + nx} '
                               f'{org + nx * ny} {org + nx * ny + 1} {org + nx * ny + nx + 1} {org + nx * ny + nx}\n')
        file.write('</DataArray>\n')
        
        # Cell Offsets
        file.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        for offset in range(8, (nx - 1) * (ny - 1) * (nz - 1) * 8 + 1, 8):
            file.write(f'{offset}\n')
        file.write('</DataArray>\n')
        
        # Cell Types
        file.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in range((nx - 1) * (ny - 1) * (nz - 1)):
            file.write('12\n')  # Hexahedral cell type
        file.write('</DataArray>\n')
        
        file.write('</Cells>\n')
        file.write('</Piece>\n')
        file.write('</UnstructuredGrid>\n')
        file.write('</VTKFile>\n')


# main driver 

# Read the data from the input file
with open(input_file_path, "r") as file:
    # Read the first line to get the number of grid points
    ngrid= int(file.readline().strip())
    # Initialize an empty list to store the data
    data = []
    # Read the rest of the lines and parse the data
    for line in file:
        # Split the line into tokens
        tokens = line.split()

        # Check if the line contains enough values
        if len(tokens) >= 3:
            # Extract x, y, z points
            x, y, z = map(float, tokens[:3])

            # Convert the rest of the tokens to floats
            values = list(map(float, tokens[3:]))

            # Append the data to the list
            data.append([x, y, z] + values)

# Convert the data to a numpy array
data = np.array(data)

# Define sorting keys (x, y, z)
sorting_keys = (data[:, 0], data[:, 1], data[:, 2])

# Perform hierarchical sorting
sorted_indices = np.lexsort(sorting_keys)

# Apply the sorted indices to reorder the original data
sorted_data = data[sorted_indices]

write_vtu_data(ngrid, sorted_data, output_file_path)

