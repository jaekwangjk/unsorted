import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
from scipy.interpolate import griddata


def read_file(filename):

    with open(filename, 'r') as file:
        lines = file.readlines()
    
    total_rows = int(lines[0].strip())  # Read the total number of rows
    data = []
    
    for line in lines[2:]:  # Skip the first line
        values = [float(v) for v in line.split()]
        data.append(values)
    
    return np.array(data), total_rows
    

def read_file_LDAD(gridfilename, stressfilename):
    # Read first file
    with open(gridfilename, 'r') as file:
        lines1 = file.readlines()
    
    # Read second file
    with open(stressfilename, 'r') as file:
        lines2 = file.readlines()
    
    # Assume total_rows from the first file
    total_rows = int(lines1[0].strip())
    
    data = []
    
    # Iterate over both files in parallel, starting from line 2
    for line1, line2 in zip(lines1[2:], lines2[2:]):
        values1 = [float(v) for v in line1.split()]
        values2 = [float(v) for v in line2.split()]
        combined_values = values1 + values2
        data.append(combined_values)
    
    return np.array(data), total_rows
    
    

def plot_stress_contour(outputfilename, data, lc, ifD=0, avgR=0,component=3):
    radius = 15 * lc 
    x = (data[:, 0] - 50 * lc) / radius
    y = (data[:, 1] - 50 * lc) / radius
    stress = data[:, component] * 160.2 / 0.1

    xi = np.linspace(np.min(x), np.max(x), 200)
    yi = np.linspace(np.min(y), np.max(y), 200)
    xi, yi = np.meshgrid(xi, yi)

    valid_mask = ~np.isnan(stress)
    zi = griddata((x[valid_mask], y[valid_mask]), stress[valid_mask], (xi, yi), method='cubic')

    # Apply mask after interpolation
    distance_squared = xi**2 + yi**2
    radius_squared = (radius + 1 * ifD + avgR)**2 / (radius**2)
    zi[distance_squared < radius_squared] = np.nan

    print(f"zi actual min = {np.nanmin(zi):.3f}, max = {np.nanmax(zi):.3f}")

    cmap = plt.colormaps.get_cmap('turbo').copy()
    
    cmap.set_bad('white')

    fig, ax = plt.subplots(figsize=(8, 6))
    
    if(component==3):
        cvmin = -1.2
        cvmax = 3.5
    else:
        cvmin = -1.0
        cvmax = 1.0
    
    
    contour = ax.contourf(xi, yi, zi, levels=200, cmap=cmap ,  vmin=cvmin, vmax=cvmax)
    
    #contour = plt.contourf(xi, yi, zi, levels=200, cmap=cmap)
    
    
    # force vmax, vmin of colorbar
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array(zi)
    mappable.set_clim(cvmin, cvmax)

    cbar= plt.colorbar(mappable, ax=ax)
    
    # if not scale colorbar 
    #cbar= plt.colorbar()
    
    if(component ==3):
        cbar.set_label(r'$\sigma_{xx} / \sigma_{\infty} $', fontsize=20)
    elif(component ==4):
        cbar.set_label(r'$\sigma_{xy} / \sigma_{\infty} $', fontsize=20)
    elif(component ==7):
        cbar.set_label(r'$\sigma_{yy} / \sigma_{\infty} $', fontsize=20)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(r'$x/R$', fontsize=20)
    plt.ylabel(r'$y/R$', fontsize=20)
    plt.tick_params(axis='both', labelsize=14)

    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    plt.savefig(outputfilename, dpi=300, bbox_inches='tight')
    plt.show()
    



models = [
#{"model_name": "SW_StillingerWeber_1985_Si__MO_405512056662_006", "lattice_constant": 5.430949784815312, "influent_distance": 7.52336},
{"model_name": "NeuralNetwork_d1__MO_000000111111_000", "lattice_constant": 5.465894, "influent_distance": 5.5},
]


for model in models:
    
    modelname = model["model_name"]
    lc = model["lattice_constant"]
    influenceD = model["influent_distance"]
    
    print(f"모델명: {modelname}")
    print(f"Lattice Constant: {lc}")
    print(f"Influent Distance: {influenceD}")
    print("-----------------------------")
    
    component = 3
    
    file1= "dEdr_NNd1_ldad_constant_stress_ref.pushedGrids"
    file2= "dEdr_NNd1_ldad_constant_stress_ref.CauchyPushed"
    
    data2, total_rows2 = read_file_LDAD(file1,file2)
    #print(data2)
    
    if(component ==3):
        compName = "xx"
    elif(component ==4):
        compName = "xy"
    elif(component ==7):
        compName = "yy"
    
    outputfilename = modelname + "_LDAD_CaucyPushed_dEdr_sig_" + compName +".png"
    avR=0
    plot_stress_contour(outputfilename, data2,lc,influenceD,avR,component)
    
