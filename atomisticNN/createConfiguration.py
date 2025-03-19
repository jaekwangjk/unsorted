import subprocess 
import sys
import math
import numpy as np 

import configLibrary

# Main execution
if __name__ == "__main__":

    print("******* Program Executed....")
    
    # Assuming that you have a input files 
    box_size, lattice_constant, D,C11,C12,C44,sigma11,modelName = configLibrary.read_input_from_stdin()
    
    positions, atom_ids = configLibrary.fill_box_with_fcc_diamond_lattice(lattice_constant, box_size, D)
    
    print("Box size :", box_size)
    print("D: ", D)
    print("Lattice Constant:", lattice_constant)
    print("C11:", C11)
    print("C12:", C12)
    print("C44:", C44)
    print("sigma11:", sigma11)
    
    # [Part A] : Creating Undeformed Configuration 
    
    print("******* PART A: Creating Undeformed Configuration ")
    
    undeformed_config_filename_In_extendXYZ = modelName+"_undeformed.xyz"
    
    # [Caution] Do not change the style of --undeformed_config_filename_In_extendXYZ--,
    # because it is linked with FORTRAN Code as well.
    # Otherwise, you should change the fortran code and compile it together.
    
    configLibrary.write_xyz_file(positions, atom_ids, undeformed_config_filename_In_extendXYZ,"Si")
    
    # [Part B] : Applying deformation using the analytic solution 
    
    print("******* PART B: Applying deformation ")
    ftr_input_file = f"./infiles/{modelName}.in"
    
    with open(ftr_input_file, "r") as infile:
        subprocess.run(["./deformPlate_lammpsOut"], stdin=infile)
    
    # [Part C] : Running lammps minimization code 
    
    print("******* PART C: Running lammps minimization code ")
    
    FileName_lammps_in = "lammps_minimization.in"
    
    FileName_ftr_out = modelName + "_forlammps_deformed.txt"
    FileName_lammps_before = modelName + "_before_minimization.lammpstrj"
    FileName_lammps_after = modelName + "_after_minimization.lammpstrj"
    
    run_lammps_command = "lmp -var modelName " + modelName \
                        + " -var inFileName " + FileName_ftr_out \
                        + " -var beforeFileName " + FileName_lammps_before \
                        + " -var afterFileName " + FileName_lammps_after \
                        + " < " + FileName_lammps_in #  + "< lammps_minimization.in"
    
    subprocess.run(run_lammps_command,shell=True)
    
    
    #[Part D]: Converting to MDstresslab file 
    
    print("******* PART D: Converting the lammps output files to MDstresslab file ")
    
    FileName_mdstresslab_config = modelName + "_config.data" 
    
    run_dump2str_command = "./dump2str " + FileName_lammps_before + " " + FileName_lammps_after + " " + FileName_mdstresslab_config \
                            + " 2 Si Si"
    
    subprocess.run(run_dump2str_command,shell=True)
     
    
