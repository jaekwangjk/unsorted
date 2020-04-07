//
//  write_binary.h
//  myproject
//
//  Created by Jaekwang Kim on 3/3/20.
//  Copyright © 2020 Jaekwang Kim. All rights reserved.
//

#include <iostream>
#include <fstream>

#ifndef write_binary_h
#define write_binary_h



template <typename T>
void SwapEnd(T& var)
{
    char* varArray = reinterpret_cast<char*>(&var);
    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
        std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}





using namespace std;



void good_write_binary_file()
{
     std::cout <<"Let us figure out this" << std::endl;
     
     double myarray[72] = {
     0.001,0.002,0,1,0,0,2,0,0,3,0,0,4,0,0,
     5,0,0,0,1,0,1,1,0,2,1,0,3,1,0,
     4,1,0,5,1,0,0,2,0,1,2,0,2,2,0,
     3,2,0,4,2,0,5,2,0,0,3,0,1,3,0,
     2,3,0,3,3,0,4,3,0,5,3,0};
     
     std::ofstream vtkstream;
     vtkstream.open("data/new_binary_test.vtk", std::ios::out | std::ios::app | std::ios::binary);
     if (vtkstream) {
     vtkstream<<"# vtk DataFile Version 2.0"<<"\n";
     vtkstream<<"Exemple"<<"\n";
     vtkstream<<"BINARY"<<"\n";
     vtkstream<<"DATASET STRUCTURED_GRID"<<std::endl;
     vtkstream<<"DIMENSIONS 6 4 1"<<std::endl;
     vtkstream<<"POINTS 24 double"<<std::endl;
     for (unsigned int i = 0; i < 72; ++i) {
     SwapEnd(myarray[i]);
     vtkstream.write((char*)&myarray[i], sizeof(double));
     }
     vtkstream.close();
     } else {
     std::cout<<"ERROR"<<std::endl;
     }
    
}

void write_binary_file()
{
    /*
    std::cout << "Write Binary function is executed" << std::endl;
    
    int n1=10;int n2=10;int n3=10;
    
    double *scalars;
    scalars=new double[n1*n2*n3];
    
    for(int i=0; i<n1*n2*n3; i++)
    {
        scalars[i]= ((double) rand() / (RAND_MAX));
    }
 
    
    //fstream wf("data/binary.vtu", ios::out | ios::binary);
    fstream wf("data/binary.vtu");
    */
    
    
    /*
    //wf << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
    wf << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
    //wf << "<Binary>" << endl;
    wf << "<ASCII>" << endl;
    wf << "<UnstructuredGrid>" << endl;
    wf << "<Piece NumberOfPoints=\" "<< n1*n2*n3<<" \" " <<
    "NumberOfCells=\" " << (n1-1)*(n2-1)*(n3-1) <<" \"> " << endl;
    wf << "<PointData Scalars=\"scalars\">" << endl;
    wf << "<DataArray type=\"Float32\" Name=\"  myname \" format=\"float\">" << endl;

    for(int i=0; i<n1*n2*n3; i++)
    {
        //wf.write(reinterpret_cast<char *>(&scalars[i]),sizeof(double));
        wf << scalars[i] <<" ";
    }*/
    
    /*
    wf << "DATASET STRUCTURED_POINTS\n";
    wf << "DIMENSIONS " << n1<< " "<< n2 <<" "<< n3 <<std::endl;
    wf << "ORIGIN " << 0 <<" " << 0 << " " << 0 <<std::endl;
    wf << "SPACING " << 0.01 << " "<< 0.01 << " "<<  0.01 <<std::endl;
    
    for(int i=0; i<n1*n2*n3; i++)
    {
        
        wf << scalars[i] <<" ";
    }
    
  
    
    wf.close();
     */
    
    
    /*
    fstream file;
    file.open("data/binary_test.vtk", std::ios::out | std::ios::app | std::ios::binary);

    file << "# vtk DataFile Version 2.0" << std::endl
    << "Comment if needed" << std::endl;
    
    file << "BINARY"<< std::endl << std::endl;
   
    file << "DATASET POLYDATA" << std::endl << "POINTS " << nb_particles << " float"  << std::endl;
    for(size_t cell_i=0;cell_i<n_cells;cell_i++)
    for(size_t i =0; i<cells[cell_i].size();++i)
    {

        double rx = cells[cell_i][field::rx][i];
        double ry = cells[cell_i][field::ry][i];
        double rz = cells[cell_i][field::rz][i];
     
        SwapEnd(rx);
        SwapEnd(ry);
        SwapEnd(rz);
     
        file.write(reinterpret_cast<char*>(&rx), sizeof(double));
        file.write(reinterpret_cast<char*>(&ry), sizeof(double));
        file.write(reinterpret_cast<char*>(&rz), sizeof(double));


     if(is_binary)
     file << std::endl;
     
     file << "POINT_DATA " << nb_particles << std::endl;
     
     if(has_id_field)
     {
     file<< "SCALARS index int 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;
     for(size_t cell_i=0;cell_i<n_cells;cell_i++)
     for(size_t i =0; i<cells[cell_i].size();++i)
     {
     if(is_binary)
     {
     uint64_t id = cells[cell_i][field::id][i];
     
     SwapEnd(id);
     
     file.write(reinterpret_cast<char*>(&id), sizeof(uint64_t));
     }
     else
     file << cells[cell_i][field::id][i] << std::endl;
     }
     if(is_binary)
     file << std::endl;
     }
     
     if(has_type_field)
     {
     file<< "SCALARS type int 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;
     for(size_t cell_i=0;cell_i<n_cells;cell_i++)
     for(size_t i =0; i<cells[cell_i].size();++i)
     {
     if(is_binary)
     {
     uint8_t type;
     
     SwapEnd(type);
     
     file.write(reinterpret_cast<char*>(&type), sizeof(uint8_t));
     }
     else
     file << static_cast<int>(cells[cell_i][field::type][i]) << std::endl;
     }
     if(is_binary)
     file << std::endl;
     }
     
     file.close();
     
     */

}






void write_binary_simple_file()
{
    /*
    # vtk DataFile Version 3.0
    vtk output
    ASCII
    DATASET POLYDATA
    POINTS 4 float
    0 0 0
    1 0 0
    1.1 1.1 0
    0 1 0
    POLYGONS 1 5
    4 0 1 2 3
    CELL_DATA 1
    POINT_DATA 4
    SCALARS nodal float
    LOOKUP_TABLE default
    0 1 2 1.1
    */
    
    std::ofstream wf;
    wf.open("data/binary.vtk", std::ios::out | std::ios::app | std::ios::binary);
    
    //std::ofstream wf("data/binary.vtk", std::ios::out | std::ios::app | std::ios::binary);
    
    /*
    wf << "# vtk DataFile Version 3.0" << endl;
    wf << "vtk output" << std::endl;
    wf << "ASCII" << std::endl;
    wf << "DATASET POLYDATA"<<endl;
    wf << "POINTS 4 float" <<endl;
    //wf << "0 0 0" <<endl;
    //wf << "1 0 0" <<endl;
    //wf << "1.1 1.1 0" <<endl;
    //wf << "0 1 0" <<endl;
    wf << "POLYGONS 1 5" <<endl;
    //wf << "4 0 1 2 3" <<endl;
    wf << "CELL_DATA 1" <<endl;
    wf << "POINT_DATA 4" <<endl;
    wf << "SCALARS nodal float" <<endl;
    wf << "LOOKUP_TABLE default" <<endl;
    //wf << "0 1 2 1.1" <<endl;
    */
    
    
    
    double line_one[3],line_two[3],line_three[3],line_four[3];
    unsigned int line_five[5];
    double line_six[4];
    
    line_one[0]=0.0; line_one[1]=0.0; line_one[2]=0.0;
    line_two[0]=1; line_two[1]=0.0; line_two[2]=0.0;
    line_three[0]=1.1; line_three[1]=1.1; line_three[2]=0.0;
    line_four[0]=0.0; line_four[1]=1; line_four[2]=0.0;
    
    line_five[0]=4; line_five[1]=0; line_five[2]=1;
    line_five[3]=2; line_five[4]=3;
    
    line_six[0]=0; line_six[1]=1; line_six[2]=2; line_six[3]=1.1;
    
    
    wf << "# vtk DataFile Version 2.0" << endl;
    wf << "vtk output" << std::endl;
    wf << "BINARY" << std::endl;
    wf << "DATASET POLYDATA"<<endl;
    wf << "POINTS 4 double" <<endl;
  
    
    for (unsigned int i = 0; i < 3; ++i) {
        std::cout <<line_one[i] << std::endl;
        SwapEnd(line_one[i]);
        wf.write((char*)&line_one[i], sizeof(double));
    }
    
     wf << endl;
    for (unsigned int i = 0; i < 3; ++i) {
        SwapEnd(line_two[i]);
        wf.write((char*)&line_two[i], sizeof(double));
    }
     wf << endl;

    for (unsigned int i = 0; i < 3; ++i) {
        SwapEnd(line_three[i]);
        wf.write((char*)&line_three[i], sizeof(double));
    }
        wf << endl;

    for (unsigned int i = 0; i < 3; ++i) {
        SwapEnd(line_four[i]);
        wf.write((char*)&line_four[i], sizeof(double));
    }
    
    wf << endl;
    wf << "POLYGONS 1 5" <<endl;
    
    /*
    for (unsigned int i = 0; i < 5; ++i) {
        SwapEnd(line_five[i]);
        wf.write(reinterpret_cast<char*>(&line_four[i]), sizeof(line_four[i]));
    }*/
    wf << "4 0 1 2 3" <<endl;
    
    wf << "CELL_DATA 1" <<endl;
    wf << "POINT_DATA 4" <<endl;
    wf << "SCALARS nodal double" <<endl;
    wf << "LOOKUP_TABLE default" <<endl;
    
    for (unsigned int i = 0; i < 4; ++i) {
        SwapEnd(line_six[i]);
        wf.write((char*)&line_six[i], sizeof(double));
    }
    wf << endl;
    //wf << "0 1 2 1.1" <<endl;
    
    wf.close();
}



#endif /* write_binary_h */
