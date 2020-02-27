#include <iostream>
#include <fstream>
#include <ios>
#include <vector>

using namespace std;

int main (); 

int main ()
{
    
    int nx=50;
    int ny=50;
    int nz=50;
    
	double x; 
	double y; 
	double z; 
	double scalar; 
	
	fstream fout; 
	
	fout.open("working.vtu",ios::out);
	fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
    fout << "<UnstructuredGrid>" << endl;
	fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " << "NumberOfCells=\" " << (nx-1)*(ny-1)*(nz-1) <<" \"> "
    << endl;
	fout << "<PointData Scalars=\"scalars\">" << endl;
    fout << "<DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">" << endl;
    
    // use for loop to outputs data (nx*ny*nz)
    // scalar function values
    for(int i=0; i<nz; i++)
    {    for(int j=0; j<ny; j++)
            {
                for(int k=0; k<nx; k++)
                {
                    double x=(1.0*k)/nx; double y=(1.0*j)/ny; double z=(1.0*i)/nz;
                    double value = x*x + y*y + z*z;
                    fout << value <<" ";
                }
                fout << endl;
            }
        
    }
        

    fout << "</DataArray>" << endl;
    fout << "</PointData>" << endl;
    
    //Provide points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {fout << k << " " << j << " " << i << endl;
        }}}
 
        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }
    
    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;
    
        // If it works, change to cell type 12
        for(int i=0; i<nz-1; i++)
        {
            for(int j=0; j<ny-1; j++)
            {
                for(int k=0; k<nx-1; k++)
                {
                    int org= i*nx*ny + j*nx +k;
                    
                    fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << " "
                    << org+nx*ny << " " << org+nx*ny+1 << " " << org+nx*ny+nx+1 << " " << org+nx*ny+nx << endl;
        }}}
        fout << "</DataArray>" << endl;
    }
    
    
    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {       for(int k=0; k<nx-1; k++)
                    {
                      fout << count*8 << " ";
                      count+=1;
                    }
                    fout << std::endl;
            }
        }
        fout << "</DataArray>" << endl;
    }
    
    // Cell Type information
    {
      fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;
        
        for(int i=0; i<nz-1; i++)
        {   for(int j=0; j<ny-1; j++)
            {for(int k=0; k<nx-1; k++)
                { fout << 12 << " ";}
                fout << std::endl;}
        }
        
        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }
    
    //VTK file format end
    
	fout.close(); 

}


