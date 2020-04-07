//
//  polycrystal.h
//  myproject
//
//  Created by Jaekwang Kim on 2/13/20.
//  Copyright © 2020 Jaekwang Kim. All rights reserved.
//
//  2D randomly distriubted Polyscrytal analysis

#ifndef polycrystal_h
#define polycrystal_h

#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ios>

using namespace std;

/*
template <typename T>
void SwapEnd(T& var)
{
    char* varArray = reinterpret_cast<char*>(&var);
    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
        std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}
*/


double periodic_distance(double x1, double x2){
    
    double diff=fabs(x1-x2);
    double dist=0.5-fabs(0.5-diff);
    
    return dist;
    
}


template<int dim>
class Polycrystal
{
public:
    Polycrystal(const unsigned int n1,
                const unsigned int n2,
                const unsigned int lcount
                );
    
    void Generate_RandomCrystal();
    void Output_vtk(string s); /*vtk output*/
    
    void import_Polycrystal(double *Xangles_pointer,
                            double *Yangles_pointer,
                            double *Zangles_pointer,
                            double *labels_pointer);
    
    void Export_2d_Polycrystal_vtu(string file_name);
     void Export_2d_Polycrystal_binary_vtu(string file_name);
    
    void Export_Orientation_info(string file_name);
    void Import_Orientation_info(string file_name);
    
    void Import_from_vtk(string file_name);
    
    void Free_Memory();
    
    void Statistics_Area(string file_name);
    
    int *labels;
    double *Xangles;
    double *Yangles;
    double *Zangles;
    
    double *area;
    
private:
    
    int lcount;
    int pcount;
    
    
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    
    double Mean_Area;
    
};


//Constructor
template<int dim>
Polycrystal<dim>::Polycrystal(const unsigned int n1,
                              const unsigned int n2,
                              const unsigned int lcount)
:nx(n1), ny(n2), lcount(lcount)
{
    labels =new int[nx*ny];
    Zangles=new double[lcount];
    area = new double[lcount];
}


template<int dim>
void Polycrystal<dim>::Generate_RandomCrystal()
{
    double seeds[2*lcount];
    double max=0;
    
    for(int l=0;l<lcount;l++){
        // Generate random seeds and angles
        seeds[2*l]= rand()/(RAND_MAX*1.0);
        seeds[2*l+1]=rand()/(RAND_MAX*1.0);
        
        Zangles[l]=(M_PI/2) * rand()/(RAND_MAX*1.0);
        
        if(Zangles[l]>max){
            max=Zangles[l];}
    }
    
    for(int i=0;i<ny;i++){
        for(int j=0;j<nx;j++)
        {
            double x=j/(nx*1.0); double y=i/(ny*1.0);
            double min=FLT_MAX;
            int minIndex=0;
            
            for(int l=0;l<lcount;l++){
                 
                double xdist=periodic_distance(x,seeds[2*l]);
                double ydist=periodic_distance(y,seeds[2*l+1]);
                double dist=xdist*xdist+ydist*ydist;
                    
                if(dist<min){
                min=dist;
                minIndex=l;}
            }
            
            labels[i*nx+j]=minIndex;
        }
    }
}


template<int dim>
void Polycrystal<dim>::Statistics_Area(string file_name)
{
  
    for(int i=0; i<lcount; i++)
    {
        //zeroing
        area[i]=0.0;
    }
    
    
    double dA= (1.0/nx) * (1.0/ny);
    
    std::cout <<"Examine Grain statistics"<< std::endl;
    std::cout <<"dA= " << dA << std::endl; 
    
    for(int j=0;j<ny;j++){
        for(int k=0;k<nx;k++)
        {
            int current_label=labels[j*nx+k];
            area[current_label]+=1.0 * dA;
        }
        
    }
    
    
    //Find alive grains and mean area
    
    int num_live_grain=0;
    double Sum_area=0.0;
    
    for(int i=0; i<lcount; i++)
    {
        std::cout << "Area of ["<<i<<"] Grain : " << area[i] <<std::endl; 
        
        if(area[i]>1e-10)
        {
            num_live_grain+=1;
            Sum_area+=area[i];
        }
    }
    
    Mean_Area = Sum_area / num_live_grain;
    
    std::cout << "Mean area : " << Mean_Area <<std::endl; 
    
    fstream dataout;
    dataout.open(file_name,ios::out);
    
    for(int i=0; i<lcount; i++)
    {
        dataout << i << " " << area[i]/Mean_Area <<std::endl;
    }
     
    dataout.close();
    //You may want to check this statistics with excel ?
}


template<int dim>
void Polycrystal<dim>::Export_Orientation_info(string file_name)
{
    cout <<"Orientatoin_Info"<<endl;
    
    fstream fout;
    fout.open(file_name,ios::out);
    
    //Angle info
    fout << lcount << std::endl;
    for(int i=0; i<lcount; i++)
    {
        fout << Zangles[i] <<std::endl;
    }
    
    //label info
    for(int j=0;j<ny;j++){
        for(int k=0;k<nx;k++)
        {
            fout << labels[j*nx+k] << " " ;
        }
        fout <<std::endl;
    }
    
    fout.close();
}


template<int dim>
void Polycrystal<dim>::Import_Orientation_info(string file_name)
{
    //Lcount
    //Zlabels[lcount]
    //Label[i*y+x]
    
    cout <<"Import ORIT. info "<<endl;
    std::string varstr;
    ifstream fin;
    fin.open(file_name);
    
    int total_n;
    fin >> total_n; getline(fin,varstr);
    
    //Angle info
    for(int i=0; i<lcount; i++)
    {
        fin >> Zangles[i];  getline(fin,varstr);
    }
    
    //label info
    for(int j=0;j<ny;j++){
        for(int k=0;k<nx;k++)
        {
            fin >> labels[j*nx+k];
        }
        getline(fin,varstr);
    }
    
    fin.close();

}

template<int dim>
void Polycrystal<dim>::Import_from_vtk(string file_name)
{
    cout <<"Import vtu File"<<endl;
}




template<int dim>
void Polycrystal<dim>::Export_2d_Polycrystal_vtu(string file_name)
{
    cout <<"Import vtu File"<<endl;
  
    nz=1;
    
    fstream fout;
    fout.open(file_name,ios::out);
  
    //Provide Data
    {
        fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
        fout << "<UnstructuredGrid>" << endl;
        fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
        "NumberOfCells=\"" << (nx-1)*(ny-1)<<"\"> " << endl;
        fout << "<PointData Scalars=\"scalars\">" << endl;
        fout << "<DataArray type=\"Float32\" Name=\"Orientation\" format=\"ascii\">" << endl;
        
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nx; k++)
            {
                //theta is magnitude of orientation respect to global rotation system
                //, which is chosen as the Identity tensor
                //You lose information of the direction of misorientation
                double theta=sqrt(Zangles[labels[j*nx+k]]*Zangles[labels[j*nx+k]]);
                
                fout << theta << " " ;
            }
            fout << endl;
            
        }
        fout << "</DataArray>" << endl;
        fout << "</PointData>" << endl;
    }
    
    
    //Provide Coordinate Points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
        }}}
        
        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }
    
    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                int org= j*nx +k;
                fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << endl;
            }
        }
        fout << "</DataArray>" << endl;
    }
    
    
    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                fout << count*4 << " ";
                count+=1;
            }
            fout << std::endl;
        }
        fout << "</DataArray>" << endl;
    }
    
    // Cell Type information
    {
        fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;
        
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                fout << 9 << " ";
            }
            fout << std::endl;
        }
        
        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }

    fout.close(); 
}


template<int dim>
void Polycrystal<dim>::Export_2d_Polycrystal_binary_vtu(string file_name)
{
    cout <<"Import vtu File"<<endl;
    
    nz=1;

     
    fstream fout;
    fout.open(file_name,ios::out| std::ios::app | std::ios::binary);

    double theta[nx*ny];
    
    for(int j=0; j<ny; j++)
    {
        for(int k=0; k<nx; k++)
        {
            theta[nx*j+k]=Zangles[labels[j*nx+k]];
        }
    }
    //Provide Data
    {
        fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
        fout << "<UnstructuredGrid>" << endl;
        fout << "<Piece NumberOfPoints=\" "<< nx*ny*nz<<" \" " <<
        "NumberOfCells=\"" << (nx-1)*(ny-1)<<"\"> " << endl;
        fout << "<PointData Scalars=\"scalars\">" << endl;
        fout << "<DataArray type=\"Float32\" Name=\"Orientation\" format=\"binary\">" << endl;
        
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nx; k++)
            {
                //theta is magnitude of orientation respect to global rotation system
                //, which is chosen as the Identity tensor
                //You lose information of the direction of misorientation
                SwapEnd(theta[nx*j+k]);
                fout.write((char*)&theta[nx*j+k], sizeof(double));
            
                //ascii : that works
                //fout << theta[nx*j+k] << " " ;
                
                
            }
            fout << endl;
            
        }
        fout << "</DataArray>" << endl;
        fout << "</PointData>" << endl;
    }
    
    
    //Provide Coordinate Points
    {
        fout << "<Points>" <<endl;
        fout << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" <<endl;
        //Increase x->y>-z
        for(int i=0; i<nz; i++)
        {for(int j=0; j<ny; j++)
        {for(int k=0; k<nx; k++)
        {fout << (1.0*k)/nx << " " << (1.0*j)/ny << " " << (1.0*i)/nz << endl;
        }}}
        
        fout << "</DataArray>" << endl;
        fout << "</Points>" <<endl;
    }
    
    // Cell Connectivity
    {
        fout << "<Cells>" <<endl;
        fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" <<endl;
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                int org= j*nx +k;
                fout << org << " " << org+1 << " " << org+1+nx << " " << org+nx << endl;
            }
        }
        fout << "</DataArray>" << endl;
    }
    
    
    // Cell Type Offset
    {
        fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" <<endl;
        int count=1;
        
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                fout << count*4 << " ";
                count+=1;
            }
            fout << std::endl;
        }
        fout << "</DataArray>" << endl;
    }
    
    // Cell Type information
    {
        fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" <<endl;
        
        for(int j=0; j<ny-1; j++)
        {
            for(int k=0; k<nx-1; k++)
            {
                fout << 9 << " ";
            }
            fout << std::endl;
        }
        
        fout << std::endl;
        fout << "</DataArray>" << endl;
        fout << "</Cells>" << endl;
        fout << "</Piece>" << endl;
        fout << "</UnstructuredGrid>" << endl;
        fout << "</VTKFile>" << endl;
    }
    
    fout.close();
}



template<int dim>
void Polycrystal<dim>::Free_Memory()
{
    delete[] labels;
   // delete[] Xangles;
   // delete[] Yangles;
    delete[] Zangles;
    delete[] area; 
}



#endif /* polycrystal_h */



