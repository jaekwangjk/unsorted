#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<numeric>
#include<array>
#include<limits>
#include "heap/heap.h"
#include "Matrix.h"
#include "Index.h"

#define MAXVAL std::numeric_limits<double>::max()

/* Movie, does this part run through shell?
 
 std::string str = "ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s " + std::to_string(nx) +
 "x" + std::to_string(ny) + " -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4";
 FILE* pipeout = popen(str.c_str(), "w");
 unsigned char* pixels = new unsigned char[pcount];
 
 */

// Create a grid
// Define a function chi on the grid
// Input: grid size, velocity field, chi describing a surface, step size
// Output: chi (t)
//

// Create a grid
// Define a function chi on the grid
// Define a surface, by initializing chi on all outer and inner grid points (w.r.t surface) to 0 and 1 respectively
// (J.K.??????)

// Tag the grid points as "alive", "close" (dead points that share a neighbor with alive points) and "far" points
// Initialize u to zero at all alive points
// loop 
//     Find the value of u at all close points by solving the quradratic equation
//     Heapify close points with back pointers to the grid
//     move the grid associated to the root to "alive"
//     identify the "neighbors of the root" that are in "far", and "move" to "close" to form the new "close"
// end
// 
// Input: grid size, velocity field, chi describing a surface, step size
// Output: chi (t)
//
int mod(const int& , const int& );

int mod(const int& i, const int& j){
    int c = i%j ;
    return (c<0)? c+j : c;
}


int initial_profile(const std::array<double,3>& x)
{
    double func = (3.*x[0]+1.)*pow(1.-3.*x[0],3)-pow((9.*pow(x[1],2)-1),2);
    return (func <= 0.) ? 1 : 0;
}

void initialize_chi(Matrix<3>& chi,const std::array<double,3>& corner, const double& h )
    {
        /* JK ..???....what dose corner represent? */
    // initialize_chi //
    std::array<double,3> x;
    int boo;

    std::ofstream chiFile;
    chiFile.open("initialChi.txt");

    for(int ix= 0; ix<chi.size(); ix++)
    {
       /* JK ..???.. data structure of Matrix<3> chi? */
        
        x[0] = corner[0]+ ix*h - 1.0;
        for(int iy= 0; iy<chi[0].size(); iy++)
        {
            x[1] = corner[1]+ iy*h - 1.0;
            for(int iz= 0; iz<chi[0][0].size(); iz++)
            {
                x[2] = corner[2]+ iz*h;
                boo = initial_profile(x);
                // Alive: 0, Dead: 1 <--Is this correct? inside the tooth is
                // alive, and it propagates outwards,
                //chi(ix,iy,iz) = (boo==1) ? 0. : MAXVAL;
                chi(ix,iy,iz) = (boo==0) ? 0. : MAXVAL; //why MAXVAL is needed here?
                // it is like...outer tooth- boo==0, 0
                // inside the tooth, boo==1, chi(x,y)= MAXVAL
                chiFile << x[0] << ", " << x[1]  << ", " << chi(ix,iy,iz) << std::endl;
            }
        }
    }
    chiFile.close();
    }

///////////////////////////////////////////////////////////////////////////////////////
int main()
{
    typedef double data_t;

    const int grayscale = 255;
    const int dim = 3;
    const int numOfIterations = 1e6;
    int nx,ny,nz,indx;
    const std::array<data_t,dim> corner{0.,0.,0.};
    //JKnote::  what does corner mean?, data_t? 
    std::array<int,dim> grid_indx;
    std::array<double,dim> x;
    data_t h;

    nx = 100;
    ny = 100;
    nz = 1;
    assert(nx==ny);
    h  = 2.0/nx;

    // total number of grid points
    int pcount = nx*ny*nz;

    // Movie
    
    std::string str = "ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt gray -s " + std::to_string(nx) + 
                      "x" + std::to_string(ny) + " -r 30 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 out.mp4";
    FILE* pipeout = popen(str.c_str(), "w");
    unsigned char* pixels = new unsigned char[pcount];
    
    
    
    Matrix<dim> chi(nx,ny,nz);
    initialize_chi(chi,corner,h);
 

    MinHeap<data_t> close(pcount);
    HeapData<data_t> hdata; // Just one element// One saves it and push this to MinHeap

    // Add all boundary points to the heap
    for(int p= 0; p<pcount; p++)
    {
        grid_indx = index(p,nx,ny,nz);
        int ix = grid_indx[0]; int iy = grid_indx[1]; int iz = grid_indx[2];
        if (chi(ix,iy,iz) == 0      && (chi(mod(ix+1,nx),iy,iz) == MAXVAL || //it has periodic!!
                                        chi(mod(ix-1,nx),iy,iz) == MAXVAL ||
                                        chi(ix,mod(iy+1,ny),iz) == MAXVAL ||
                                        chi(ix,mod(iy-1,ny),iz) == MAXVAL ||
                                        chi(ix,iy,mod(iz+1,nz)) == MAXVAL ||
                                        chi(ix,iy,mod(iz-1,nz)) == MAXVAL))
        {
            // form the heap data k
            hdata.set_val(0.);
            hdata.set_bckPtr(p);
            // insert the key into the heap
            close.insertKey(hdata); //insert key increase heap size 1 and
            // change pixel
            pixels[p] = grayscale;
        }
    }
    // At this point the heap contains all the boundary points

    const std::array<int,4> xgrid = {-1, 1, 0, 0};
    const std::array<int,4> ygrid = { 0, 0, 1,-1};
    double chi_east,chi_west,chi_north,chi_south,chi_top,chi_bottom;
    double chi_x,chi_y,chi_z;
    double val;

    clock_t b,e;
    b=clock();
    // Loop begins here
    //     extract the root, i.e remove it from heap
    //     Find its neighbors that are in the interior
    //     Compute chi at all interior neighbors of the extracted root and add them to the heap
    // end loop
    HeapData<data_t> root;
    int it;
    for(it=0;it<numOfIterations && close.heap_size!=0 ;it++)
    {   //Previous Declation MinHeap<data_t> close(pcount);
        
        // Extract the root
        root = close.extractMin(); /// aren't they all 0???, which is minimum then?
        int bckPtr = root.get_bckPtr();
        grid_indx = index(bckPtr,nx,ny,nz);


        pixels[bckPtr] = grayscale;
        fwrite(pixels, 1, pcount, pipeout);

        int ix = grid_indx[0];
        int iy = grid_indx[1];
        int iz = grid_indx[2];

        chi(ix,iy,iz) = root.get_val();

        // indices of neighbours of (ix,iy,iz)
        int ixn,iyn,izn;
        data_t a,b,c,discriminant;

        // Look for its neighbors that are in the interior
        for (int k = 0; k<4; k++)
        {
            // Solve the quadratic for each of the 
            // interior neighbors (ixn,iyn,izn)
            ixn= mod(ix + xgrid[k],nx);
            iyn= mod(iy + ygrid[k],ny);
            izn= iz;
            indx = ixn*ny*nz+iyn*nz+izn;

            // check if (ixn,iyn,izn) is in the interior
            if (chi(ixn,iyn,izn) == MAXVAL || close.locations[indx] != -1)
            {
                chi_east = chi(mod(ixn+1,nx),iyn,izn);
                chi_west = chi(mod(ixn-1,nx),iyn,izn);
                chi_x = fmin(chi_east,chi_west);

                chi_north= chi(ixn,mod(iyn+1,ny),izn);
                chi_south= chi(ixn,mod(iyn-1,ny),izn);
                chi_y = fmin(chi_north,chi_south);

                chi_top= chi(ixn,iyn,mod(izn+1,nz));
                chi_bottom= chi(ixn,iyn,mod(izn-1,nz));
                chi_z = fmin(chi_top,chi_bottom);
                
                // This is choosing upwind ??
                /*h^2 = (chi-chi_x)^2 + (chi-chi_y)^2 + (chi-chi_z)^2
                2d: 0 = 2 chi^2 - 2*chi*(chi_x + chi_y) + (chi_x^2+chi_y^2-h^2)
                   => a = 2, b = - 2*(chi_x + chi_y), c = chi_x^2+chi_y^2-h^2
                   if chi_x = max, then a = 1, b = -2*chi_y, c = chi_y^2-h^2
                   if chi_y = max, then a = 1, b = -2*chi_x, c = chi_x^2-h^2*/

                a= 2.;
                b= -2.*(chi_x+chi_y);
                c= pow(chi_x,2)+pow(chi_y,2)-pow(h,2);
                discriminant = pow(b,2)-4.*a*c;

                if (discriminant>0) //Real solution exists...
                {
                    val = -b+sqrt(pow(b,2)-4.*a*c);
                    val = val/(2.*a);
                }
                else if (chi_x < chi_y)
                {
                    val = chi_x + h; //why val = chi_x +h ? meaning??
                }
                else
                {
                    val = chi_y + h; 
                }

                // if chi_x and chi_y = MAXVAL, then that means their common neighbor ixn,iyn,izn = MAXVAL
                // IMPOSSIBLE
                assert(!(chi_x == MAXVAL && chi_y == MAXVAL));

                // is the neighnor already in the heap? If yes
                if (close.locations[indx] != -1)
                {
                    // update the heap node
                    close.harr[close.locations[indx]].set_val(val);
                    assert(close.harr[close.locations[indx]].get_bckPtr() == indx);
                }
                // else add it to the heap
                else
                {
                    hdata.set_val(val);
                    hdata.set_bckPtr(indx);
                    close.insertKey(hdata);
                }
            }
        }
    }
    e=clock();
    std::cout << (e-b)/(CLOCKS_PER_SEC*1.0) << std::endl;;


    
    std::ofstream closeFile;
    closeFile.open("closeChi.txt");
    
    std::cout << "Number of iterations: " << it << std::endl;

    for(int ix= 0; ix<chi.size(); ix++)
    {
        x[0] = corner[0]+ ix*h - 1.0;
        for(int iy= 0; iy<chi[0].size(); iy++)
        {
            x[1] = corner[1]+ iy*h - 1.0;
            for(int iz= 0; iz<chi[0][0].size(); iz++)
            {
                x[2] = corner[2]+ iz*h;
                indx = ix*ny*nz+iy*nz+iz;
                closeFile << x[0] << ", " << x[1]  << ", " << chi(ix,iy,iz) << std::endl;
            }
        }
    }

     fflush(pipeout);
     pclose(pipeout);
     delete[] pixels;
    
    //std::cout <<"MAX value is "  << MAXVAL << std::endl;
    return 0;
}


