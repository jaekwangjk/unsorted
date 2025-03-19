#include "typedef.h"
#include "SpatialHash.h"
#include <iostream>
#include <fstream>
#include <iomanip>

int main()
{
	int numberOfPoints;

	std::ifstream fileA; 
	fileA.open("NN_full_grid_eps.dat",std::ios::in); 

	if ( fileA.is_open() ) { // always check whether the file is open
	  fileA >> numberOfPoints; // pipe file's content into stream
	  std::cout << "Tgrid Number is " << numberOfPoints <<std::endl;
	}else{
	  std::cout <<"Problem in reading the grid file. Exit the present program"<<std::endl; 
	  exit(1); 
	}

	MatrixXd configuration(numberOfPoints,DIM);
//	configuration= MatrixXd::Random(numberOfPoints,DIM);

	for(int i=0; i<numberOfPoints; i++){
		fileA >> configuration(i,0); 
		fileA >> configuration(i,1);
		fileA >> configuration(i,2);  
	}
/*
	for(int i=0; i<numberOfPoints; i++){
		std::cout << configuration(i,0) <<" " << configuration(i,1 ) << " " <<
		configuration(i,2) << std::endl; 
	}
*/
	fileA.close();

	Vector3d origin,step;
	origin.setConstant(0);
	step.setConstant(40);

	ConstSpatialHash hash(origin,step,configuration);	

	//to count the number of each grid in the part grid 
	//Segmentation Fault will happen if maxBinNumber is too small.
	int maxBinNumber =500; 
	int numOfpointInPgrid[maxBinNumber]; 
	for(int i=0; i<maxBinNumber; i++){
		numOfpointInPgrid[i]=0; 
	}

	int binNumber= 0;
	
	for(auto pair : hash.hashTable)
	{
		for (auto i_particle : pair.second) // loop of all atoms in side the box
		{// pair.first is triplet
		 // pair.second is a vector of integer
		 	numOfpointInPgrid[binNumber]+=1; 
			 // file << std::setw(10) << binNumber << std::setw(15) << i_particle << std::setw(10) << pair.first << std::setw(15) << configuration.row(i_particle) << std::endl;
		}
		binNumber++;
	}



	binNumber=0; 
	for(auto pair : hash.hashTable)
	{
		std::string subgridFileName="subGridEPS_"; 

		char file_number[3];
		snprintf(file_number,3,"%03d",binNumber);
		subgridFileName+=file_number; 
		subgridFileName+=".dat"; 

		std::fstream fout;
		fout.open(subgridFileName,std::ios::out);

		fout << numOfpointInPgrid[binNumber] << std::endl; 
		for (auto i_particle : pair.second) // loop of all atoms in side the box
		{
		 fout <<  configuration(i_particle,0) << std::endl;
		 fout <<  configuration(i_particle,1) << std::endl;
		 fout <<  configuration(i_particle,2) << std::endl;
		}
		binNumber++;

	}



}
