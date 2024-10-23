//Write a function that read grid and split grid 
void splitGrid(string inFileName, string outFileName, int nOfGridFiles)
{
	int Tgrid;  
	fstream fin;
	fin.open(inFileName,std::ios::in);

	if ( fin.is_open() ) { // always check whether the file is open
	  fin >> Tgrid; // pipe file's content into stream
	  std::cout << "Tgrid Number is " << Tgrid <<std::endl; // pipe stream's content to standard output
	}else{
	  std::cout <<"Problem occured in splitGrid function"<<std::endl; 
	  exit(1); 
	}

	int nEachFile = Tgrid/nOfGridFiles; 
    int nLeftPoints = Tgrid%nOfGridFiles; 

	string currentFileName; 
	for(int i=0; i<nOfGridFiles; i++) {

		string currentFileName=outFileName; 
		char file_number[3];
		snprintf(file_number,3,"%02d",i);
  		currentFileName+=file_number;
  		currentFileName+=".dat";
		cout << "currentFileName: " << currentFileName << std::endl; 
    
		fstream fout;
		fout.open(currentFileName,std::ios::out);

		if(i< nOfGridFiles-1){

			fout << nEachFile << endl; 
			for(int j=0; j<nEachFile; j++){
				double x,y,z; 
				fin >> x >> y >> z;
				fout << x << std::endl; 
				fout << y << std::endl; 
				fout << z << std::endl; 
			}
		}else{

			fout << nEachFile + nLeftPoints << endl; 
			for(int j=0; j< nEachFile + (nLeftPoints) ; j++){
				double x,y,z; 
				fin >> x >> y >> z;
				fout << x << std::endl; 
				fout << y << std::endl; 
				fout << z << std::endl; 
			
			}
		}

		fout.close(); 
    }
	
	fin.close();

}
