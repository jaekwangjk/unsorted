// ==========================================================    
class Oldroyd_FluidClass
{
  
public:

  double lam1; // Solution Relaxation Time
  double etap; // Polymer Viscosity
  double etas; // Solvent Viscosity
  double lam2; // Solution Retardation Time


  Oldroyd_FluidClass() // default constructor
    {
      std::string varstr;
      
      std::ifstream inputfile("fluid_prms.dat");
      inputfile >> lam1;       getline(inputfile,varstr);
      inputfile >> etap;       getline(inputfile,varstr);
      inputfile >> etas;       getline(inputfile,varstr);
      lam2 = (etas)/(etap+etas) * lam1;
        
    std::cout
    << "lam1:         " << lam1     << "\n"
    << "lam2:         " << lam2     << "\n"
    << "etap:         " << etap     << "\n"
    << "etas:         " << etas     << "\n"
	<< std::endl;
        
        

    }
    
};

