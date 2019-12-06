
#ifndef bc_function_h
#define bc_function_h

using namespace dealii;


template <int dim>
class Uinlet : public Function<dim>
{
public:
    double Ux;
    double h;
    
    Uinlet    () : Function<dim>(dim+1+dim*dim) {
        std::string varstr;
        
        std::ifstream inputfile("flow_prms.dat");
        inputfile >> Ux; getline(inputfile,varstr);
        inputfile >> h; getline(inputfile,varstr);
        
    }
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
    
};

template <int dim>
double
Uinlet<dim>::value (const Point<dim>  &p,
                    const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
    
    double y=p[1];
    
    if (component == 0)
        return Ux * y/h ; //exact solution for incoming flow
    else
        return 0;
}


template <int dim>
void
Uinlet<dim>::vector_value (const Point<dim> &p,
                           Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = Uinlet<dim>::value (p, c);
}


template <int dim>
class AdvectionBoundaryValues : public Function<dim>
{
  public:
    //AdvectionBoundaryValues () : Function<dim>(dim+1+dim*dim) {}
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    double   Ux,h;
    double   Txy,Txx,Tyy;
    
    AdvectionBoundaryValues  () : Function<dim>(dim+1+dim*dim) {
      std::string varstr;

    double  lam1, etap, etas,lam2;
       
      std::ifstream inputfile_flow("flow_prms.dat");
      inputfile_flow >> Ux; getline(inputfile_flow,varstr);
      
      inputfile_flow >> h; getline(inputfile_flow,varstr);
      std::cout
	  << "Ux:        " << Ux     << "\n"
	  << "h:         " << h     << "\n"
	  << std::endl;

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
        
      //Steady state analytic solution
        
       Txy = (etap)*Ux/h;  //T is only Polymeric Part of the stress
       Tyy = 0.;
       Txx = 2 * (etas+etap) *(lam1-lam2) * (Ux/h) * (Ux/h);
    
       std::cout << "boundary values at inflow is" << std::endl;
        
	   std::cout<< "Txy:         " << Txy     << "\n"
	   << "Tyy:         " << Tyy     << "\n"
       << "Txx:         " << Txx     << "\n"
	   << std::endl;
       
       std::cout << "Dimensionless Number" << std::endl;
        
       std::cout<< "Wi= lam1*(Ux/h)         " << lam1 * Ux     << "\n"
        << "alpha= etas/(etap+etas):         " << etas/(etap+etas)   << "\n"
        << std::endl;
}
};

    
template <int dim>
void AdvectionBoundaryValues<dim>::vector_value (const Point<dim> &p,
					      Vector<double>   &values) const
{
    values(0) = 0;
    values(1) = 0;
    values(2) = 0;
    values(3) = Txx;
    values(4) = Txy;
    values(5) = Txy;
    values(6) = Tyy;
}

template <int dim>
double AdvectionBoundaryValues<dim>::value (const Point<dim>   &p,
				       const unsigned int component ) const
{
    double value;
    
    if (component == 4 || component ==5)
      value = Txy;
    else if (component == 6)
      value = Tyy;
    else if (component == 3)
      value = Txx;
    else
      value = 0;

    std::cout << component << "..";
    
    return value;
}




template <int dim>
class ConstantU : public Function<dim>
{
   public:
   double U;
  
   ConstantU  () : Function<dim>(dim+1+dim*dim) {
   std::string varstr;
    
   std::ifstream inputfile("flow_prms.dat");
   inputfile >> U; getline(inputfile,varstr);
   std::cout << "Moving plate U:         " << U     << "\n"
	      << std::endl;
    
  }
        
  virtual double value (const Point<dim>   &p,
			const unsigned int component = 0) const;
        
  virtual void vector_value (const Point<dim> &p,
			     Vector<double>   &value) const;
  
};
    
template <int dim>
double
ConstantU<dim>::value (const Point<dim>  &p,
			const unsigned int component) const
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));
                
  if (component == 0)
    return U;
  else
    return 0;
}


template <int dim>
void
ConstantU<dim>::vector_value (const Point<dim> &p,
			       Vector<double>   &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = ConstantU<dim>::value (p, c);
}



#endif /* header_h */
