//  bc_and_rhs.h
//  Created by Kim, Jaekwang on 7/23/17.

#ifndef bc_and_rhs_h
#define bc_and_rhs_h

#include <stdio.h>

using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide () : Function<dim>(dim+2) {}
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
    double f_x =0; double f_z=0;
    
    
    values(0) = f_x; //Force in r-direction
    values(1) = f_z; //Force in z-direction
    values(2) = 0; //RHS for continuity equation
}


// I should pass k_a from fluid class
template <int dim>
class RightHandSide_trpt : public Function<dim>
{
public:
    RightHandSide_trpt () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    void read_parameter(double fluid_k_a);
    double k_a;
};

template <int dim>
void RightHandSide_trpt<dim>::read_parameter (double fluid_k_a)
{
    k_a = fluid_k_a;
}

template <int dim>
double RightHandSide_trpt<dim>::value (const Point<dim>   &p,
                                const unsigned int  component) const
{
    Assert (component == 0, ExcIndexRange (component, 0, 1));
    
    double rhs = k_a;
    return rhs;
}


//Boundary condition of transport equation 
template <int dim>
class StructureBoundaryValues : public Function<dim>
{
public:
    StructureBoundaryValues () : Function<dim>(1) {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


template <int dim>
double
StructureBoundaryValues<dim>::value (const Point<dim> &p,
                                     const unsigned int /*component */) const
{
    double value = 1.0;
    return value;
}
////////


#endif /* header_h */
