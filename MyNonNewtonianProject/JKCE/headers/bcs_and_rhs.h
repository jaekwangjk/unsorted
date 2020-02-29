#ifndef bc_and_rhs_h
#define bc_and_rhs_h
#include <stdio.h>

using namespace dealii;


//Declaration
template <int dim>
class ConstantUx : public Function<dim>
{
public:
    double U_inflow;
    ConstantUx  (double U) : Function<dim>(dim+2) {
        U_inflow = U;
    }
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};

template <int dim>
class StructureBoundaryValues : public Function<dim>
{
public:
    StructureBoundaryValues () : Function<dim>(1) {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};

template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide () : Function<dim>(dim+2) {}
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
    
};

template <int dim>
class RightHandSide_trpt : public Function<dim>
{
public:
    double k_a;
    
    RightHandSide_trpt (double parameter_k_a) : Function<dim>()
    {
        k_a=parameter_k_a;
    }
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


//Implementation

//ConstantUx<dim>Class
template <int dim>
double
ConstantUx<dim>::value (const Point<dim>  &p,
                        const unsigned int component) const
{
    if (component == 0)
        return U_inflow;
    else
        return 0;
}

template <int dim>
void
ConstantUx<dim>::vector_value (const Point<dim> &p,
                               Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = ConstantUx<dim>::value (p, c);
}


//StructuralBoundaryValues<dim>Class
template <int dim>
double
StructureBoundaryValues<dim>::value (const Point<dim> &p,
                                     const unsigned int /*component */) const
{
    double value = 1;
    return value;
}


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


template <int dim>
double
RightHandSide_trpt<dim>::value (const Point<dim>   &p,
                                const unsigned int  component) const
{
    Assert (component == 0, ExcIndexRange (component, 0, 1));
    double rhs = k_a ;
    return rhs;
    
}

#endif /* header_h */
