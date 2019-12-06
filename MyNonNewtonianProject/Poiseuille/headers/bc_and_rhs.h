//  bc_and_rhs.h
//  Created by Kim, Jaekwang on 7/23/17.

#ifndef bc_and_rhs_h
#define bc_and_rhs_h

#include <stdio.h>

using namespace dealii;

double U_inflow=1;


template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
    
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                           const unsigned int /*component*/) const
{
    return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = RightHandSide<dim>::value (p, c);
}


#endif /* header_h */
