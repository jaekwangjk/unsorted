//
//  bc_and_rhs.h
//  
//
//  Created by Kim, Jaekwang on 7/23/17.
//
//

#ifndef bc_and_rhs_h
#define bc_and_rhs_h

#include <stdio.h>

using namespace dealii;

double U_inflow=10.0;


//Declaration
template <int dim>
class ConstantUx : public Function<dim>
{
public:
    double U_inflow;
    ConstantUx  (double U) : Function<dim>(dim+1) {
        U_inflow = U;
    }
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
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


//defines Newtonian Freespace soluton

template <int dim>
class NT_Freespace_Solution : public Function<dim>
{
public:
    NT_Freespace_Solution  () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
    
    
};

template <int dim>
double
NT_Freespace_Solution<dim>::value (const Point<dim>  &p,
                                const unsigned int component) const
{
    
    double U=U_inflow;
    double obj_rad=0.2; //Object Radius
    double r_sph = sqrt( pow(p[0],2) + pow(p[1],2) );
    
    double phi = 2*3.141592653589793238462643 - acos(p[1]/r_sph);
    
    
    //Newtonian Exact Velocity
    double q_r_sph=U*cos(phi)*(1+ 0.5*pow(obj_rad,3)/pow(r_sph,3)-1.5*obj_rad/r_sph );
    double q_phi= -U*sin(phi)*(1 -0.25 *pow(obj_rad,3)/pow(r_sph,3)-0.75*obj_rad/r_sph);
    
    
    double q_x = q_r_sph * (-sin(phi)) + q_phi * (-cos(phi));
    double q_z = q_r_sph * cos(phi) + q_phi * (-sin(phi));
    
    double pres= -(3/2)*1*(U)* obj_rad *cos(phi) / pow(r_sph,3); //questionable
    
    double stream  = 0.5 *U* ( std::pow(r_sph,2) + std::pow(obj_rad,3)/(2*r_sph) - 3*obj_rad*r_sph/2) * std::pow(sin(phi),2);
    
    
    if (component == 0)
        return q_x;
    if (component == 1)
        return q_z;
    if (component == 2)
        return pres;  // it might be pressure value then?
    if (component == 3)
        return r_sph * phi;  //Manufactured Solution 1 for lambda  = r * phi
    else
        return 0;
    
}


template <int dim>
void
NT_Freespace_Solution<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = NT_Freespace_Solution<dim>::value (p, c);
}

/////////Boundary Surface Force Function- This will be replaced by non-linear term


///////

template <int dim>
class VelBndValues : public Function<dim>
{
public:
    VelBndValues () : Function<dim>(dim+1) {}
    
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
VelBndValues<dim>::value (const Point<dim>  &p,
                          const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
    
    if (component == 0)
        return U_inflow;
    else
        return 0;
}



template <int dim>
void
VelBndValues<dim>::vector_value (const Point<dim> &p,
                                 Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = VelBndValues<dim>::value (p, c);
}

////////////////


template <int dim>
class Slip_Condition : public Function<dim>
{
public:
    Slip_Condition (BlockVector<double> &solution,
                    FESystem<dim>        &fe,
                    MappingQ<dim>        &mapping,
                    DoFHandler<dim>      &dof_handler
                    ) : Function<dim>(dim+1), solution_obj(&solution),fe_obj(&fe),mapping_obj(&mapping),
                        dof_obj(&dof_handler) {}
    
    BlockVector<double> *solution_obj;  //should I use smart pointer?
    SmartPointer<const FESystem<dim>> fe_obj;
    SmartPointer<const Mapping<dim>> mapping_obj;
    SmartPointer<const DoFHandler<dim> > dof_obj;
    
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int component = 0) const;
    
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
};


template <int dim>
double
Slip_Condition<dim>::value (const Point<dim>  &p,
                          const unsigned int component) const
{
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
   
    //const FEValuesExtractors::Vector velocities (0);
    //const FEValuesExtractors::Scalar radialvel (0);
    //const FEValuesExtractors::Scalar pressure (dim);
    
    
    //const ComponentSelectFunction<dim> pressure_mask (dim, dim+1);
    //double pressure;
    Vector<double> podi2;
    
    Point<2> P2;
    P2[0]=0.3;
    P2[1]=0.3;
    podi2.reinit(3);
    
    VectorTools::point_difference (*mapping_obj,
                                   *dof_obj,
                              *solution_obj,
                              ZeroFunction<dim>(dim+1),
                              podi2,
                              P2);
    
    VectorTools::point_value (*mapping_obj,
                                   *dof_obj,
                                   *solution_obj,
                                   p,
                                   podi2
                                   );
    
                              
    // and then point graident and then bla bla bla

    if (component == 0)
        return 0;
    else
        return 0;
}



template <int dim>
void
Slip_Condition<dim>::vector_value (const Point<dim> &p,
                                 Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = Slip_Condition<dim>::value (p, c);
}



////////////////

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
