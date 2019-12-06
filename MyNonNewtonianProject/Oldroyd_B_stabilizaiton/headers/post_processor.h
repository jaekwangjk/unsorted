//
//  post_processor.h
//  
//
//  Created by Kim, Jaekwang on 7/23/17.
//
//

#ifndef post_processor_h
#define post_processor_h

#include <stdio.h>

using namespace dealii;

//header file tester :
/*
  int square (int a);

  int square (int a)
  {
  return a*a;
  }
*/
//end of header file tester


// Post processor 1 ; The second-invariant-tensor of Shear-rate field Calculator
template <int dim>
class ShearRate : public DataPostprocessorScalar<dim>
{
    
 public:
  ShearRate ();
  virtual
    void compute_derived_quantities_vector (const std::vector<Vector<double> >               &uh,
                                            const std::vector<std::vector<Tensor<1, dim> > > &duh,
                                            const std::vector<std::vector<Tensor<2, dim> > > &dduh,
                                            const std::vector<Point<dim> >                   &normals,
                                            const std::vector<Point<dim> >                   &evaluation_points,
                                            std::vector<Vector<double> >                     &computed_quantities) const;
};


template <int dim>

ShearRate<dim>::ShearRate ():DataPostprocessorScalar<dim> ("ShearRate", update_values | update_gradients | update_quadrature_points)

{}

template <int dim>
void ShearRate<dim>::compute_derived_quantities_vector (
                                                        const std::vector<Vector<double> >                 & uh,
                                                        const std::vector<std::vector<Tensor<1, dim> > >   & duh,
                                                        const std::vector<std::vector<Tensor<2, dim> > >   & /*dduh*/,
                                                        const std::vector<Point<dim> >                     & /*normals*/,
                                                        const std::vector<Point<dim> >                     & evaluation_points,
                                                        std::vector<Vector<double> >                       & computed_quantities) const
{
  for (unsigned int i=0; i<computed_quantities.size(); i++)
        
    { Assert(computed_quantities[i].size() == 1,
             ExcDimensionMismatch (computed_quantities[i].size(), 1));
      Assert(uh[i].size() == 7, ExcDimensionMismatch (uh[i].size(), 7));
        
      Tensor<2,dim> grad_u;
      for (unsigned int d=0; d<dim; ++d)
	grad_u[d] = duh[i][d]; // extract the first derivative of velocity component only (from full solution: origianl FEM system (v,p,strucutre))
        
      const SymmetricTensor<2,dim> shear_rate_tensor = symmetrize (grad_u);

      computed_quantities[i](0) = std::sqrt(2*(shear_rate_tensor*shear_rate_tensor));
    }
        
        
}


#endif /* header_h */
