// Developed by Jaekwang Kim
// Welcome to the world of viscoelastic fluids
// Coding from 18th May ~
// Viscoelastic Models 
// Oldroyd B model 
// upper convective maxwel polymer + a Newtonian solvent  
// Mixed finite element (also called viscous formulation) will be considered
// (a) Starting from large viscous contribution from a Newtonian Solvent 
// (b) Low Wissenberg Number
// (c) Stabilization Method 
// (d) Becareful on the LBB condition; Continuous Discretization of T_tensor?
//what is the convention of grad_u_[0][1] ?


// Check List, should I consider a periodic solution ?

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h> // CG used for viscous flow
#include <deal.II/lac/solver_gmres.h> //Matrix solver for transport equation, unsymmetric matrix
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h> //For the periodic boundary condition
#include <iostream>
#include <fstream>
#include <sstream>

#include "headers/post_processor.h"
#include "headers/precondition.h"
#include "headers/bcfunctions.h"
#include "headers/fluidclass.h"

namespace Viscoelastic
{
  using namespace dealii;
    
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem (const unsigned int degree);
    void run ();
        
  private:
        
    void readMesh ();
    void setup_dofs ();
    void setup_bc_constraints ();
    void assemble_stokes (); // Flow system
    void solve_stokes ();      // Solve with CG
    void assemble_polymer_system (); // UCM advection
    void solve_polymer_system ();           // Solve with GMRES

    void refine_mesh (const unsigned int refinement_cycle);
 
    void output_results (const unsigned int refinement_cycle);
    const unsigned int        degree;

    Triangulation<dim>        triangulation;
    FESystem<dim>             fe;
    DoFHandler<dim>           dof_handler;
    ConstraintMatrix          constraints;
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
        
    // Previous solution for Piccard iterations
    BlockVector<double>      previous_solution;
        
    AdvectionBoundaryValues<dim>   advection_boundary_values;
    ConstantU<dim>                 Uplate; // BC Function

    Oldroyd_FluidClass                  fluid;
    double lam1 = fluid.lam1;
    double lam2 = fluid.lam2;
    double etap = fluid.etap;
    double etas = fluid.etas;
    
    std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
};


//Class constructor
template <int dim> StokesProblem<dim>::StokesProblem (const unsigned int degree)
:degree(degree),
 triangulation(Triangulation<dim>::allow_anisotropic_smoothing),
 fe (FE_Q<dim>(degree+1),dim,    // u
 FE_Q<dim>(degree),1,    // p
 FE_Q<dim>(degree),dim*dim), //tau_p
 dof_handler (triangulation){}
// Discretization Space: I consider the space of tau_p is same with the graident space of u
  
  
template <int dim>
void StokesProblem<dim>::setup_bc_constraints ()
{
  constraints.clear ();
    
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar xvel(0);
  FEValuesExtractors::Scalar yvel(1);
            
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);

   
  // INLET (Full inlet analytical profile is provided)
  VectorTools::interpolate_boundary_values (dof_handler,
                                            4,
                                            Uinlet<dim>(),
                                            constraints,
                                            fe.component_mask(velocities));
    
  
  // ENDS - parallel flow
  VectorTools::interpolate_boundary_values (dof_handler,
                                            2,
                                            ZeroFunction<dim>(dim+1+dim*dim),
                                            constraints,
                                            fe.component_mask(yvel));
  
  
  //Periodicity
  /*
  std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>> make_pair;
      
  GridTools::collect_periodic_faces(dof_handler,4,2
                                   ,0 //direction
                                  ,make_pair);

    
  DoFTools::make_periodicity_constraints(dof_handler,2,4
                                         ,0,constraints,
                                         fe.component_mask(velocities));
  */
  
  // Top - Uniform wall motion:
  
  VectorTools::interpolate_boundary_values (dof_handler,
                                            3,
                                            Uplate,
                                            constraints,
                                            fe.component_mask(velocities));
  
  
  // Bottom - Fixed wall
  VectorTools::interpolate_boundary_values (dof_handler,
                                            1,
                                            ZeroFunction<dim>(dim+1+dim*dim),
                                            constraints,
                                            fe.component_mask(velocities));
    constraints.close ();
}

template <int dim>
void StokesProblem<dim>::setup_dofs (){

  std::cout << "   SET_UP_DOF -- START" <<std::endl;
    
  A_preconditioner.reset ();
  system_matrix.clear ();
    
  dof_handler.distribute_dofs (fe);
  DoFRenumbering::Cuthill_McKee (dof_handler);
        
  std::vector<unsigned int> block_component (dim+1+dim*dim,0); // initial 0 for u
  block_component[dim]   = 1;      // pressure
  for (unsigned int i=0;i<dim*dim;i++)
    block_component[dim+1+i] = 2; //block for polymer_T variable...
   
  DoFRenumbering::component_wise (dof_handler, block_component);
  std::vector<types::global_dof_index> dofs_per_block (3);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);

  const unsigned int n_u = dofs_per_block[0];
  const unsigned int n_p = dofs_per_block[1];
  const unsigned int n_s = dofs_per_block[2];
        
  std::cout << "   n_u: " << n_u << "   n_p: " << n_p << "   n_s: " << n_s << std::endl;

  // APPLY BC CONSTRAINTS
  setup_bc_constraints();

  {
    BlockDynamicSparsityPattern dsp (3,3);
            
    dsp.block(0,0).reinit (n_u, n_u);
    dsp.block(1,0).reinit (n_p, n_u);
    dsp.block(0,1).reinit (n_u, n_p);
    dsp.block(1,1).reinit (n_p, n_p);
    dsp.block(1,2).reinit (n_p, n_s);
    dsp.block(2,1).reinit (n_s, n_p);
    dsp.block(2,0).reinit (n_s, n_u);
    dsp.block(0,2).reinit (n_u, n_s);
    dsp.block(2,2).reinit (n_s, n_s);
        
    dsp.collect_sizes();
            
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, true);
    sparsity_pattern.copy_from (dsp);
  }
        
  system_matrix.reinit (sparsity_pattern);
        
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.block(2).reinit (n_s);
  solution.collect_sizes ();
        
  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.block(2).reinit (n_s);
  system_rhs.collect_sizes ();
        
  std::cout << " Active cells: "<< triangulation.n_active_cells() << std::endl
    << " DoFs: " << dof_handler.n_dofs() << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
    << std::endl;
}

  
// Assemble Momentum and Continuity Equation
template <int dim>
void StokesProblem<dim>::assemble_stokes () {

  system_matrix=0;
  system_rhs=0;
  const MappingQ<dim> mapping (degree);
  QGauss<dim>   quadrature_formula(2*degree+1);
        
  FEValues<dim> fe_values (mapping, fe, quadrature_formula,
			     update_values    |
			     update_quadrature_points  |
			     update_JxW_values |
			     update_gradients);
        
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();
        
  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);
        
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  std::vector<Vector<double>> rhs_values (n_q_points, Vector<double>(dim+1+dim*dim));
        
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
    
  const FEValuesExtractors::Scalar Txx (dim+1);
  const FEValuesExtractors::Scalar Txy (dim+2);
  const FEValuesExtractors::Scalar Tyx (dim+3);
  const FEValuesExtractors::Scalar Tyy (dim+4);
      
  std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
  std::vector<Tensor<2,dim> >          grad_phi_u   (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
  std::vector<double>                  phi_ur      (dofs_per_cell);

  std::vector<double>  phi_i_s_xx(dofs_per_cell);
  std::vector<double>  phi_i_s_xy(dofs_per_cell);
  std::vector<double>  phi_i_s_yx(dofs_per_cell);
  std::vector<double>  phi_i_s_yy(dofs_per_cell);
      
      
  //Previous Solution
  std::vector<SymmetricTensor<2,dim> > local_symgrad_u (n_q_points);
  std::vector<SymmetricTensor<2,dim> > local_T (n_q_points);
  std::vector<double>                  local_Txx (n_q_points);
  std::vector<double>                  local_Txy (n_q_points);
  std::vector<double>                  local_Tyx (n_q_points);
  std::vector<double>                  local_Tyy (n_q_points);
      
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell) {
  
	  fe_values.reinit (cell);
	  local_matrix = 0;
	  local_rhs = 0;
            
	  fe_values[velocities].get_function_symmetric_gradients(previous_solution, local_symgrad_u);
	  fe_values[Txx].get_function_values(previous_solution,local_Txx);
    fe_values[Txy].get_function_values(previous_solution,local_Txy);
	  fe_values[Tyx].get_function_values(previous_solution,local_Tyx);
    fe_values[Tyy].get_function_values(previous_solution,local_Tyy);
    
	  for (unsigned int q=0; q<n_q_points; ++q) {
    
	    for (unsigned int k=0; k<dofs_per_cell; ++k) {
        symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
        grad_phi_u[k]    = fe_values[velocities].gradient (k, q);
        div_phi_u[k]     = fe_values[velocities].divergence (k, q);
        phi_p[k]         = fe_values[pressure].value (k, q);
      }
                
	    
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
		    
        for (unsigned int j=0; j<=i; ++j) {
		   
          local_matrix(i,j) +=(
                      //(a) Solvent Contribution
                       2.*etas*(symgrad_phi_u[i]*symgrad_phi_u[j])
                       //(b) presure term and continuity term
                       - div_phi_u[i]*phi_p[j]
                       - phi_p[i]*div_phi_u[j]
                       + phi_p[i]*phi_p[j])  * fe_values.JxW(q);
        }
        
          // Look for your bilinear form and understand why it is
          local_rhs(i) -=( grad_phi_u[i][0][0]*(local_Txx[q])
                          +grad_phi_u[i][1][0]*(local_Txy[q])
                          +grad_phi_u[i][0][1]*(local_Tyx[q])
                          +grad_phi_u[i][1][1]*(local_Tyy[q])
                         ) *fe_values.JxW(q);
    
      } // end of dof per cell
      
	  }//end of quadrature point
  
  //symmetrize
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	      local_matrix(i,j) = local_matrix(j,i);
            
	  cell->get_dof_indices (local_dof_indices);
	
	  constraints.distribute_local_to_global(local_matrix,
      local_rhs,local_dof_indices,system_matrix, system_rhs);
 
  }
    A_preconditioner= std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
    A_preconditioner->initialize (system_matrix.block(0,0),typename InnerPreconditioner<dim>::type::AdditionalData());
}
    


template <int dim>
void StokesProblem<dim>::assemble_polymer_system ()
{
  std::cout << "   ASSEMBLE POLYMER CONSTITUTIVE" << std::endl;
        
  const MappingQ<dim> mapping (degree);
  QGauss<dim>   quadrature_formula(2*degree+1);
  QGauss<dim-1> face_quadrature_formula(2*degree+1);
        
  FEValues<dim> fe_values (mapping, fe, quadrature_formula,
			     update_values |
			     update_quadrature_points |
			     update_JxW_values |
			     update_gradients);
       
  FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature_formula,
				      update_values |
				      update_normal_vectors |
				      update_gradients |
				      update_quadrature_points |
				      update_JxW_values);
    
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();
  const unsigned int n_face_q_points = face_quadrature_formula.size();
    
  std::vector<Vector<double> > solution_values_face(n_face_q_points, Vector<double>(dim+1+dim*dim));
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        
  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);
   
  //Flow solution;
  std::vector<SymmetricTensor<2,dim> >  local_symgrad_phi_u (n_q_points);
  std::vector<Tensor<2,dim> >           local_grad_phi_u (n_q_points);
  std::vector<Tensor<1,dim> >           local_phi_u (n_q_points);
  std::vector<Tensor<1,1> >              local_pres (n_q_points);
        
  const FEValuesExtractors::Vector  velocities (0);
  const FEValuesExtractors::Scalar  Txx (dim+1);
  const FEValuesExtractors::Scalar  Txy (dim+2);
  const FEValuesExtractors::Scalar  Tyx (dim+3);
  const FEValuesExtractors::Scalar  Tyy (dim+4);
  const FEValuesExtractors::Scalar  pres (dim);
    
  std::vector<Tensor<1,dim>>  face_advection_directions (n_face_q_points);
  std::vector<Vector<double>> face_boundary_values (n_face_q_points,Vector<double>(dim+1+dim*dim));
  Tensor<1,dim> normal;
        
  double  vel_magnitude;
  double  phi_i_s_xx, phi_i_s_xy, phi_i_s_yy, phi_i_s_yx;
  double  phi_j_s_xx, phi_j_s_xy, phi_j_s_yy, phi_j_s_yx;

  Tensor<1,dim> grad_phi_i_s_xx, grad_phi_i_s_xy, grad_phi_i_s_yy, grad_phi_i_s_yx;
  Tensor<1,dim> grad_phi_j_s_xx, grad_phi_j_s_xy, grad_phi_j_s_yy, grad_phi_j_s_yx;
        
  typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
  
  for(; cell!=endc; ++cell) {
  
	fe_values.reinit (cell);
	local_matrix = 0;
	local_rhs = 0;
                        
	// Use of local solution 
	fe_values[velocities].get_function_symmetric_gradients (solution, local_symgrad_phi_u);
	fe_values[velocities].get_function_gradients (solution, local_grad_phi_u);
	fe_values[velocities].get_function_values (solution, local_phi_u);
    
	// STABILIZATION
	double delta = 0.1  * cell->diameter ();
          
	for (unsigned int q=0; q<n_q_points; ++q) {
  
    vel_magnitude = std::sqrt(local_phi_u[q]*local_phi_u[q]);

    for (unsigned int i=0; i<dofs_per_cell; ++i) {
		 
      phi_i_s_xx = fe_values[Txx].value(i,q);
      phi_i_s_xy = fe_values[Txy].value(i,q);
      phi_i_s_yx = fe_values[Tyx].value(i,q);
      phi_i_s_yy = fe_values[Tyy].value(i,q);
              
		  grad_phi_i_s_xx = fe_values[Txx].gradient (i, q);
		  grad_phi_i_s_xy = fe_values[Txy].gradient (i, q);
		  grad_phi_i_s_yx = fe_values[Tyx].gradient (i, q);
      grad_phi_i_s_yy = fe_values[Tyy].gradient (i, q);
              
		  for (unsigned int j=0; j<dofs_per_cell; ++j) {
		   
        phi_j_s_xx = fe_values[Txx].value(j,q);
        phi_j_s_xy = fe_values[Txy].value(j,q);
        phi_j_s_yx = fe_values[Tyx].value(j,q);
        phi_j_s_yy = fe_values[Tyy].value(j,q);
              
        grad_phi_j_s_xx = fe_values[Txx].gradient (j, q);
        grad_phi_j_s_xy = fe_values[Txy].gradient (j, q);
        grad_phi_j_s_yx = fe_values[Tyx].gradient (j, q);
        grad_phi_j_s_yy = fe_values[Tyy].gradient (j, q);
              
        local_matrix(i,j) +=(
          //self term
          + phi_i_s_xx * phi_j_s_xx
          + phi_i_s_xy * phi_j_s_xy
          + phi_i_s_yx * phi_j_s_yx
          + phi_i_s_yy * phi_j_s_yy
               
          + lam1*(
          //advection
          + phi_i_s_xx * (local_phi_u[q]*grad_phi_j_s_xx ) //u \cdot \grad \tau_xx
          + phi_i_s_xy * (local_phi_u[q]*grad_phi_j_s_xy )
          + phi_i_s_yx * (local_phi_u[q]*grad_phi_j_s_yx )
          + phi_i_s_yy * (local_phi_u[q]*grad_phi_j_s_yy )
			       
         // rotation term //
         // Gradient Tensor convention (dux/dy)=local_grad_phi_u[q][0][1])
          - phi_i_s_xx * ( 2.*local_grad_phi_u[q][0][0]*phi_j_s_xx +
                           local_grad_phi_u[q][0][1]*phi_j_s_yx +
                           local_grad_phi_u[q][0][1]*phi_j_s_xy )
          //I am skeptical about this
          - phi_i_s_xy * ( local_grad_phi_u[q][1][0]*phi_j_s_xx +
                           local_grad_phi_u[q][0][0]*phi_j_s_xy +  //(a)
                           local_grad_phi_u[q][1][1]*phi_j_s_xy +  //(b)
                           local_grad_phi_u[q][0][1]*phi_j_s_yy )
                    
          - phi_i_s_yx * ( local_grad_phi_u[q][1][0]*phi_j_s_xx +
                           local_grad_phi_u[q][0][0]*phi_j_s_yx +  //(c)
                           local_grad_phi_u[q][1][1]*phi_j_s_yx +  //(d)
                           local_grad_phi_u[q][0][1]*phi_j_s_yy  )
                            
          //actually (a+b+c+d) will be zero effectively
                            
          - phi_i_s_yy * ( 2.*local_grad_phi_u[q][1][1]*phi_j_s_yy +
                           local_grad_phi_u[q][1][0]*phi_j_s_yx +
                           local_grad_phi_u[q][1][0]*phi_j_s_xy )
                       
                            //std::cout << "local_"
          )
          // stabilization term: not proportional to lambda.
          // Check whether really it is in Crochet!!
          + delta/vel_magnitude *(
              local_phi_u[q]*grad_phi_i_s_xx //(u \cdot \nabla S)
             *local_phi_u[q]*grad_phi_j_s_xx
				     +
				      local_phi_u[q]*grad_phi_i_s_xy
				     *local_phi_u[q]*grad_phi_j_s_xy
				     +
             local_phi_u[q]*grad_phi_i_s_yx
				     *local_phi_u[q]*grad_phi_j_s_yx
             +
             local_phi_u[q]*grad_phi_i_s_yy
             *local_phi_u[q]*grad_phi_j_s_yy
				     )
        )*fe_values.JxW(q);
		  } //close 'j' cycle
      //I should have stabilization herm also in RHS
      //No, using the Crochet and Legat (1992), you do not need to consider
      //this for RHS.
      local_rhs(i) +=etap*(
                    + phi_i_s_xx * 2.*local_grad_phi_u[q][0][0]
                    + phi_i_s_xy * (local_grad_phi_u[q][1][0] + local_grad_phi_u[q][0][1])
                    + phi_i_s_yy * 2.*local_grad_phi_u[q][1][1]
                    + phi_i_s_yx * (local_grad_phi_u[q][1][0] + local_grad_phi_u[q][0][1])
                    )*fe_values.JxW(q);
    } //close i-lopop
  } // q

	  //inflow boundary condition
    //In this problem, we apply boundary condition weakly
    //(tau,S)_on boundary = (B.C, S)_on boundary
    
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no) {
	    if (cell->face(face_no)->boundary_id()==1) { // strong inflow boundary condition
        
        // inflow face
        fe_face_values.reinit (cell, face_no);
        fe_face_values.get_function_values (solution,solution_values_face);
                    
        advection_boundary_values.vector_value_list
          (fe_face_values.get_quadrature_points(),face_boundary_values);
                    
        for (unsigned int q=0; q<n_face_q_points; ++q) {
		      Tensor<1,dim> present_u_face;
                        
		      for (unsigned int d=0; d<dim; ++d)
		        present_u_face[d] = solution_values_face[q](d);
              
		      //if-inflow condition
		      if (fe_face_values.normal_vector(q) * present_u_face < 0) {
                     
            for (unsigned int i=0; i<dofs_per_cell; ++i) {
            
			        phi_i_s_xx = fe_face_values[Txx].value(i,q);
			        phi_i_s_xy = fe_face_values[Txy].value(i,q);
			        phi_i_s_yy = fe_face_values[Tyy].value(i,q);
			        phi_i_s_yx = fe_face_values[Tyx].value(i,q);

			        for (unsigned int j=0; j<dofs_per_cell; ++j) {
              
			          phi_j_s_xx = fe_face_values[Txx].value(j,q);
			          phi_j_s_xy = fe_face_values[Txy].value(j,q);
			          phi_j_s_yy = fe_face_values[Tyy].value(j,q);
			          phi_j_s_yx = fe_face_values[Tyx].value(j,q);
                    
			          local_matrix(i,j) -= (
						       phi_i_s_xx*phi_j_s_xx
						       +phi_i_s_xy*phi_j_s_xy
						       +phi_i_s_yy*phi_j_s_yy
						       +phi_i_s_yx*phi_j_s_yx
						    )*fe_face_values.JxW(q);     
              } //close cycle j
			  
			      local_rhs(i) -=
              (
			          + phi_i_s_xx * face_boundary_values[q][3]
                + phi_i_s_xy * face_boundary_values[q][4]
                + phi_i_s_yx * face_boundary_values[q][5]
                + phi_i_s_yy * face_boundary_values[q][6]
              )* fe_face_values.JxW(q);
            
			    } // i
        } // if-inflow
		  } // q
    } // cell
  } // face

	cell->get_dof_indices (local_dof_indices);
	constraints.distribute_local_to_global
    (local_matrix, local_rhs,local_dof_indices, system_matrix, system_rhs);
  }
}
    
    
template <int dim>
void StokesProblem<dim>::solve_stokes () {

  const InverseMatrix<SparseMatrix<double>,
  typename InnerPreconditioner<dim>::type>
  A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  
  {
    Vector<double> schur_rhs (solution.block(1).size());
    A_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);
            
    SchurComplement<typename InnerPreconditioner<dim>::type>
	  schur_complement (system_matrix, A_inverse);
            
    SolverControl solver_control ( 100 * solution.block(1).size(),
				     1e-6*schur_rhs.l2_norm()); // solver_control here
    
    SolverCG<>    cg (solver_control);
            
    SparseILU<double> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1),
    SparseILU<double>::AdditionalData());
            
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
	  m_inverse (system_matrix.block(1,1), preconditioner);
            
    cg.solve (schur_complement, solution.block(1),
              schur_rhs, m_inverse);
            
    constraints.distribute (solution);
            
    std::cout << "FLOW SOLVED:" << solver_control.last_step()
		<< " outer CG Schur complement iterations for pressure"
		<< std::endl;
  }
    
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
      
    A_inverse.vmult (solution.block(0), tmp);
    constraints.distribute (solution);
  }
  
}
    
    
template <int dim>
void StokesProblem<dim>::solve_polymer_system () {
  std::cout << "   SOLVE ADVECTION" << std::endl;
        
  SolverControl solver_control (std::pow(10,8), std::pow(10,-4));
  unsigned int restart = 500;
  SolverGMRES< Vector<double> >::AdditionalData gmres_additional_data(restart+2);
  SolverGMRES< Vector<double> > solver(solver_control, gmres_additional_data);
  //'gmres_additional_data' means how much temporary vector you will going to use for grmres solver
  //the more the faster, but the more memory will be consummed.
  //make preconditioner
  SparseILU<double>::AdditionalData additional_data(0,1000); // (0 , additional diagonal terms)
  SparseILU<double> preconditioner;
  preconditioner.initialize (system_matrix.block(2,2), additional_data);
  solver.solve (system_matrix.block(2,2), solution.block(2), system_rhs.block(2), preconditioner);
  constraints.distribute (solution);
}
 
template <int dim>
void StokesProblem<dim>::output_results (const unsigned int refinement_cycle)
{
  std::vector<std::string> solution_names(dim,"Velocity");
  solution_names.push_back ("Pressure");
  solution_names.push_back ("Txx");
  solution_names.push_back ("Txy");
  solution_names.push_back ("Tyx");
  solution_names.push_back ("Tyy");
    
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  
  //push back T
  for (unsigned int i=0; i<dim*dim; i++) {
      data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  }

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, solution_names,
  DataOut<dim>::type_dof_data,
  data_component_interpretation);
  data_out.build_patches ();
        
  std::ostringstream filenameeps;
  filenameeps << "solution-"<< Utilities::int_to_string (refinement_cycle, 2)<< ".vtk";
        
  std::ofstream output (filenameeps.str().c_str());
  data_out.write_vtk (output);
  //data_out.write_tecplot (output);
      
}
    
template <int dim>
void StokesProblem<dim>::readMesh()  //Generate mesh and designate boundary
{
  // Parallel plate
  const Point<2> end (5.0,1.0);
  const Point<2> start (-5.0,0.0);
  GridGenerator::hyper_rectangle (triangulation, start, end, false);
  
  for(typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
       cell!=triangulation.end(); ++cell) {
    
    for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
     // You need this line to ensure that this doesn't affect anything other than bondary faces!
      if (cell->face(f)->at_boundary() == true){
        //This is an internal faces!
        //You never want to set this to anything other than its default value.
        // The above check should make sure that we never could execute this part of the code.
        if (std::abs(cell->face(f)->center()[0]-0.5)< 1e-12 ){
          cell->face(f)->set_boundary_id(10);
        }
                  
        // Boundary faces
        double tol = 1e-12;
        
        if ( std::abs(cell->face(f)->center()[0]+5.0)< tol) {
          cell->face(f)->set_boundary_id(4);
        }
        
        if ( std::abs(cell->face(f)->center()[0]-5.0)< tol) {
          cell->face(f)->set_boundary_id(2);
        }
        
        if ( std::abs(cell->face(f)->center()[1]-0.0)< tol) {
          cell->face(f)->set_boundary_id(1);
        }
        
        if ( std::abs(cell->face(f)->center()[1]-1.0)< tol) {
           cell->face(f)->set_boundary_id(3);
        }
      
      }//end of 'cell->face(f)->at_boundary() == true)'
    }//end of facecs_per_cell
  }// end of cell iterator
      
  triangulation.refine_global (3);

}


template <int dim>
void StokesProblem<dim>::run () {

  readMesh(); //To read mesh file from outside
      
  for (unsigned int refinement_cycle = 0; refinement_cycle<1; ++refinement_cycle) {
  
    std::cout << "Refinement cycle " << refinement_cycle << std::endl;
                  
	  setup_dofs ();
	  
    if (refinement_cycle < 1 ) {
      previous_solution = solution;
	    previous_solution = 0.;
	  }       
	 
    //first solve
    assemble_stokes ();
	  solve_stokes ();
    assemble_polymer_system ();
	  solve_polymer_system ();

    BlockVector<double> difference;
	  int iteration_number=0;
	  previous_solution=solution;
    
          
	  do {
	    iteration_number +=1;
	    assemble_stokes ();
	    solve_stokes ();
	    assemble_polymer_system ();
	    solve_polymer_system ();
                
	    difference = solution;
	    difference -= previous_solution;
	    previous_solution=solution;
                
	    std::cout << "    Iteration Number: " << iteration_number
      << "    Difference Norm : " << difference.l2_norm()
      << std::endl << std::flush;
                
	  } while (difference.l2_norm()>pow(10,-6)* dof_handler.n_dofs());
    
    output_results (refinement_cycle);
   // refine_mesh (refinement_cycle) ;
     
  }//end of refinement cycle
}// end of run driver

}//end of Viscoelastic namespace

int main ()
{
  try{
  
    using namespace dealii;
    using namespace Viscoelastic;
    
    StokesProblem<2> flow_problem1(1);
    flow_problem1.run ();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
    std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
        
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
    std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
    return 1;
  }
    
  return 0;
}
