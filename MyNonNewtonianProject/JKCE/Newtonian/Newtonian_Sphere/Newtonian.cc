#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
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
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "headers/bc_and_rhs.h"



const unsigned int max_cycle=5;

namespace MyStokes
{
  using namespace dealii;
  template <int dim>
  struct InnerPreconditioner;

  template <>
  struct InnerPreconditioner<2>
  {
    typedef SparseDirectUMFPACK type;
  };

  template <>
  struct InnerPreconditioner<3>
  {
    typedef SparseILU<double> type;
  };

  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem (const unsigned int degree);
    void run ();
    
    double stress_rphi (const Point<dim>   &p);
      
  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void compute_drag ();
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();
      
    const unsigned int   degree;
      
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    MappingQ<dim>        mapping;
    DoFHandler<dim>      dof_handler;
    ConstraintMatrix     constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;
      
    std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
      
  };

  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Subscriptor
  {
  public:
    InverseMatrix (const Matrix         &m,
                   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const SmartPointer<const Preconditioner> preconditioner;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
                                                       const Preconditioner &preconditioner)
    :
    matrix (&m),
    preconditioner (&preconditioner)
  {}


  template <class Matrix, class Preconditioner>
  void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
                                                    const Vector<double> &src) const
  {
    SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
    SolverCG<>    cg (solver_control);

    dst = 0;

    cg.solve (*matrix, dst, src, *preconditioner);
  }

  template <class Preconditioner>
  class SchurComplement : public Subscriptor
  {
  public:
    SchurComplement (const BlockSparseMatrix<double> &system_matrix,
                     const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > A_inverse;

    mutable Vector<double> tmp1, tmp2;
  };



  template <class Preconditioner>
  SchurComplement<Preconditioner>::
  SchurComplement (const BlockSparseMatrix<double> &system_matrix,
                   const InverseMatrix<SparseMatrix<double>,Preconditioner> &A_inverse)
    :
    system_matrix (&system_matrix),
    A_inverse (&A_inverse),
    tmp1 (system_matrix.block(0,0).m()),
    tmp2 (system_matrix.block(0,0).m())
  {}


  template <class Preconditioner>
  void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
                                               const Vector<double> &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    A_inverse->vmult (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
  }

    
  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree)
    :
    degree (degree),
    triangulation (Triangulation<dim>::maximum_smoothing),
    fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1),
    mapping (degree),
    dof_handler (triangulation)
  {}

  template <int dim>
  void StokesProblem<dim>::setup_dofs ()
  {
    A_preconditioner.reset ();
    system_matrix.clear ();

    dof_handler.distribute_dofs (fe);
    DoFRenumbering::Cuthill_McKee (dof_handler);

    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);

    
    {
      constraints.clear ();

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar radialvel(1);
      
      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints);
   
            std::cout << "   No-Slip Boundary condition implied" << std::endl;
           
            //axis
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      15,
                                                      ZeroFunction<dim>(dim+1),
                                                      constraints,
                                                      fe.component_mask(radialvel));
            
            //inlet and outlet
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      12,
                                                      ConstantUx<dim>(U_inflow),  // Boundary Condition for the case sphere falling in tube
                                                      constraints,
                                                      fe.component_mask(velocities));
            
            //inlet and outlet
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      13,
                                                      ConstantUx<dim>(U_inflow),  // Boundary Condition for the case sphere falling in tube
                                                      constraints,
                                                      fe.component_mask(velocities));
            
            //wall
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      14,
                                                      ConstantUx<dim>(U_inflow),
                                                      constraints,
                                                      fe.component_mask(velocities));
            
            //Sphere-1
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      10,
                                                      ZeroFunction<dim>(dim+1),
                                                      constraints,
                                                      fe.component_mask(velocities));
        
    }

    constraints.close ();

    std::vector<types::global_dof_index> dofs_per_block (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
    const unsigned int n_u = dofs_per_block[0],
                       n_p = dofs_per_block[1];

    std::cout << "   Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << n_u << '+' << n_p << ')'
              << std::endl;

    {
      BlockDynamicSparsityPattern dsp (2,2);

      dsp.block(0,0).reinit (n_u, n_u);
      dsp.block(1,0).reinit (n_p, n_u);
      dsp.block(0,1).reinit (n_u, n_p);
      dsp.block(1,1).reinit (n_p, n_p);

      dsp.collect_sizes();

      DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
      sparsity_pattern.copy_from (dsp);
    }

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (2);
    solution.block(0).reinit (n_u);
    solution.block(1).reinit (n_p);
    solution.collect_sizes ();

    system_rhs.reinit (2);
    system_rhs.block(0).reinit (n_u);
    system_rhs.block(1).reinit (n_p);
    system_rhs.collect_sizes ();
  }


  template <int dim>
  void StokesProblem<dim>::assemble_system () // no-slip
  {
    const long double pi = 3.141592653589793238462643;
      
    system_matrix=0;
    system_rhs=0;

    QGauss<dim>   quadrature_formula(2*degree+2);
    QGaussLobatto<dim-1> face_quadrature_formula(2*degree+3);

    FEValues<dim> fe_values (mapping,fe, quadrature_formula,
                             update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();
      
      
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double> >      rhs_values (n_q_points,
                                                  Vector<double>(dim+1));

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar radialvel (1);
    const FEValuesExtractors::Scalar pressure (dim);

   
    std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);
    std::vector<double>                  phi_ur      (dofs_per_cell);
    std::vector<double>                  phi_uz      (dofs_per_cell);
    
    double                               r_point;
      
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        local_matrix = 0;
        local_rhs = 0;

        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
                r_point = fe_values.quadrature_point (q)[1]; // radial location
			
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                phi_p[k]         = fe_values[pressure].value (k, q);
                phi_ur[k]        = fe_values[radialvel].value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  {
                    local_matrix(i,j) += 2* pi * (2 * r_point * (symgrad_phi_u[i] * symgrad_phi_u[j])
                                                    + 2 * phi_ur[i] * phi_ur[j] / r_point
                                          - (r_point * div_phi_u[i] + phi_ur[i]) * phi_p[j]
                                          - phi_p[i] * (r_point * div_phi_u[j] + phi_ur[j])
                                          + r_point * phi_p[i] * phi_p[j])
                                         * fe_values.JxW(q);

                  }


                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                local_rhs(i) += 2*pi* r_point* pi* fe_values.shape_value(i,q) *
                                rhs_values[q](component_i) *
                                fe_values.JxW(q);
              }
          }
          
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }

      std::cout << "   Computing preconditioner..." << std::endl << std::flush;

      A_preconditioner
        = std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
      A_preconditioner->initialize (system_matrix.block(0,0),
                                  typename InnerPreconditioner<dim>::type::AdditionalData());

  }


  template <int dim>
  void StokesProblem<dim>::compute_drag ()
  {
    const long double pi = 3.141592653589793238462643;
    
    QGauss<dim-1>   quadrature_formula_face(degree+2);

    double viscous_drag=0;
    double pressure_drag=0;
    double total_drag =0;

    FEFaceValues<dim> fe_face_values (mapping,fe, quadrature_formula_face,
				      update_JxW_values |
				      update_quadrature_points |
				      update_gradients |
				      update_values |
				      update_normal_vectors);

    const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;
    const unsigned int   n_q_face_points = fe_face_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar radialvel (1);
    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<double>                  local_pressure_values (n_q_face_points);
    std::vector<SymmetricTensor<2,dim>>  local_sym_vel_gradient (n_q_face_points);
    Tensor<1,dim>                        normal;
    
    double                               r_point;
    double                               drag=0;
		

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {

	for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
	  if (cell->face(face_no)->boundary_id()==10)
	  {
	    
	    fe_face_values.reinit (cell, face_no);	    	   	    
	    fe_face_values[pressure].get_function_values (solution,local_pressure_values);
	    fe_face_values[velocities].get_function_symmetric_gradients (solution,local_sym_vel_gradient);
	    for (unsigned int q=0; q<n_q_face_points; ++q)
	      {
              
                normal = fe_face_values.normal_vector (q);
                r_point = fe_face_values.quadrature_point (q)[1];
              
                pressure_drag += 2.* pi * r_point  * (normal[0]*local_pressure_values[q]) * fe_face_values.JxW (q) ;
                viscous_drag += 2.* pi * r_point * (-2.*normal[0]*local_sym_vel_gradient[q][0][0]
                                                    -2.*normal[1]*local_sym_vel_gradient[q][0][1]
                                                    )*fe_face_values.JxW (q);
              
                drag += 2.*pi*r_point*(normal[0]*local_pressure_values[q] -
                                       2.*normal[0]*local_sym_vel_gradient[q][0][0] -
                                       2.*normal[1]*local_sym_vel_gradient[q][0][1]
                                       )*fe_face_values.JxW (q);

	      }
	  }
      
      }


    const double a = 0.025;
    const double Drat = a/0.1;
    const double K1 = 1./(+1.
			  -2.10443*pow(Drat,1)
			  +2.08877*pow(Drat,3)
			  -0.94813*pow(Drat,5)
			  -1.372  *pow(Drat,6)
			  +3.87   *pow(Drat,8)
			  -4.19   *pow(Drat,10)
			  );
      
    const double drag_exact = K1* 6.*pi*a*U_inflow;
      
      std::cout << std::fixed << std::setprecision(10) <<
      "   DRAG = " << drag << "    " <<
      "   Exact Solution Drag = " << drag_exact << std::endl;
      /*
      std::cout << std::fixed << std::setprecision(10) <<
      "   pressure drag = " << pressure_drag<< "    " <<
      "   viscous drag = " << viscous_drag <<
      std::endl;
       */
  }
    
   
  
  template <int dim>
  void StokesProblem<dim>::solve ()
  {
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

      SolverControl solver_control (solution.block(1).size(),
                                    1e-6*schur_rhs.l2_norm());
      SolverCG<>    cg (solver_control);

      SparseILU<double> preconditioner;
      preconditioner.initialize (system_matrix.block(1,1),
                                 SparseILU<double>::AdditionalData());

      InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (system_matrix.block(1,1), preconditioner);

      cg.solve (schur_complement, solution.block(1), schur_rhs,
                m_inverse);

      constraints.distribute (solution);

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
  void
  StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
  {
      std::vector<std::string> solution_names (dim, "velocity");
      solution_names.push_back ("pressure");
      
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);
      
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names,
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << Utilities::int_to_string (refinement_cycle, 3)
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }


  template <int dim>
  void
  StokesProblem<dim>::refine_mesh ()
  {
      
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar radialvel (0);
    const FEValuesExtractors::Scalar pressure (dim);

    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(degree+1),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell,
                                        fe.component_mask(pressure));

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.0);
    triangulation.execute_coarsening_and_refinement ();
      
    
  }


  template <int dim>
  void StokesProblem<dim>::run ()
  {
    {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 4;

      GridIn<dim> grid_in;
      grid_in.attach_triangulation (triangulation);
      std::ifstream input_file("one_half_sphere.msh");

      Assert (dim==2, ExcInternalError());
      
      grid_in.read_msh (input_file);

      const Point<2> center (0,0);
      // spherical minifold cannot compute no-normal flux alone.
      // you need to add boundary description
      static const SphericalManifold<dim> manifold_description;
      static HyperShellBoundary<dim> boundary;
        
      triangulation.set_manifold (10, manifold_description);
      triangulation.set_manifold (10, boundary);
        
        
    }

    triangulation.refine_global (0);

    for (unsigned int refinement_cycle = 0; refinement_cycle<max_cycle;
         ++refinement_cycle)
      {
        std::cout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
           // triangulation.refine_global();
            refine_mesh ();
        
        std::cout << "   setting up system..." << std::endl << std::flush;
        setup_dofs ();
        std::cout << "   Assembling..." << std::endl << std::flush;
        assemble_system ();
        std::cout << "   Solving..." << std::endl  << std::flush;
        solve ();
        std::cout << "   Computing drag..." << std::endl  << std::flush;
        compute_drag ();
        std::cout << "   Refinement..." << std::endl  << std::flush;
        output_results (refinement_cycle);
        std::cout << std::endl;
        
          
      } //end of refinement cycle.
      
      
 
  }
}

int main ()
{
  try
    {
      using namespace dealii;
      using namespace MyStokes;

      StokesProblem<2> flow_problem(2);
      flow_problem.run ();
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



//Newtonain Stress Solution


/*
double  tau_analytic;
{//analytic solution
    double r_sph = sqrt( pow(elements[track_id][0] ,2) + pow(elements[track_id][1],2) );
    double phi = 2*3.141592653589793238462643 - acos(elements[track_id][1]/r_sph);
    double obj_rad=0.2;
    double W=1;
    
    double q_r_sph=W*cos(phi)*(1+ 0.5*pow(obj_rad,3)/pow(r_sph,3)-1.5*obj_rad/r_sph );
    double q_phi= -W*sin(phi)*(1 -0.25 *pow(obj_rad,3)/pow(r_sph,3)-0.75*obj_rad/r_sph);
    
    double strain_rr = 2.0*W*cos(phi)*( (3.0*obj_rad)/(2.0*pow(r_sph,2))-(3.0*pow(obj_rad,3))/(2.0*pow(r_sph,4)))   ;
    double strain_rtheta = (-1.5) * W * pow(obj_rad,3) / pow(r_sph,4) * sin(phi);
    double strain_tt =1.5 * W * cos(phi) * (pow(obj_rad,3)/pow(r_sph,4)-obj_rad/pow(r_sph,2));
    double strain_phiphi = 2.0 * (q_r_sph/r_sph +  q_phi/r_sph/tan(phi));
    
    tau_analytic= sqrt( (pow(strain_rr,2)
                         + 2*pow(strain_rtheta,2)
                         + pow(strain_tt,2)
                         + pow(strain_phiphi,2))/2
                       );
    
}*/
