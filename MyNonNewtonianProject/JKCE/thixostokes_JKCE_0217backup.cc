#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
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
#include <deal.II/grid/grid_out.h>
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
#include <iostream>
#include <fstream>
#include <sstream>

#include <deal.II/numerics/derivative_approximation.h>

#include "headers/post_processor.h"
#include "headers/precondition.h"
#include "headers/miscellaneous.h"

double U_inflow;

double G_drag_force_total;
double G_drag_force_pressure;
double G_drag_froce_viscous;

double k_d = 0.5; //destruction
double k_a = 1.0; //reconstruction
double c = 0.1;

//unsigned int max_refinement_cycle=7;
unsigned int max_refinement_cycle=1;


namespace MyStokes
{
    using namespace dealii;
    
    template <int dim>
    class StokesProblem
    {
    public:
        StokesProblem (const unsigned int degree);
        void run ();
        void second_run ();
        
    private:
        
        void setup_dofs ();
        void assemble_system ();
        void solve_flow ();
        void solve_transport (const unsigned int refinement_cycle);
        void assemble_transport_system ();
    
        void refine_mesh (const unsigned int refinement_cycle);
        void set_anisotropic_flags ();
        void compute_drag (const unsigned int refinement_cycle); // not written yet
        void output_results (const unsigned int refinement_cycle, const char diff_switch);
        const unsigned int   degree;
        
        double get_shear_rate (const SymmetricTensor<2,dim> symgrad_u, const double ur, const double r_point);
        void post_processing ();
        Vector<double> cellwise_shear_rate;
        Vector<double> cellwise_stress_field;
        Vector<double> cellwise_viscosity;
        
        Triangulation<dim>   triangulation;
        FESystem<dim>        fe;
        DoFHandler<dim>      dof_handler;
        ConstraintMatrix     constraints;
        BlockSparsityPattern      sparsity_pattern;
        BlockSparseMatrix<double> system_matrix;
        BlockVector<double> solution;
        BlockVector<double> system_rhs;
        
        BlockVector<double> exact_solution;
        
        BlockVector<double> previous_solution;
        
        
        Vector<double> cellwise_errors ;
        Vector<double> previous_str; //for comparison plot
    
        double eta_lambda= 100.0; //depends on stress
        double eta_0=1.0 ; //purely viscous
        
        std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
    };
    
    

    
    double viscosity_lambda (const double lambda, const double eta_0, const double eta_lambda)  // (lambda, model parameter #1, model parameter #2 ... )
    {
        double eta = eta_0 + eta_lambda * lambda;
        
        if (eta<eta_0)
            eta=eta_0;
    
        return eta;
    }
    
    
    
    template <int dim>
    class ConstantUz : public Function<dim>
    {
    public:
        ConstantUz  () : Function<dim>(dim+2) {}
        
        virtual double value (const Point<dim>   &p,
                              const unsigned int component = 0) const;
        
        virtual void vector_value (const Point<dim> &p,
                                   Vector<double>   &value) const;
        
        
    };
    
    template <int dim>
    double
    ConstantUz<dim>::value (const Point<dim>  &p,
                            const unsigned int component) const
    {
        if (component == 1)
            return U_inflow;
        else
            return 0;
    }
    
    
    template <int dim>
    void
    ConstantUz<dim>::vector_value (const Point<dim> &p,
                                   Vector<double>   &values) const
    {
        for (unsigned int c=0; c<this->n_components; ++c)
            values(c) = ConstantUz<dim>::value (p, c);
    }
    
    
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
        
        double value = 1;
    
        return value;
    }
    
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
    

    template <int dim>
    class RightHandSide_trpt : public Function<dim>
    {
    public:
        RightHandSide_trpt () : Function<dim>() {}
        
        virtual double value (const Point<dim>   &p,
                              const unsigned int  component = 0) const;
    };
    
    
    template <int dim>
    double
    RightHandSide_trpt<dim>::value (const Point<dim>   &p,
                               const unsigned int  component) const
    {
        Assert (component == 0, ExcIndexRange (component, 0, 1));
        
        double rhs = k_a ;

        return rhs;
        
        
    }
    
    
    template <int dim>
    StokesProblem<dim>::StokesProblem (const unsigned int degree)
    :
    degree (degree),
    triangulation (Triangulation<dim>::allow_anisotropic_smoothing),
        fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1,
        FE_Q<dim>(degree+1), 1),
    dof_handler (triangulation)
    {}
    
    
    template<int dim>
    double StokesProblem<dim>::get_shear_rate(const SymmetricTensor<2,dim> symgrad_u, const double ur, const double r_point)
    {
        return std::sqrt( 2.*(symgrad_u*symgrad_u + ur*ur/(r_point*r_point) ) );
    }
    
    template <int dim>
    void StokesProblem<dim>::setup_dofs ()
    {
        std::cout << "   -set_up_dof; start" <<std::endl;
    
        A_preconditioner.reset ();
        system_matrix.clear ();
        dof_handler.distribute_dofs (fe);
        DoFRenumbering::Cuthill_McKee (dof_handler);
        
        std::vector<unsigned int> block_component (dim+2,0);   // dim+1 int components, initialized with 0
        block_component[dim] = 1;
        block_component[dim+1] = 2; //Block for structure variable...
        
        DoFRenumbering::component_wise (dof_handler, block_component);
        
        std::vector<types::global_dof_index> dofs_per_block (3);
        DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);

        const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1], n_s=dofs_per_block[2];
        
        std::cout << "   -n_u: " << n_u << "   -n_p: " << n_p <<"   -n_s: " << n_s << std::endl;
        
        
        {
            constraints.clear ();
            
            FEValuesExtractors::Vector velocities(0);
            FEValuesExtractors::Scalar radialvel(0);
            
            DoFTools::make_hanging_node_constraints (dof_handler,
                                                     constraints);
            
            
            
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      2,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(radialvel));
          
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      4,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(radialvel));
           
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      3,
                                                      ConstantUz<dim>(),  // Boundary Condition for the case sphere falling in tube
                                                      constraints,
                                                      fe.component_mask(velocities));
            
            
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      5,
                                                      ConstantUz<dim>(),  // Boundary Condition for the case sphere falling in tube
                                                      constraints,
                                                      fe.component_mask(velocities));
            
           
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      1,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(velocities));
            
        }
        
        constraints.close ();
        
        {
            BlockDynamicSparsityPattern dsp (3,3);
            
            dsp.block(0,0).reinit (n_u, n_u);
            dsp.block(1,0).reinit (n_p, n_u);
            dsp.block(2,0).reinit (n_s, n_u);
            dsp.block(0,1).reinit (n_u, n_p);
            dsp.block(1,1).reinit (n_p, n_p);
            dsp.block(2,1).reinit (n_s, n_p);
            dsp.block(0,2).reinit (n_u, n_s);
            dsp.block(1,2).reinit (n_p, n_s);
            dsp.block(2,2).reinit (n_s, n_s);
            
            dsp.collect_sizes();
            
            DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, true);
            sparsity_pattern.copy_from (dsp);
            sparsity_pattern.compress(); // what does this do? meaningful?
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
        
        std::cout << "   Number of active cells: "
        << triangulation.n_active_cells()
        << std::endl
        << "   Number of degrees of freedom: "
        << dof_handler.n_dofs()
        << " (" << n_u << '+' << n_p << '+'<<n_s <<')'
        << std::endl;

        
        std::cout << "   -set_up_dof; end" <<std::endl;
       
    }
    
    
    template <int dim>
    void StokesProblem<dim>::assemble_system ()
    {
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
        
        const RightHandSide<dim>          right_hand_side;
        std::vector<Vector<double> >      rhs_values (n_q_points,
                                                      Vector<double>(dim+2));
        
        const FEValuesExtractors::Vector velocities (0);
        const FEValuesExtractors::Scalar radialvel (0);
        const FEValuesExtractors::Scalar pressure (dim);
        const FEValuesExtractors::Scalar structure (dim+1);
        
        std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
        std::vector<double>                  div_phi_u   (dofs_per_cell);
        std::vector<double>                  phi_p       (dofs_per_cell);
        std::vector<double>                  phi_ur      (dofs_per_cell);
        
        std::vector<double>                  local_previous_solution_structure (n_q_points);  //local lambda value
        
        double                               r_point;
        double                               lambda;
        
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
            fe_values[structure].get_function_values (previous_solution, local_previous_solution_structure);
            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                
                r_point = fe_values.quadrature_point (q)[0]; // radial location
             
                double structure = local_previous_solution_structure[q];
                
                double viscosity = viscosity_lambda (structure, eta_0, eta_lambda);
             
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                    symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                    div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                    phi_p[k]         = fe_values[pressure].value (k, q);
                    phi_ur[k]        = fe_values[radialvel].value (k, q);
                }
                
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    for (unsigned int j=0; j<=i; ++j)
                    {
                        local_matrix(i,j) += (2 * r_point * viscosity * (symgrad_phi_u[i] * symgrad_phi_u[j])
                                              + 2 * viscosity * phi_ur[i] * phi_ur[j] / r_point
                                              - (r_point * div_phi_u[i] + phi_ur[i]) * phi_p[j]
                                              - phi_p[i] * (r_point * div_phi_u[j] + phi_ur[j])
                                              + r_point * phi_p[i] * phi_p[j]
                                              )
                        * fe_values.JxW(q);
                        
                    }
                    
                    
                    const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                    local_rhs(i) += fe_values.shape_value(i,q) *
                    rhs_values[q](component_i) * r_point *
                    fe_values.JxW(q);
                }
            }
            
            for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=i+1; j<dofs_per_cell; ++j)
                    local_matrix(i,j) = local_matrix(j,i);
            
            cell->get_dof_indices (local_dof_indices);
            
            constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                    local_dof_indices,
                                                    system_matrix, system_rhs);
            
            
        }
        
        std::map<types::global_dof_index,double> boundary_values;
        
        A_preconditioner
        = std_cxx11::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
        A_preconditioner->initialize (system_matrix.block(0,0),
                                      typename InnerPreconditioner<dim>::type::AdditionalData());
        
     
    }
    
    template <int dim>
    void StokesProblem<dim>::assemble_transport_system ()
    {
        std::cout << "   -assemble transport equation" << std::endl;
        
        const MappingQ<dim> mapping (degree);
        QGauss<dim>   quadrature_formula(2*degree+2);
        QGauss<dim-1> face_quadrature_formula(2*degree+2);
        
        FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                                 update_values    |
                                 update_quadrature_points  |
                                 update_JxW_values |
                                 update_gradients);
       
        FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature_formula,
                                          update_values    | update_normal_vectors | update_gradients |
                                          update_quadrature_points  | update_JxW_values);
    
        const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
        const unsigned int   n_q_points      = quadrature_formula.size();
        const unsigned int   n_face_q_points = face_quadrature_formula.size();
    
        std::vector<Vector<double> > solution_values(n_q_points, Vector<double>(dim+2));
        std::vector<Vector<double> > solution_values_face(n_face_q_points, Vector<double>(dim+2));
        
        std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        
        FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
        Vector<double>       local_rhs (dofs_per_cell);
   
        const RightHandSide_trpt<dim>      right_hand_side_trpt;
        std::vector<double>         rhs_values_trpt (n_q_points);
        
        std::vector<SymmetricTensor<2,dim> > local_symgrad_phi_u (n_q_points);
        std::vector<double>  local_solution (n_q_points);
        std::vector<double>  local_phi_ur (n_q_points);
        
        const FEValuesExtractors::Vector velocities (0);
        const FEValuesExtractors::Scalar structure (dim+1);
        const FEValuesExtractors::Scalar radialvel (0);
        
        std::vector<double>                  structure_lambda      (dofs_per_cell);
        std::vector<Tensor<1,dim>>           structure_lambda_grad     (dofs_per_cell);
        
        StructureBoundaryValues<dim>         structure_boundary_values;
        std::vector<Tensor<1,dim>>           face_advection_directions (n_face_q_points);
        std::vector<double>                  face_boundary_values (n_face_q_points);
        Tensor<1,dim>                        normal;
        
        double r_point;
        double shear_rate;
        
        
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
            fe_values.reinit (cell);
            local_matrix = 0;
            local_rhs = 0;
            
            fe_values.get_function_values (solution, solution_values); //use result from flow solution
            right_hand_side_trpt.value_list(fe_values.get_quadrature_points(),rhs_values_trpt);
            
            fe_values[velocities].get_function_symmetric_gradients (solution, local_symgrad_phi_u);
            fe_values[radialvel].get_function_values (solution, local_phi_ur);
            
            
            double h = cell -> diameter();
            double delta = c  * pow(h,degree);  // pow(h,degree+1.5) may work more accurately, but I found that then the system is less likely solved to by GRES solver, (stabilization term decays too fast...)
            
            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                
                Tensor<1,dim> advection_field;
                
                for (unsigned int d=0; d<dim; ++d) //extract velocity solution
                {advection_field[d] = solution_values[q](d);}
                
                double vel_magnitude= std::sqrt(advection_field*advection_field);
                
                r_point = fe_values.quadrature_point (q)[0]; // radial location
                
                shear_rate = get_shear_rate(local_symgrad_phi_u[q], local_phi_ur[q], r_point);
                
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    const double phi_i_s = fe_values[structure].value(i,q);
                    
                    const Tensor<1,dim> grad_phi_i_s = fe_values[structure].gradient (i, q);
                    
                    
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                        const double phi_j_s = fe_values[structure].value(j,q);
                        const Tensor<1,dim> grad_phi_j_s = fe_values[structure].gradient (j, q);
                        

                        local_matrix(i,j) += (phi_i_s * ( advection_field * grad_phi_j_s  )
                                             + (k_d * shear_rate + k_a) *phi_i_s * phi_j_s
                                            
                                             + delta * (advection_field * grad_phi_i_s)
                                              * ( advection_field * grad_phi_j_s + (k_d * shear_rate + k_a) * phi_j_s)
                                              
                                              )
                                              *  r_point * fe_values.JxW(q);

                    } //close 'j' cycle
                    
                        local_rhs(i) += phi_i_s *
                                        rhs_values_trpt[q] * r_point *
                                        fe_values.JxW(q);
                   
                } //close 'i' cycle
            } //close 'q' cycle
            
            for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            {
                if (cell->face(face_no)->boundary_id()==3) // strongly designated  boundary condition for the sphere falling in tube
                {
                    
                        fe_face_values.reinit (cell, face_no);
                        fe_face_values.get_function_values (solution,solution_values_face); //advection field coming from flow solution
                    
                        structure_boundary_values.value_list (fe_face_values.get_quadrature_points(),face_boundary_values);
                    
                    for (unsigned int q=0; q<n_face_q_points; ++q)
                    {
                        r_point = fe_face_values.quadrature_point (q)[0]; // radial location
                        
                        double z_point = fe_face_values.quadrature_point (q)[1]; // vertical location
                        
                       
                        Tensor<1,dim> present_u_face; // Present u at face
                        
                        for (unsigned int d=0; d<dim; ++d)
                            present_u_face[d] = solution_values_face[q](d);
                        
                        if (fe_face_values.normal_vector(q) * present_u_face < 0){
                     
                            for (unsigned int i=0; i<dofs_per_cell; ++i)
                                
                            {
            
                                const double phi_i_s = fe_face_values[structure].value(i,q);
                                const Tensor<1,dim> grad_phi_i_s = fe_face_values[structure].gradient (i, q);
                            
                                
                                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                                    {
                                    
                                            const double phi_j_s = fe_face_values[structure].value(j,q);
                                            const Tensor<1,dim> grad_phi_j_s = fe_face_values[structure].gradient (j, q);
                                                        local_matrix(i,j) -= (present_u_face *
                                                      fe_face_values.normal_vector(q) *
                                                      phi_i_s *
                                                      phi_j_s *
                                                      r_point *
                                                      fe_face_values.JxW(q));
                             
                                    } //close cycle j
                                
                                local_rhs(i) -= (present_u_face *
                                                 fe_face_values.normal_vector(q) *
                                                 face_boundary_values[q]         *
                                                 phi_i_s *
                                                 r_point *
                                                 fe_face_values.JxW(q));
                                
                                
                            } //close cycle i
                        }//close if-inflow condition
                    } //close quadrature for loop
                }
                
           
             
                
            }

            cell->get_dof_indices (local_dof_indices);
         
            constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                    local_dof_indices,
                                                    system_matrix, system_rhs);
        }
     
    }
    
    
    template <int dim>
    void StokesProblem<dim>::solve_flow ()
    {
        std::cout << "   -solve stokes flow begins"<< std::endl;
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
            
            SolverControl solver_control ( 10 * solution.block(1).size(),
                                          1e-10*schur_rhs.l2_norm());
            SolverCG<>    cg (solver_control);
            
            SparseILU<double> preconditioner;
            preconditioner.initialize (system_matrix.block(1,1),
                                       SparseILU<double>::AdditionalData());
            
            InverseMatrix<SparseMatrix<double>,SparseILU<double> >
            m_inverse (system_matrix.block(1,1), preconditioner);
            
            cg.solve (schur_complement, solution.block(1), schur_rhs,
                      m_inverse);
            
            constraints.distribute (solution);
            
            std::cout << "   -"
            << solver_control.last_step()
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
    void StokesProblem<dim>::solve_transport (const unsigned int refinement_cycle)
    {
        std::cout << "   -solve transport begins"<< std::endl;
        SolverControl           solver_control (std::pow(10,6), system_rhs.block(2).l2_norm() * pow(10,-3));
        
        unsigned int restart = 500;
        SolverGMRES< Vector<double> >::AdditionalData gmres_additional_data(restart+2);
        SolverGMRES< Vector<double> > solver(solver_control, gmres_additional_data);
        
        SparseILU<double>::AdditionalData additional_data(0,1000); // (0 , additional diagonal terms)
        SparseILU<double> preconditioner;
        preconditioner.initialize (system_matrix.block(2,2), additional_data);
        
        solver.solve (system_matrix.block(2,2), solution.block(2), system_rhs.block(2), preconditioner);
        constraints.distribute (solution);
    }
    
    template <int dim>
    void
    StokesProblem<dim>::set_anisotropic_flags ()
    {
        typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active()
        ,endc=dof_handler.end();
        
        for (; cell!=endc; ++cell)
        {
            if (cell->refine_flag_set())
            {
                cell->set_refine_flag(RefinementCase<dim>::cut_y);
            }
        }
        
        
    }
    
    
    template <int dim>
    void
    StokesProblem<dim>::refine_mesh (const unsigned int refinement_cycle) //AMR based on VELOCITY GRADIENT
    {
        
        
        if(refinement_cycle <= 2)
        {
            GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                             cellwise_shear_rate,
                                                            0.9, 0.0);
        }
        
        else if(refinement_cycle>2 && refinement_cycle<5)
        {
            GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                             cellwise_shear_rate,
                                                             0.9, 0.0);
            set_anisotropic_flags();
        }
        
        else
        {
            GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                             cellwise_shear_rate,
                                                             0.05, 0.0);
            set_anisotropic_flags();
        }
     
        SolutionTransfer<dim,BlockVector<double>> solution_transfer(dof_handler);
        triangulation.prepare_coarsening_and_refinement ();
        solution_transfer.prepare_for_coarsening_and_refinement(solution);
        triangulation.execute_coarsening_and_refinement ();
        
        
        dof_handler.distribute_dofs (fe); //whats the role of this? I should understand this..
        DoFRenumbering::Cuthill_McKee (dof_handler);
        
        std::vector<unsigned int> block_component (dim+2,0);   // dim+1 int components, initialized with 0
        block_component[dim] = 1;
        block_component[dim+1] = 2; //Block for structure variable...
        DoFRenumbering::component_wise (dof_handler, block_component);
        
        std::vector<types::global_dof_index> dofs_per_block (3);
        DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
        const unsigned int n_u = dofs_per_block[0];
        const unsigned int n_p = dofs_per_block[1];
        const unsigned int n_s = dofs_per_block[2];
        
        previous_solution.reinit (3);
        previous_solution.block(0).reinit (n_u);
        previous_solution.block(1).reinit (n_p);
        previous_solution.block(2).reinit (n_s);
        previous_solution.collect_sizes ();
        
        solution_transfer.interpolate(solution, previous_solution);
        
    }
    
 
    template <int dim>
    void
    StokesProblem<dim>::output_results (const unsigned int refinement_cycle,const char diff_switch)
    {
	    if ( refinement_cycle == max_refinement_cycle -1 ) 
        {
            
    
        std::vector<std::string> solution_names (dim, "Velocity");
        solution_names.push_back ("Pressure");
        
        if(diff_switch=='n')
        {
            solution_names.push_back ("Structure");
            
        }else
        {
            solution_names.push_back ("Str_Diff");
            //solution.block(2) -= previous_str;
            
            
            std::vector<unsigned int> block_component (dim+2,0);   // dim+1 int components, initialized with 0
            block_component[dim] = 1;
            block_component[dim+1] = 2; //Block for structure variable...
            
            DoFRenumbering::component_wise (dof_handler, block_component);
            
            std::vector<types::global_dof_index> dofs_per_block (3);
            DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
            
            
            for(unsigned int i=0; i <dofs_per_block[2]; i++)
            {
                solution.block(2)[i]= std::abs(solution.block(2)[i]);
            }
            
            
        }
      
        
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation
        (dim, DataComponentInterpretation::component_is_part_of_vector);
        data_component_interpretation
        .push_back (DataComponentInterpretation::component_is_scalar);
        
        data_component_interpretation
        .push_back (DataComponentInterpretation::component_is_scalar);
        
        //ShearRate<dim> shear_rate; // This should be declared before Dataout Class
        //Viscosity<dim> viscosity_field;
        
        DataOut<dim> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (solution, solution_names,
                                  DataOut<dim>::type_dof_data,
                                  data_component_interpretation);
       
            
        //data_out.add_data_vector (solution, shear_rate); //shear rate is calculated from post_processor
        //data_out.add_data_vector (solution, viscosity_field); //shear rate is calculated from post_processor
            
        data_out.build_patches ();

        
        
        std::ostringstream filenameeps;
        filenameeps << "a_mean_kd_"<<k_d <<"ka" <<k_a << "U"<< U_inflow << ".vtk";
        std::ofstream output (filenameeps.str().c_str());
        data_out.write_vtk (output);
        
        }
        
    }
    
    
    template <int dim>
    void
    StokesProblem<dim>::compute_drag (const unsigned int refinement_cycle)
    {
        std::cout << "compute drag" << std::endl;
        
        const long double pi = 3.141592653589793238462643;
        
        const MappingQ<dim> mapping (degree);
        
        double r_point;
        double z_point;
        
        double theta_wise_pressure_drag;
        double theta_wise_viscous_drag;
        
        double viscous_drag=0;
        double pressure_drag=0;
        double total_drag =0;
        
        QGauss<dim-1>   quadrature_formula_face(2*degree+1); //check weather this is enough
        
        FEFaceValues<dim> fe_face_values (mapping, fe, quadrature_formula_face,
                                          update_JxW_values |
                                          update_quadrature_points |
                                          update_gradients |
                                          update_values |
                                          update_normal_vectors);
        
        const FEValuesExtractors::Vector velocities (0);
        const FEValuesExtractors::Scalar radialvel (0);
        const FEValuesExtractors::Scalar pressure (dim);
        const FEValuesExtractors::Scalar structure (dim+1);
        
        const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;
        const unsigned int   n_q_face_points = fe_face_values.n_quadrature_points;
        
        std::vector<double>                  local_pressure_values (n_q_face_points);
        std::vector<double>                  local_ur_values (n_q_face_points);
        std::vector<double>                  local_structure_values (n_q_face_points); // extract structure value to calculate local viscosity
        
        std::vector<SymmetricTensor<2,dim>>  local_sym_vel_gradient (n_q_face_points);
        Tensor<1,dim>                        normal;
        
        double lambda_s;
        double viscosity;
        double shear_rate;
        
        

    
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
            
            for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
                if (cell->face(face_no)->boundary_id()==1)
                {
                    
                    fe_face_values.reinit (cell, face_no);
                    fe_face_values[pressure].get_function_values (solution,local_pressure_values);
                    fe_face_values[velocities].get_function_symmetric_gradients (solution,local_sym_vel_gradient);
                    fe_face_values[radialvel].get_function_values (solution, local_ur_values);
                    fe_face_values[structure].get_function_values (solution, local_structure_values);
                   

                    for (unsigned int q=0; q<n_q_face_points; ++q)
                    {
                        lambda_s = local_structure_values [q];
                        
                        viscosity = viscosity_lambda (lambda_s, eta_0, eta_lambda);
                        
                        r_point = fe_face_values.quadrature_point (q)[0];
                        z_point = fe_face_values.quadrature_point (q)[1];
                        
                        
                        
                        normal = fe_face_values.normal_vector (q);
                        
                   
                        shear_rate = get_shear_rate(local_sym_vel_gradient[q], local_ur_values[q], r_point);
                        
                        pressure_drag += 2.* pi * r_point  * (normal[1]*local_pressure_values[q]) * fe_face_values.JxW (q) ;
                
                        viscous_drag += 2.* pi * r_point * (-2.*viscosity*normal[0]*local_sym_vel_gradient[q][0][1] +
                                                            -2.*viscosity*normal[1]*local_sym_vel_gradient[q][1][1]
                                                            )*fe_face_values.JxW (q);
                        
                        
                        
                                               
                    }

                    
                }
            
        }
        
        
        total_drag= pressure_drag + viscous_drag;
        
        
        
        std::cout << std::fixed << std::setprecision(10) <<
        "   DRAG = " << total_drag << "    " << "Pressure Drag = " << pressure_drag << "   " <<"Viscous Drag = " << viscous_drag <<"   " << std::endl ;
        
        G_drag_force_total=total_drag;
        G_drag_force_pressure=pressure_drag;
        G_drag_froce_viscous=viscous_drag;
        
    }
    
    
    template <int dim>
    void StokesProblem<dim>::run ()
    {
        
        {   //Read Mesh
            std::vector<unsigned int> subdivisions (dim, 1);
            subdivisions[0] = 4;
            
            GridIn<dim> grid_in;
            grid_in.attach_triangulation (triangulation);
            std::ifstream input_file("two_sphere.msh");
            Assert (dim==2, ExcInternalError());
            grid_in.read_msh (input_file);
            
            static const SphericalManifold<dim> manifold_description_1;
            triangulation.set_manifold (10, manifold_description_1);
            
            static const SphericalManifold<dim> manifold_description_2;
            triangulation.set_manifold (11, manifold_description_2);
            
            triangulation.refine_global (1);
            print_mesh_info (triangulation, "yourMesh.eps");
        }
        
        /*
        for (unsigned int refinement_cycle = 0; refinement_cycle<max_refinement_cycle ; ++refinement_cycle)
        {
            std::cout << "Refinement cycle " << refinement_cycle << std::endl;
            
            if (refinement_cycle > 0)
                refine_mesh (refinement_cycle);
            
            setup_dofs ();
            
            
            if (refinement_cycle == 0 )  //initialize previous solution at first cycle
            {
                previous_solution = solution;
                previous_solution = 0.;
            }
            
            assemble_system ();
            solve_flow ();
            assemble_transport_system ();
            solve_transport (refinement_cycle);
            
            BlockVector<double> difference;
            int iteration_number=0 ;
            
            previous_solution =solution ;
            
            do{
                iteration_number +=1;
                
                
                assemble_system ();
                solve_flow ();
                assemble_transport_system ();
                solve_transport (refinement_cycle);
                
                difference = solution;
                difference -= previous_solution;
                previous_solution=solution;
                
                std::cout << "   Iteration Number showing : " << iteration_number << "     Difference Norm : " << difference.l2_norm() << std::endl << std::flush;
                
            }while (difference.l2_norm()> 5* pow(10,-5)* dof_handler.n_dofs());  //defines iteration tolerance
            
           post_processing (); //calculate cellwise-shear-rate
           compute_drag (refinement_cycle);
           output_results (refinement_cycle,'n');
            
        }
         */
         
    } // end of StokesProblem<dim>::run ()
    
    
    template <int dim>
    void StokesProblem<dim>::second_run ()
    {
        previous_str = solution.block(2);
        
        setup_dofs ();
        assemble_system ();
        solve_flow ();
        assemble_transport_system ();
        solve_transport (max_refinement_cycle);
        
        BlockVector<double> difference;
        int iteration_number=0 ;
        
        previous_solution =solution ;
        
        do{
            iteration_number +=1;
            
            assemble_system ();
            solve_flow ();
            assemble_transport_system ();
            solve_transport (max_refinement_cycle);
            
            difference = solution;
            difference -= previous_solution;
            previous_solution=solution;
            
            std::cout << "   Iteration Number showing : " << iteration_number << "     Difference Norm : " << difference.l2_norm() << std::endl << std::flush;
            
        }while (difference.l2_norm()> 5* pow(10,-5)* dof_handler.n_dofs());  //defines iteration tolerance

         compute_drag (max_refinement_cycle-1);
         output_results (max_refinement_cycle-1,'y');
    }
}

int main ()
{
    try
    {
        using namespace dealii;
        using namespace MyStokes;
    
        U_inflow=2.5;
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
