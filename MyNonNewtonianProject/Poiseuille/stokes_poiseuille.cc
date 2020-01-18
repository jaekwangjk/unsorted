#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
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
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "headers/bc_and_rhs.h"

//Mesh type info
// 1: Rectangular
// 2: One cylinder (Gmsh)
// 3: ...

int mesh_type = 0;

const unsigned int max_cycle=1;
double c = 0.1; //Stabilization constant

// ==========================================================
class MooreThixotropy_Fluid
{
    
public:
    
    double eta_str; // Solution Relaxation Time
    double eta_infty; // Polymer Viscosity
    double k_d; // Solvent Viscosity
    double k_a; // Solution Retardation Time
    
    double lambda;
    double viscosity;
    
    MooreThixotropy_Fluid() // default constructor
    {
        std::cout << "Moore Thixotropy Fluid is constructed" <<std::endl;
        std::string varstr;
        std::ifstream inputfile("fluid_prms.dat");
        inputfile >> eta_str;       getline(inputfile,varstr);
        inputfile >> eta_infty;       getline(inputfile,varstr);
        inputfile >> k_d;       getline(inputfile,varstr);
        inputfile >> k_a;       getline(inputfile,varstr);
        std::cout
        << "eta_str:         " << eta_str     << "\n"
        << "eta_infty:         " << eta_infty    << "\n"
        << "k_d:         " << k_d     << "\n"
        << "k_a:         " << k_a     << "\n"
        << std::endl;
    }
    double compute_viscosity (double lambda);
};

double MooreThixotropy_Fluid::compute_viscosity (double lambda)
{
    double visosity = lambda * eta_str + eta_infty;
    
    //regularization
    if (viscosity < eta_infty)
        viscosity = eta_infty;
    
    return viscosity;
}

MooreThixotropy_Fluid fluid;

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
        
    private:
        
        void read_mesh (int mesh_type_num);
        void setup_dofs (int mesh_type_num);
        void assemble_stokes_system ();
        void solve_flow ();
        
        void assemble_transport_system ();
        void solve_transport ();
       
        void output_results (const unsigned int refinement_cycle);
        void refine_mesh ();
        
        void measure_pressure_profile_in_x (); //at the center line
        void measure_velocity_profile_in_y (); //at the outlet
        
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
        BlockVector<double> previous_solution;
        
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
    triangulation (Triangulation<dim>::allow_anisotropic_smoothing),
    fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1,
        FE_Q<dim>(degree+1), 1),
    mapping (degree),
    dof_handler (triangulation)
    {}
    
    
    template <int dim>
    void StokesProblem<dim>::setup_dofs (int mesh_type_num)
    {
        A_preconditioner.reset ();
        system_matrix.clear ();
        dof_handler.distribute_dofs (fe);
        DoFRenumbering::Cuthill_McKee (dof_handler);
        
        std::vector<unsigned int> block_component (dim+2,0);
        block_component[dim] = 1;
        block_component[dim+1] = 2;
        DoFRenumbering::component_wise (dof_handler, block_component);
        
        std::vector<types::global_dof_index> dofs_per_block (3);
        DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
        
        const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1], n_s=dofs_per_block[2];
    
        {
            constraints.clear ();
            
            FEValuesExtractors::Vector velocities(0);
            FEValuesExtractors::Scalar xvel(0);
            FEValuesExtractors::Scalar yvel(1);
            
            DoFTools::make_hanging_node_constraints (dof_handler,
                                                     constraints);
            
            //boundary_id =1,3 bottom and top wall
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      1,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(velocities));
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      3,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(velocities));
            if(mesh_type_num == 1)
            {
                //boundary_id =12 :: Cylinder
                VectorTools::interpolate_boundary_values (dof_handler,
                                                      12,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(velocities));
            }
            
            //inlet(4) and outlet(2) -
            /*
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      4,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(yvel));
            
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      2,
                                                      ZeroFunction<dim>(dim+2),
                                                      constraints,
                                                      fe.component_mask(yvel));
           */
            
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
            
            DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
            sparsity_pattern.copy_from (dsp);
            sparsity_pattern.compress();
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
    void StokesProblem<dim>::assemble_stokes_system ()
    {
        std::cout << "   -assemble stokes equation" << std::endl;
        
        system_matrix=0;
        system_rhs=0;
        const MappingQ<dim> mapping (degree);
        QGauss<dim>   quadrature_formula(degree+2);
        QGaussLobatto<dim-1> face_quadrature_formula(2*degree+1);
        
        FEValues<dim> fe_values (mapping,fe, quadrature_formula,
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
        
        std::vector<double> local_previous_solution_structure (n_q_points);
        
        double lambda;
        
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
            
            fe_values[structure].get_function_values (previous_solution,
                                                      local_previous_solution_structure);
            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
               
                double lambda = local_previous_solution_structure[q];
                double viscosity = fluid.compute_viscosity(lambda);
                
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                    symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                    div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                    phi_p[k]         = fe_values[pressure].value (k, q);
                }
                
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                        local_matrix(i,j) +=  (2.0  * viscosity *
                                              (symgrad_phi_u[i] * symgrad_phi_u[j])
                                               - (div_phi_u[i]) * phi_p[j]
                                               - phi_p[i] * (div_phi_u[j])
                                               + phi_p[i] * phi_p[j]) * fe_values.JxW(q);
                    }
                    
                    const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                    local_rhs(i) += fe_values.shape_value(i,q) *
                                    rhs_values[q](component_i) *
                                    fe_values.JxW(q);
                }
            }// close q
            
            
            for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            {
                if (cell->face(face_no)->boundary_id()==4)//inlet
                {
                  
                    fe_face_values.reinit (cell, face_no);
                    
                    for (unsigned int q=0; q<n_face_q_points; ++q)
                    {
                        for (unsigned int i=0; i<dofs_per_cell; ++i)
                        {
                            const Tensor<1, dim> phi_i_u =
                            fe_face_values[velocities].value(i, q);
                        
                            double pressure_in =20.0;
                    
                            local_rhs(i) +=
                            -(phi_i_u * fe_face_values.normal_vector(q) *
                              //note that it is (vector dot vector)
                              pressure_in * fe_face_values.JxW(q));
                        }
                    }
                }
            }
            
            cell->get_dof_indices (local_dof_indices);
            constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                    local_dof_indices,
                                                    system_matrix, system_rhs);
        }
        
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
        
        RightHandSide_trpt<dim>      right_hand_side_trpt;
        right_hand_side_trpt.read_parameter(fluid.k_a);
        
        std::vector<double>         rhs_values_trpt (n_q_points);
        
        std::vector<SymmetricTensor<2,dim> > local_symgrad_phi_u (n_q_points);
        std::vector<double>  local_solution (n_q_points);
        std::vector<double>  local_phi_ur (n_q_points);
        
        const FEValuesExtractors::Vector velocities (0);
        const FEValuesExtractors::Scalar structure (dim+1);
        const FEValuesExtractors::Scalar radialvel (0);
        
        std::vector<double>                  structure_lambda      (dofs_per_cell);
        std::vector<Tensor<1,dim>>           structure_lambda_grad     (dofs_per_cell);
        
        
        //can this just be ? StructureBoundaryValues<1>? It is actually a oneD
        //you need to consider when cleaning up the code
        StructureBoundaryValues<dim>         structure_boundary_values;
        
        
        std::vector<Tensor<1,dim>>           face_advection_directions (n_face_q_points);
        std::vector<double>                  face_boundary_values (n_face_q_points);
        Tensor<1,dim>                        normal;
        
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
            
            double h = cell -> diameter();
            double delta = c  * pow(h,degree);  // pow(h,degree+1.5) may work more accurately, but I found that then the system is less likely solved to by GRES solver, (stabilization term decays too fast...)
            
            for (unsigned int q=0; q<n_q_points; ++q)
            {
                Tensor<1,dim> advection_field;
                
                for (unsigned int d=0; d<dim; ++d) //extract velocity solution
                {advection_field[d] = solution_values[q](d);}
                
                double vel_magnitude= std::sqrt(advection_field*advection_field);
                double shear_rate = std::sqrt(local_symgrad_phi_u[q]* local_symgrad_phi_u[q] *2);
                
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                    const double phi_i_s = fe_values[structure].value(i,q);
                    const Tensor<1,dim> grad_phi_i_s = fe_values[structure].gradient (i, q);
                    
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                        const double phi_j_s = fe_values[structure].value(j,q);
                        const Tensor<1,dim> grad_phi_j_s = fe_values[structure].gradient (j, q);
                        
                        local_matrix(i,j) += (phi_i_s * ( advection_field * grad_phi_j_s  )
                                              + (fluid.k_d * shear_rate + fluid.k_a) *phi_i_s * phi_j_s
                                              //original test function
                                              
                                              + delta * (advection_field * grad_phi_i_s)
                                              * (advection_field * grad_phi_j_s + (fluid.k_d * shear_rate + fluid.k_a) * phi_j_s)
                                              //stabilized term (step-9 tutorial)
                                              )* fe_values.JxW(q);
                    }
                    local_rhs(i) += phi_i_s * rhs_values_trpt[q] * fe_values.JxW(q);
                }
            }// end of q cycle
            
            for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            {
                if (cell->face(face_no)->boundary_id()==4)
                {
                    fe_face_values.reinit (cell, face_no);
                    fe_face_values.get_function_values (solution,solution_values_face);
                    //advection field coming from flow solution
                    
                    structure_boundary_values.value_list (fe_face_values.get_quadrature_points(),face_boundary_values);
                    
                    for (unsigned int q=0; q<n_face_q_points; ++q){
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
                                                          fe_face_values.JxW(q));
                                    
                                } //close cycle j
                                
                                local_rhs(i) -= (present_u_face *
                                                 fe_face_values.normal_vector(q) *
                                                 face_boundary_values[q]         *
                                                 phi_i_s *
                                                 fe_face_values.JxW(q));
                            }//close cycle i
                        }//close if-inflow condition
                    }//close face_quadrature_loop
                } //end of cell-face with boundary_id() =4
            }//end of face cycle
            
            cell->get_dof_indices (local_dof_indices);
            
            constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                    local_dof_indices,
                                                    system_matrix, system_rhs);
            
        }
        
        
    }
    template <int dim>
    void StokesProblem<dim>::solve_flow ()
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
            
            SolverControl solver_control (solution.block(1).size()*100,
                                          1e-6*schur_rhs.l2_norm());
            SolverCG<>    cg (solver_control);
            
            SparseILU<double> preconditioner;
            preconditioner.initialize (system_matrix.block(1,1),
                                       SparseILU<double>::AdditionalData());
            
            InverseMatrix<SparseMatrix<double>,SparseILU<double> >
            m_inverse (system_matrix.block(1,1), preconditioner);
            
            cg.solve (schur_complement, solution.block(1),
                      schur_rhs, m_inverse);
            
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
    void StokesProblem<dim>::solve_transport ()
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
    StokesProblem<dim>::output_results (const unsigned int refinement_cycle)
    {
        
        std::vector<std::string> solution_names (dim, "velocity");
        solution_names.push_back ("pressure");
        solution_names.push_back ("structure");
        
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation
        (dim, DataComponentInterpretation::component_is_part_of_vector);
        data_component_interpretation
        .push_back (DataComponentInterpretation::component_is_scalar);
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
        << Utilities::int_to_string (refinement_cycle, 2)
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
                                                         0.3, 0.1);
        triangulation.execute_coarsening_and_refinement ();
    }
    
    template <int dim>
    void StokesProblem<dim>::read_mesh (int mesh_type_num)
    {
        
        if(mesh_type_num==0)
        {
            
            double Lx=10;
            double Ly=1;
            std::vector< unsigned int > repetitions_a(2); //number of discretization
            repetitions_a[0]=4;  //3;
            repetitions_a[1]=1;   //2;
            GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions_a,
                                                       Point<2>(0.0,0.0),
                                                       Point<2>(Lx,Ly));
            
            triangulation.refine_global(3);
            //To designate boundray condition--loop over all faces.
            
            for (typename Triangulation<dim>::active_cell_iterator
                 cell=triangulation.begin_active();
                 cell!=triangulation.end(); ++cell)
            {
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                    
                    if (cell->face(f)->at_boundary() == true)
                    {
                        
                        static const double tol = 1e-12;
                        
                        // Boundary faces
                        if ( std::abs(cell->face(f)->center()[0]-0.0)< tol )
                        {
                            //Inlet
                            cell->face(f)->set_boundary_id(4);
                           // std::cout << "boundary 4 is created" << std::endl;
                        }
                        if ( std::abs(cell->face(f)->center()[1]-0.0)< tol )
                        {
                            //Bottom wall
                            cell->face(f)->set_boundary_id(1);
                           // std::cout << "boundary 1 is created" << std::endl;
                        }
                        if ( std::abs(cell->face(f)->center()[0]-Lx)< tol )
                        {
                            //Outlet
                            cell->face(f)->set_boundary_id(2);
                           // std::cout << "boundary 2 is created" << std::endl;
                        }
                        if ( std::abs(cell->face(f)->center()[1]-Ly)< tol )
                        {
                            //Topwall
                            cell->face(f)->set_boundary_id(3);
                          //  std::cout << "boundary 3 is created" << std::endl;
                        }
                    }
                }
            }
            
            
        }else if(mesh_type_num==1) //one cylinder
        {
            GridIn<dim> gridin;
            gridin.attach_triangulation(triangulation);
            std::ifstream f("one_cylinder.msh");
            gridin.read_msh(f);
        
            const Point<2> center (2.5,2.5);
            static const SphericalManifold<dim> manifold_description_1;
            triangulation.set_manifold (12, manifold_description_1);
        }
        
    }
    
    template <int dim>
    void StokesProblem<dim>::run ()
    {
        
        read_mesh (mesh_type);
        triangulation.refine_global(2);
        setup_dofs (mesh_type);
        
        // previoius solution initialization
        previous_solution = solution;
        previous_solution = 0.;
        
        assemble_stokes_system ();
        solve_flow ();
        assemble_transport_system ();
        solve_transport ();
       
        BlockVector<double> difference;
        int iteration_number=0 ;
        previous_solution =solution ;
        
        do{
            iteration_number +=1;
            
            assemble_stokes_system ();
            solve_flow ();
            assemble_transport_system ();
            solve_transport ();
            
            difference = solution;
            difference -= previous_solution;
            previous_solution=solution;
            
            std::cout << "   Iteration Number showing : " << iteration_number << "     Difference Norm : " << difference.l2_norm() << std::endl << std::flush;
            
        }while (difference.l2_norm()> 5* pow(10,-5)* dof_handler.n_dofs());  //defines iteration tolerance

        output_results (1);
        
        // triangulation.refine_global(2);
        // output_results (1);
        
    
        /*
        for (unsigned int refinement_cycle = 0; refinement_cycle<max_cycle;
             ++refinement_cycle)
        {
        }
         */
        
    }
    
    
    
    template <int dim>
    void StokesProblem<dim>::measure_velocity_profile_in_y ()
    {
        
        std::ostringstream filename;
        filename << "velocity-profile.txt";
        std::ofstream output (filename.str().c_str());
       
    
        const MappingQ<dim> mapping (degree);
        QGauss<dim-1>   quadrature_formula_face(2*degree+1);
        
        FEFaceValues<dim> fe_face_values (mapping, fe, quadrature_formula_face,
                                          update_JxW_values |
                                          update_quadrature_points |
                                          update_gradients |
                                          update_values |
                                          update_normal_vectors);
        
        const FEValuesExtractors::Scalar xvel (0);
        
        const unsigned int   faces_per_cell  = GeometryInfo<dim>::faces_per_cell;
        const unsigned int   n_q_face_points = fe_face_values.n_quadrature_points;
        
        std::vector<double>  local_ux_values (n_q_face_points);
        
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
        {
            for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
            {
                if (cell->face(face_no)->boundary_id()==2) //outlet
                {
                    fe_face_values.reinit (cell, face_no);
                    fe_face_values[xvel].get_function_values (solution, local_ux_values);
                    
                    for (unsigned int q=0; q<n_q_face_points; ++q)
                    {
                       double x_point = fe_face_values.quadrature_point (q)[0];
                       double y_point = fe_face_values.quadrature_point (q)[1];
                       double u_x = local_ux_values[q];
                       //std::cout << "(" << x_point<<","<<y_point<<"), u_x =" << u_x <<std::endl;
                       output << u_x << " " << y_point << std::endl;
                    }
                }
            }
        }
            
        
         output.close();
        
        
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

