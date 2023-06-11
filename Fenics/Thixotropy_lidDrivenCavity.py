#
#
# Non dimensionallized...
# Expand function space...
# Make it time... ...
#
# 
# Impulsively started thixotropy... time dependent flow...

from fenics import *

#hyper parameters 
rho_0 = 1
mu = 0.01
u_x = 1

# ##### Mesh parameters ##### #
mesh = UnitSquareMesh(40, 40, diagonal="crossed")

# ##### Define function spaces ##### #
# Talyer-Hood element
V_ele = VectorElement("CG", mesh.ufl_cell(), 2)
Q_ele = FiniteElement("CG", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([V_ele, Q_ele, Q_ele]))


# ##### Boundary conditions ##### #
lid_location = "near(x[1],  1.)"
fixed_wall_locations = "near(x[0], 0.) | near(x[0], 1.) | near(x[1], 0.)"
noslip = DirichletBC(W.sub(0), (0, 0), fixed_wall_locations)

u_x = 1
velocity_x = Constant((u_x, 0.0))
top_velocity = DirichletBC(W.sub(0), velocity_x, lid_location)
bcu = [top_velocity, noslip]

# ##### Define variational parameters ##### #
u_, p_ , xi_ = TestFunctions(W)
upxi = Function(W)
u, p, xi = split(upxi)

f = Constant((0, 0)) # force 

#Initial Newtonian Solution 

F = (
     rho_0*inner(grad(u)*u, u_)*dx + mu*inner(grad(u), grad(u_))*dx
     - inner(p, div(u_))*dx + inner(p_, div(u))*dx
     + inner(f, u_)*dx
)

F += inner(xi,xi_)* dx 

J = derivative(F, upxi , TrialFunction(W)) # derivative of F w.r.t up in W space 

problem = NonlinearVariationalProblem(F, upxi, bcu, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0

solver.solve()

(u,p,xi) = upxi.split(True);


ofile_u = File("thixo_u.pvd")
ofile_u << u

'''
# Time step
t=0
dt=5
t_end=10




# execute the time steps....

while t < t_end:
    print("\n *** Time: %.2f / %.2f *** Time step: *** \n"%((t),(t_end)) )
    
    F = (
         inner( (u-u0)/dt, u_) *dx  # time dependent term 
         +rho_0*inner(grad(u)*u, u_)*dx + mu*inner(grad(u), grad(u_))*dx
         - inner(p, div(u_))*dx + inner(p_, div(u))*dx
         + inner(f, u_)*dx
         )     
    J = derivative(F, upxi , TrialFunction(W)) # derivative of F w.r.t up in W space 
    
    t=t+dt
    
    problem=NonlinearVariationalProblem(F,upxi,bcs,J)
    solver=NonlinearVariationalSolver(problem)
    
    #prm = solver.parameters
    #prm['newton_solver']['absolute_tolerance'] = 1E-8
    #prm['newton_solver']['relative_tolerance'] = 1E-7
    #prm['newton_solver']['maximum_iterations'] = 25
    #prm['newton_solver']['relaxation_parameter'] = 1.0
    
    solver.solve()
    
    (u,p,xi) = upxi.split(True);
    # assign previous
    u0=assign(u)
    p0=assign(p)
    xi0=assign(xi)
    
    
    # todo...
    # subs u0,p0,xi0


'''

#check solution...
#ofile_xi = File("thixo_xi.pvd")
#ofile_xi << xi
#ofile_p = File("thixo_p.pvd")
#ofile_p << p


# Initialize structures 

#xi_init = Constant(0.0)


#xi = (interpolate(xi_init, Q_ele))
#


#xi.assign(interpolate(xi_init, Q_ele))

#ci.assign(interpolate(cinit,L)) # L is the Lagrange Function Space 




print("Code is sucessfully finished...")
