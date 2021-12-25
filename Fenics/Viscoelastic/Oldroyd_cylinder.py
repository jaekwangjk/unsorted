from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from mshr import *

'''
DE = lambda * V / a, or Wi
a=1
'''

# Vavg and mu0 is fixed 1 
# and if it only plays with lambda!
# Did any paper specified this?
mu0=1
beta =0.41
mu_P = mu0 * beta
mu_N = mu0 * (1-beta)
lamda= 0.5
Vavg = 1.0 
radius = 1

De = lamda * Vavg  / radius

# mesh generation

channel = Rectangle ( Point ( -12 , -2) ,Point (16 , 2))
cylinder = Circle ( Point (0, 0) , radius, 100) 
domain = channel - cylinder
mesh = generate_mesh (domain , 150) # it was 105 initially

# Refinement
zone_1 = MeshFunction ('bool', mesh , 2)
for i in cells ( mesh ):
    if (i. midpoint ().x() > -3 and i. midpoint ().x() < 3):
        zone_1 [i] = True
        
mesh=refine(mesh,zone_1)

# Construct facet markers
bndry = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
for f in facets(mesh):
    mp = f.midpoint()
    if near(mp[0], -12): # inflow
        bndry[f] = 1
    elif near(mp[0], 16): # outflow
        bndry[f] = 2
    elif near(mp[1], -2.0) or near(mp[1], 2.0): # walls
        bndry[f] = 3
    elif mp.distance(Point(0,0)) <= 1: # cylinder of radius 1 
        bndry[f] = 5
        
#List boundary condition
inlet_profile = Expression(("V* (3./2.0)*(1. - (x[1]*x[1])/4.0)","0"),V=Vavg,degree=2)
noslip = Expression(("0","0"),degree=2)


P = FiniteElement ("CG", triangle , 1)
V = VectorElement ("CG", triangle , 2)
Tau = TensorElement ("CG", triangle , 2)

E = MixedElement([P, V, Tau ])
W = FunctionSpace (mesh ,E)
w = Function(W)


inlet_velocity = DirichletBC (W.sub (1) ,inlet_profile, bndry,1 )
noslip_cylinder = DirichletBC (W.sub (1) ,noslip , bndry,5 )
noslip_walls = DirichletBC (W.sub (1) ,noslip , bndry,3)


bcs=[inlet_velocity, noslip_cylinder, noslip_walls]


# geometric dimension
I=Identity(2)

(p,v, tau ) = split (w)
(p_ ,v_ , tau_ ) = TestFunctions (W)


def L(v):
    return grad(v)

def D(v):
    return ( grad (v)+ grad (v).T) /2.0


## Stabilization term 
h = CellDiameter(mesh)

def eta(v):
    eta1 = 1.0* h / Vavg
    #eta1 = h / Vavg
    return eta1

T = -p*I + 2*mu_N* D(v) + tau

n = FacetNormal(mesh)
# link boundary information to 'ds' 
ds = Measure("ds", subdomain_data=bndry)


w10 = div(v)*p_*dx
w11 = inner ((-p*I +2.* mu_N *D(v)+tau),grad(v_))*dx - dot (dot(T,n),v_)*ds
w13 = inner (( tau+ lamda *( dot (v, nabla_grad( tau ))-dot (L(v),tau )-dot (tau ,L(v).T)) -2.* mu_P *D(v)),tau_\
+ dot ( eta(v)*v,nabla_grad ( tau_ )))*dx

weak_form = w10 + w11 + w13

derivative_of_weak_form = derivative (weak_form ,w)


print("start solving the variational problem.")
my_problem = NonlinearVariationalProblem(weak_form ,w,bcs ,derivative_of_weak_form )
my_solver = NonlinearVariationalSolver (my_problem)


parameters['linear_algebra_backend'] = 'PETSc'
prm = my_solver.parameters
prm['newton_solver']['absolute_tolerance'] = 6E-2
prm['newton_solver']['relative_tolerance'] = 3E-3
prm['newton_solver']['maximum_iterations'] = 100
prm['newton_solver']['relaxation_parameter'] = 1.

my_solver.solve()




#create output vtk
p_,u_,tau_ = w.split()

# Define stress tensor
def sigma(u, p,mu_N):
    return 2*mu_N*D(u) - p*Identity(len(u))

force = dot(-p_*I + 2.0*mu_N*sym(grad(u_)), n) + dot(tau_,n)

D=(force[0])* ds(5)
L=(force[1])* ds(5)

Drag = -assemble(D)
Dcof = Drag /(mu0*Vavg)

print("Wi=", De)
print("Drag=", Drag)
print("Dcof=", Dcof)

pfile_pvd = File("pressure_1224.pvd")
pfile_pvd << p_

ufile_pvd = File("velocity_1224.pvd")
ufile_pvd << u_

taufile_pvd = File("tau_1224.pvd")
taufile_pvd << tau_

