#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
import sympy as sym

#Original scale where parameters have defined
L0=100.

alpha20=0.001
beta20= 0.0002
s0=0.001
fe0=0.0005
mu0=(1e-6)**.5
b_phi0=1.

# Grid scale
L=1
scale=L/L0 # scale =0.01

# Parameter scaling
alpha2=alpha20 * scale
beta2= beta20 * scale

fe=fe0/scale
s=s0
b_phi=b_phi0/scale

mu=mu0/scale

#Mob1 = Mob_mu/beta20**.5
mvmax = 1./b_phi
mvmin = 1.e-5*mvmax

def b_theta(phi):
    return 1./(mvmin+(1.-phi**3*(10.-15.*phi+6.*phi**2))*(mvmax-mvmin))

def ffun(phi):
    return fe*(phi-1)*(phi-1)

# e=2.1 * 10^(-5) for L=100

def gfun(phi):
    #return -2.*(ln(1.02-phi)+phi)
    return -1.*ln(1.02-phi)


## Cruzial: An smooth definition of p(grad) and p'(grad)

def p(gradu):
    return sqrt(inner(gradu,gradu)+mu*mu)
#
def dp(gradu):
    return gradu/p(gradu)

## Periodic
class PeriodicBoundary(SubDomain):
  
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
    
    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1] - 1.0 #top and bottom


pbc = PeriodicBoundary()


mesh = UnitSquareMesh.create( 80, 80, CellType.Type.quadrilateral)
x = mesh.coordinates()
scaling_factor = L
x[:, :] *= scaling_factor

element_esc = FiniteElement('Lagrange',mesh.ufl_cell(),1) #element order 
element_Z = FiniteElement('Lagrange', mesh.ufl_cell(),1)
element_mixed = MixedElement([element_esc,element_Z])


#V = FunctionSpace(mesh, element_mixed)
V = FunctionSpace(mesh, element_mixed, constrained_domain=pbc)

# L=1
f2 = open('crystal_unitscale.txt','r')
lines2=f2.readlines()
j=0
points=[]
angle=[]
for i in lines2:
    j=+1
    spl=i.split()
    points.append([L*float(spl[0]),L*float(spl[1]) ])
    angle.append(float(spl[3]))
angle=np.array(angle)*2.
f2.close()
points=np.array(points)


class InitialCondition(UserExpression):
    def eval_cell(self, value, x, ufc_cell):
        dcamin=4*L
        i=0
        crystal=0
        for point in points:            
            dca=np.linalg.norm(x-point) 
   
            if(dca<=dcamin):
                dcamin=dca
                crystal=i   
            i+=1

        value[0]=1.
        value[1]=angle[crystal]

    def value_shape(self):
        return (2,)
    
class boundary_grains(UserExpression):
    def eval_cell(self, value, x, ufc_cell):
        dcamin=4*L
        i=0
        crystal=0
        for point in points:            
            dca=np.linalg.norm(x-point) 
   
            if(dca<=dcamin):
                dcamin=dca
                crystal=i   
            i+=1
        value[0]=angle[crystal]


expression_GB=boundary_grains(degree=2)
expression_un=InitialCondition(degree=2)

element_esc = FiniteElement('Lagrange',mesh.ufl_cell(),1)
element_Z = FiniteElement('Lagrange', mesh.ufl_cell(),1)
element_mixed = MixedElement([element_esc,element_Z])

V = FunctionSpace(mesh, element_mixed)
q_degree = 6 
dx = dx(metadata={'quadrature_degree': q_degree})

def boundary(x, on_boundary):
    return on_boundary

# bound_phi = Constant('1')
# bc_phi = DirichletBC(V.sub(0), bound_phi, boundary)  # BC for initial field
# bc_theta = DirichletBC(V.sub(1), expression_GB, boundary)

bc = []


du = TestFunction(V)
dphi,dtheta = split(du)
u = Function(V) #solution field in t+dt
u_n = Function(V) #solution field in t
u_n2 = Function(V)
utrial = TrialFunction(V)
phi,theta= split(u) # Labels for the two coupled fields

u_n = interpolate(expression_un,V)
phi_n,theta_n = split(u_n)


# Define variational problem
# Note that in F might enter python functions of u and expressions of x

num_plot=500 # number of frames saved

# This is to force a given integration rule, here the exact one is too large
#q_degree = 6
#dx = dx(metadata={'quadrature_degree': q_degree})

# Variational non-linear form, find u such F(u,w)=0
# Note that implicit backward Euler is used for time integration
t=0.
dt=1e-3*b_theta(1)  # b_theta - mvmax (10^5) when L0=100
tend=2*b_phi*1e6*scale


Energy_functional = .5*alpha2*dot(grad(phi),grad(phi))*dx+\
                    ffun(phi)*dx+s*gfun(phi)*p(grad(theta))*dx+\
                    .5*beta2*( dot(grad(theta),grad(theta) ) )*dx

def vari_shape(dt):
    return derivative(Energy_functional, u, du) \
      +(b_phi/dt)*dphi*phi*dx -(b_phi/dt)*dphi*phi_n*dx \
      +(b_theta(phi_n)/dt)*dot(dtheta,theta)*dx-(b_theta(phi_n)/dt)*dot(dtheta,theta_n)*dx
  

      
F = vari_shape(dt)    
J  = derivative(F, u, utrial) 
problem = NonlinearVariationalProblem(F, u, bc, J)
solver  = NonlinearVariationalSolver(problem)

parameters['linear_algebra_backend'] = 'PETSc'
prm = solver.parameters

#prm['newton_solver']['absolute_tolerance'] = 1E-8
#prm['newton_solver']['relative_tolerance'] = 1E-7
#prm['newton_solver']['maximum_iterations'] = 15
#prm['newton_solver']['relaxation_parameter'] = 1.0
#:q
prm["newton_solver"]["linear_solver"] = "gmres"
#prm["newton_solver"]["preconditioner"] = "sor"
prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-14
prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-10
prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 300
prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = True
prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
prm['newton_solver']['absolute_tolerance'] = 1E-9
prm['newton_solver']['relative_tolerance'] = 5E-7
prm['newton_solver']['maximum_iterations'] = 10
prm['newton_solver']['relaxation_parameter'] = 1. 

# Create VTK files for visualization output
vtkfile_phi= File('2D_KWC/phi.pvd')
vtkfile_theta= File('2D_KWC/theta.pvd')

u.assign(u_n)

phiv,thetav =u.split()
vtkfile_phi << (phiv,t)
vtkfile_theta << (thetav,t)

Energy_evol=[]
tlist=[]
tlist.append(0.)
Energy_evol.append(assemble(Energy_functional))

inc=0
itertot=0
cut=0
while t<tend :
    t+=dt
    a=False    
    while(a == False):          
        try:            
            a=solver.solve()            
        except:
             cut+=1
             dt=dt/2.
             print('No convergency, decrease dt',dt)
             F = vari_shape(dt)    
             J  = derivative(F, u, utrial) 
             problem = NonlinearVariationalProblem(F, u, bc, J)
             solver  = NonlinearVariationalSolver(problem)
          
             prm = solver.parameters
#             prm['newton_solver']['absolute_tolerance'] = 1E-8
#             prm['newton_solver']['relative_tolerance'] = 1E-7
#             prm['newton_solver']['maximum_iterations'] = 15
#             prm['newton_solver']['relaxation_parameter'] = 1.0
#             #
             prm["newton_solver"]["linear_solver"] = "gmres"
             #prm["newton_solver"]["preconditioner"] = "sor"
             prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-14
             prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-10
             prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 300
             prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = True
             prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
             prm['newton_solver']['absolute_tolerance'] = 1E-9
             prm['newton_solver']['relative_tolerance'] = 5E-7
             prm['newton_solver']['maximum_iterations'] = 10
             prm['newton_solver']['relaxation_parameter'] = 1. 

               
             u.assign(u_n)
             a=False  
    inc+=1  
    itertot+=a[0]
#    u_n2.assign(u_n)
    u_n.assign(u)
    dt_old=dt
    tlist.append(t)
    Energy_evol.append(assemble(Energy_functional))
#    
    
    print('Converged with n=',a[0],t,t/tend,inc,itertot,dt)
    
    if( int(num_plot*(t/tend)) != int(num_plot*((t-dt)/tend))): # to store solution only
        print('writting plot file')
        phiv,thetav =u.split()
        vtkfile_phi << (phiv,t)
        vtkfile_theta << (thetav,t) 
    
    if(a[0]<=3):          
        print('Increase dt')
        dt=dt*2.
        F = vari_shape(dt)    
        J  = derivative(F, u, utrial) 
        problem = NonlinearVariationalProblem(F, u, bc, J)
        solver  = NonlinearVariationalSolver(problem)
       
        prm = solver.parameters
#        prm['newton_solver']['absolute_tolerance'] = 1E-8
#        prm['newton_solver']['relative_tolerance'] = 1E-7
#        prm['newton_solver']['maximum_iterations'] = 15
#        prm['newton_solver']['relaxation_parameter'] = 1.0
        #
        prm["newton_solver"]["linear_solver"] = "gmres"
        #prm["newton_solver"]["preconditioner"] = "sor"
        prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-14
        prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-10
        prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 300
        prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = True
        prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
        prm['newton_solver']['absolute_tolerance'] = 1E-9
        prm['newton_solver']['relative_tolerance'] = 5E-7
        prm['newton_solver']['maximum_iterations'] = 10
        prm['newton_solver']['relaxation_parameter'] = 1.


