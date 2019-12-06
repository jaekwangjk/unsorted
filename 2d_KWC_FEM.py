
"""
Created on Wed Jan 24 21:34:49 2018

@author: javierseguradoescudero

@ Jaekwang's memo
It seems that the codes is dealing with shirinking boundaries

This does larger time stepping at later times... why?

"""


# 3D-KWC model implementation using quaternions
# See notes for the details of implementation

from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
import math

# Some orientation manipulation function
# From rotation to quaternion (to define initial BCs)

#Constants
alpha2=0.00530
epsilon2= 0.00021333
s=0.0017
e=0.0021
gamma=250
b_phi=1E3
b_theta=1E3
L=10.
mu=1E-5

# f and g functions in KWC model
phi=sym.symbols('phi')
f=e*(phi-1)**2 # Symbolic f
fprim=sym.diff(f,phi) # Symbolic df/dphi
f_code=sym.printing.ccode(f) # Expression f
fprim=sym.lambdify(phi,fprim) # Python function df/dphi

g=phi**2 # Symbolic g
gprim=sym.diff(g,phi) # Symbolic dg/dphi
g_code=sym.printing.ccode(g) # Expression g
gprim = sym.lambdify(phi,gprim)# Python function dg/dphi
g=sym.lambdify(phi,g)
#

## Cruzial: An smooth definition of p(grad) and p'(grad)


def p(gradu):
    return sqrt(dot(gradu,gradu)+mu)
#
def dp(gradu):
    return gradu/(sqrt(dot(gradu,gradu)+mu))


mesh = RectangleMesh(Point(-L,-L), Point(L, L), 20, 20, "right")

# Function Spaces
P2 = FiniteElement('CG',mesh.ufl_cell(), 2) 
element = MixedElement([P2, P2])
V = FunctionSpace(mesh, element) 

# Define boundary condition
theta1=0.
theta2=np.pi/12.


def boundary(x, on_boundary):
    return on_boundary
#smooth function of orientation change near center
bound = Expression(('1.','(x[0]/L)*(x[0]/L)+(x[1]/L)*(x[1]/L)<=0.25  ? theta2 : theta1'), degree=2,L=L,theta1=theta1,theta2=theta2)
# Define funcions in the FE space
w1,w2 = TestFunctions(V) 
u = Function(V) #solution field in t+dt
u_n = Function(V) #solution field in t

phi,theta= split(u) # Labels for the two coupled fields

u_n = interpolate(bound,V) # Initial field value u(t=0) 
phi_n,theta_n = split(u_n)
bc = DirichletBC(V, bound, boundary)  # BC for initial field

# Define variational problem
# Note that in F might enter python functions of u and expressions of x

dt=.001*b_phi

num_plot=100 # number of frames saved

# This is to force a given integration rule, here the exact one is too large
q_degree = 3
dx = dx(metadata={'quadrature_degree': q_degree})
#

# Variational non-linear form, find u such F(u,w)=0
# Note that implicit backward Euler is used for time integration
F=  (b_phi/dt)*w1*phi*dx  \
+ alpha2 * dot( grad(w1) , grad(phi) )*dx \
+ w1*fprim(phi)*dx \
+ s*w1*gprim(phi)*p(grad(theta))*dx \
+ -(b_phi/dt)*w1*phi_n*dx  \
+ (b_theta/dt)*w2*theta*dx \
+ epsilon2*dot( grad(w2),grad(theta) )*dx \
+ s*g(phi)*dot( dp(grad(theta)) , grad(w2) )*dx \
+ -(b_theta/dt)*w2*theta_n*dx

##Full implicit

# Create VTK files for visualization output
vtkfile_phi= File('2D_KWC/phi.pvd')
vtkfile_theta= File('2D_KWC/theta.pvd')

# Create progress bar
#progress = Progress('Time-stepping')
#set_log_level(PROGRESS)

## Time intrementation 

print('first step')
inc=0
dt=.001*b_phi
t=0
tend=1E6*dt
inc=0.
while t<tend:
    F=  (b_phi/dt)*w1*phi*dx  \
    + alpha2 * dot( grad(w1) , grad(phi) )*dx \
    + w1*fprim(phi)*dx \
    + s*w1*gprim(phi)*p(grad(theta))*dx \
    + -(b_phi/dt)*w1*phi_n*dx  \
    + (b_theta/dt)*w2*theta*dx \
    + epsilon2*dot( grad(w2),grad(theta) )*dx \
    + s*g(phi)*dot( dp(grad(theta)) , grad(w2) )*dx \
    + -(b_theta/dt)*w2*theta_n*dx
    inc=1
    while (inc/100.) != int(inc/100): ##Do this until inc becomes a hundreds
        inc+=1
        t+=dt
        solve(F==0,u,bc)
        u_n.assign(u)
        if( int(num_plot*(t/tend)) != int(num_plot*((t-dt)/tend))): # to store solution only 
            #print t
            phi_v , theta_v=u.split()
## Save solution to file (VTK)
            vtkfile_phi << (phi_v,t)
            vtkfile_theta << (theta_v, t)
    #progress.update(t / tend)
    dt=dt*10
    print('new dt',dt)


