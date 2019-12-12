from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
import math

#I need to fix scaling: does it convergent to all orders of parameters?

#Constants
alpha2=0.00530
epsilon2= 0.00021333
epsilon = math.sqrt(epsilon2)
s=0.0017
e=0.0021
gamma=250
b_phi=1E3
b_theta=1E3
L=10.
mu=1E-5

# f and g functions in KWC model

def b_theta(phi):
    #return 1./(mvmin+(1.-phi**3*(10.-15.*phi+6.*phi**2))*(mvmax-mvmin))
    return 1

def ffun(phi):
    return 0.5*(1/epsilon2)*(phi-1)*(phi-1)

def gfun(phi):
    return -2.*(ln(1.02-phi)+phi)


def p(gradu):
    return sqrt(dot(gradu,gradu)+mu*mu)

def dp(gradu):
    return gradu/p(gradu)



mesh = UnitSquareMesh.create( 80, 80, CellType.Type.quadrilateral)

element_esc = FiniteElement('Lagrange',mesh.ufl_cell(),1) #element order
element_Z = FiniteElement('Lagrange', mesh.ufl_cell(),1) #linear element
element_mixed = MixedElement([element_esc,element_Z])

V = FunctionSpace(mesh, element_mixed)
q_degree = 6
dx = dx(metadata={'quadrature_degree': q_degree})



# Define Initial Condition boundary condition
theta1=0.
theta2=np.pi/6.
theta3=np.pi/3.

class InitialCondition(UserExpression):
    def eval(self, value, x):

        value[0]=1 #bc for phi
        
        if (x[0] <0.2):
            temp1= (1.0 - 0.5)/(0.0-0.2) * (x[0]-0.2) + 0.5
            temp2= -(1.0 - 0.5)/(0.0-0.2) * (x[0]-0.2) + 0.5
            
            if (x[1]>temp1):
                value[1]=theta3
            
            elif(x[1]<temp1 and x[1]>temp2):
                value[1]=theta2
            
            else:
                value[1]=theta1

        else:
            if (x[1]>0.5):
                value[1]=theta3
            else:
                value[1]=theta1

    def value_shape(self):
        return (2,)


def boundary(x, on_boundary):
    return on_boundary


class boundary_grains(UserExpression):
    def eval(self, value, x):
        
        if (x[0] <0.2):
        
            temp1= (1.0 - 0.5)/(0.0-0.2) * (x[0]-0.2) + 0.5
            temp2= -(1.0 - 0.5)/(0.0-0.2) * (x[0]-0.2) + 0.5
        
            if (x[1]>temp1):
                value[0]=theta3
            
            elif(x[1]<temp1 and x[1]>temp2):
                value[0]=theta2
        
            else:
                value[0]=theta1
        
        else:
            if (x[1]>0.5):
                value[0]=theta3
            else:
                value[0]=theta1

            value[0]=value[0]

    def value_shape(self):
        return ()  #return scalar



expression_un=InitialCondition(degree=2)
expression_GB=boundary_grains(degree=2)

bound_phi = Constant('1')
bc_phi = DirichletBC(V.sub(0), bound_phi, boundary)  # BC for initial field
bc_theta = DirichletBC(V.sub(1), expression_GB, boundary)

bc = [bc_theta] ## It seems that only bc of theta is implemented?

du = TestFunction(V)
dphi,dtheta = split(du)
u = Function(V) #solution field in t+dt
u_n = Function(V) #solution field in t
u_n2 = Function(V)
utrial = TrialFunction(V)
phi,theta= split(u) # Labels for the two coupled fields

u_n = interpolate(expression_un,V)
phi_n,theta_n = split(u_n)


num_plot=2 # number of frames saved


# Variational non-linear form, find u such F(u,w)=0
# Note that implicit backward Euler is used for time integration
t=0.
dt=1e-3
tend=dt * 3

# Create progress bar, does not work with my fenics version
# progress = Progress('Time-stepping')
# set_log_level(PROGRESS)


Energy_functional = .5*alpha2*dot(grad(phi),grad(phi))*dx+\
    ffun(phi)*dx+s*gfun(phi)*p(grad(theta))*dx+\
    .5*epsilon2*( dot(grad(theta),grad(theta) ) )*dx

def vari_shape(dt):
    return derivative(Energy_functional, u, du) \
        +(b_phi/dt)*dphi*phi*dx -(b_phi/dt)*dphi*phi_n*dx \
        +(b_theta(phi_n)/dt)*dot(dtheta,theta)*dx \
        -(b_theta(phi_n)/dt)*dot(dtheta,theta_n)*dx


F = vari_shape(dt)
J  = derivative(F, u, utrial)
problem = NonlinearVariationalProblem(F, u, bc, J)
solver  = NonlinearVariationalSolver(problem)

parameters['linear_algebra_backend'] = 'PETSc'
prm = solver.parameters

prm["newton_solver"]["linear_solver"] = "gmres"
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
            prm["newton_solver"]["linear_solver"] = "gmres"
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
        prm["newton_solver"]["linear_solver"] = "gmres"
        prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-14
        prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-10
        prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 300
        prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = True
        prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True
        prm['newton_solver']['absolute_tolerance'] = 1E-9
        prm['newton_solver']['relative_tolerance'] = 5E-7
        prm['newton_solver']['maximum_iterations'] = 10
        prm['newton_solver']['relaxation_parameter'] = 1.


    plt.plot(tlist,Energy_evol)

    print(inc)
    print(itertot)
    print(cut)



print("Progam finished successfully")




## Time intrementation
'''
print('first step')
inc=0
dt=.001*b_phi
t=0
#tend=1E6*dt
tend=1
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


'''



'''
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
+ epsilon2*dot( grad(w2),grad(theta) )*dx
#\ + s*g(phi)*dot( dp(grad(theta)) , grad(w2) )*dx \
+ -(b_theta/dt)*w2*theta_n*dx

##Full implicit

# Create VTK files for visualization output
vtkfile_phi= File('2D_KWC/phi.pvd')
vtkfile_theta= File('2D_KWC/theta.pvd')






'''
'''
    phi=sym.symbols('phi')
    f=e*(phi-1)**2 # Symbolic f
    fprim=sym.diff(f,phi) # Symbolic df/dphi
    f_code=sym.printing.ccode(f) # Expression f
    fprim=sym.lambdify(phi,fprim) # Python function df/dphi
    
    g=phi**2 # Symbolic g
    
    #g= -ln(1-phi)
    gprim=sym.diff(g,phi) # Symbolic dg/dphi
    g_code=sym.printing.ccode(g) # Expression g
    gprim = sym.lambdify(phi,gprim)# Python function dg/dphi
    g=sym.lambdify(phi,g)
    
'''
