program deformPlate
!
!  This is the main driver for the molecular dynamics program.
!
use Plate
use Continuum_Stress
implicit none

double precision,parameter :: PI = 3.141592653589793D0
integer, parameter :: DIM=3         ! change DIM to switch between 2D and 3D !
integer :: N=0
integer, dimension(:), allocatable :: atom_ids
double precision lattice_param
double precision diameter
double precision, dimension(DIM) :: pos, new_pos,Mass_Center
double precision, dimension(:,:), allocatable    :: positions
character*100     :: modelName        !  Name of file containing input sample
character*150     :: SampIn        !  Name of file containing input sample
character*150     :: SampOut       !  Name of file to contain output sample
character*80			  ::   garbage_char
integer :: i,lis
integer :: x_periodic, y_periodic, z_periodic
double precision box_strain 
double precision C11,C12,C44,garLambda,garAlpha,garMu

double precision          ::      xo,yo,zo,garbage
! Read E, nu, G, rhole, psi
! system size 
read(*,*,end=200,err=300) xo !Lx 
read(*,*,end=200,err=300) yo !Ly 
read(*,*,end=200,err=300) zo !Lz 
read(*,*,end=200,err=300) lattice_param !Lattice Constant
read(*,*,end=200,err=300) diameter ! size of hole
! correct length scale from [times] to [Angstrom]

xo = xo * lattice_param 
yo = yo * lattice_param
zo = zo * lattice_param 
diameter = diameter * lattice_param

rhole = diameter * 0.5

print*, "Since the last update, the current code take C11, C12, and C44, instead of E,nu,G"
print*, "make sure your current inputs are in the above form..."



read(*,*,end=200,err=300) C11
read(*,*,end=200,err=300) C12
read(*,*,end=200,err=300) C44
read(*,*,end=200,err=300) force

! angle of stress is restricted to 0 
! read(*,*,end=200,err=300) psi ! do not read this line 

psi = 0 
psi = psi * PI/180.d0                   ! Convert to radians


!periodic condition is restricted to 0,0,1
!read(*,*,end=200,err=300) x_periodic, y_periodic, z_periodic ! do not read this line
x_periodic = 0
y_periodic = 0 
z_periodic = 1

print*, 'periodicity from input ', x_periodic, y_periodic, z_periodic
!!!!
!!!!
!!!!
!!!

garLambda=C12 
garMu=(C11-C12)*0.5
garAlpha=garMu/garLambda

! not going to use it ! 
nu = 1.0/(2.0 * (garAlpha +1.0))
E = 2*garMu*(1+nu)
G = C44

! Read modelName 


read(*,'(a)',end=200,err=300) modelName 

SampIn  = trim(modelName) // "_undeformed.xyz"
SampOut = trim(modelName) // "_forlammps_deformed.txt"


! You should define sample in and sample out 
! (Caution) Must be set same as that of the python code !1
! read(*,'(a)',end=200,err=300) SampIn                       ! "
! read(*,'(a)',end=200,err=300) SampOut                      ! "

print*, 'rhole in Angstrom: ', rhole
print*, 'E in bars: ',E
print*, 'nu: ',nu 
print*, 'G in bars: ',G 
print*, 'force in eV/Angstrom: ', force
print*, 'psi in degree: ', psi
print*, 'inputfilename: ', SampIn
print*, 'outputfilename: ', SampOut

call Prep_constants

! open the input file
lis = len_trim(SampIn)
open(unit=1,file=SampIn(1:lis),status='old', &
     action='read',err=700)   
! open the output file
lis = len_trim(SampOut)
open(unit=2,file=SampOut(1:lis),status='unknown', &
     action='write',err=600)   

	 ! unit 1 - input file 
	 ! unit 2 - output file 
	 ! end and err

! First, read lammps dump file  
! Read the number of atoms
read(1, *) N
! Skip line "Atoms in FCC lattice"
read(1, *)

print*, "N:",N

allocate(positions(DIM,N))
allocate(atom_ids(N))

do i=1,N
   read(1,*,end=800,err=900) garbage_char, atom_ids(i), pos
   positions(:,i) = pos
enddo

! Check Data Read
!do i=1,N
!	print*, atom_ids(i), positions(:,i)
!enddo

!  compute center of mass coordinates
!
Mass_Center = sum( positions , dim=2 ) / N


print*, "Mass_Center", Mass_Center

! Start writing output config file that can be read from 

box_strain = (1-2*nu)*(1+nu)/(1-nu)*force/E ! strain solution without without a hole

print*, "box strain:",box_strain


!! Begins File writing 
!! Working from here!!
write(2,*,err=1000) "LAMMPS Description" 
write(2,*) 
write(2,*,err=1000) N, "atoms"
write(2,*,err=1000) 2, "atom types"
write(2,*) 
! write box reference and final 

write(2,*,err=1000) -1000.0, 1000.0, "xlo ", "xhi"
write(2,*,err=1000) -1000.0, 1000.0, "ylo ", "yhi"
write(2,*,err=1000) 0.0, zo, "zlo ", "zhi"
write(2,*) 


write(2,*) "Masses"
write(2,*)
write(2,*) 1, 28.0855
write(2,*) 2, 28.0855
write(2,*)

write(2,*) "Atoms"
write(2,*) 
!!!!

!! Just add atom id
do i=1,N
   pos = positions(:,i) - Mass_Center
   call Anal_Sol(pos(1),pos(2))
   new_pos(1) = positions(1,i)+u
   new_pos(2) = positions(2,i)+v
   new_pos(3) = positions(3,i)
   ! print*, 'sigX', sigX

   write(2, *) i, atom_ids(i), new_pos(1), new_pos(2), new_pos(3)
enddo


stop 

close(unit=1)
close(unit=2)
return
!
200 continue
   print*,'Read_Input: FATAL: premature end-of-file in standard input'
   stop
300 continue
   print*,'Read_Input: FATAL: read error in standard input'
   stop
600 continue
   print*,'Read_Sample: FATAL: ',SampOut(1:lis),' not found.'
   stop
700 continue
   print*,'Read_Sample: FATAL: ',SampIn(1:lis),' not found.'
   stop
800 continue
   print*,'Read_Sample: FATAL: premature end-of-file at atom ',i
   close(unit=1)
   stop   
900 continue
   print*,'Read_Sample: FATAL: read error in ',SampIn(1:lis)
   close(unit=1)
   stop
1000 continue
   print*,'Write_Sample: WARNING: error when writing ',SampOut(1:lis)
   close(unit=2)
   return

end program deformPlate



