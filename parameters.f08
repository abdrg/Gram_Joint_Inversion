program parameters
! This program setup the parameters for the joint inversion of gravity and magnetic data.
! Author: Abraham Del Razo, IPICyT. Last update: August 2023
! Mail: abraham.delrazo@ipicyt.edu.mx
!***********************************************************************************************************************************

use iso_fortran_env , dp => real64
!use gram_joint_inversion 
use gram_joint_inversion_host

implicit none

real, allocatable, dimension(:) :: d_gv, stdevd_gv, m_gv, mapr_gv, stdevm_gv
real, allocatable, dimension(:) :: d_mg, stdevd_mg, m_mg, mapr_mg, stdevm_mg
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real(dp), allocatable :: A_gv(:,:), A_mg(:,:)
real :: step
integer :: nx, ny, nz, m, n, num_iters, TYPE
real :: H, Fi, Fd, err0
real :: reg(10), dip, strike

real :: start0, finish0

!***********************************************************************************************************************************
! INVERSION PARAMETERS

! Magnetic Field
H = 40000 !nT
Fi = 45 !degree
Fd = 45 !degree

!*************************
! REGULARIZATION PARAMETERS 
reg(1) = 1e2 !alpha_gv /smoothness
reg(2) = 1e2 !alpha_mg
reg(3) = 2e0 !beta_gv /a-priori information
reg(4) = 9e-1 !beta_mg
reg(5) = 1e5 !gamma_gv /dip and strike
reg(6) = 1e5 !gamma_mg
reg(7) = 0 !eta_gv /verticalily
reg(8) = 0 !eta_mg
reg(9) = 1e3 !mu_gv /structural coupling
reg(10) = 1e3 !mu_mg
dip = -45 !degree / only apply with gamma different if zero  !Positive from east to north.
strike = 0 !degree  !Positive from north to west.
!*************************

! STOP CONDITIONS
num_iters = 30
err0 = 0.1 !% of model total variation for each iteration
!*************************

!***********************************************************************************************************************************
! READ DATA

! ANOMALIES 
call readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
m = size(d_gv)
print*, "Total number of data: ", size(d_gv)

! overwrite default values of DATA STANDARD DEVIATION
stdevd_gv = 0.02*( maxval(d_gv)-minval(d_gv) ) !2% flexibility to fit the data
stdevd_mg = 0.02*( maxval(d_mg)-minval(d_mg) )

! MODEL DIMENSIONS
call readmodel(xcell,ycell,zcell,m_gv,m_mg)
n = size(m_gv)
print*, "Total number of parameters: ", n

! overwrite default values of INITIAL MODEL 
m_gv = 0 !Homogeneous distribution of initial parameters
m_mg = 0

! MODEL STANDARD DEVIATION
nx = 20 ; ny = 20 ; nz = 10 ; step = 50 !meters
allocate( stdevm_gv(n), stdevm_mg(n) )
call model1D(nx,ny,nz,1,stdevm_gv)
call model1D(nx,ny,nz,1,stdevm_mg)

! A-PRIORI MODEL 
allocate( mapr_gv(n), mapr_mg(n) )
mapr_gv = 0
mapr_mg = 0
call model1D(nx,ny,nz,2,mapr_gv)
!call model1D(nx,ny,nz,2,mapr_mg)

! READ LINEAR OPERATOR MATRICES (binary format)
allocate( A_gv(m,n), A_mg(m,n) )
OPEN(UNIT=1, FILE="input_a_gv.dat", ACTION="read", FORM="unformatted")
READ(1) A_gv
CLOSE(UNIT=1)

OPEN(UNIT=2, FILE="input_a_mg.dat", ACTION="read", FORM="unformatted")
READ(2) A_mg
CLOSE(UNIT=2)
!***********************************************************************************************************************************


!***********************************************************************************************************************************
! MODELING BY JOINT INVERSION
TYPE = 1

!Code version of inversion is defined for variable TYPE
! = 1: CPU version.
! = 2: GPU version. /Only available using gram_joint_inversion.cuf module
! = 3: CPU/GPU version. /Only available using gram_joint_inversion.cuf module
! = 0: no inversion, only input data read test.
!***********************************************************************************************************************************
call cpu_time(start0)

if (TYPE == 1) then
      call jointGramCPU(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                      d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TYPE == 2) then
      !call jointGramGPU(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
      !                  d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TYPE == 3) then
      !call jointGramCPUGPU(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
      !                  d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TYPE == 0) then
      PRINT*, "Test input data"
      print *, "d_gv", d_gv(1:5)
      print *, "m_gv", m_gv(1:5)
      print *, "A_gv:", A_gv(1:5,1:5)

else 
      PRINT*, 'The inversion type was wrongly selected, please choose an appropriate TYPE'
end if

call cpu_time(finish0)

!***********************************************************************************************************************************
! Archivo salida de parametros
OPEN(unit=13,file='output_parameters_inv.txt',status='unknown')
write(13,*) 'Number of data:', m
write(13,*) 'Volume size/Number of cells:', n
write(13,*) 'nx =',nx
write(13,*) 'nx =',ny
write(13,*) 'nz =',nz
write(13,*) 'Cell size =',step,'mts'
write(13,*) ' '
write(13,*) 'Magnetic Field'
write(13,*) 'H =',H,', Fi =',Fi,', Fd =',Fd
write(13,*) ' '
write(13,*) 'Regularization Parameters'
write(13,*) 'alpha_gv =', reg(1)
write(13,*) 'alpha_mg =', reg(2)
write(13,*) 'beta_gv =', reg(3)
write(13,*) 'beta_mg =', reg(4)
write(13,*) 'gamma_gv =', reg(5)
write(13,*) 'gamma_mg =', reg(6)
write(13,*) 'eta_gv =', reg(7)
write(13,*) 'eta_mg =', reg(8)
write(13,*) 'mu_gv =', reg(9)
write(13,*) 'mu_mg =', reg(10)
write(13,*) 'dip =', dip
write(13,*) 'strike =', strike
write(13,*) ' '
write(13,*) 'StdDeviation gv data = 2% rango de anomalia'
write(13,*) 'StdDeviation mg data = 2% rango de anomalia'
write(13,*) 'StdDeviation gv model = 1 todo modelo y 0.3 para cima del dominio gr/cm3 *MANUALLY setup in model1D subroutine'
write(13,*) 'StdDeviation mg model = 1 todo modelo y 0.3 para cima del dominio A/m  *MANUALLY setup in model1D subroutine'
write(13,*) ' '
write(13,*) 'Max Iterations Joint Inversion =',num_iters
write(13,*) 'Stop Condition, % Convergence less than', err0
write(13,*) ' '
write(13,*) 'Total execution time for TYPE:', TYPE
write(13,*) (finish0-start0),'seg =',(finish0-start0)/60, 'min'
close(unit=13)
!***********************************************************************************************************************************

end program parameters
