program parameters
! Programa que utiliza una serie de subrutinas relacionados con modelado e inversion de métodos potenciales.
! Autor: Abraham Del Razo, IPICyT. Inicio de desarrollo: May-2020. Ultima actualización: Nov-2022.
! Mail: abraham.delrazo@ipicyt.edu.mx
!***********************************************************************************************************************************

use iso_fortran_env , dp => real64
use gram_joint_inversion 

implicit none

real, allocatable, dimension(:) :: d_gv, stdevd_gv, m_gv, mapr_gv, stdevm_gv
real, allocatable, dimension(:) :: d_mg, stdevd_mg, m_mg, mapr_mg, stdevm_mg
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real(dp), allocatable :: A_gv(:,:), A_mg(:,:)
real :: step
integer :: nx, ny, nz, m, n, num_iters, TIPO_INVERSION
real :: H, Fi, Fd, err0
real :: alpha_gv, alpha_mg, beta_gv, beta_mg, gamma_gv, gamma_mg, eta_gv, eta_mg, mu_gv, mu_mg, reg(10), dip, strike

real :: start0, finish0

!***********************************************************************************************************************************
! PARAMETROS DE INVERSION

! Campo magnetico
H = 40000 !nT
Fi = 45 !grados
Fd = 45 !grados

!*************************
! PARAMETROS DE REGULARIZACION
alpha_gv = 1e1  !curvatura
alpha_mg = 1e1
beta_gv = 5e-2!1e-1  !apriori
beta_mg = 5e-2!1e-1
gamma_gv = 5e3  !dip and strike
gamma_mg = 5e3
eta_gv = 0!-1e3  !verticalidad
eta_mg = 0!-1e3
mu_gv = 1  !acoplamiento estructural
mu_mg = 1
dip = -45 !grados !solo aplica con eta diferente de cero.
strike = 0 !grados	

!Arreglo de parametros de regularizacion 
reg(1) = alpha_gv
reg(2) = alpha_mg
reg(3) = beta_gv
reg(4) = beta_mg
reg(5) = gamma_gv
reg(6) = gamma_mg
reg(7) = eta_gv
reg(8) = eta_mg
reg(9) = mu_gv
reg(10) = mu_mg
!*************************

! CONDICIONES DE PARO
num_iters = 30
err0 = 0.1 !% variacion total de modelos en cada iteracion
!*************************

!***********************************************************************************************************************************
! LECTURA DE DATOS

! ANOMALIAS
call readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)

m = size(d_gv) !numero de observaciones
print*, "Total number of data: ", size(d_gv)

! sobreescribir valores default DESVIACION STD DE DATOS
stdevd_gv = 0.2*( maxval(d_gv)-minval(d_gv) ) !20% del rango de anomalia de flexibilidad para ajustar datos
stdevd_mg = 0.2*( maxval(d_mg)-minval(d_mg) )

! DIMENSIONES DE MODELO
call readmodel(xcell,ycell,zcell,m_gv,m_mg)

n = size(m_gv) !numero de parametros
print*, "Total number of parameters: ", n

! DESVIACION STD DE MODELOS
nx = 20 ; ny = 20 ; nz = 10 ; step = 50 !mts
allocate( stdevm_gv(n), stdevm_mg(n) )
call model1D(nx,ny,nz,1,stdevm_gv)
call model1D(nx,ny,nz,1,stdevm_mg)

! MODELO INICIAL para INVERSION CONJUNTA 
! Distribucion homogenea de parametros en m0 = 0) 
m_gv = 0
m_mg = 0

! MODELO APRIORI 
allocate( mapr_gv(n), mapr_mg(n) )
!mapr_gv = 0
mapr_mg = 0
call model1D(nx,ny,nz,2,mapr_gv)
!call model1D(nx,ny,nz,2,mapr_mg)

! ARCHIVO BINARIO DE MATRICES A
allocate( A_gv(m,n), A_mg(m,n) )
OPEN(UNIT=1, FILE="input_a_gv.dat", ACTION="read", FORM="unformatted")
READ(1) A_gv
CLOSE(UNIT=1)

OPEN(UNIT=2, FILE="input_a_mg.dat", ACTION="read", FORM="unformatted")
READ(2) A_mg
CLOSE(UNIT=2)
!***********************************************************************************************************************************


!***********************************************************************************************************************************
! MODELADO MEDIANTE INVERSION CONJUNTA
TIPO_INVERSION = 0

!Se define el tipo de inversion a realizar a partir de la variable TIPO_INVERSION
! = 0: Inversion separada version en CPU
! = 1: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores Version en CPU
! = 2: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores version en GPU
! = 3: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores version en CPU/GPU
! = -1: no hacer inversion, solo prueba de datos de entrada
!***********************************************************************************************************************************
call cpu_time(start0)

if (TIPO_INVERSION == 0) then
      call separateInv(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                        d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TIPO_INVERSION == 1) then
      call jointGram5(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                      d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TIPO_INVERSION == 2) then
      call jointGramPL5(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                        d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TIPO_INVERSION == 3) then
      call jointGramPL6(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                        d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,reg,dip,strike,num_iters,err0)

elseif (TIPO_INVERSION == -1) then
      PRINT*, "Solo para prueba de datos de entrada"
else 
      PRINT*, 'El tipo de inversión fue seleccionado erroneamente, porfavor elige un TIPO_INVERSION apropiado'
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
write(13,*) 'H =',H,', Fi =',Fi,', Fd =',Fd
write(13,*) ' '
write(13,*) 'StdDeviation gv data = 2% rango de anomalia'
write(13,*) 'StdDeviation mg data = 2% rango de anomalia'
write(13,*) 'StdDeviation gv model = 1 todo modelo y 0.2 para cima del dominio gr/cm3 *MANUAL'
write(13,*) 'StdDeviation mg model = 1 todo modelo y 0.2 para cima del dominio A/m  *MANUAL'
write(13,*) ' '
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
write(13,*) ' '
write(13,*) 'Max Iterations Joint Inversion =',num_iters
write(13,*) 'Stop Condition, % Convergence less than', err0
write(13,*) ' '
write(13,*) 'Total execution time for TIPO_INVERSION:', TIPO_INVERSION
write(13,*) (finish0-start0),'seg =',(finish0-start0)/60, 'min'
close(unit=13)
!print*, '~(^-^)~  Process Finished  ~(^-^)~'
!***********************************************************************************************************************************

end program parameters
