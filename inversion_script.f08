program parameters
! Programa que utiliza una serie de subrutinas relacionados con modelado e inversion de métodos potenciales.
! Autor: Abraham Del Razo, IPICyT. Inicio de desarrollo: May-2020. Ultima actualización: Nov-2022.
! Mail: abraham.delrazo@ipicyt.edu.mx
!***********************************************************************************************************************************

use,intrinsic :: ieee_arithmetic
use iso_fortran_env , dp => real64
use auxiliar, only: readdata, std_dev_data, model_apr, readmodel, std_dev_model, mesh2D, coord_model, coord_model2, readtopo
use forward_modeling, only: A_gen1, A_gen2, A_gen3, forward_mod
!use inversion_serial
use inversion_parallel

implicit none

real, allocatable, dimension(:) :: d_gv, stdevd_gv, m_gv, mapr_gv, stdevm_gv
real, allocatable, dimension(:) :: d_mg, stdevd_mg, m_mg, mapr_mg, stdevm_mg
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real(dp), allocatable :: A_gv(:,:), A_mg(:,:)
real, allocatable, dimension(:) :: x, y, z
real :: xmin, xmax, ymin, ymax, zmax, fact_incr, step
integer :: nx, ny, nz, m, n, num_iters, TIPO_INVERSION, i, ii, jj, kk, num_dat_terr, count, count2D, zdom
real, allocatable, dimension(:) :: topo_model_3D
real :: H, Fi, Fd, err0
real :: param_reg(10), dip, strike
real :: start0, finish0

!***********************************************************************************************************************************
! PARAMETROS DE INVERSION

!Cuerpos, step 50, 30x20x15
!H = 40000
!Fi = 45
!Fd = 15

!Popo, step 1000
!H = 42200.2
!Fi = 90!47.15
!Fd = 0!5.76

!Sur India 01/09/2017, step 6500
!H = 29360.3 !nT
!Fi = 0 !grados Reducido al Ecuador
!Fd = 0  !grados

!Aeromag 01/04/1998
!H = 44465.0
!Fi = 51.216
!Fd = 7.333

!Terrestre 08/12/2022, step 100
H = 41922.8 !nT
Fi = 90 !Reducido al Polo
Fd = 0
!Fi = 50.86 !grados
!Fd = 4.509  !grados
!*************************

! REGULARIZATION PARAMETERS 
param_reg(1) = 1e7 !alpha_gv /smoothness
param_reg(2) = 1e7 !alpha_mg
param_reg(3) = 1!1e0 !beta_gv /a-priori information
param_reg(4) = 1!1e0 !beta_mg
param_reg(5) = 0!1e4 !gamma_gv /dip and strike
param_reg(6) = 0!1e7 !gamma_mg
param_reg(7) = 1e4!1e0 !eta_gv /verticalily
param_reg(8) = 1e7!1e0 !eta_mg
param_reg(9) = 1!1e7 !mu_gv /structural coupling
param_reg(10) = 1!1e7 !mu_mg
strike = 0 !degree  !Positive from north to west. *only apply with gamma different if zero
dip = -90 !degree  !Positive from east to north.

!*************************
!Condiciones de paro, método iterativo (aplica para inv. conjunta).
num_iters = 25
err0 = -0.1 !% variacion total de modelos en cada iteracion *negativo para desactivar restriccion
!*************************


!***********************************************************************************************************************************
! MISMOS nombres de variables se deben utilizar en subrutinas que usen estos arreglos, para poder utilizar es mismo arreglo en la declaracion de variables como allocable y no tener que usar otro espacio de memoria del mismo tamano.
!allocate( d_gv(rows), stdevd_gv(rows), d_mg(rows), stdevd_mg(rows), xobs(rows), yobs(rows), zobs_gv(rows), zobs_mg(rows) )
!allocate( xcell(nx+1), ycell(ny+1), zcell(nz+1) )


! Registrar tiempos de ejecución
OPEN(unit=0,file='output_time.txt',status='unknown')
! Evaluar eficiencia.
call cpu_time(start0)


!***********************************************************************************************************************************
! INPUT DATOS
!subroutine mesh2D(xmin,xmax,ymin,ymax,cell_size)
!call mesh2D(xmin,xmax,ymin,ymax,step)
!call mesh2D(366461,370461,2530721,2534721,200)

!subroutine readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
call readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)

m = size(d_gv) !numero de observaciones
print*, "Total number of data: ", size(d_gv)

xmin = minval(xobs)
xmax = maxval(xobs)
ymin = minval(yobs)
ymax = maxval(yobs)

!***********************************************************************************************************************************
! CONSTRUCCION DE MALLA 3D o dominio de celdas cubicas
step = 60 !mts
!zmax = 1300 !m.s.n.m elevacion plana
zmax = maxval(zobs_gv) !m.s.n.m >0 sobre el mar
!zmax = minval(zobs_gv) !m.s.n.m <0 debajo del mar
nz = 13 !numero de celdas en z
!fact_incr = 1.5 !para cuando se pueda celda variable en z

! Para modelos que usan las coordenadas de los datos como base
!subroutine coord_model(xmin,xmax,ymin,ymax,cell_size,zmax,nz,xcell,ycell,zcell,nx,ny,x,y,z) !Para celdas simétricas (cubos)
call coord_model(xmin,xmax,ymin,ymax,step,zmax,nz,xcell,ycell,zcell,nx,ny,x,y,z)

!subroutine coord_model2(xmin,xmax,ymin,ymax,cell_size,zmax,nz,fact_incr,xcell,ycell,zcell,nx,ny) !Para celdas variables en z (prismas rectangulares)
!call coord_model2(xmin,xmax,ymin,ymax,step,zmax,nz,fact_incr,xcell,ycell,zcell,nx,ny) !para coordenadas celdas rectangulares en vez de cubicas

!Dimensiones de malla (o número de parametros). (num celdas = num vertices-1)
print*, "nx, ny, nz: ", nx, ny, nz
n = nx*ny*nz !numero de parametros
print*, "Total number of parameters: ", n

!***********************************************************************************************************************************
! CONSTRUCCION DE OPERADORES A_gv Y A_mg
allocate( A_gv(m,n), m_gv(n), mapr_gv(n), stdevm_gv(n) )
allocate( A_mg(m,n), m_mg(n), mapr_mg(n), stdevm_mg(n) )

!subroutine A_gen3(xobs,yobs,zobs_gv,zobs_mg,x,y,z,Fi,Fd,A_gv,A_mg)
call A_gen3(xobs,yobs,zobs_gv,zobs_mg,x,y,z,Fi,Fd,A_gv,A_mg)
deallocate (x,y)!,z)

!print*, A_mg
!do ii = 1, m
!  do jj = 1, n
!      if (IEEE_IS_NAN(A_gv(ii,jj))) A_gv(ii,jj) = 0
!      if (IEEE_IS_NAN(A_mg(ii,jj))) A_mg(ii,jj) = 0
!  end do
!end do
!print*, A_mg


!***********************************************************************************************************************************
! DESVIACION STD DE DATOS
!porcentaje arbitrario a todos los datos.
!stdevd_gv = 0.05*( maxval(d_gv)-minval(d_gv) ) !5% del rango de anomalia de flexibilidad para ajustar datos
!stdevd_mg = 0.05*( maxval(d_mg)-minval(d_mg) )

!separando porcentajes para datos terrestres y sat/aero
num_dat_terr = 122
stdevd_gv(1:num_dat_terr) = 0.05*( maxval(d_gv)-minval(d_gv) ) !terr
stdevd_mg(1:num_dat_terr) = 0.05*( maxval(d_mg)-minval(d_mg) )

stdevd_gv(num_dat_terr+1:) = 0.10*( maxval(d_gv)-minval(d_gv) ) !sat/aero
stdevd_mg(num_dat_terr+1:) = 0.10*( maxval(d_mg)-minval(d_mg) )


!subroutine std_dev_data(nx,ny,d_type,stdev_data)
!call std_dev_data(nx,ny,1,stdev_gv)
!call std_dev_data(nx,ny,2,stdev_mg)
!print*, stdevd_gv

!***********************************************************************************************************************************
! MODELO INICIAL para INVERSION CONJUNTA 
! distribucion homogenea de parametros en m0 = 0) 
m_gv = 0
m_mg = 0

!subroutine readmodel(m_gv,m_mg)
!call readmodel(m_gv,m_mg)

!***********************************************************************************************************************************
! MODELO APRIORI 
mapr_gv = 0
mapr_mg = 0

!subroutine model_apr(nx,ny,nz,type,mapr)
!call model_apr(nx,ny,nz,3,mapr_gv)
!call model_apr(nx,ny,nz,3,mapr_mg)

!subroutine readmodel(m_gv,m_mg)
!call readmodel(mapr_gv,mapr_mg)

!***********************************************************************************************************************************
! DESVIACION STD DE MODELOS
!valor arbitrario amplio para no dar prioridad a ninguna celda
!stdevm_gv = 1e8
!stdevm_mg = 1e8

!Desviacion estandar sobre topografia muy cercana a cero para peso sea casi nulo
call readtopo(topo_model_3D)
count = 0
count2D = 0
do ii = 1,nx,1
      do jj = 1,ny,1
            count2D = count2D + 1
            zdom = topo_model_3D(count2D)
            do kk = 1,nz,1
                  count = count + 1
                  if ( (z(kk)) > (zdom) ) then
                        stdevm_gv(count) = 1e-5!0.00001
                        stdevm_mg(count) = 1e-5
                  else
                        stdevm_gv(count) = 1e8
                        stdevm_mg(count) = 1e8
                  end if
            end do
      end do
end do
deallocate (z,topo_model_3D)

!generar archivo con desviaciones estandar de modelo 3D
OPEN(unit=10,file='output_stdevm_gv.csv',status='unknown') !crear y escribir en archivo
! Header
WRITE(10,*) 'stdevm_gv[gr/cm3]'
! Body
count = 0
do ii = 1,nx,1
      do jj = 1,ny,1
            do kk = 1,nz,1
                  count = count + 1
                  WRITE(10,'(E9.1)') stdevm_gv(count)
            end do
      end do
end do
CLOSE(unit=10)


!subroutine std_dev_model(nx,ny,nz,d_type,stdev_model)
!call std_dev_model(nx,ny,nz,1,stdevm_gv)
!call std_dev_model(nx,ny,nz,2,stdevm_mg)
!print*, stdevm_gv

!***********************************************************************************************************************************
!Se define el tipo de inversion a realizar a partir de la variable TIPO_INVERSION
! = 0: no hacer inversion, solo prueba de datos de entrada
! = 1: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores Version en CPU
! = 2: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores version en GPU
! = 3: Inversion conjunta mediante enfoque Gramiano con 5 funcionales regularizadores version en CPU/GPU
! = 4: Inversion conjunta enfoque Gramiano con 5 funcionales regularizadores. Ver. Automatizada p/ experimentos multiples

TIPO_INVERSION = 1
!***********************************************************************************************************************************

if (TIPO_INVERSION == 1) then
      ! MODELADO MEDIANTE INVERSION CONJUNTA
      !subroutine jointGram(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
      !                     d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
      call jointGram5(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                      d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,dip,strike,num_iters,err0,1)

elseif (TIPO_INVERSION == 2) then
      ! MODELADO MEDIANTE INVERSION CONJUNTA
      !subroutine jointGramPL2(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
      !                  d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
      call jointGramPL5(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                        d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,dip,strike,num_iters,err0,1)

elseif (TIPO_INVERSION == 3) then
      ! MODELADO MEDIANTE INVERSION CONJUNTA
      !subroutine jointGramPL2(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
      !                  d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
      call jointGramPL6(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                        d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,dip,strike,num_iters,err0,1)

elseif (TIPO_INVERSION == 4) then
      ! Modo ciclos para hacer varios experimentos de corrido
      do i = 2, 8, 2

            param_reg(1) = 1e2 !alpha_gv /smoothness
            param_reg(2) = 1e2 !alpha_mg
            param_reg(3) = 1 !beta_gv /a-priori information
            param_reg(4) = 1 !beta_mg
            param_reg(5) = 1e1 !gamma_gv /dip and strike
            param_reg(6) = 1e1 !gamma_mg
            param_reg(7) = 0 !eta_gv /verticalily
            param_reg(8) = 0 !eta_mg
            param_reg(9) = 0!1e3 !mu_gv /structural coupling
            param_reg(10) = 0!1e3 !mu_mg
            dip = 0 !degree / only apply with gamma different if zero  !Positive from east to north.
            strike = 0 !degree  !Positive from north to west.

            !call jointGramPL4(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
            !            d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,dip,strike,num_iters,err0,1)
      enddo

elseif (TIPO_INVERSION == 0) then
      PRINT*, "Solo prueba de datos de entrada"
else 
      PRINT*, 'El tipo de inversión fue seleccionado erroneamente, porfavor elige un TIPO_INVERSION apropiado'
end if

!***********************************************************************************************************************************
! Archivo salida de parametros
OPEN(unit=13,file='output_parameters_inv.txt',status='unknown')
write(13,*) 'Number of data:', m
write(13,*) 'Volume size/Number of cells:', n
write(13,*) 'nx =',nx
write(13,*) 'nx =',ny
write(13,*) 'nz =',nz
write(13,*) 'Cell size =',step,'mts'
write(13,*) 'zmax =',zmax,'mts'
write(13,*) 'H =',H,', Fi =',Fi,', Fd =',Fd
write(13,*) ' '
write(13,*) 'StdDeviation gv data = 5% rango de anomalia para terrestre y 10% para satelital'
write(13,*) 'StdDeviation mg data = 5% rango de anomalia para terrestre y 10% para aero'
write(13,*) 'StdDeviation gv model = 1e8 todo modelo y 0.2 para fondo del dominio gr/cm3 *MANUAL'
write(13,*) 'StdDeviation mg model = 1e8 todo modelo y 0.2 para fondo del dominio A/m  *MANUAL'
write(13,*) ' '
write(13,*) 'alpha_gv =', param_reg(1)
write(13,*) 'alpha_mg =', param_reg(2)
write(13,*) 'beta_gv =', param_reg(3)
write(13,*) 'beta_mg =', param_reg(4)
write(13,*) 'gamma_gv =', param_reg(5)
write(13,*) 'gamma_mg =', param_reg(6)
write(13,*) 'eta_gv =', param_reg(7)
write(13,*) 'eta_mg =', param_reg(8)
write(13,*) 'mu_gv =', param_reg(9)
write(13,*) 'mu_mg =', param_reg(10)
write(13,*) ' '
write(13,*) 'Max Iterations Joint Inversion =',num_iters
write(13,*) 'Stop Condition, % Convergence less than', err0
close(unit=13)

!***********************************************************************************************************************************

! Evaluar eficiencia.
call cpu_time(finish0)
WRITE(0,*) 'Total execution time:', (finish0-start0), 'seg =', (finish0-start0)/60, 'min'
WRITE(0,*) '~(^-^)~  Process Finished  ~(^-^)~'
CLOSE(unit=0)


!deallocate(A_gv,d_gv,m_gv,mapr_gv,stdevd_gv,stdevm_gv,zobs_gv)
!deallocate(A_mg,d_mg,m_mg,mapr_mg,stdevd_mg,stdevm_mg,zobs_mg)
!deallocate(xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell)


end program parameters
!***********************************************************************************************************************************
! Para evitar escribir compilar cada modulo correr el script bash desde directorio de proyecto en terminal:
!./run_modelado_host.sh
!./run_modelado_device.sh
