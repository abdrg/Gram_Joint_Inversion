module auxiliar
! Contiene una serie de subrutinas auxiliares para el funcionamiento del programa principal
!
! Autor: Abraham Del Razo, IPICyT. Abr-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
contains

!*******************************************************************************
subroutine readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
! Subrutina que genera los arreglos d_gv y d_mg observados, ademas de stdev a partir de dos archivos tipo .csv

use iso_fortran_env, dp => real64

implicit none

real, dimension(:), allocatable :: d_gv, d_mg, stdevd_gv, stdevd_mg
real, dimension(:), allocatable :: xobs, yobs, zobs_gv, zobs_mg !coordenadas de cada punto anomalia 2d, no se usarán para nada pero es necesario leerlas.
character*256 :: CTMP
integer :: i = 0, IERR = 0, rows = 0

!OPEN(unit=3,file='input_d_gv_test.csv',status='old',access='sequential', form='formatted', action='read')
!OPEN(unit=5,file='input_d_mg_test.csv',status='old',access='sequential', form='formatted', action='read')
OPEN(unit=3,file='input_d_gv.csv',status='old',access='sequential', form='formatted', action='read')
OPEN(unit=5,file='input_d_mg.csv',status='old',access='sequential', form='formatted', action='read')


! Get number of lines
READ(3, * )   !! skip the header
READ(5, * )

do while (IERR == 0)
      rows = rows + 1
      READ(3,*,iostat=IERR) CTMP
      !READ(5,*,iostat=IERR) CTMP !solo necesario estimar tamaño de 1 archivo

end do
rows = rows - 1
!print*, "Total number of layers: ", rows

allocate( d_gv(rows), d_mg(rows), stdevd_gv(rows), stdevd_mg(rows), &
          xobs(rows), yobs(rows), zobs_gv(rows), zobs_mg(rows) )

! Read the file content
rewind(3)
rewind(5)
READ(3, * )   !! skip the header
READ(5, * )

do i = 1,rows,1
  READ(3,*) xobs(i), yobs(i), zobs_gv(i), d_gv(i), stdevd_gv(i)
  READ(5,*) xobs(i), yobs(i), zobs_mg(i), d_mg(i), stdevd_mg(i)
end do

CLOSE(unit=3)
CLOSE(unit=5)

!debug
!print*, 'shape=', shape(d_gv)!, 'd_gv=',d_gv
!print*, 'shape=', shape(stdevd_mg)!, 'stdevd_mg=',stdevd_mg

!deallocate(xobs,yobs,zobs_gv,zobs_mg) !el resto se desocupan despues del llamado del modulo principal

end subroutine readdata
!*******************************************************************************


subroutine readmodel(xcell,ycell,zcell,m_gv,m_mg)
! Subrutina que genera arreglos de parametros modelo para gv y mg a partir de dos archivos tipo .csv

use iso_fortran_env, dp=>real64

implicit none

real, dimension(:), allocatable, intent(out) :: m_gv, m_mg
real, dimension(:), allocatable, intent(out) :: xcell, ycell, zcell !coordenadas de cada punto anomalia 2d, no se usarán para nada pero es necesario leerlas.
character*256 :: CTMP
integer :: i = 0, IERR = 0, rows = 0

OPEN(unit=1,file='input_m_gv.csv',status='old',access='sequential', form='formatted', action='read')
OPEN(unit=2,file='input_m_mg.csv',status='old',access='sequential', form='formatted', action='read')

! Get number of lines
READ(1, * )   !! skip the header
READ(2, * )

do while (IERR == 0)
      rows = rows + 1
      READ(1,*,iostat=IERR) CTMP

end do
rows = rows - 1
!print(*,'(A,I0)') "Total number of layers: ", rows

allocate( m_gv(rows), m_mg(rows), xcell(rows), ycell(rows), zcell(rows) )

! Read the file content
rewind(1)
rewind(2)
READ(1, * )   !! skip the header
READ(2, * )

do i = 1,rows,1
  READ(1,*) xcell(i), ycell(i), zcell(i), m_gv(i)
  READ(2,*) xcell(i), ycell(i), zcell(i), m_mg(i)
end do

!debug
!print*, 'm_gv=',m_gv

CLOSE(unit=1)
CLOSE(unit=2)

!deallocate(xcell,ycell,zcell) !el resto se desocupan despues del llamado del modulo principal

end subroutine readmodel
!*******************************************************************************


subroutine readtopo(topo_model_3D)
! Subrutina que genera los arreglos d_gv y d_mg observados, ademas de stdev a partir de dos archivos tipo .csv

use iso_fortran_env, dp => real64

implicit none

real, dimension(:), allocatable, intent(out) :: topo_model_3D
character*256 :: CTMP
integer :: i = 0, IERR = 0, rows = 0

OPEN(unit=6,file='input_dominio3D_topo.csv',status='old',access='sequential', form='formatted', action='read')

! Get number of lines
READ(6, * )   !! skip the header

do while (IERR == 0)
      rows = rows + 1
      READ(6,*,iostat=IERR) CTMP

end do
rows = rows - 1
!print*, "Total number of layers: ", rows

allocate( topo_model_3D(rows) )

! Read the file content
rewind(6)
READ(6, * )   !! skip the header

do i = 1,rows,1
  READ(6,*) topo_model_3D(i)
end do

CLOSE(unit=6)

end subroutine readtopo
!*******************************************************************************

subroutine model_apr(nx,ny,nz,type,mapr)
! Subrutina que genera un vector con los valores de un modelo apriori

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz, type
real, intent(out) :: mapr(nx*ny*nz)
integer :: i, j, k, count
real, allocatable, dimension(:,:,:) :: mapr3D

allocate( mapr3D(nx,ny,nz) )

if (type == 1) then
      mapr3D = 0
      mapr3D(7:15,4:9,1:5) = -1

else if (type == 2) then
      mapr3D = 0
      mapr3D(7:15,4:9,1:5) = 1

else if (type == 3) then
      mapr3D = 0
      mapr3D(4:7,4:8,3:6) = 1 !cubo
      do i = 1,11
            mapr3D(24-i:27-i,3:7,1+i) = 1 !dique
      end do
      mapr3D(3:16,13:17,6:8) = 1 !loza
      mapr3D(24:27,13:17,3:6) = -1 !cavidad

end if

!conversión a arreglo de 1D
count = 0
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      mapr(count) = mapr3D(i,j,k)
    end do
  end do
end do

deallocate (mapr3D)

return
end subroutine model_apr
!*******************************************************************************


subroutine std_dev_model(nx,ny,nz,d_type,stdev_model)
! Subrutina que genera un vector con las desviaciones estandard por sectores del modelo.
! ***Pendiente conversion para pedir info desde archivo o teclado.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz, d_type
real, intent(out) :: stdev_model(nx*ny*nz)
real :: stdev_model_3D(nx,ny,nz)
integer :: i, j, k, count

! stdevm = 1 significa sin pesos por ser neutro multiplicativo.
! stdevm no puede ser cero (se indetermina calculo de pesos), pero puede ser muy cercano a cero.

if (d_type == 1) then
  stdev_model_3D = 1e8 !2e-1 !todas las celdas. 20% del contraste entre promedio y max o min del modelo apriori (+-1 para m0 del 1 al 4).
  !stdev_model_3D(1:nx,1:ny,1) = 0.2 !celdas limite superior. Para evitar propagacion anomalia desde superficie.
  !stdev_model_3D(1:nx,1:ny,nz) = 0.1 !celdas limite inferior. Para evitar propagacion de anomalia hasta fondo.
  !stdev_model_3D(1,1:ny,1:nz) = 0.001 !celdas limites laterales
  !stdev_model_3D(nx,1:ny,1:nz) = 0.001
  !stdev_model_3D(1:nx,1,1:nz) = 0.001
  !stdev_model_3D(1:nx,ny,1:nz) = 0.001
  !determinacion de topografía en modelo3D

else if (d_type == 2) then
  stdev_model_3D = 1e8
  !stdev_model_3D(1:nx,1:ny,1) = 0.2
  !stdev_model_3D(1:nx,1:ny,nz) = 0.1
  !stdev_model_3D(1,1:ny,1:nz) = 0.001
  !stdev_model_3D(nx,1:ny,1:nz) = 0.001
  !stdev_model_3D(1:nx,1,1:nz) = 0.001
  !stdev_model_3D(1:nx,ny,1:nz) = 0.001

end if

!conversión a arreglo de 1D
count = 0
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      stdev_model(count) = stdev_model_3D(i,j,k)
    end do
  end do
end do

return
end subroutine std_dev_model
!*******************************************************************************


subroutine std_dev_data(nx,ny,d_type,stdev_data)
! Subrutina que genera un vector con las desviaciones estandard para cada dato observado
! escogidas arbitrariamente, idealmente relacionadas con la precision (en mGal o nT) de la adquisición.
! SUMA DE TODAS DESVIACIONES DEBERIA SER SIEMPRE IGUAL A 0.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in):: nx, ny, d_type
real, intent(out) :: stdev_data(nx*ny)
integer :: i, j, count
real :: v

if (d_type == 1) then
  !v = 0.002  !Gal, basado en la desviacion std del cg5 +-1microGal. !Valores en correspodencia (magnitud) similara anomalia generada, si no residual do-dc divergen en lugar de converger a cero.
  v = 0.06 !experimento oldenburg

else if (d_type == 2) then
  !v = 1 !nT, basado en la precision de magnetometro g-857.
  v = 5 !experimento oldenburg
end if

!desviacion std para cada dato default
count = 0
do i = 1,nx
  do j =1,ny
    count = count + 1
    stdev_data(count) = v
  end do
end do

return

end subroutine std_dev_data
!*******************************************************************************


subroutine cov(num,type,std_dev,covar)
!Subrutina que calcula la matriz de covarianza INVERSA para los datos o los parametros del modelo.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: num, type
real, intent(in) :: std_dev(num)
real, intent(out) :: covar(num,num)
integer :: i, j

!construir la matriz diagonal partir de valores propios.
do i = 1,num
  do j = 1,num
    if ((i == j).and.(type == 1)) then
      covar(i,j) = 1/(std_dev(i)**2) !autovalores son inverso de varianza. Resulta matriz de covarianza C^-1.
    else if ((i == j).and.(type == 2)) then
      covar(i,j) = 1/std_dev(i) !autovalores son inverso de desviacion estandar. Resulta matriz de pesos W.
    else
      covar(i,j) = 0
    end if
  end do
end do
!print *, shape(covar)

return
end subroutine cov
!*******************************************************************************


subroutine mesh2D(xmin,xmax,ymin,ymax,cell_size)
! Subrutina que genera una malla regular de puntos 2D 

use iso_fortran_env, dp=>real64

implicit none

real, intent(in) :: xmin,xmax,ymin,ymax,cell_size
integer :: i, j, count, nx, ny
real :: divx, divy
real, allocatable, dimension(:) :: x, y

divx = (xmax-xmin)/cell_size
divy = (ymax-ymin)/cell_size

nx = divx+3
ny = divy+3

allocate( x(nx), y(ny) )

! Centro celdas en x
x(1) = xmin-(cell_size/2)
do i = 2,nx,1
      x(i) = x(i-1)+cell_size
end do

! Centro celdas en y
y(1) = ymin-(cell_size/2)
do j = 2,ny,1
      y(j) = y(j-1)+cell_size
end do


OPEN(unit=7,file='output_mesh2d.csv',status='unknown')
WRITE(7,*) 'x[m],y[m]'
!conversión a arreglo de 1D
count = 0
do i = 1,nx
  do j = 1,ny
      WRITE(7,'(F10.3,",",F11.3)') x(i), y(j)
  end do
end do

close(unit=7)

return
end subroutine mesh2D
!*******************************************************************************


subroutine coord_model(xmin,xmax,ymin,ymax,cell_size,zmax,nz,xcell,ycell,zcell,nx,ny,x,y,z)
! Subrutina que genera valores para CENTRO DE PRISMAS RECTOS SIMETRICOS basado en xobs_min, yobs_min. nx, ny, nz son el numero de celdas y
! por lo tanto el numero de parametros en el modelo, pero cada celda se define por 4 esquinas por lo que se necesita
! un vector de esquinas con nx+1, ny+1, nz+1 entradas (para poder definir la ultima celda). Siempre num esquinas+1 = num celdas

use iso_fortran_env, dp=>real64

implicit none

real, intent(in) :: xmin, xmax, ymin, ymax, cell_size, zmax
integer, intent(in) :: nz
real, allocatable, intent(out), dimension(:) :: xcell, ycell, zcell
integer, intent(out) :: nx, ny
real, allocatable, intent(out), dimension(:) :: x, y, z
real  :: divx, divy
integer :: i, j, k, count

divx = ( (xmax-xmin)/cell_size )+1
divy = ( (ymax-ymin)/cell_size )+1

nx = divx+2 !una celda extra al inicio y al fin del dominio
ny = divy+2

allocate( x(nx), y(ny), z(nz) )
! Centro celdas en x
x(1) = xmin-(cell_size)!/2)
do i = 2,nx,1
      x(i) = x(i-1)+cell_size
end do

! Centro celdas en y
y(1) = ymin-(cell_size)!/2)
do j = 2,ny,1
      y(j) = y(j-1)+cell_size
end do

! Centro celdas en z
z(1) = zmax-(cell_size/2)-1
do k = 2,nz+1,1
      z(k) = z(k-1)-cell_size

end do


allocate( xcell(nx*ny*nz), ycell(nx*ny*nz), zcell(nx*ny*nz) )
count= 0 
do i = 1,nx,1
      do j = 1,ny,1
            do k = 1,nz,1    
                  count = count + 1              
                  xcell(count) = x(i) 
                  ycell(count) = y(j) 
                  zcell(count) = z(k) 
            end do
      end do
end do

!deallocate (x,y,z)

end subroutine coord_model
!*******************************************************************************


subroutine coord_model2(xmin,xmax,ymin,ymax,cell_size,zmax,nz,fact_incr,xcell,ycell,zcell,nx,ny) !**POR TERMINAR DE PROGRAMAR
! Subrutina que genera valores para CENTROS DE PRISMAS RECTANGULARES DE TAMAÑO VARIABLE basado en xobs_min, yobs_min. nx, ny, nz son el numero de celdas y
! por lo tanto el numero de parametros en el modelo, pero cada celda se define por 4 esquinas por lo que se necesita
! un vector de esquinas con nx+1, ny+1, nz+1 entradas (para poder definir la ultima celda). Siempre num esquinas+1 = num celdas

use iso_fortran_env, dp=>real64

implicit none

real, intent(in) :: xmin, xmax, ymin, ymax, cell_size, zmax, fact_incr
integer, intent(in) :: nz
real, allocatable, dimension(:) :: xcell, ycell, zcell
integer, intent(out) :: nx, ny
real, allocatable, dimension(:) :: x, y, z
integer :: divx, divy, i, j, k, count

divx = (xmax-xmin)/cell_size
divy = (ymax-ymin)/cell_size

nx = divx+1+4
ny = divy+1+4

allocate(x(divx+1+4), y(divy+1+4), z(nz) )

! Centro celdas en x
x(1) = xmin-(cell_size*4)
x(2) = xmin-(cell_size*2)
x(3) = xmin
do i = 4,divx+1,1
      x(i) = x(i-1)+cell_size
end do
x(divx+2) = xmax+(cell_size*2)
x(divx+3) = xmax+(cell_size*4)

! Centro celdas en y
y(1) = ymin-(cell_size*4)
y(2) = ymin-(cell_size*2)
y(3) = ymin
do j = 4,divy+1,1
      y(i) = y(i-1)+cell_size
end do
y(divy+2) = ymax+(cell_size*2)
y(divy+3) = ymax+(cell_size*4)

! Centro celdas en z
z(1) = zmax-(cell_size/4)
z(2) = zmax-(cell_size/2)
z(3) = zmax-(cell_size)
do k = 4,nz,1
      z(k) = z(k-1)-(z(k-1)*fact_incr)
end do


allocate( xcell(nx*ny*nz), ycell(nx*ny*nz), zcell(nx*ny*nz) )

count = 0 
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      xcell(count) = x(i)
      ycell(count) = y(j)
      zcell(count) = z(k)
    end do
  end do
end do

deallocate (x,y,z)

end subroutine coord_model2
!*******************************************************************************


subroutine coord_voxel(nx,ny,nz,block_size,xcell,ycell,zcell)
! Subrutina que asigna valores en km para esquinas de las celdas cubicas (VOXELS) donde nx, ny, nz son el numero de celdas y
! por lo tanto el numero de parametros en el modelo, pero cada celda se define por 4 esquinas por lo que se necesita
! un vector de esquinas con nx+1, ny+1, nz+1 entradas (para poder definir la ultima celda). Siempre num esquinas+1 = num celdas

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(in) :: block_size !tamaño de celda
real :: x(nx+1), y(ny+1), z(nz+1)
real, intent(out) :: xcell(nx), ycell(ny), zcell(nz)
integer :: i, j, k

! Celdas en x
do i = 1,nx+1
      x(i) = (i-1)*block_size
end do

! Celdas en y
do j = 1,ny+1
      y(j) = (j-1)*block_size
end do

! Celdas en z
do k = 1,nz+1
      z(k) = -(k-1)*block_size
end do


! Crear puntos de obseración a la mitad de cada celda.
do i = 1,nx,1
  xcell(i) = x(i) + ( abs(x(i)-x(i+1))/2 )
end do
do j = 1,ny,1
  ycell(j) = y(j) + ( abs(y(j)-y(j+1))/2 )
end do
do k = 1,nz,1
  zcell(k) = z(k) - ( abs(z(k)-z(k+1))/2 )
end do


end subroutine coord_voxel
!*******************************************************************************

subroutine gaussian_noise(num_ele,percentage,data,noised_data)
! Subrutina que calcula una cantidad (igual a el numero de mediciones) de numeros aleatorios
! gaussianos (promedio cero y desviacion estandar 1) mediante ele método BoxMuller, posteriomente
! normaliza entre -1 y 1 y finalmente los escala respecto a un porcentaje del valor max y min de las mediociones.

use iso_fortran_env, dp=>real64

implicit none

real, parameter :: pi = 3.1415927

integer, intent(in) :: num_ele
real, intent(in) :: percentage
real, intent(in) :: data(num_ele)
real :: u1, u2
integer :: i!, sd=1, mean=0
real :: gauss_noise(num_ele), regul_noise(num_ele), escaled_noise(num_ele)
real, intent(out) :: noised_data(num_ele)

!no jaló box-muller porque se ven bandas cada 2 espacios en y al graficar.
!call random_number(u1)
!call random_number(u2)

!if ( (u1==0).or.(u2==0) ) then
!  call random_number(u1)
!  call random_number(u2)

!else
  !do i=2,num_ele,2 
  !    gauss_noise(i) = ( sqrt(-2*log(u1)) * cos(2*pi*u2) * sd ) + mean
  !    gauss_noise(i-1) = ( sqrt(-2*log(u1)) * sin(2*pi*u2) * sd ) + mean
  !end do
!end if

do i = 1,num_ele,1
      call random_number(u1)
      gauss_noise(i) = u1
end do

! normalizacion y escalado segun maximo y minimo de datos
do i=1,num_ele,1
  regul_noise(i) = ( ( (gauss_noise(i)-minval(gauss_noise)) / (maxval(gauss_noise)-minval(gauss_noise)) )*2)-1
  escaled_noise(i) = regul_noise(i) * (maxval(data)-minval(data))*(percentage/100)

  noised_data(i) = data(i) + escaled_noise(i)
end do

end subroutine gaussian_noise
!*******************************************************************************


subroutine partialder(nx,ny,nz,step,dir,pd)
! Subrutina que construye una matriz de derivada parcial en x, y o z de un dominio en 3D.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(in) :: step
character(len=1), intent(in) :: dir
real(dp), intent(out) :: pd(nx*ny*nz,nx*ny*nz) !Tamaño de matriz derivada parcial
integer :: i, j, k, count
real :: fac

fac = 1/(step)
!matriz de ceros
pd(:,:) = 0

count = 0
do i=1,nx,1
  do j=1,ny,1
    do  k=1,nz,1
      count = count + 1 

      !1 caso para para centro del cubo (sin considerar los nodos frontera)
      if ( (i /= 1) .and. (i /= nx) .and. (j /= 1) .and. (j /= ny) .and. &
         (k /= 1) .and. (k /= nz) ) then

         if (dir == 'z') then
           pd(count,count-1) = -1 * 0.5*fac  !nodos en z
           pd(count,count+1) = 1 * 0.5*fac
         else if (dir == 'y') then
           pd(count,count-nz) = -1 * 0.5*fac  !nodos en y
           pd(count,count+nz) = 1 * 0.5*fac
         else if (dir == 'x') then
           pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
           pd(count,count+(ny*nz)) = 1 * 0.5*fac
         end if

      !8 casos diferentes para vertices del cubo
      else if ( (i==1) .and. (j==1) .and. (k==1) ) then !vertice luf
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
        end if

      else if ( (i==nx) .and. (j==1) .and. (k==1) ) then !vertice ruf
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      else if ( (i==1) .and. (j==1) .and. (k==nz) ) then !vertice ldf
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
        end if

      else if ( (i==nx) .and. (j==1) .and. (k==nz) ) then !vertice rdf
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      else if ( (i==1) .and. (j==ny) .and. (k==1) ) then !vertice lub
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
        end if

      else if ( (i==nx) .and. (j==ny) .and. (k==1) ) then !vertice rub
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      else if ( (i==1) .and. (j==ny) .and. (k==nz) ) then !vertice ldb
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
        end if

      else if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then !vertice rdb
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      !12 casos diferentes para aristas del cubo
    else if ( (i==1) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista lf
      !pd(count,count) = -1 * fac  !nodo a derivar
      if (dir == 'z') then
        pd(count,count-1) = -1 * 0.5*fac  !nodos en z
        pd(count,count+1) = 1 * 0.5*fac
      else if (dir == 'y') then
        pd(count,count) = -1 * fac  !nodo a derivar
        pd(count,count+nz) = 1 * fac  !nodos en y
      else if (dir == 'x') then
        pd(count,count) = -1 * fac  !nodo a derivar
        pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
      end if

      else if ( (i==nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista rf
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==1) ) then !arista uf
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==nz) ) then !arista df
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+nz) = 1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i==1) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista lb
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac  !nodos en x
        end if

      else if ( (i==nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista rb
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac  !nodos en x
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==1) ) then !arista ub
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==nz) ) then !arista db
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista lu
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = -1 * fac !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac !nodos en x
        end if

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista ld
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = -1 * fac !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac !nodos en x
        end if

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista ru
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = 1 * fac !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac !nodos en x
        end if

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista rd
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = 1 * fac !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac !nodos en x
        end if

      !6 casos diferentes para caras del cubo
    else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !cara f
      if (dir == 'z') then
        pd(count,count-1) = -1 * 0.5*fac  !nodos en z
        pd(count,count+1) = 1 * 0.5*fac
      else if (dir == 'y') then
        pd(count,count) = -1 * fac  !nodo a derivar
        pd(count,count+nz) = 1 * fac  !nodos en y
      else if (dir == 'x') then
        pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
        pd(count,count+(ny*nz)) = 1 * 0.5*fac
      end if

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !cara b
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-nz) = -1 * fac  !nodos en y
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara l
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = -1 * fac !nodo a derivar
          pd(count,count+(ny*nz)) = 1 * fac !nodos en x
        end if

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara r
        if (dir == 'z') then
          pd(count,count-1) = -1 * 0.5*fac  !nodos en z
          pd(count,count+1) = 1 * 0.5*fac
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count) = 1 * fac !nodo a derivar
          pd(count,count-(ny*nz)) = -1 * fac !nodos en x
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !cara u
        if (dir == 'z') then
          pd(count,count) = -1 * fac  !nodo a derivar
          pd(count,count+1) = 1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !cara d
        if (dir == 'z') then
          pd(count,count) = 1 * fac  !nodo a derivar
          pd(count,count-1) = -1 * fac  !nodos en z
        else if (dir == 'y') then
          pd(count,count-nz) = -1 * 0.5*fac !nodos en y
          pd(count,count+nz) = 1 * 0.5*fac
        else if (dir == 'x') then
          pd(count,count-(ny*nz)) = -1 * 0.5*fac  !nodos en x
          pd(count,count+(ny*nz)) = 1 * 0.5*fac
        end if

      end if

    end do
  end do
end do

end subroutine partialder
!*******************************************************************************


subroutine laplaciano(nx,ny,nz,step,L)
! Contruye el laplaciano o matriz de suavidad de una matriz de datos 3D reorganizada como un vector (triple for anidado en orden ijk),
! a partir de segundas derivadas numericas calculadas adelantadas, atrasadas o centradas dependiendo de
! la posición del nodo en la matriz 3d original, considerando deltax=deltay=deltaz como el tamaño de paso.
! l=left, r=right, f=front, b=back, u=up, d=down

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(in) :: step
real(dp), intent(out) :: L(nx*ny*nz,nx*ny*nz) !Tamaño de matriz laplaciano
integer :: i, j, k, count
real :: fac

fac = 1/(step**2)
!matriz de ceros
L(:,:) = 0

count = 0
do i=1,nx,1
  do j=1,ny,1
    do  k=1,nz,1
      count = count + 1  !este sistema de conteo no funcionará en paralelizacion

      !1 caso para para centro del cubo (sin considerar los nodos frontera)
      if ( (i /= 1) .and. (i /= nx) .and. (j /= 1) .and. (j /= ny) .and. &
         (k /= 1) .and. (k /= nz) ) then

         L(count,count) = -6 * fac  !nodo a derivar
         L(count,count-1) = 1 * fac  !nodos en z
         L(count,count+1) = 1 * fac
         L(count,count-nz) = 1 * fac  !nodos en y
         L(count,count+nz) = 1 * fac
         L(count,count-(ny*nz)) = 1 * fac  !nodos en x
         L(count,count+(ny*nz)) = 1 * fac

      !8 casos diferentes para vertices del cubo
      else if ( (i==1) .and. (j==1) .and. (k==1) ) then !vertice luf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k==1) ) then !vertice ruf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==1) .and. (k==nz) ) then !vertice ldf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k==nz) ) then !vertice rdf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k==1) ) then !vertice lub
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k==1) ) then !vertice rub
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k==nz) ) then !vertice ldb
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then !vertice rdb
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      !12 casos diferentes para aristas del cubo
    else if ( (i==1) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista lf
        !L(count,count) = 0 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista rf
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==1) ) then !arista uf
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==nz) ) then !arista df
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista lb
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista rb
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==1) ) then !arista ub
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==nz) ) then !arista db
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista lu
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista ld
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista ru
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista rd
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      !6 casos diferentes para caras del cubo
    else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !cara f
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !cara b
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara l
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara r
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !cara u
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !cara d
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      end if

    end do
  end do
end do

end subroutine laplaciano
!*******************************************************************************


subroutine test_L(nx,ny,nz,step,d_type,m_vec,m_smooth)
! subrutina escrita para poder graficar el resultado de aplicar el operador laplaciano
! en modelo m codificado como un vector 1D.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz, d_type
real, intent(in) :: step
real, intent(in) :: m_vec(nx*ny*nz)
real, intent(out) :: m_smooth(nx*ny*nz)
real(dp) :: L_mat(nx*ny*nz,nx*ny*nz)
real :: x(nx+1), y(ny+1), z(nz+1), xcell(nx), ycell(ny), zcell(nz)
integer :: i, j, k, count, n

call laplaciano(nx,ny,nz,step,L_mat)

n = nx*ny*nz
!function DGEMV(TRANSA,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) !Operacion BLAS Matriz por Vector.
!call SGEMV('N',n,n,1.0,L_mat,n,m_vec,1,0.0,m_smooth,1)
m_smooth = MATMUL(L_mat,m_vec)

! Crear puntos centro de celda.
call coord_voxel(nx,ny,nz,step,x,y,z)
do i = 1,nx,1
  xcell(i) = x(i) + ( abs(x(i+1)-x(i))/2 )
end do

do j = 1,ny,1
  ycell(j) = y(j) + ( abs(y(j+1)-y(j))/2 )
end do

do k = 1,nz,1
  zcell(k) = z(k) - ( abs(z(k+1)-z(k))/2 )
end do

! Escribir la posicion nodos en un texto formato x,y,z,parametro.
if (d_type == 1) then
  OPEN(unit=9,file='output_gv_m_smooth.csv',status='unknown') !crear y escribir en archivo
  write(9,*) 'x [km],y [km],z [km],Density Contrast [g/cm**3]'
else if (d_type == 2) then
  OPEN(unit=9,file='output_mg_m_smooth.csv',status='unknown') !crear y escribir en archivo
  write(9,*) 'x [km],y [km],z [km],Magnetization Contrast [A/m]'
end if

! Body
count = 0
do i = 1,nx,1
  do j = 1,ny,1
    do k = 1,nz,1
      count = count + 1
      write(9,*) xcell(i),',',ycell(j),',',zcell(k),',',m_smooth(count)
    end do
  end do
end do
close(unit=9)

end subroutine test_L
!*******************************************************************************


subroutine xgrad(nx,ny,nz,step,m1,m2,xgi,xgj,xgk)
! Subrutina que calcula el gradiente cruzado de dos distribuciones de propiedades
! en 3D haciendo uso de la rutina partialder para obtener los componentes del cada gradiente
! y posteriormente hacer el producto vectorial de estos dos gradientes.

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(in) :: step
real, intent(in) :: m1(nx*ny*nz), m2(nx*ny*nz)
real, intent(out) :: xgi(nx*ny*nz), xgj(nx*ny*nz), xgk(nx*ny*nz)
real(dp) :: dx(nx*ny*nz,nx*ny*nz), dy(nx*ny*nz,nx*ny*nz), dz(nx*ny*nz,nx*ny*nz)
real :: dxm1(nx*ny*nz), dym1(nx*ny*nz), dzm1(nx*ny*nz), dxm2(nx*ny*nz), &
                    dym2(nx*ny*nz), dzm2(nx*ny*nz)
real :: x(nx+1), y(ny+1), z(nz+1), xcell(nx), ycell(ny), zcell(nz)
integer :: i, j, k, count, n
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

! construir matrices de derivadas
call partialder(nx,ny,nz,step,'x',dx)
call partialder(nx,ny,nz,step,'y',dy)
call partialder(nx,ny,nz,step,'z',dz)

! calcular derivadas en x, y, z para ambas distribuciones
n = nx*ny*nz
!function DGEMV(TRANSA,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
dxm1 = MATMUL(dx,m1)
dym1 = MATMUL(dy,m1)
dzm1 = MATMUL(dz,m1)
dxm2 = MATMUL(dx,m2)
dym2 = MATMUL(dy,m2)
dzm2 = MATMUL(dz,m2)

! calcular producto vectorial grad_m1 x grad_m2
do i = 1,n,1

  xgi(i) = ( dym1(i)*dzm2(i) ) - ( dym2(i)*dzm1(i) )
  xgj(i) = ( dxm2(i)*dzm1(i) ) - ( dxm1(i)*dzm2(i) )
  xgk(i) = ( dxm1(i)*dym2(i) ) - ( dxm2(i)*dym1(i) )

end do

call coord_voxel(nx,ny,nz,step,x,y,z)
!call quasi_inf(xmax,ymax,zmax,x,y,z)

! Crear puntos centro de celda.
do i = 1,nx,1
  xcell(i) = x(i) + ( abs(x(i+1)-x(i))/2 )
end do

do j = 1,ny,1
  ycell(j) = y(j) + ( abs(y(j+1)-y(j))/2 )
end do

do k = 1,nz,1
  zcell(k) = z(k) - ( abs(z(k+1)-z(k))/2 )
end do

! Escribir la posicion y valor de componentes del vector en un texto formato x,y,z,parametro.
OPEN(unit=7,file='output_xgrad.csv',status='unknown') !crear y escribir en archivo
write(7,*) 'x [km],y [km],z [km], xgrad_i, xgrad_j, xgrad_k' !Encabezado

! Cuerpo
count = 0
do i = 1,nx,1
  do j = 1,ny,1
    do k = 1,nz,1
      count = count + 1
      write(7,*) xcell(i),',',ycell(j),',',zcell(k),',',xgi(count),',',xgj(count),',',xgk(count)
    end do
  end do
end do

close(unit=7)

! Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: xgrag =', (finish-start), 'seg =', (finish-start)/60, 'min'

end subroutine xgrad
!*******************************************************************************


subroutine quasi_inf(nx,ny,nz,x,y,z)
! Subrutina que agrega 6 celdas por lado intentando simular un espacio quasi-infinito.
! tambien agrega una ultima celda al final de los arreglos para poder cuadrar al considerar pares de datos para los calculos.

use iso_fortran_env , dp => real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(out) :: x(nx+1), y(ny+1), z(nz+1)
integer :: i, j, k

! Celdas en x agregando espacio quasi-infinito al inicio y final.
do i = 1,nx+1,1
  if (i == 1) then
    x(i) = -500
  else if (i == 2) then
    x(i) = -100
  else if (i == 3) then
    x(i) = -50
  else if (i == 4) then
    x(i) = -20
  else if (i == 5) then
    x(i) = -5
  else if (i == 6) then
    x(i) = -2
  else if ( (i >= 7).and.(i <= nx-6) ) then !valores del perfil
    x(i) = i-7  !i-7 porque -6 de las celdas anteriores -1 de que el contador inicia en cero
  else if (i == nx-5) then
    x(i) = (nx-13)+2 !nx-13 porque -12 celdas agregadas -1 de que el contador inicia en cero
  else if (i == nx-4) then
    x(i) = (nx-13)+5
  else if (i == nx-3) then
    x(i) = (nx-13)+20
  else if (i == nx-2) then
    x(i) = (nx-13)+50
  else if (i == nx-1) then
    x(i) = (nx-13)+100
  else if (i == nx) then
    x(i) = (nx-13)+500
  else if ( i == nx+1) then !celda extra para cuadrar calculos en pares de datos.
    x(i) = (nx-13)+501
  end if
end do

! Celdas en y agregando espacio quasi-infinito al inicio y final.
do j = 1,ny+1,1
  if (j == 1) then
    y(j) = -500
  else if (j == 2) then
    y(j) = -100
  else if (j == 3) then
    y(j) = -50
  else if (j == 4) then
    y(j) = -20
  else if (j == 5) then
    y(j) = -5
  else if (j == 6) then
    y(j) = -2
  else if ( (j >= 7).and.(j <= ny-6) ) then  !valores del perfil
    y(j) = j-7
  else if (j == ny-5) then
    y(j) = (ny-13)+2
  else if (j == ny-4) then
    y(j) = (ny-13)+5
  else if (j == ny-3) then
    y(j) = (ny-13)+20
  else if (j == ny-2) then
    y(j) = (ny-13)+50
  else if (j == ny-1) then
    y(j) = (ny-13)+100
  else if (j == ny) then
    y(j) = (ny-13)+500
  else if ( j == ny+1) then !celda extra para cuadrar calculos en pares de datos.
    y(j) = (ny-13)+501
  end if
end do

! Celdas en z. Pendiente crearlas de forma log aumentando en profundidad.
!z = (/ (k-1,k=1,nz+1,1) /)  !escrito como do implicito
do k = 1,nz+1,1
  z(k) = -(k-1)
end do

end subroutine quasi_inf
!*******************************************************************************



subroutine laplaciano2(nx,ny,nz,step,L)
! Contruye el laplaciano o matriz de suavidad de una matriz de datos 3D reorganizada como un vector (triple for anidado en orden ijk),
! a partir de segundas derivadas numericas calculadas adelantadas, atrasadas o centradas dependiendo de
! la posición del nodo en la matriz 3d original, considerando deltax=deltay=deltaz como el tamaño de paso.
! l=left, r=right, f=front, b=back, u=up, d=down

use iso_fortran_env, dp=>real64

implicit none

integer, intent(in) :: nx, ny, nz
real, intent(in) :: step
real(dp), intent(out) :: L(nx*ny*nz,nx*ny*nz) !Tamaño de matriz laplaciano
integer :: i, j, k, count
real :: fac

fac = 1/(step**2)
!matriz de ceros
L(:,:) = 0

count = 0
do i=1,nx,1
  do j=1,ny,1
    do  k=1,nz,1
      count = count + 1  !este sistema de conteo no funcionará en paralelizacion

      !1 caso para para centro del cubo (sin considerar los nodos frontera)
      if ( (i /= 1) .and. (i /= nx) .and. (j /= 1) .and. (j /= ny) .and. &
         (k /= 1) .and. (k /= nz) ) then

         L(count,count) = -6 * fac  !nodo a derivar
         L(count,count-1) = 1 * fac  !nodos en z
         L(count,count+1) = 1 * fac
         L(count,count-nz) = 1 * fac  !nodos en y
         L(count,count+nz) = 1 * fac
         L(count,count-(ny*nz)) = 1 * fac  !nodos en x
         L(count,count+(ny*nz)) = 1 * fac

      !8 casos diferentes para vertices del cubo
      else if ( (i==1) .and. (j==1) .and. (k==1) ) then !vertice luf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k==1) ) then !vertice ruf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==1) .and. (k==nz) ) then !vertice ldf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k==nz) ) then !vertice rdf
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k==1) ) then !vertice lub
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k==1) ) then !vertice rub
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k==nz) ) then !vertice ldb
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k==nz) ) then !vertice rdb
        L(count,count) = 3 * fac  !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      !12 casos diferentes para aristas del cubo
    else if ( (i==1) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista lf
        !L(count,count) = 0 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !arista rf
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==1) ) then !arista uf
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k==nz) ) then !arista df
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista lb
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !arista rb
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==1) ) then !arista ub
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k==nz) ) then !arista db
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista lu
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista ld
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !arista ru
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !arista rd
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      !6 casos diferentes para caras del cubo
    else if ( (i/=1) .and. (i/=nx) .and. (j==1) .and. (k/=1) .and. (k/=nz) ) then !cara f
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count+nz) = -2 * fac  !nodos en y
        L(count,count+(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j==ny) .and. (k/=1) .and. (k/=nz) ) then !cara b
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = -2 * fac  !nodos en y
        L(count,count-(2*nz)) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i==1) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara l
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count+(ny*nz)) = -2 * fac  !nodos en x
        L(count,count+(2*ny*nz)) = 1 * fac

      else if ( (i==nx) .and. (j/=1) .and. (j/=ny) .and. (k/=1) .and. (k/=nz) ) then !cara r
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = 1 * fac  !nodos en z
        L(count,count+1) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = -2 * fac  !nodos en x
        L(count,count-(2*ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==1) ) then !cara u
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count+1) = -2 * fac  !nodos en z
        L(count,count+2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      else if ( (i/=1) .and. (i/=nx) .and. (j/=1) .and. (j/=ny) .and. (k==nz) ) then !cara d
        L(count,count) = -3 * fac !nodo a derivar
        L(count,count-1) = -2 * fac  !nodos en z
        L(count,count-2) = 1 * fac
        L(count,count-nz) = 1 * fac  !nodos en y
        L(count,count+nz) = 1 * fac
        L(count,count-(ny*nz)) = 1 * fac  !nodos en x
        L(count,count+(ny*nz)) = 1 * fac

      end if

    end do
  end do
end do

end subroutine laplaciano2
!*******************************************************************************



end module auxiliar
!*******************************************************************************