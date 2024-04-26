module forward_modeling
! Contiene subrutinas para hacer modelado directo a partir de modelo sintético creado por model_ini.mod
! para el caso gravimétricos forward_gv y magnéticos forward_mg o
! a partir de parametros de modelo invertidos m y matriz A  para el caso de forward_post_inv
!
! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
contains


!*******************************************************************************
subroutine forward_gv(nx,ny,nz,step,zobs_gv,m_type,lvnoise,A_gv,d_gv,m_gv)
! Programa que calcula una anomalía gravimétrica a partir de modelado directo de un modelo sintético,
! mediante un grid 3D de prismas regulares con distribucion de densidades g/cm**3.

use iso_fortran_env, dp => real64
use operators
use model_ini
use auxiliar, only: gaussian_noise, std_dev_model, coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, m_type
real, intent(in) :: step, lvnoise
real, allocatable, intent(in) :: zobs_gv(:)
real(dp), allocatable :: A_gv(:,:)
real, allocatable, dimension(:) :: d_gv, m_gv
integer :: a, b, i, j, k, count_M, count_N, count
real(dp) :: gz
real :: x(nx+1), y(ny+1), z(nz+1), xobs(nx), yobs(ny), z2(nz) !z2 variable sin uso pero necesaria para subrutina voxel
real :: data_gv(nx*ny)
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

!subroutine m0_gv(nx,ny,nz,step,m_type,delta_rho_1d)
call m_ini_gv(nx,ny,nz,step,m_type,m_gv)

! Crear puntos de obseración a la mitad de cada celda.
!subroutine coord_voxel(nx,ny,nz,block_size,xcell,ycell,zcell)
call coord_voxel(nx,ny,nz,step,xobs,yobs,z2)

! Celdas en x
do i = 1,nx+1
      x(i) = (i-1)*step
end do
! Celdas en y
do j = 1,ny+1
      y(j) = (j-1)*step
end do
! Celdas en z
do k = 1,nz+1
      z(k) = -(k-1)*step
end do

! Puntos de observacion podrían se mayores a nx o ny lo que cambiaria el tamaño de la matriz de datos,
! si esta es > matriz de parametros m hay que cambiar dimensiones de matrices en operaciones LAPACK porque este programa esta pensado para un problema subdeterminado.

! Calculo del operador considera la posición de los parametros en el centro de la celda, aunque por practicidad estan codificados en el vector como si estubieran en el vertice superior izq.
count_M = 0
do a = 1,nx,1 !x del punto de observación
  do b = 1,ny,1 !y del punto de observación
      count_M = count_M + 1

      gz = 0
      count_N = 0
      do i = 1,nx,1
        do j = 1,ny,1
          do k = 1,nz,1
            count_N = count_N + 1

            !subroutine op_gv(x0,y0,z0,x1,y1,z1,x2,y2,z2,operator)
            call op_gv(xobs(a),yobs(b),zobs_gv(count_M),x(i),y(j),z(k),x(i+1),y(j+1),z(k+1),A_gv(count_M,count_N))
            
            !calcular datos a partir de operador directo y parametros modelo
            gz = gz + ( A_gv(count_M,count_N)*m_gv(count_N) )

          end do
        end do
      end do

      !constriur vector columna para modelado inverso
      data_gv(count_M) = gz

  end do
end do

! Archivo salida en binario de opepador A
OPEN(UNIT=2, FILE="output_a_gv.dat", ACTION="write", STATUS="replace", FORM="unformatted")
WRITE(2) A_gv
CLOSE(UNIT=2)

! Escribir archivo con datos para graficar anomalia calculada.
OPEN(unit=3,file='output_d_gv.csv',status='unknown') !crear y escribir en archivo
WRITE(3,*) 'x[m],y[m],z[m],gv[mGal],Std_dev[mGal]' !encabezado

! Datos con ruido
call gaussian_noise(nx*ny,lvnoise,data_gv,d_gv)
OPEN(unit=4,file='output_d_gv_noised.csv',status='unknown')
WRITE(4,*) 'x[m],y[m],z[m],gv[mGal],Std_dev[mGal]'

! Cuerpo del archivo de salida.
count = 0
do a = 1,nx,1
  do b = 1,ny,1
    count = count + 1
      WRITE(3,'(F10.3,",",F11.3,",",F10.3,",",E11.4,",",F7.3)') xobs(a),yobs(b),zobs_gv(count),data_gv(count),0.1
      WRITE(4,'(F10.3,",",F11.3,",",F10.3,",",E11.4,",",F7.3)') xobs(a),yobs(b),zobs_gv(count),d_gv(count),0.1
  end do
end do

CLOSE(unit=3)
CLOSE(unit=4)


!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
WRITE(0,*) 'Execution time: forward_gv_m0 =', (finish-start), 'seg =', (finish-start)/60, 'min'
!CLOSE(unit=0)

return
end subroutine forward_gv

!*******************************************************************************
subroutine forward_mg(nx,ny,nz,step,zobs_mg,m_type,lvnoise,H,Fi,Fd,A_mg,d_mg,m_mg)
! Subrutina que calcula una anomalía magnética a partir de modelado directo de un modelo sintético,
! mediante un grid 3D de prismas regulares con distribucion de magnetización A/m.

use iso_fortran_env, dp => real64
use model_ini
use operators
use auxiliar, only: gaussian_noise, std_dev_model, coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, m_type
real, intent(in) :: step, lvnoise, H, Fi, Fd
real, allocatable, intent(in) :: zobs_mg(:)
real(dp), allocatable :: A_mg(:,:)
real, allocatable, dimension(:) :: d_mg, m_mg
integer :: a, b, i, j, k, count_M, count_N, count
real(dp) :: T
!real(dp) :: sensibility_mg_t, sensibility_mg_b, T_top, T_bottom, T
real :: x(nx+1), y(ny+1), z(nz+1), xobs(nx), yobs(ny), z2(nz) !z2 variable sin uso
real :: data_mg(nx*ny)
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

!subroutine m0_mg(nx,ny,nz,step,m_type,H,delta_J)
call m_ini_mg(nx,ny,nz,step,m_type,H,m_mg)

! Crear puntos de obseración a la mitad de cada celda.
!coord_voxel(nx,ny,nz,block_size,xcell,ycell,zcell)
call coord_voxel(nx,ny,nz,step,xobs,yobs,z2)

! Celdas en x
do i = 1,nx+1
      x(i) = (i-1)*step
end do
! Celdas en y
do j = 1,ny+1
      y(j) = (j-1)*step
end do
! Celdas en z
do k = 1,nz+1
      z(k) = -(k-1)*step
end do

! Puntos de observacion podrían se mayores a nx o ny, lo que cambiaria el tamaño de la matriz de datos,
! si esta es > matriz de parametros m hay que cambiar dimensiones de matrices en operaciones LAPACK porque este programa esta pensado para un problema subdeterminado.

! Calculo del operador considera la posición de los parametros en el centro de la celda, aunque por practicidad estan codificados en el vector como si estubieran en el vertice superior izq.
count_M = 0
do a = 1,nx,1 !x del punto de observación
  do b = 1,ny,1 !y del punto de observación
      count_M = count_M + 1

      T = 0
      count_N = 0
      do i = 1,nx,1
        do j = 1,ny,1
          do k = 1,nz,1
            count_N = count_N + 1

            !subroutine op_mg(x0,y0,z0,x1,y1,z1,x2,y2,z2,Fi,Fd,operator)
            call op_mg(xobs(a),yobs(b),zobs_mg(count_M),x(i),y(j),z(k),x(i+1),y(j+1),z(k+1),Fi,Fd,A_mg(count_M,count_N))

            !calcular datos a partir de operador directo y parametros modelo
            T = T + ( A_mg(count_M,count_N)*m_mg(count_N) )

          ! Calcular anomalia magnetica por Battacharya.
            !subroutine op_mg2(x0,y0,z0,x1,y1,x2,y2,z,Fi,Fd,operator)
          !  call op_mg2(xobs(a),yobs(b),c,x(i),y(j),x(i+1),y(j+1),z(k),Fi,Fd,sensibility_mg_t) !en ztop
          !  T_top = sensibility_mg_t * m_mg(count_N)

          !  call op_mg2(xobs(a),yobs(b),c,x(i),y(j),x(i+1),y(j+1),z(k+1),Fi,Fd,sensibility_mg_b) !en zbottom
          !  T_bottom = sensibility_mg_b * (-1*m_mg(count_N))

            !constriur matriz ixj y vectores columna para modelado inverso
          !  A_mg(count_M,count_N) = sensibility_mg_t - sensibility_mg_b
            !calcular datos a partir de operador directo y parametros modelo
          !  T = T + (T_top-T_bottom)

          end do
        end do
      end do

      !constriur vector columna para modelado inverso
      data_mg(count_M) = T

  end do
end do

! Archivo salida en binario de opepador A
OPEN(UNIT=2, FILE="output_a_mg.dat", ACTION="write", STATUS="replace", FORM="unformatted")
WRITE(2) A_mg
CLOSE(UNIT=2)

! Escribir archivo con datos para graficar anomalia calculada.
OPEN(unit=5,file='output_d_mg.csv',status='unknown') !crear y escribir en archivo
WRITE(5,*) 'x[m],y[m],z[m],mg[nT],Std_dev[nT]' !encabezado

! Datos con ruido
call gaussian_noise(nx*ny,lvnoise,data_mg,d_mg)
OPEN(unit=6,file='output_d_mg_noised.csv',status='unknown')
WRITE(6,*) 'x[m],y[m],z[m],mg[nT],Std_dev[nT]'

! Cuerpo del archivo de salida.
count = 0
do a = 1,nx,1
  do b = 1,ny,1
    count = count + 1

      WRITE(5,'(F10.3,",",F11.3,",",F10.3,",",E11.4,",",F7.3)') xobs(a),yobs(b),zobs_mg(count),data_mg(count),0.1
      WRITE(6,'(F10.3,",",F11.3,",",F10.3,",",E11.4,",",F7.3)') xobs(a),yobs(b),zobs_mg(count),d_mg(count),0.1
  end do
end do

CLOSE(unit=5)
CLOSE(unit=6)


!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: forward_mg_m0 =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

return
end subroutine forward_mg


!*******************************************************************************
subroutine forward_mod(nx,ny,nz,step,zobs,d_type,A_mat,m_vec)
! Programa que calcula una anomalía gravimétrica/magnética a partir de una matriz u operador de sensibility y un vector de parámetros modelo,
! mediante un grid 3D de prismas regulares con distribucion de densidad gr/cm**3 / magnetización A/m.

use iso_fortran_env, dp => real64
use auxiliar, only: coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, d_type
real, intent(in) :: step, zobs(nx*ny)
real, intent(in) :: m_vec(nx*ny*nz)
real(dp), intent(in) :: A_mat(nx*ny,nx*ny*nz)
real :: d_pred(nx*ny)
integer :: a, b, i, j, k, count
real :: xobs(nx), yobs(ny), z2(nz) !z2 variable sin uso
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

!generar coordenadas de celdas
call coord_voxel(nx,ny,nz,step,xobs,yobs,z2)

! Escribir archivo con datos para graficar anomalia calculada.
if (d_type == 1) then
      OPEN(unit=7,file='output_anomaly_gv.csv',status='unknown') !crear y escribir en archivo
      WRITE(7,*) 'x[m],y[m],z[m],gv[mGal]'
else if (d_type == 2) then
      OPEN(unit=7,file='output_anomaly_mg.csv',status='unknown') !crear y escribir en archivo
      WRITE(7,*) 'x[m],y[m],z[m],mg[nT]'
end if

d_pred = MATMUL(A_mat,m_vec)

count = 0
do a = 1,nx,1 !x del punto de observación.
      do b = 1,ny,1 !y del punto de observación.
            count = count + 1
            !Cuerpo del archivo de salida.
            WRITE(7,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(a),yobs(b),zobs(count),d_pred(count)
      end do
end do

CLOSE(unit=7)


!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: forward_post_inv =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

return
end subroutine forward_mod
!*******************************************************************************


subroutine A_gen1(nx,ny,nz,step,elev,d_type,Fi,Fd,A_op)
!Subrutina que genera el operador A para datos gravimétricos o magneticos a partir de las dimensiones del dominio de inversión.

use iso_fortran_env, dp => real64
use operators
use auxiliar, only: coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, d_type
real, intent(in) :: step, Fi, Fd,  elev(nx*ny)
integer :: a, b, i, j, k, count_M, count_N, count
real :: xobs(nx), yobs(ny), z2(nz) !z2 variable sin uso
real :: x(nx+1), y(ny+1), z(nz+1)
real(dp), intent(out) :: A_op(nx*ny,nx*ny*nz)
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

!generar coordenadas de celdas
call coord_voxel(nx,ny,nz,step,xobs,yobs,z2)

! Crear esquinas de cada celda.
! Celdas en x
do i = 1,nx+1
      x(i) = (i-1)*step
end do
! Celdas en y
do j = 1,ny+1
      y(j) = (j-1)*step
end do
! Celdas en z
do k = 1,nz+1
      z(k) = -(k-1)*step
end do


! Calculo del operador considera la posición de los parametros en el centro de la celda, aunque por practicidad estan codificados en el vector como si estubieran en el vertice superior izq.
count_M = 0
do a = 1,nx,1 !x del punto de observación
      do b = 1,ny,1 !y del punto de observación
            count_M = count_M + 1

            count_N = 0
            do i = 1,nx,1
                  do j = 1,ny,1
                        do k = 1,nz,1
                              count_N = count_N + 1

                              if (d_type == 1) then
                                    !subroutine op_gv(x0,y0,z0,x1,y1,z1,x2,y2,z2,operator)
                                    call op_gv(xobs(a),yobs(b),elev(count_M),x(i),y(j),z(k), &
                                               x(i+1),y(j+1),z(k+1),A_op(count_M,count_N))
                              else if (d_type == 2) then
                                    !subroutine op_mg(x0,y0,z0,x1,y1,z1,x2,y2,z2,Fi,Fd,operator)
                                    call op_mg(xobs(a),yobs(b),elev(count_M),x(i),y(j),z(k), &
                                               x(i+1),y(j+1),z(k+1),Fi,Fd,A_op(count_M,count_N))
                              end if

                        end do
                  end do
            end do

      end do
end do

!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: Agenerator =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

end subroutine A_gen1
!*******************************************************************************


subroutine A_gen2(xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,Fi,Fd,A_gv,A_mg)
!Subrutina que genera el operador A para datos gravimétricos o magneticos a partir de datos observados y modelos leidos de archivo.

use iso_fortran_env, dp => real64
use operators

implicit none

real, intent(in) :: Fi, Fd
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real(dp), allocatable, dimension(:,:) :: A_gv, A_mg
real, allocatable, dimension(:) :: x, y, z
integer :: i, a, b, c, count, size_x, size_y, size_z!,n1, n2
real :: x1, y1, z1, x2, y2, z2
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)


! Convertir arreglo de celdas en coordenadas de esquinas de celdas.
x(1) = xcell(1)
count = 1
do i = 2,size(xcell),1
      if (xcell(i) /= x(count)) then
            count = count + 1
            x(count) = xcell(i) 
      end if
end do
!INCOMPLETO, falta recostruir arreglos x,y,z (centros de celda) a partir de arreglos xcell,ycell,zcell (centros de celda en ordenados para graficar)

! Calculo del operador considera la posición de los parametros en el centro de la celda, aunque por practicidad estan codificados en el vector como si estubieran en el vertice superior izq.
do i = 1,size(xobs),1

      count = 0
      do a = 1,size(x),1
            if (a < size(x) ) then
                  size_x = x(a+1)-x(a)
            elseif (a == size(x) ) then
                  size_x = x(a)-x(a-1)
            end if
            
            do b = 1,size(y),1
                  if (b < size(y) ) then
                        size_y = y(b-1)-y(b)
                  elseif (b == size(y) ) then
                        size_y = y(b)-y(b-1)
                  end if
                  
                  do c = 1, size(z),1
                        if (c < size(z) ) then
                              size_z = z(c+1)-z(c)
                        elseif (c == size(z) ) then
                              size_z = z(c)-z(c-1)
                        end if

                        !Calcular esquinas a partir de centro de celda
                        x1 = x(a)-(size_x/2)
                        x2 = x(a)+(size_x/2)
                        y1 = y(b)-(size_y/2)
                        y2 = y(b)+(size_y/2)
                        z1 = z(c)-(size_z/2)
                        z2 = z(c)+(size_z/2)
                        !print*, z1,z2
                        
                        count = count + 1
                        !subroutine op_gv(x0,y0,z0,x1,y1,z1,x2,y2,z2,operator)
                        call op_gv(xobs(i),yobs(i),zobs_gv(i),x1,y1,z1,x2,y2,z2,A_gv(i,count))
                        !subroutine op_mg(x0,y0,z0,x1,y1,z1,x2,y2,z2,Fi,Fd,operator)
                        call op_mg(xobs(i),yobs(i),zobs_gv(i),x1,y1,z1,x2,y2,z2,Fi,Fd,A_mg(i,count))

                  end do
            end do
      end do
      
end do

!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: Agenerator =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

end subroutine A_gen2
!*******************************************************************************


subroutine A_gen3(xobs,yobs,zobs_gv,zobs_mg,x,y,z,Fi,Fd,A_gv,A_mg)
!Subrutina que genera el operador A para datos gravimétricos o magneticos a partir de datos observados y modelo generado a partir de dominio de datos es superficie.

use iso_fortran_env, dp => real64
use operators

implicit none

real, intent(in) :: Fi, Fd
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg, x, y, z
real(dp), allocatable, dimension(:,:) :: A_gv, A_mg
integer :: i, a, b, c, count, size_x, size_y, size_z
real :: x1, y1, z1, x2, y2, z2
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)

! Calculo del operador considera la posición de los parametros en el centro de la celda, aunque por practicidad estan codificados en el vector como si estubieran en el vertice superior izq.
do i = 1,size(xobs),1

      count = 0
      do a = 1,size(x),1
            if (a < size(x) ) then
                  size_x = x(a+1)-x(a)
            elseif (a == size(x) ) then
                  size_x = x(a)-x(a-1)
            end if
            
            do b = 1,size(y),1
                  if (b < size(y) ) then
                        size_y = y(b+1)-y(b)
                  elseif (b == size(y) ) then
                        size_y = y(b)-y(b-1)
                  end if
                  
                  do c = 1, size(z),1
                        if (c < size(z) ) then
                              size_z = z(c)-z(c+1)
                        elseif (c == size(z) ) then
                              size_z = z(c-1)-z(c)
                        end if

                        !Calcular esquinas a partir de centro de celda
                        x1 = x(a)-(size_x/2)
                        x2 = x(a)+(size_x/2)
                        y1 = y(b)-(size_y/2)
                        y2 = y(b)+(size_y/2)
                        z1 = z(c)+(size_z/2)
                        z2 = z(c)-(size_z/2)
                        !print*, z1,z2
                        
                        count = count + 1
                        !subroutine op_gv(x0,y0,z0,x1,y1,z1,x2,y2,z2,operator)
                        call op_gv(xobs(i),yobs(i),zobs_gv(i),x1,y1,z1,x2,y2,z2,A_gv(i,count))
                        !subroutine op_mg(x0,y0,z0,x1,y1,z1,x2,y2,z2,Fi,Fd,operator)
                        call op_mg(xobs(i),yobs(i),zobs_gv(i),x1,y1,z1,x2,y2,z2,Fi,Fd,A_mg(i,count))

                  end do
            end do
      end do
      
end do

!Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
write(0,*) 'Execution time: Agenerator =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

end subroutine A_gen3
!*******************************************************************************

end module forward_modeling
!*******************************************************************************