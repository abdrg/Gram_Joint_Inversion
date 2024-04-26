module model_ini
! Contiene subrutinas para generar un modelo inicial m0_gv para el caso gravimétrico
! m0_mg para el caso magnético con caracteristicas simples para modelado directo.
! ref_method = 1: ref para calcular contrastes es cte, ref_method = 2: promedio de todas las celdas.
! m_type = 0: homogenea, m_type = 1: cubo centrado, m_type = 2: cubos descentrados -1 y -1, m_type = 3: cubos descentrados +1 y +1, m_type = 4 cubos descentrados -1 y -1, m_type >= 5 figura irregular en profundidad.
! casos m_type 1-5 unicamente programados para num vertices (16+12)*(16+12)*16. Ademas en mg, H = 40,000.

! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
contains

!*******************************************************************************
subroutine m_ini_gv(nx,ny,nz,step,m_type,drho_1d)
! Generar una malla 3D de tamaño nx*ny*nz con valores de densidad rho (g/cm^3).
! Para el nx y ny, 12 celdas son utilizadas para generar espacio pseudo-infinito por la función quasi_inf.
! por lo que la parte de interés del modelo tiene num de celdas: nx-12, ny-12, nz.

use iso_fortran_env , dp => real64
use auxiliar, only: coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, m_type
real, intent(in) :: step
real, intent(out) :: drho_1d(nx*ny*nz)
real :: delta_rho(nx,ny,nz)
real :: xcell(nx), ycell(ny), zcell(nz)
integer :: i, j, k, count

! Establecer valores de densidad rho dentro del volumen.
delta_rho = 0
! Experimento 0, sin constrastes de densidad
if (m_type == 0) then
  delta_rho = 0
  !delta_rho(1:nx,1:ny,1:nz) = 1 !slice en lugar de ciclos for anidados (rendimiento es el mismo)
  !do i = 1,nx,1
  !  do j = 1,ny,1
  !    do k = 1,nz,1
  !      delta_rho(i,j,k) = 2.670
  !    end do
  !  end do
  !end do

! Experimentos tesis Emilia Fregoso, 2010. Modelo de 16x16x16 y paso de 0.005km=5mts.
! Experimento 1, cubo centrado constraste +1 g/cm**3 (valor esquinas de celdas)
else if (m_type == 1) then
  delta_rho(7:10,7:10,7:10) = 1

! Experimento 2, cubos decentrados constraste +1 y -1 g/cm**3
else if (m_type == 2) then
  delta_rho(7:10,3:6,7:10) = 1
  delta_rho(7:10,11:14,7:10) = -1

! Experimento 3, cubos decentrados constraste +1 y +1 g/cm**3
else if (m_type == 3) then
  delta_rho(7:10,3:6,7:10) = 1
  delta_rho(7:10,11:14,7:10) = 1

! Experimento 4, cubos decentrados constraste -1 y -1 g/cm**3
else if (m_type == 4) then
  delta_rho(7:10,3:6,7:10) = -1
  delta_rho(7:10,11:14,7:10) = -1

! Experimento paper Li y Oldenburg, 1998. Modelo de 20x20x10 con paso de 0.05km=50mts.
! Experimento 5, figura irregular en profundidad con constraste +1 k/cm**3
else if (m_type == 5) then
  do i = 1,7
    delta_rho(13-i:17-i,8:13,1+i:1+i) = 1 !Oldenburg
    !delta_rho(8:13,13-i:17-i,1+i:1+i) = 1 !Fregoso
  end do

! Modelo realista. Modelo de 30x20x15 y paso de 2km.
! Experimento 6, figura de
else if (m_type == 6) then
  delta_rho(1:nx,1:ny,1:nz) = 2.6 !riolita
  delta_rho(1:nx,1:ny,nz-8:nz-4) = 2.8 !caliza
  do i = 1,3
    delta_rho(22-i:nx,1:ny,3+i) = 2.8 !caliza bloque superior
  end do
  delta_rho(1:nx,1:ny,nz-3:nz) = 2.4 !arenisca
  do i = 1,11
    delta_rho(23-i:25-i,1:ny,i) = 2.2 !zona falla
  end do
  delta_rho(13:15,6:8,1:3) = 3.0 !basalto
  delta_rho(8:11,13:15,1:3) = 3.0 !basalto

! Ejemplo cuerpos más complejo. Modelo 30x20x15
else if (m_type == 7) then
  delta_rho(4:7,4:8,3:6) = 1 !cubo
  do i = 1,11
    delta_rho(24-i:27-i,3:7,1+i) = 1 !dique
  end do
  delta_rho(3:16,13:17,6:8) = 1 !loza
  delta_rho(24:27,13:17,3:6) = -1 !cavidad


! Modelos para m_apriori (valor en centro de celdas)
else if (m_type == -1) then !mapr cubo
  delta_rho(7:10,7:10,7:10) = 1 !mismo cubo
  !delta_rho(6:11,6:11,6:11) = 0.7 !cubo más grande
  !delta_rho(8:9,8:9,5:12) = 0.7 !dique

else if (m_type == -2) then !mapr cubos
  delta_rho(5:12,1:8,5:12) = 0.9
  delta_rho(5:12,9:ny,5:12) = -0.9
  !delta_rho(5:12,9:ny,5:12) = ref_rho +  0.9 !para test cubos gv y mg caracteristicas contrarias.

else if (m_type == -6) then !mapr realista
  delta_rho(1:nx,1:ny,8:nz) = 0.5 !2da capa
  do i = 1,nz,1
    delta_rho(22:24,1:ny,i) = -0.3 !falla o dique vertical
  end do

else if (m_type == -7) then !mapr dos capas Joya LC
  delta_rho(1:nx,1:ny,1:nz) = 0.8
  delta_rho(20:nx,1:ny,1:4) = 0
  do i = 1,4
    delta_rho(1:19,5+i:14-i,i) = 0
  end do

else if (m_type == -10) then
  delta_rho = 1


end if

! PROBAR UTILIZANDO RESHAPE en lugar de ciclos anidados****
!conversión a arreglo de 1D
count = 0
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      drho_1d(count) = delta_rho(i,j,k)
    end do
  end do
end do

! Generar coordenadas centro de celdas !el signo - en z para que 0 sea superficie ya esta conteplado desde coord_voxel
call coord_voxel(nx,ny,nz,step,xcell,ycell,zcell)


if (m_type == -1 .or. m_type == -2 .or. m_type == -3 .or. m_type == -6 .or. m_type == -7 .or. m_type == -10) then !modelo apriori centro de celdas.
  ! Escribir la posicion nodos en un texto formato x,y,z,rho combinando la subrutina quasi_inf y los arreglos creados previamente.
  OPEN(unit=1,file='output_mapr_gv.csv',status='unknown') !crear y escribir en archivo
  ! Header
  WRITE(1,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'
  ! Body
  count = 0
  do i = 1,nx,1
    do j = 1,ny,1
      do k = 1,nz,1
          count = count + 1
          WRITE(1,'(F10.3,",",F11.3,",",F10.3,",",F8.4)') xcell(i),ycell(j),zcell(k),drho_1d(count)
      end do
    end do
  end do
  CLOSE(unit=1)

else  !modelo inicial centro de celdas.
  ! Escribir la posicion nodos en un texto formato x,y,z,rho combinando la subrutina quasi_inf y los arreglos creados previamente.
  OPEN(unit=1,file='output_m_gv_ini.csv',status='unknown') !crear y escribir en archivo
  ! Header
  WRITE(1,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'
  ! Body
  count = 0
  do i = 1,nx,1
    do j = 1,ny,1
      do k = 1,nz,1
          count = count + 1
          WRITE(1,'(F10.3,",",F11.3,",",F10.3,",",F8.4)') xcell(i),ycell(j),zcell(k),drho_1d(count)
      end do
    end do
  end do
  CLOSE(unit=1)

end if

return
end subroutine m_ini_gv
!*******************************************************************************

subroutine m_ini_mg(nx,ny,nz,step,m_type,H,delta_J)
! Generar una malla 3D de tamaño nx*ny*nz con valores de susceptibilidad magnética Xm (S.I.).
! archivo de salida se guarda en Magnetizacion J (A/m) mediante la relación J=Xm*H.

use iso_fortran_env , dp => real64
use auxiliar, only: coord_voxel

implicit none

integer, intent(in) :: nx, ny, nz, m_type
real, intent(in) :: step, H
real :: delta_Xm(nx,ny,nz),  dXm_1d(nx*ny*nz)
real, intent(out) :: delta_J(nx*ny*nz)
real :: xcell(nx), ycell(ny), zcell(nz)
integer :: i, j, k, count

! Establecer valores de susceptibilidad dentro del volumen.
delta_Xm = 0
! Experimento 0, sin constrastes de susceptibilidad magnetica (S.I)
if (m_type == 0) then
  delta_Xm = 0

! Experimentos tesis Emilia Fregoso, 2010.
! Experimento 1, cubo centrado constraste +1 A/m (considerando H=40,000) (valor esquinas de celdas)
else if (m_type == 1) then
  delta_Xm(7:10,7:10,7:10) = 2.5e-5

! Experimento 2, cubos decentrados constraste +1 y -1 A/m (considerando H=40,000)
else if (m_type == 2) then
  delta_Xm(7:10,3:6,7:10) = 2.5e-5
  delta_Xm(7:10,11:14,7:10) = -2.5e-5

! Experimento 3, cubos decentrados constraste +1 y +1 A/m (considerando H=40,000)
else if (m_type == 3) then
  delta_Xm(7:10,3:6,7:10) = 2.5e-5
  delta_Xm(7:10,11:14,7:10) = 2.5e-5

! Experimento 4, cubos decentrados constraste -1 y -1 A/m (considerando H=40,000)
else if (m_type == 4) then
  delta_Xm(7:10,3:6,7:10) = -2.5e-5
  delta_Xm(7:10,11:14,7:10) = -2.5e-5

! Experimento paper Li y Oldenburg, 1998.
! Experimento 5, figura irregular en profundidad con constraste +1 A/m (considerando H=40,000)
else if (m_type == 5) then
  do i = 1,7
    delta_Xm(13-i:17-i,8:13,1+i:1+i) = 2.5e-5 !Oldemburg
    !delta_Xm(8:13,13-i:17-i,1+i:1+i) = 2.5e-5 !Fregoso
  end do

! Modelo realista. Modelo de 30x20x15 y paso de 2km.
! Experimento 6, figura de
else if (m_type == 6) then
  delta_Xm(1:nx,1:ny,1:nz) = 519e-6 !riolita
  delta_Xm(1:nx,1:ny,nz-8:nz-4) = 250e-6 !caliza
  do i = 1,3
    delta_Xm(22-i:nx,1:ny,3+i) = 250e-6 !caliza bloque superior
  end do
  delta_Xm(1:nx,1:ny,nz-3:nz) = 133e-6 !arenisca
  do i = 1,11
    delta_Xm(23-i:25-i,1:ny,i) = 650e-6 !zona falla
  end do
  delta_Xm(13:15,6:8,1:3) = 773e-6 !basalto
  delta_Xm(8:11,13:15,1:3) = 773e-6 !basalto

! Ejemplo cuerpos más complejo. Modelo 30x20x15
else if (m_type == 7) then
  delta_Xm(4:7,4:8,3:6) = 2.5e-5 !cubo
  do i = 1,11
    delta_Xm(24-i:27-i,3:7,1+i) = 2.5e-5 !dique
  end do
  delta_Xm(3:16,13:17,6:8) = 2.5e-5 !loza
  delta_Xm(24:27,13:17,3:6) = -2.5e-5 !cavidad

  
! Modelos para m_apriori (valor en centro de celdas)
else if (m_type == -1) then !mapr cubo
  delta_Xm(7:10,7:10,7:10) = 2.5e-5 !mismo cubo
  !delta_Xm(6:11,6:11,6:11) = 2e-5 !cubo mas grande
  !delta_Xm(8:9,8:9,5:12) = 2e-5 !dique

else if (m_type == -2) then !mapr cubos
  delta_Xm(5:12,1:8,5:12) = 2.25e-5
  delta_Xm(5:12,9:ny,5:12) = -2.25e-5

else if (m_type == -6) then !mapr realista
  delta_Xm(1:nx,1:ny,8:nz) = -300e-6 !2da capa
  do i = 1,nz,1
    delta_Xm(22:24,1:ny,i) = 200e-6 !falla o dique vertical
  end do

else if (m_type == -7) then !mapr dos capas Joya LC
  delta_Xm(1:nx,1:ny,1:nz) = -4.5e-6 !*44089=-0.2
  delta_Xm(20:nx,1:ny,1:4) = 0
  do i = 1,4
    delta_Xm(1:19,5+i:14-i,i) = 0
  end do

else if (m_type == -10) then
  delta_Xm = 1


end if


!conversión a arreglo de 1D
count = 0
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      dXm_1d(count) = delta_Xm(i,j,k)
    end do
  end do
end do

!conversión a contrastes de magnetización (J=Xm*H)
do i = 1,nx*ny*nz
      delta_J(i) = dXm_1d(i)*H
end do

! Generar coordenadas centro de celdas !el signo - en z para que 0 sea superficie ya esta conteplado desde coord_space
call coord_voxel(nx,ny,nz,step,xcell,ycell,zcell)


if (m_type == -1 .or. m_type == -2 .or. m_type == -3 .or. m_type == -6 .or. m_type == -7.or. m_type == -10) then  !modelo apriori centro de celdas.
  ! Escribir la posicion nodos en un texto formato x,y,z,Xm,Mi,Md combinando la subrutina quasi_inf y los arreglos creados previamente.
  OPEN(unit=2,file='output_mapr_mg.csv',status='unknown') !crear y escribir en archivo
  ! Header
  WRITE(2,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'!,Inclinacion [Deg],Declinacion [Deg]'
  ! Body
  count = 0
  do i = 1,nx,1
    do j = 1,ny,1
      do k = 1,nz,1
          count = count + 1
          WRITE(2,'(F10.3,",",F11.3,",",F10.3,",",F11.4)') xcell(i),ycell(j),zcell(k),delta_J(count)
      end do
    end do
  end do
  CLOSE(unit=2)

else  !modelo inicial centro de celdas.
  ! Escribir la posicion nodos en un texto formato x,y,z,Xm,Mi,Md combinando la subrutina quasi_inf y los arreglos creados previamente.
  OPEN(unit=2,file='output_m_mg_ini.csv',status='unknown') !crear y escribir en archivo
  ! Header
  WRITE(2,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'!,Inclinacion [Deg],Declinacion [Deg]'
  ! Body
  count = 0
  do i = 1,nx,1
    do j = 1,ny,1
      do k = 1,nz,1
          count = count + 1
          WRITE(2,'(F10.3,",",F11.3,",",F10.3,",",F11.4)') xcell(i),ycell(j),zcell(k),delta_J(count)
      end do
    end do
  end do
  CLOSE(unit=2)

end if

return
end subroutine m_ini_mg
!*******************************************************************************


end module model_ini
!*******************************************************************************