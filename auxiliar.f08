module auxiliar
! Contains auxiliar subroutines to perform joint inversion of gravimetric and magnetic data in parallel using Fortran.
! Autor: Abraham Del Razo, IPICyT. Last Update: August 2023.
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************

contains

!*******************************************************************************
subroutine readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
! Stores in memory the file content from two files (input_d_gv.csv and input_d_mg.csv)

use iso_fortran_env, dp => real64
implicit none

real, allocatable, dimension(:) :: d_gv, d_mg, stdevd_gv, stdevd_mg
real, allocatable, dimension(:) :: xobs, yobs, zobs_gv, zobs_mg
character*256 :: CTMP
integer :: i = 0, IERR = 0, rows = 0


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

!deallocate(xobs,yobs,zobs_gv,zobs_mg)

end subroutine readdata
!*******************************************************************************

subroutine readmodel(xcell,ycell,zcell,m_gv,m_mg)
! Stores in memory the file content from two files (input_m_gv.csv and input_m_mg.csv)

use iso_fortran_env, dp=>real64
implicit none

real, allocatable, intent(out), dimension(:) :: m_gv, m_mg
real, allocatable, intent(out), dimension(:) :: xcell, ycell, zcell !coordenadas de cada punto anomalia 2d, no se usarán para nada pero es necesario leerlas.
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

allocate( xcell(rows), ycell(rows), zcell(rows), m_gv(rows), m_mg(rows) )

! Read the file content
rewind(1)
rewind(2)
READ(1, * )   !! skip the header
READ(2, * )

do i = 1,rows,1
  READ(1,*) xcell(i), ycell(i), zcell(i), m_gv(i)
  READ(2,*) xcell(i), ycell(i), zcell(i), m_mg(i)
end do

CLOSE(unit=1)
CLOSE(unit=2)

!deallocate(xcell,ycell,zcell)

end subroutine readmodel
!*******************************************************************************

subroutine model1D(nx,ny,nz,type,model)
! Generates a linearized 3-D property distribution.

use iso_fortran_env, dp=>real64
implicit none

integer, intent(in) :: nx, ny, nz, type
real, intent(out) :: model(nx*ny*nz)
real :: model_3D(nx,ny,nz)
integer :: i, j, k, count

! Standard Deviation
if (type == 1) then
      model_3D = 1 !all cells with same value
      model_3D(1:nx,1:ny,1) = 0.3 !top cells /To avoid anomaly propagation from surface.
      !model_3D(1:nx,1:ny,nz) = 0.1 !bottom cells /To avoid anomaly propagation to bottom.
      !model_3D(1,1:ny,1:nz) = 0.0001 !lateral cells
      !model_3D(nx,1:ny,1:nz) = 0.0001
      !model_3D(1:nx,1,1:nz) = 0.001
      !model_3D(1:nx,ny,1:nz) = 0.001

! Physical property contrast
else if (type == 2) then
      model_3D = 0
      model_3D(14,10,2:4) = 1 !well_1
      model_3D(8,10,6:8) = 1 !well_2

      !do i = 1,7
            !model_3D(13-i:17-i,8:13,1+i:1+i) = 1 original synthetic model
      !enddo
end if

!linearize 3D distribution
count = 0
do i = 1,nx
  do j = 1,ny
    do k = 1,nz
      count = count + 1
      model(count) = model_3D(i,j,k)
    end do
  end do
end do


return
end subroutine model1D
!*******************************************************************************

end module auxiliar
