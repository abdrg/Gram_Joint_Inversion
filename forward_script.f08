program for_modeling
! Programa que utiliza una serie de subrutinas para hacer modelado directo sintético de métodos potenciales.
! Autor: Abraham Del Razo, IPICyT. May-2023.
! Mail: abraham.delrazo@ipicyt.edu.mx
!***********************************************************************************************************************************

use iso_fortran_env , dp => real64
use auxiliar
use model_ini
use forward_modeling

implicit none

real, allocatable, dimension(:) :: d_gv, m_gv, zobs_gv, d_mg, m_mg, zobs_mg
real, allocatable, dimension(:) :: xcell, ycell, zcell, xobs, yobs, stdevd_gv, stdevd_mg, x, y, z
real(dp), allocatable :: A_gv(:,:), A_mg(:,:)
integer :: m, n, nx, ny, nz, m_ini_type, mapr_type, d_type, a
real :: H, Fi, Fd, step, lvnoise, TIPO_MODELADO
real :: start0, finish0
real(dp), allocatable, dimension(:,:) :: dx, dy, dz, opL

!***********************************************************************************************************************************
! Registrar tiempos de ejecución
OPEN(unit=0,file='output_time.txt',status='unknown')
! Evaluar eficiencia.
call cpu_time(start0)

!***********************************************************************************************************************************
! PARAMETROS DE MODELO DIRECTO SINTETICO

!Dimensiones de malla (o número de parametros). (num celdas = num vertices-1)
!*************************

!Reproducir experimentos Emilia Fregoso, 2010. m_ini_type = 1,2,3,4. 
!nx = 16!+12 !incremento lineal c/extremos logaritmicos
!ny = 16!+12
!nz = 16 !incremento lineal, !***PENDIENTE contemplar incremento log
!step = 5 !mts.  !16 celdas necesarias para modelos de cubos

!Campo Magnético Ambiental en grados.
!H = 40000  !nT necesarios para obtener magnetización=1 en cubos (J=Xm*H).
!Fi = 51 !Emilia Fregoso 2010.
!Fd = 0

!Tipo de modelo inicial para modelado directo sintético.
!m_ini_type = 1
!lvnoise = 2!%
!*************************

!Reproducir experimentos Li y Oldenburg, 1998. m_ini_type = 5. 
!nx = 20
!ny = 20
!nz = 10
!step = 50 !mts

!H = 40000 !nT necesarios para obtener magnetización=1 en cubos (J=Xm*H).
!Fi = 45 !Li y Oldenburg, 1998.
!Fd = 45

!Tipo de modelo inicial para modelado directo sintético.
!m_ini_type = 5
!lvnoise = 2!%
!*************************

!Para experimento 'modelo realista'. m_ini_type=6.
!nx = 30
!ny = 20
!nz = 15
!step = 2000 !2km, similar a resolucion datos satelitales=1.8km.

!H = 44734  !nT  !Datos campo en SLP para experimento 'realista'
!Fi = 51.88 !grados
!Fd = 7.83  !grados

!Tipo de modelo inicial para modelado directo sintético.
!m_ini_type = 6
!lvnoise = 5!%
!*************************

!Para datos satelitales Joja LC, SLP. m_ini_type=6.
!nx = 25
!ny = 25
!nz = 101mts
!step = 200 

!H = 44089  !nT  !Datos campo en Joya LC, SLP
!Fi = 51.13!1995 !50.99 19/10/22 !grados
!Fd = 7.20  !7.31 19/10/22 ! grados
!*************************

!Para pruebas de desempeño.
nx = 30
ny = 20
nz = 15
step = 50 !mts

H = 40000  !nT
Fi = 45 !grados
Fd = 15  !grados

!Tipo de modelo inicial para modelado directo sintético.
m_ini_type = 7
lvnoise = 5!%
!*************************

!Tipo de dato geofísico: d_type = 1(gv), 2(mg).
TIPO_MODELADO = 1 !1,2,3

!***********************************************************************************************************************************
! MODELADO DIRECTO

if (TIPO_MODELADO == 1) then  !Modelado directo de modelo sintetico preprogramado
      
      m = nx*ny
      n = nx*ny*nz
      ! MISMOS nombres de variables se deben utilizar en subrutinas que usen estos arreglos, para poder utilizar es mismo arreglo en la declaracion de variables como allocable y no tener que usar otro espacio de memoria del mismo tamano.
      allocate( A_gv(m,n), d_gv(m), m_gv(n), zobs_gv(m) )
      allocate( A_mg(m,n), d_mg(m), m_mg(n), zobs_mg(m) )

      !Elevación fija de datos. Valido para adquisición aérea con vuelo de altura fija. Vuelo de contorno o terrestres necesario indicar elevacion para cada estación.
      zobs_gv = 0.01 !para gravimetricos a 0mts
      zobs_mg = 0.01 !para aeromagneticos a 300mts

      !GENERANDO DATOS SINTETICOS
      !subroutine forward_gv(nx,ny,nz,step,zobs_gv,m_type,lvnoise,A_gv,d_gv,m_gv)
      !subroutine forward_mg(nx,ny,nz,step,zobs_mg,m_type,lvnoise,H,Fi,Fd,A_mg,d_mg,m_mg)
      call forward_gv(nx,ny,nz,step,zobs_gv,m_ini_type,lvnoise,A_gv,d_gv,m_gv)
      call forward_mg(nx,ny,nz,step,zobs_mg,m_ini_type,lvnoise,H,Fi,Fd,A_mg,d_mg,m_mg)

      a = 0 !instruccion para guardar binarios de operadores A_gv, A_mg, dx, dy, dz, opL
      if (a == 1) then
            ! Archivo salida en binario de opepador A
            !A_gv = 0 !test de que funciona la lectura
            !OPEN(UNIT=2, FILE="output_a_gv.dat", ACTION="read", FORM="unformatted")
            !READ(2) A_gv
            !CLOSE(UNIT=2)
            !print*, 'A_gv=', A_gv(1:3,1:3)

            allocate( dx(n,n),dy(n,n),dz(n,n),opL(n,n) )
            !subroutine partialder(nx,ny,nz,step,dir,pd)
            call partialder(nx,ny,nz,step,'x',dx)
            call partialder(nx,ny,nz,step,'y',dy)
            call partialder(nx,ny,nz,step,'z',dz)
            
            OPEN(UNIT=3, FILE="output_dx.dat", ACTION="write", STATUS="replace", FORM="unformatted")
            WRITE(3) dx
            CLOSE(UNIT=3)
            
            OPEN(UNIT=4, FILE="output_dy.dat", ACTION="write", STATUS="replace", FORM="unformatted")
            WRITE(4) dy
            CLOSE(UNIT=4)
            
            OPEN(UNIT=5, FILE="output_dz.dat", ACTION="write", STATUS="replace", FORM="unformatted")
            WRITE(5) dz
            CLOSE(UNIT=5)

            !subroutine laplaciano(nx,ny,nz,step,L)
            call laplaciano(nx,ny,nz,step,opL)
      
            OPEN(UNIT=6, FILE="output_Lapl.dat", ACTION="write", STATUS="replace", FORM="unformatted")
            WRITE(6) opL
            CLOSE(UNIT=6)
      end if

elseif (TIPO_MODELADO == 2) then !Modelado directo de modelo sintetico leido de archivo
      
      m = nx*ny
      n = nx*ny*nz
      ! MISMOS nombres de variables se deben utilizar en subrutinas que usen estos arreglos, para poder utilizar es mismo arreglo en la declaracion de variables como allocable y no tener que usar otro espacio de memoria del mismo tamano.
      allocate( A_gv(m,n), d_gv(m), m_gv(n), zobs_gv(m) )
      allocate( A_mg(m,n), d_mg(m), m_mg(n), zobs_mg(m) )

      !Elevación fija de datos. Valido para adquisición aérea con vuelo de altura fija. Vuelo de contorno o terrestres necesario indicar elevacion para cada estación.
      zobs_gv = 0.01 !para gravimetricos a 0mts
      zobs_mg = 0.01 !para aeromagneticos a 300mts

      ! INPUT DATOS
      !subroutine readmodel(xcell,ycell,zcell,m_gv,m_mg)
      call readmodel(xcell,ycell,zcell,m_gv,m_mg)

      ! CONSTRUCCION DE OPERADORES A_gv Y A_mg
      !subroutine A_gen1(nx,ny,nz,step,elev,d_type,Fi,Fd,A_op)
      call A_gen1(nx,ny,nz,step,zobs_gv,1,0.,0.,A_gv)
      call A_gen1(nx,ny,nz,step,zobs_mg,2,Fi,Fd,A_mg)

      ! MODELADO DIRECTO
      !subroutine forward_post_inv(nx,ny,nz,step,zobs,d_type,A_mat,m_vec)
      call forward_mod(nx,ny,nz,step,zobs_gv,1,A_gv,m_gv)
      call forward_mod(nx,ny,nz,step,zobs_mg,2,A_mg,m_mg)

elseif (TIPO_MODELADO == 3) then !INCOMPLEO Modelado directo de modelo sintetico a distribución en superficie leidos de archivos
      
      ! INPUT DATOS
      !subroutine readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
      call readdata(xobs,yobs,zobs_gv,zobs_mg,d_gv,d_mg,stdevd_gv,stdevd_mg)
      !subroutine readmodel(xcell,ycell,zcell,m_gv,m_mg)
      call readmodel(xcell,ycell,zcell,m_gv,m_mg)

      m = size(d_gv)
      n = size(m_gv)
      allocate( A_gv(m,n), A_mg(m,n) )

      ! CONSTRUCCION DE OPERADORES A_gv Y A_mg
      allocate ( x(16), y(16), z(16) )
      x = (/2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5/)
      y = (/2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5/)
      z = (/-2.5,-7.5,-12.5,-17.5,-22.5,-27.5,-32.5,-37.5,-42.5,-47.5,-52.5,-57.5,-62.5,-67.5,-72.5,-77.5/)

      !subroutine A_gen3(xobs,yobs,zobs_gv,zobs_mg,x,y,z,Fi,Fd,A_gv,A_mg)
      call A_gen3(xobs,yobs,zobs_gv,zobs_mg,x,y,z,Fi,Fd,A_gv,A_mg)

      ! MODELADO DIRECTO
      ! Escribir archivo con datos para graficar anomalia calculada.
      OPEN(unit=7,file='output_gv_anomaly.csv',status='unknown') !crear y escribir en archivo
      WRITE(7,*) 'x[m],y[m],z[m],gv[mGal]'
      OPEN(unit=8,file='output_mg_anomaly.csv',status='unknown') !crear y escribir en archivo
      WRITE(8,*) 'x[m],y[m],z[m],mg[nT]'

      d_gv = MATMUL(A_gv,m_gv)
      d_mg = MATMUL(A_mg,m_mg)

      do a = 1,size(d_gv),1
            !Cuerpo del archivo de salida.
            WRITE(7,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(a),yobs(a),zobs_gv(a),d_gv(a)*1000
            WRITE(8,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(a),yobs(a),zobs_mg(a),d_mg(a)
      end do
      CLOSE(unit=7)
      CLOSE(unit=8)

else
      PRINT*, 'El tipo de modelado directo fue seleccionado erroneamente, porfavor elige un TIPO_MODELADO apropiado'
end if
!***********************************************************************************************************************************

! Evaluar eficiencia.
call cpu_time(finish0)
WRITE(0,*) 'Total execution time:', (finish0-start0), 'seg =', (finish0-start0)/60, 'min'
WRITE(0,*) '~(^-^)~  Its over dude'
CLOSE(unit=0)

!deallocate(A_gv,d_gv,m_gv,zobs_gv)!,stdevd_gv)
!deallocate(A_mg,d_mg,m_mg,zobs_mg)!,stdevd_mg)
!deallocate(xcell,ycell,zcell,xobs,yobs)

end program for_modeling
!***********************************************************************************************************************************
! Para evitar escribir compilar cada modulo correr el script bash desde directorio de proyecto en terminal:
!./run_modelado_host.sh

