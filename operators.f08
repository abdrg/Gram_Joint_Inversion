module operators
! Contiene subrutinas para calcular forward lineal operator para datos gravimétricos op_gv.
! y de datos magnéticos op_mg.

! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
contains

!*******************************************************************************
subroutine op_gv(x0,y0,z0,x1,y1,z1,x2,y2,z2,operator)
!Subrutina que calcula la Atraccion Gravitacional Vertical de un prisma rectangular.
!Los lados del prisma son paralelos a los ejes x,y,z y el eje z es positivo hacia abajo.

!Utilizando el desarrollo de Plouff (1976).

!Parametros de entrada:
  !El punto de observacion es  (x0,y0,z0).
  !El prisma se extiende desde x1 a x2, desde y1 a y2 y desde z1 a z2.
  !Todas las distancias estan en mts, la densidad del prisma es rho-rho_base, rho en g/cm^3.
!Parametros de salida:
  !Operador lineal o funcion de Green para un par punto de observacion-prisma rectangular.
  !Para obtener atraccion gravitacional vertical falta multiplicar por parametro de modelo rho.

! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************

use iso_fortran_env , dp => real64

implicit none

real, intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2!, rho, rho_ref
real(dp), intent (out) :: operator
integer :: i, j, k
real, dimension(2) :: x, y, z, signo
real(dp) :: Gamma, SI2mGal, kgm2grcm, arg1, arg2, arg3, r_ijk, mu_ijk, suma!, km2m!, delta_rho

signo = (/-1,1/) !declara un arreglo de 2 elementos de forma implicita
Gamma = 6.67e-11 !cte de proporcionalidad en m^3/(kg*s^2)
SI2mGal = 1e5 !conversion de unidades 1mGal=1e-5 m/s^2
kgm2grcm = 1e3 !conversion de kg/m^3 a gr/cm^3
!km2m = 1e3 !conversion de unidades (entradaba en km y para mGal usa mts)

!delta_rho = rho-rho_ref  !densidad anómala, por eso delta calcula el contraste

x(1) = x0-x1
y(1) = y0-y1
z(1) = z0-z1
x(2) = x0-x2
y(2) = y0-y2
z(2) = z0-z2
!x(1) = abs(x0-x1)
!y(1) = abs(y0-y1)
!z(1) = abs(z0-z1)
!x(2) = abs(x0-x2)
!y(2) = abs(y0-y2)
!z(2) = abs(z0-z2)

suma = 0

do i = 1,2
  do j = 1,2
    do k = 1,2

      r_ijk = sqrt( x(i)**2 + y(j)**2 + z(k)**2 )
      mu_ijk = signo(i)*signo(j)*signo(k)

      arg1 = (x(i)*y(j)) / (z(k)*r_ijk) !+1e-6
      arg2 = r_ijk + y(j) !+1e-6
      arg3 = r_ijk + x(i) !+1e-6

      suma = suma + mu_ijk*( z(k)*atan(arg1) - x(i)*log(arg2) - y(j)*log(arg3) )
      !funciones intrinsecas alternatibas de fortran
      !suma = suma + mu_ijk*( z(k)*atan2(x(i)*y(j),z(k)*r_ijk) - x(i)*log(arg2) - y(j)*log(arg3) )

    end do
  end do
end do

operator = suma*Gamma*SI2mGal*kgm2grcm!*km2m!*delta_rho

return
end subroutine op_gv
!*******************************************************************************

subroutine op_mg(x0,y0,z0,x1,y1,z1,x2,y2,z2,Fi,Fd,operator)
!Subrutina que calcula el Campo Total de un prisma rectangular.
!Los lados del prisma son paralelos a los ejes x,y,z y el eje z es positivo hacia abajo y se extiende hasta el infinito.

!Utilizando el desarrollo de Gallardo (1997).
!Requiere subrutina cos_dir.

!Parametros de entrada:
  !Punto de observacion es (x0,y0,z0).
  !El prima se extiende desde x1 a x2, desde y1 a y2 y desde z1 al z2.
  !Magnetizacion definida por su inclinacion Ji, declinacion Jd e intensidad J que es producto de susceptibilidad en S.I (suscept)*intensidad de campo en nT (H).
  !Campo Ambiente definido por su intensidad H, inclinacion Fi y declinacion Fd. Eje x tiene declinacion theta respecto al norte verdadero. Las unidades de distancia son irrelevantes pero deben ser consistentes.
  !Angulos estan de grados, con inclinación (0-90) positiva bajo la horizontal y declinacion (0-360) positiva hacia el Este. Magnetizacion en A/m.
  !Todas las distancias estan en mts. Para simplificar se supone que no hay magnetización remanente y que la magnetización en el material es paralela al campo aplicado, por lo que Ji=Fi y Jd=Fd.

!Parametros de salida:
  !Operador lineal o funcion de Green para un par punto de observacion-prisma rectangular.
  !Para obtener Anomalia de Campo Total T, en nT hace falta multiplicar por parametro de modelo Magnetización (susceptibilidad*Intensidad de campo)

! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx 
!*******************************************************************************

use iso_fortran_env , dp => real64

implicit none

real, intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, Fi, Fd!, Ji, Jd!, H, Xm, J_ref
real(dp), intent(out) :: operator
integer :: i, j, k
real, dimension(2) :: x, y, z, signo
real(dp) :: Cm, T2nT, m2km, a, b, g, arg1, arg2, arg3, arg4, arg5, arg6, &
                    arg1b, arg2b, arg3b, arg4b, arg5b, arg6b, suma, r_ijk, mu_ijk!, J, delta_J

signo = (/-1,1/) !declarar un arreglo de 2 elementos
Cm = 1e-7 !cte de proporcionalidad Mu_0/4*pi
T2nT = 1e9 !conversion de unidades

!J = Xm*H  !Magnetización A/m. Susceptibilidad en S.I e Intensidad de Campo Magnético en nT
!delta_J = J - J_ref  !magnetizacion anómala, por eso delta calcula el contraste.

!x(1) = x0-x1
!y(1) = y0-y1
!z(1) = z0-z1
!x(2) = x0-x2
!y(2) = y0-y2
!z(2) = z0-z2
x(1) = x1-x0
y(1) = y1-y0
z(1) = z1-z0
x(2) = x2-x0
y(2) = y2-y0
z(2) = z2-z0

call cos_dir(Fi,Fd,a,b,g)
!call cos_dir(Mi,Md,a2,b2,g2) !Se considera incl y decl de magnetización en cuerpos paralelos a la campo ambiental.

arg1 = a**2
arg2 = b**2
arg3 = g**2
arg4 = 2*a*b
arg5 = 2*a*g
arg6 = 2*b*g

suma = 0

do i = 1,2
  do j = 1,2
    do k = 1,2

      r_ijk = sqrt( x(i)**2 + y(j)**2 + z(k)**2 )
      mu_ijk = signo(i)*signo(j)*signo(k)

      arg1b = (y(j)*z(k)) / (x(i)*r_ijk)
      arg2b = (x(i)*z(k)) / (y(j)*r_ijk)
      arg3b = (x(i)*y(j)) / (z(k)*r_ijk)
      arg4b = r_ijk + z(k)
      arg5b = r_ijk + y(j)
      arg6b = r_ijk + x(i)

      suma = suma + mu_ijk*( -arg1*atan(arg1b) - arg2*atan(arg2b) - arg3*atan(arg3b) &
                              + arg4*log(arg4b) + arg5*log(arg5b) + arg6*log(arg6b) )
    end do
  end do
end do

operator = suma*Cm*T2nT!*delta_J

return
end subroutine op_mg
!*******************************************************************************

subroutine op_mg2(x0,y0,z0,x1,y1,x2,y2,z1,Fi,Fd,operator)
!Subrutina que calcula el Campo Total T de un prisma rectangular.
!Los lados del prisma son paralelos a los ejes x,y,z y el eje z es positivo hacia abajo y se extiende hasta el infinito.

!Utilizando el desarrollo de Bhattacharyya (1964).
!Requiere subrutina cos_dir.

!Parametros de entrada:
  !Punto de observacion es (x0,y0,z0).
  !El prima se extiende desde x1 a x2, desde y1 a y2 y desde z1 al infinito.
  !Magnetizacion definida por su inclinacion Ji, declinacion Jd e intensidad J que es producto de susceptibilidad en S.I (suscept)*intensidad de campo en nT (H).
  !Campo Ambiente definido por su intensidad H, inclinacion Fi y declinacion Fd. Eje x tiene declinacion theta respecto al norte verdadero. Las unidades de distancia son irrelevantes pero deben ser consistentes.
  !Angulos estan de grados, con inclinación (0-90) positiva bajo la horizontal y declinacion (0-360) positiva hacia el Este. Magnetizacion en A/m.
  !Todas las distancias estan en mts.
  !Para simplificar se supone que no hay magnetización remanente y que la magnetización en el material es paralela al campo aplicado, por lo que Ji=Fi y Jd=Fd.

!Parametros de salida:
  !Operador lineal o funcion de Green para un par punto de observacion-prisma rectangular.
  !Para obtener Anomalia de Campo Total t, en nT hace falta multiplicar por parametro de modelo Magnetización (susceptibilidad*Intensidad de campo)

! Autor: Abraham Del Razo, IPICyT. Feb-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
use iso_fortran_env , dp => real64

implicit none

real, intent(in) :: x0, y0, z0, x1, y1, x2, y2, z1, Fi, Fd!, Ji, Jd!, H, Xm, J_ref
real(dp), intent(out) :: operator
integer :: i, j
real, dimension(2) :: x, y, signo
real(dp) :: Cm, T2nT, Mx, My, Mz, Fx, Fy, Fz, arg1, arg2, arg3, arg4, arg5, arg6, z, suma, mu_ij, &
                    arg1b, arg2b, arg3b, arg4b, arg5b, arg6b, zz, xx, rr, r, xy!, J, delta_J

Cm = 1e-7 !cte de proporcionalidad Mu_0/4*pi
T2nT = 1e9 !conversion de unidades
signo = (/-1,1/) !declarar un arreglo de 2 elementos

!J = Xm*H  !Magnetización A/m. Susceptibilidad en S.I e Intensidad de Campo Magnético en nT
!delta_J = J - J_ref  !magnetizacion anómala, por eso delta calcula el contraste

x(1) = x1-x0
y(1) = y1-y0
x(2) = x2-x0
y(2) = y2-y0
z = z1-z0

call cos_dir(Fi,Fd,Fx,Fy,Fz)
call cos_dir(Fi,Fd,Mx,My,Mz) !considerar que la orientación del material magnético es la misma que la del campo aplicado.
!call esfe2cart(Fi,Fd,Fx,Fy,Fz)
!call esfe2cart(Fi,Fd,Mx,My,Mz)

arg1 = (Mx*Fy + My*Fx) * 0.5
arg2 = (Mx*Fz + Mz*Fx) * 0.5
arg3 = (My*Fz + Mz*Fy) * 0.5
arg4 = Mx*Fx
arg5 = My*Fy
arg6 = Mz*Fz

suma = 0

do i = 1,2
    do j = 1,2
    mu_ij = signo(i)*signo(j) !para evitar el if de a continuación
    !signo = 1
    !if (i /= j) then
    !  sign = -1
    !end if

    zz = z**2
    xx = x(i)**2
    rr = xx + y(j)**2 + z**2
    r = sqrt(rr)
    xy = x(i)*y(j)

    arg1b = (r-x(i)) / (r+x(i))
    arg2b = (r-y(j)) / (r+y(j))
    arg3b = r + z
    arg4b = xy / (xx + (r*z) + zz)
    arg5b = xy / (rr + (r*z) + xx)
    arg6b = xy / (r*z)

    suma = suma + mu_ij*( arg1*log(arg1b) + arg2*log(arg2b) - arg3*log(arg3b) &
                          - arg4*atan(arg4b) - arg5*atan(arg5b) + arg6*atan(arg6b) )

  end do
end do

operator = suma*Cm*T2nT!*delta_J

return
end subroutine op_mg2
!*******************************************************************************


subroutine cos_dir(incl,decl,a,b,g)
!Subrutina que calcula los cosenos directores o las contribuciones
!de un vector con inclinacion y declinacion a sus componentes unitarias de x, y, z.

!Parametros de entrada:
  !incl: inclinacion en grados, positiva bajo la horizoltal.
  !decl: declinacion en grados, positiva hacia el Este desde Norte geografico.
  !azim: azimut en grados, positivo hacia el Este desde el Norte.	

!Parametros de salida:
  !x,y,z: coseno director para cada eje cartesiano.

use iso_fortran_env , dp => real64

implicit none

real, intent(in) :: incl, decl
real(dp), intent(out) :: a, b, g !alpha, beta, gamma
real, parameter :: azim=270 !permite rotar el sistema en -x grados. !270grados para colocar declinacion=0 en direccion al norte.
real(dp) :: deg2rad, incl_rad, decl_rad, azim_rad

deg2rad = 0.017453293  !1deg*pi/180
incl_rad=incl*deg2rad
decl_rad=decl*deg2rad
azim_rad=azim*deg2rad

!a=cos(incl_rad)*cos(decl_rad-azim_rad)
!b=cos(incl_rad)*sin(decl_rad-azim_rad)
a=cos(incl_rad)*cos(azim_rad-decl_rad)
b=cos(incl_rad)*sin(azim_rad-decl_rad)
g=sin(incl_rad)

return
end subroutine cos_dir
!*******************************************************************************

subroutine esfe2cart(incl,decl,x_hat,y_hat,z_hat)!,inten
!Subrutina que calcula la transformación los datos de intensidad, inclinación y declinación
!del campo (coordenadas esféricas) a sus componentes cartesianas x, y, z.

!Parametros de entrada:
  !inten: intensidad del campo en nT.
  !incl: inclinacion en grados, positiva bajo la horizoltal.
  !decl: declinacion en grados, positiva hacia el Este.

!Parametros de salida:
  !x,y,z: las 3 direcciones en sistema cartesiano.

use iso_fortran_env , dp => real64

implicit none

real, intent(in) :: incl, decl!, inten
real(dp), intent(out) :: x_hat, y_hat, z_hat!, x, y ,z
!real(dp), intent(out) :: x, y ,z

!Componentes x, y, z del vector expresado en coord esfericas.
!x = inten * sin(incl)*cos(decl)
!y = inten * sin(incl)*sin(decl)
!z = inten * cos(incl)

!Componetes unitarias (componente direccion / magnitud de vector).
x_hat = sin(incl)*cos(decl)
y_hat = sin(incl)*sin(decl)
z_hat = cos(incl)

return
end subroutine esfe2cart
!*******************************************************************************

end module operators
!*******************************************************************************
!Compilar en consola: gfortran -c operators.f08
!Al no tener funcion main la bandera -c indica que se salte esa advertencia.