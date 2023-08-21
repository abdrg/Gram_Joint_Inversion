module gram_joint_inversion_host
! Contains subroutines to perform joint inversion of gravimetric and magnetic data in parallel using Fortran.
! Autor: Abraham Del Razo, IPICyT. Last Update: August 2023.
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************

contains

!*******************************************************************************
subroutine jointGramCPU(m,n,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                     d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,dip,strike,num_iters,err0)
! Subrutina que calcula inversion conjunta de dos modelos gravimétrico y magnético
! considerando similitud estructural evaluada por gradientes cruzados implementada mediante el enfoque Gramiano.

use iso_fortran_env , dp => real64
implicit none

integer, value :: m, n 
integer, intent(in) :: num_iters
real, intent(in) :: err0, dip, strike
real, dimension(10) :: param_reg
!real, intent(in), dimension(nx*ny) :: zobs_gv, zobs_mg
!real, intent(in), dimension(nx*ny) :: d_gv, d_mg, stdevd_gv, stdevd_mg
!real, intent(in), dimension(nx*ny*nz) :: m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
!real(dp), intent(in), dimension(nx*ny,nx*ny*nz) :: A_gv, A_mg
real, allocatable, dimension (:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real, allocatable, dimension(:) :: d_gv, d_mg, stdevd_gv, stdevd_mg, m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
real(dp), allocatable, dimension(:,:) :: A_gv, A_mg

real(dp), allocatable, dimension(:,:) :: dx, dy, dz
real(dp), allocatable, dimension(:,:) :: opL, opM, opL2, opM2, opN2

real, allocatable, dimension(:,:) :: aL2, gM2, eN2
real, allocatable, dimension(:,:) :: bWm2, Wd2, ATWd2, ATWd2A
real(dp), dimension(n,n) :: Ia_gv, Ia_mg
real(dp), dimension(2*n) :: Ib

real, allocatable, dimension(:) :: dxm1, dym1, dzm1, dxm2, dym2, dzm2, diva, divb, divc, dive, &
                                   ter1, ter2, ter3, dum6, dum7, dum8, ter1J, ter2J, &
                                   dxI1, dyI1, dzI1, dxI2, dyI2, dzI2

real(dp), allocatable, dimension(:) :: mIg, mJg, I_new, I, Iconj_new, Iconj, kn_2a
real(dp) :: kn_1, kn_2b, kn_2c, kn
real :: Am(2*m), dr(2*m), m_n(2*n), m_new(2*n)

real :: conv_gv, conv_mg, rms_gv, rms_mg
real :: Wd_gv(m), Wd_mg(m)
integer :: ii, jj, kk
character(len=42) :: filename00, filename1, filename2, filename3, filename4


!*****************************************
allocate ( dx(n,n),dy(n,n),dz(n,n),opL(n,n),opM(n,n),opL2(n,n),opM2(n,n),opN2(n,n) )

!matrices para derivadas
OPEN(UNIT=3, FILE="input_dx.dat", ACTION="read", FORM="unformatted")
READ(3) dx
CLOSE(UNIT=3)

OPEN(UNIT=4, FILE="input_dy.dat", ACTION="read", FORM="unformatted")
READ(4) dy
CLOSE(UNIT=4)

OPEN(UNIT=5, FILE="input_dz.dat", ACTION="read", FORM="unformatted")
READ(5) dz
CLOSE(UNIT=5)

!Operador Suavidad
OPEN(UNIT=6, FILE="input_Lapl.dat", ACTION="read", FORM="unformatted")
READ(6) opL
CLOSE(UNIT=6)
opL2 = MATMUL(transpose(opL),opL)

!Operador variacion horizontal total
opM = dx*cos(dip)*cos(strike) + dy*cos(dip)*sin(strike) + dz*sin(dip)
opM2 = MATMUL(transpose(opM),opM)

!Operador variacion vertical
!opN = dz
opN2 = MATMUL(transpose(dz),dz)

deallocate(opL,opM)


!***********************
!para gv
allocate ( aL2(n,n),gM2(n,n),eN2(n,n),bWm2(n,n),Wd2(m,m),ATWd2(n,m),ATWd2A(n,n) )

!Calculando matrices de pesos W
do ii = 1, n
  do jj = 1, n
    if (ii == jj) then
      bWm2(ii,jj) = param_reg(3) * (1/(stdevm_gv(ii)**2))
    else
      bWm2(ii,jj) = 0
    end if
  end do
end do

do ii = 1, m
  do jj = 1, m
    if (ii == jj) then
      Wd2(ii,jj) = 1/(stdevd_gv(ii)**2)
    else
      Wd2(ii,jj) = 0
    end if
  end do
end do
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_gv = 1/stdevd_gv

!operadores al cuadrado, pesados por parametros regularizadores
aL2 = param_reg(1) * opL2
gM2 = param_reg(5) * opM2
eN2 = param_reg(7) * opN2

ATWd2 = MATMUL(transpose(A_gv),Wd2)
ATWd2A = MATMUL(ATWd2,A_gv)

!Para calcular direccion de descenso I(m)
Ia_gv = ATWd2A + aL2 + bWm2 + gM2 + eN2          
Ib(1:n) = MATMUL(ATWd2,d_gv) + MATMUL(bWm2,mapr_gv)

m_n(1:n) = m_gv
deallocate(stdevm_gv,stdevd_gv,mapr_gv,m_gv)

!**********************
!para mg
aL2=0; gM2=0; eN2=0; bWm2=0; Wd2=0; ATWd2=0; ATWd2A=0

!Calculando matrices de pesos W
do ii = 1, n
  do jj = 1, n
    if (ii == jj) then
      bWm2(ii,jj) = param_reg(4) * (1/(stdevm_mg(ii)**2))
    else
      bWm2(ii,jj) = 0
    end if
  end do
end do

do ii = 1, m
  do jj = 1, m
    if (ii == jj) then
      Wd2(ii,jj) = 1/(stdevd_mg(ii)**2)
    else
      Wd2(ii,jj) = 0
    end if
  end do
end do
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_mg = 1/stdevd_mg

!operadores al cuadrado, pesados por parametros regularizadores
aL2 = param_reg(2) * opL2
gM2 = param_reg(6) * opM2
eN2 = param_reg(8) * opN2

ATWd2 = MATMUL(transpose(A_mg),Wd2)
ATWd2A = MATMUL(ATWd2,A_mg)

!Para calcular direccion de descenso I(m)
Ia_mg = ATWd2A + aL2 + bWm2 + gM2 + eN2            
Ib(n+1:2*n) = MATMUL(ATWd2,d_mg) + MATMUL(bWm2,mapr_mg)

m_n(n+1:2*n) = m_mg
deallocate(stdevm_mg,stdevd_mg,mapr_mg,m_mg)

deallocate (opL2,opM2,opN2,aL2,gM2,eN2,bWm2,Wd2,ATWd2,ATWd2A)


!**********************
allocate ( dxm1(n),dym1(n),dzm1(n),dxm2(n),dym2(n),dzm2(n),diva(n),divb(n),divc(n),dive(n), &
           ter1(n),ter2(n),ter3(n),dum6(n),dum7(n),dum8(n),ter1J(n),ter2J(n), &
           dxI1(n),dyI1(n),dzI1(n),dxI2(n),dyI2(n),dzI2(n) )
allocate ( mIg(2*n), mJg(2*n),I_new(2*n), I(2*n),Iconj(2*n),Iconj_new(2*n),kn_2a(2*n) )


WRITE (filename00, '("output_misfit_jointGram.txt")' )
OPEN(unit=16,file=filename00,status='unknown')
WRITE(16,*) 'iteracion,convergencia%_gv,convergencia%_mg,RMS_gv,RMS_mg'

!Método iterativo Gradiente Conjugado Regularizado RCG.
conv_gv = 1e6 !error inicial (valor arbitrario grande)
conv_mg = 1e6
do ii = 1,num_iters,1
      if ((conv_gv > err0).and.(conv_mg > err0)) then

            !Preparativos para calcular I(m). Enfoque, menos operaciones
            !matrices de derivadas
            dxm1 = MATMUL(dx,m_n(1:n))
            dym1 = MATMUL(dy,m_n(1:n))
            dzm1 = MATMUL(dz,m_n(1:n))
            dxm2 = MATMUL(dx,m_n(n+1:2*n))
            dym2 = MATMUL(dy,m_n(n+1:2*n))
            dzm2 = MATMUL(dz,m_n(n+1:2*n))
      
            !Ig
            ter1 = (dxm1*dxm2)+(dym1*dym2)+(dzm1*dzm2) !producto interno nabla_m1,nabla_m2
            dum6=0; dum7=0; dum8=0
            dum6 = ter1*dxm2 !array operation
            dum7 = ter1*dym2
            dum8 = ter1*dzm2
            diva = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            dum6=0; dum7=0; dum8=0
            dum6 = ter1*dxm1
            dum7 = ter1*dym1
            dum8 = ter1*dzm1
            divc = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            ter2 = (dxm2**2) + (dym2**2) + (dzm2**2) !modulo cuadrado de m2
            dum6=0; dum7=0; dum8=0
            dum6 = ter2*dxm1 
            dum7 = ter2*dym1
            dum8 = ter2*dzm1
            divb = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            ter3 = (dxm1**2) + (dym1**2) + (dzm1**2) !modulo cuadrado de m1
            dum6=0; dum7=0; dum8=0
            dum6 = ter3*dxm2
            dum7 = ter3*dym2
            dum8 = ter3*dzm2
            dive = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            mIg(1:n) = param_reg(9) * (diva - divb)
            mIg(n+1:2*n) = param_reg(10) * (divc - dive)

            !Calculando direccion de ascenso I(m) por descenso más pronunciado.
            I_new(1:n) = MATMUL(Ia_gv,m_n(1:n)) - Ib(1:n) + mIg(1:n)
            I_new(n+1:2*n) = MATMUL(Ia_mg,m_n(n+1:2*n)) - Ib(n+1:2*n) + mIg(n+1:2*n)


            !Reset variables temporales.
            diva=0; divb=0; divc=0; dive=0

            !Jg
            dxI1 = MATMUL(dx,I_new(1:n))
            dyI1 = MATMUL(dy,I_new(1:n))
            dzI1 = MATMUL(dz,I_new(1:n))
            dxI2 = MATMUL(dx,I_new(n+1:2*n))
            dyI2 = MATMUL(dy,I_new(n+1:2*n))
            dzI2 = MATMUL(dz,I_new(n+1:2*n))

            ter1J = (dxm2*dxI1)+(dym2*dyI1)+(dzm2*dzI1) !producto interno nabla_m2,nabla_I
            dum6=0; dum7=0; dum8=0
            dum6 = ter1J*dxm2
            dum7 = ter1J*dym2
            dum8 = ter1J*dzm2
            diva = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            !ter2 = (dxm2**2) + (dym2**2) + (dzm2**2) !modulo cuadrado de m2
            dum6=0; dum7=0; dum8=0            
            dum6 = ter2*dxI1
            dum7 = ter2*dyI1
            dum8 = ter2*dzI1
            divb = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            ter2J = (dxm1*dxI2)+(dym1*dyI2)+(dzm1*dzI2) !producto interno nabla_m1,nabla_I
            dum6=0; dum7=0; dum8=0
            dum6 = ter2J*dxm1
            dum7 = ter2J*dym1
            dum8 = ter2J*dzm1
            divc = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            !ter3 = (dxm1**2) + (dym1**2) + (dzm1**2) !modulo cuadrado de m1
            dum6=0; dum7=0; dum8=0
            dum6 = ter3*dxI2
            dum7 = ter3*dyI2
            dum8 = ter3*dzI2
            dive = MATMUL(dx,dum6) + MATMUL(dy,dum7) + MATMUL(dz,dum8)

            mJg(1:n) = param_reg(9) * (diva - divb)
            mJg(n+1:2*n) = param_reg(10) * (divc - dive)


            !Direccion de ascenso mediante gradientes conjugados.
            if (ii == 1) then
                  Iconj_new = I_new
            else
                  Iconj_new = I_new + ( (DOT_PRODUCT(I_new,I_new)/DOT_PRODUCT(I,I)) * Iconj )
            end if


            ! Tamaño de paso kn en la dirección I_n.
            kn_1 = DOT_PRODUCT(Iconj_new,I_new)

            kn_2a(1:n) = MATMUL(Ia_gv,Iconj_new(1:n))
            kn_2a(n+1:2*n) = MATMUL(Ia_mg,Iconj_new(n+1:2*n))

            kn_2b = DOT_PRODUCT(Iconj_new,kn_2a)
            kn_2c = DOT_PRODUCT(Iconj_new,mJg)

            kn = kn_1 / (kn_2b + kn_2c)
            
            !Calculo m_new
            m_new = m_n - (kn*Iconj_new) !array operation


            ! Condicion de paro
            !Convergencia
            conv_gv = 100 * SQRT( SUM( ((m_new(1:n)-m_n(1:n))**2)/((m_n(1:n)**2)+1e-10) )/n )
            conv_mg = 100 * SQRT( SUM( ((m_new(n+1:2*n)-m_n(n+1:2*n))**2)/((m_n(n+1:2*n)**2)+1e-10) )/n )
            
            !Actualizando variables para la siguiente iteración
            m_n = m_new
            I = I_new
            Iconj = Iconj_new

            !Obteniendo datos calculados
            Am(1:m) = MATMUL(A_gv,m_new(1:n))
            Am(m+1:2*m) = MATMUL(A_mg,m_new(n+1:2*n))

            !RMS
            dr(1:m) = Am(1:m) - d_gv !array operation
            dr(m+1:2*m) = Am(m+1:2*m) - d_mg
            dr(1:m) = Wd_gv * dr(1:m) !array operation
            dr(m+1:2*m) = Wd_mg * dr(m+1:2*m)

            rms_gv = DOT_PRODUCT(dr(1:m),dr(1:m))
            rms_gv = SQRT(rms_gv/m)

            rms_mg = DOT_PRODUCT(dr(m+1:2*m),dr(m+1:2*m))
            rms_mg = SQRT(rms_mg/m)
            
            WRITE(16,'(I3,",",F16.4,",",F16.4,",",F12.4,",",F12.4)') ii,conv_gv,conv_mg,rms_gv,rms_mg


            !Escribir resultados
            ! Escribir archivos de modelos invertidos
            WRITE (filename1, '("output_jointGram_gv_minv_iter",I2,".csv")' ) ii
            OPEN(unit=9,file=filename1,status='unknown') !crear y escribir en archivo
            WRITE(9,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'

            WRITE (filename2, '("output_jointGram_mg_minv_iter",I2,".csv")' ) ii
            OPEN(unit=10,file=filename2,status='unknown') !crear y escribir en archivo
            WRITE(10,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'

            do jj = 1,n,1
                  WRITE(9,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(jj)
                  WRITE(10,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(n+jj)
            end do

            CLOSE(unit=9)
            CLOSE(unit=10)

            ! Escribir archivos de datos calculados a partir de modelos invertidos
            WRITE (filename3, '("output_jointGram_gv_anomaly_iter",I2,".csv")' ) ii
            OPEN(unit=11,file=filename3,status='unknown') !crear y escribir en archivo
            WRITE(11,*) 'x[m],y[m],z[m],Anomaly[mGal]'

            WRITE (filename4, '("output_jointGram_mg_anomaly_iter",I2,".csv")' ) ii
            OPEN(unit=12,file=filename4,status='unknown') !crear y escribir en archivo
            WRITE(12,*) 'x[m],y[m],z[m],Anomaly[nT]'

            do jj = 1,m,1
                  WRITE(11,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_gv(jj),Am(jj)
                  WRITE(12,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_mg(jj),Am(m+jj)
            end do

            CLOSE(unit=11)
            CLOSE(unit=12)

      end if
end do

CLOSE(unit=16)

deallocate ( dxm1,dym1,dzm1,dxm2,dym2,dzm2,diva,divb,divc,dive, &
             ter1,ter2,ter3,dum6,dum7,dum8,ter1J,ter2J, &
             dxI1,dyI1,dzI1,dxI2,dyI2,dzI2 )
deallocate ( dx,dy,dz,mIg,mJg,I_new,I,Iconj,Iconj_new,kn_2a)
deallocate ( A_gv,A_mg,d_gv,d_mg, xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell )

return
end subroutine jointGramCPU
!*******************************************************************************

end module gram_joint_inversion_host