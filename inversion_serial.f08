module inversion_serial
! Contiene subrutinas para hacer modelado inverso de datos gravimétricos inverse_gv
! y magnéticos inverse_mg.

! Autor: Abraham Del Razo, IPICyT. May-2020
! Mail: abraham.delrazo@ipicyt.edu.mx
!*******************************************************************************
contains

subroutine jointGram3(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                     d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
! Subrutina que calcula inversion conjunta de dos modelos gravimétrico y magnético
! considerando similitud estructural evaluada por gradientes cruzados implementada mediante el enfoque Gramiano.

use iso_fortran_env , dp => real64
use auxiliar, only: cov, partialder, laplaciano, coord_voxel
!use model_ini, only: quasi_inf

implicit none

integer, value :: m, n 
integer, intent(in) :: num_iters, nx, ny, nz, ciclo
real, intent(in) :: step, err0, strike
real, dimension(10) :: param_reg
!real, intent(in), dimension(nx*ny) :: zobs_gv, zobs_mg
!real, intent(in), dimension(nx*ny) :: d_gv, d_mg, stdevd_gv, stdevd_mg
!real, intent(in), dimension(nx*ny*nz) :: m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
!real(dp), intent(in), dimension(nx*ny,nx*ny*nz) :: A_gv, A_mg
real, allocatable, dimension (:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real, allocatable, dimension(:) :: d_gv, d_mg, stdevd_gv, stdevd_mg, m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
real(dp), allocatable, dimension(:,:) :: A_gv, A_mg

integer :: ii, jj, kk
real :: conv_gv, conv_mg, rms_gv, rms_mg

real :: m_n(2*n), Am(2*m), dr(2*m)
real :: Wd_gv(m), Wd_mg(m)

!real, dimension(nx*ny,nx*ny) :: Wd2_gv, Wd2_mg
!real, dimension(nx*ny*nz,nx*ny*nz) :: Wm2_gv, Wm2_mg
real, allocatable, dimension(:,:) :: Wd2_gv, Wd2_mg
real, allocatable, dimension(:,:) :: Wm2_gv, Wm2_mg
!real, dimension(nx*ny*nz,nx*ny*nz) :: bWm2
real, allocatable :: bWm2(:,:)

real(dp), dimension(n,n) :: dx, dy, dz!, opM, opN, opL, opM2, opN2, opL2
real(dp), allocatable, dimension(:,:) :: opM, opN, opL, opM2, opN2, opL2
!real, dimension(nx*ny*nz,nx*ny*nz) :: aL2, gM2, eN2
real, allocatable, dimension(:,:) :: aL2, gM2, eN2
!real(dp), dimension(nx*ny*nz,nx*ny) :: ATWd2_gv, ATWd2_mg
real(dp), allocatable, dimension(:,:) :: ATWd2_gv, ATWd2_mg
real(dp), allocatable, dimension(:,:) :: ATWd2A_gv, ATWd2A_mg
real(dp), dimension(n,n) :: Ia_gv, Ia_mg!, ATWd2A_gv, ATWd2A_mg

real, dimension(n) :: dxm1, dym1, dzm1, dxm2, dym2, dzm2, ter1, ter2, ter3, &
                             diva, divb, divc, dive, ax, ay, az, bx, by, bz, cx, cy, cz, ex, ey, ez
real(dp), dimension(n) :: dxI1, dyI1, dzI1, dxI2, dyI2, dzI2, ter1J, ter2J!, &
                                 !axJ, ayJ, azJ, bxJ, byJ, bzJ, cxJ, cyJ, czJ, exJ, eyJ, ezJ, divaJ, divbJ, divcJ, diveJ

real(dp), dimension(2*n) :: I_new, I, Iconj_new, Iconj, kn_2a, mIg, mJg, Ib
real(dp) :: kn_1, kn_2b, kn_2c, kn
real, dimension(2*n) :: m_new

character(len=42) :: filename00, filename1, filename2, filename3, filename4, filename5

real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)
PRINT*, 'Inicia inversion conjunta en serie'

!*****************************************
!reservar memoria para arreglos temporales
!allocate ( d_gv(m),d_mg(m),stdevd_gv(m),stdevd_mg(m),m_gv(n),m_mg(n),mapr_gv(n),mapr_mg(n),stdevm_gv(n),stdevm_mg(n), &
!          A_gv(m,n),A_mg(m,n),xcell(n), ycell(n), zcell(n), xobs(m), yobs(m), zobs_gv(m), zobs_mg(m) )

allocate ( Wd2_gv(m,m),Wd2_mg(m,m),Wm2_gv(n,n),Wm2_mg(n,n),bWm2(n,n), &
opM(n,n),opN(n,n),opL(n,n),opM2(n,n),opN2(n,n),opL2(n,n),aL2(n,n),gM2(n,n),eN2(n,n), &
ATWd2_gv(n,m),ATWd2A_gv(n,n),ATWd2_mg(n,m),ATWd2A_mg(n,n) )

!matrices para derivadas
call partialder(nx,ny,nz,step,'x',dx)
call partialder(nx,ny,nz,step,'y',dy)
call partialder(nx,ny,nz,step,'z',dz)

!Operador Suavidad
call laplaciano(nx,ny,nz,step,opL)
opL2 = MATMUL(transpose(opL),opL)
!Operador variacion horizontal total
opM = dx*cos(strike) + dy*sin(strike)
opM2 = MATMUL(transpose(opM),opM)
!Operador variacion vertical
!opN = dz
opN2 = MATMUL(transpose(dz),dz)

!para gv
m_n(1:n) = m_gv
!Obteniendo datos calculados
Am(1:m) = MATMUL(A_gv,m_gv)

!***********************
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_gv = 1/stdevd_gv
!Calculando matrices de pesos
!subroutine cov(num,type,std_dev,covar)
call cov(m,1,stdevd_gv,Wd2_gv)
call cov(n,1,stdevm_gv,Wm2_gv)
!operadores al cuadrado, pesados por parametros regularizadores
bWm2 = param_reg(3) * Wm2_gv
aL2 = param_reg(1) * opL2
gM2 = param_reg(5) * opM2
eN2 = param_reg(7) * opN2

!**********************
!Para calcular direccion de descenso I(m)
ATWd2_gv = MATMUL(transpose(A_gv),Wd2_gv)
ATWd2A_gv = MATMUL(ATWd2_gv,A_gv)
Ia_gv = ATWd2A_gv + aL2 + bWm2 + gM2 + eN2            
Ib(1:n) = MATMUL(ATWd2_gv,d_gv) + MATMUL(bWm2,mapr_gv)

!para mg
m_n(n+1:2*n) = m_mg
!Obteniendo datos calculados
Am(m+1:2*m) = MATMUL(A_mg,m_mg)

!***********************
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_mg = 1/stdevd_mg
!Calculando matrices de pesos
!subroutine cov(num,type,std_dev,covar)
call cov(m,1,stdevd_mg,Wd2_mg)
call cov(n,1,stdevm_mg,Wm2_mg)
!operadores al cuadrado, pesados por parametros regularizadores
bWm2 = param_reg(4) * Wm2_mg
aL2 = param_reg(2) * opL2
gM2 = param_reg(6) * opM2
eN2 = param_reg(8) * opN2

!**********************
!Para calcular direccion de descenso I(m)
ATWd2_mg = MATMUL(transpose(A_mg),Wd2_mg)
ATWd2A_mg = MATMUL(ATWd2_mg,A_mg)
Ia_mg = ATWd2A_mg + aL2 + bWm2 + gM2 + eN2            
Ib(n+1:2*n) = MATMUL(ATWd2_mg,d_mg) + MATMUL(bWm2,mapr_mg)


!**********************
!liberar memoria para arreglos temporales
deallocate (m_gv,m_mg,mapr_gv,mapr_mg,stdevd_gv,stdevd_mg,stdevm_gv,stdevm_mg)!d_gv,d_mg,A_gv,A_mg)
deallocate (Wd2_gv,Wd2_mg,Wm2_gv,Wm2_mg,bWm2,opL,opL2,opM,opM2,opN,opN2,aL2,gM2,eN2,ATWd2_gv,ATWd2A_gv,ATWd2_mg,ATWd2A_mg)

!**********************
!Método iterativo Gradiente Conjugado Regularizado RCG.

WRITE (filename00, '("output_misfit_jointGram_c",I1,".txt")' ) ciclo
OPEN(unit=16,file=filename00,status='unknown')
WRITE(16,*) 'iteracion,convergencia%_gv,convergencia%_mg,RMS_gv,RMS_mg'

conv_gv = 1e6 !error inicial (valor arbitrario grande)
conv_mg = 1e6
do ii = 1,num_iters,1
      if ((conv_gv > err0).and.(conv_mg > err0)) then
          
            !Preparativos para calcular I(m).
            !matrices de derivadas
            dxm1 = MATMUL(dx,m_n(1:n))
            dym1 = MATMUL(dy,m_n(1:n))
            dzm1 = MATMUL(dz,m_n(1:n))
            dxm2 = MATMUL(dx,m_n(n+1:2*n))
            dym2 = MATMUL(dy,m_n(n+1:2*n))
            dzm2 = MATMUL(dz,m_n(n+1:2*n))
      
            !Ig
            ter1 = (dxm1*dxm2)+(dym1*dym2)+(dzm1*dzm2) !producto interno nabla_m1,nabla_m2
            ax = ter1*dxm2 !array operation
            ay = ter1*dym2
            az = ter1*dzm2
            
            ter2 = (dxm2**2) + (dym2**2) + (dzm2**2) !modulo cuadrado de m2
            bx = ter2*dxm1 
            by = ter2*dym1
            bz = ter2*dzm1
            
            cx = ter1*dxm1
            cy = ter1*dym1
            cz = ter1*dzm1

            ter3 = (dxm1**2) + (dym1**2) + (dzm1**2) !modulo cuadrado de m1
            ex = ter3*dxm2
            ey = ter3*dym2
            ez = ter3*dzm2
            
            diva = MATMUL(dx,ax) + MATMUL(dy,ay) + MATMUL(dz,az)
            divb = MATMUL(dx,bx) + MATMUL(dy,by) + MATMUL(dz,bz)
            divc = MATMUL(dx,cx) + MATMUL(dy,cy) + MATMUL(dz,cz)
            dive = MATMUL(dx,ex) + MATMUL(dy,ey) + MATMUL(dz,ez)

            mIg(1:n) = param_reg(9) * (diva - divb)
            mIg(n+1:2*n) = param_reg(10) * (divc - dive)


            !Calculando direccion de ascenso I(m) por descenso más pronunciado.
            I_new(1:n) = MATMUL(Ia_gv,m_n(1:n)) - Ib(1:n) + mIg(1:n)
            I_new(n+1:2*n) = MATMUL(Ia_mg,m_n(n+1:2*n)) - Ib(n+1:2*n) + mIg(n+1:2*n)
            
            !Reset variables temporales.
            ax=0; ay=0; az=0; bx=0; by=0; bz=0; cx=0; cy=0; cz=0; ex=0; ey=0; ez=0
            diva=0; divb=0; divc=0; dive=0

            !Jg
            dxI1 = MATMUL(dx,I_new(1:n))
            dyI1 = MATMUL(dy,I_new(1:n))
            dzI1 = MATMUL(dz,I_new(1:n))
            dxI2 = MATMUL(dx,I_new(n+1:2*n))
            dyI2 = MATMUL(dy,I_new(n+1:2*n))
            dzI2 = MATMUL(dz,I_new(n+1:2*n))

            ter1J = (dxm2*dxI1)+(dym2*dyI1)+(dzm2*dzI1) !producto interno nabla_m2,nabla_I
            ax = ter1J*dxm2
            ay = ter1J*dym2
            az = ter1J*dzm2
            
            bx = ter2*dxI1
            by = ter2*dyI1
            bz = ter2*dzI1

            ter2J = (dxm1*dxI2)+(dym1*dyI2)+(dzm1*dzI2) !producto interno nabla_m1,nabla_I
            cx = ter2J*dxm1
            cy = ter2J*dym1
            cz = ter2J*dzm1

            ex = ter3*dxI2
            ey = ter3*dyI2
            ez = ter3*dzI2
            
            diva = MATMUL(dx,ax) + MATMUL(dy,ay) + MATMUL(dz,az)
            divb = MATMUL(dx,bx) + MATMUL(dy,by) + MATMUL(dz,bz)
            divc = MATMUL(dx,cx) + MATMUL(dy,cy) + MATMUL(dz,cz)
            dive = MATMUL(dx,ex) + MATMUL(dy,ey) + MATMUL(dz,ez)

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

            kn = kn_1 / (kn_2b + kn_2c +1e-10) !para evitar divisiones por cero

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
            
            WRITE(16,'(I3,",",F14.4,",",F14.4,",",F10.4,",",F10.4)') ii,conv_gv,conv_mg,rms_gv,rms_mg


            !Escribir resultados
            ! Escribir archivos de modelos invertidos
            WRITE (filename1, '("output_jointGram_gv_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=9,file=filename1,status='unknown') !crear y escribir en archivo
            WRITE(9,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'

            WRITE (filename2, '("output_jointGram_mg_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=10,file=filename2,status='unknown') !crear y escribir en archivo
            WRITE(10,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'

            do jj = 1,n,1
                  WRITE(9,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(jj)
                  WRITE(10,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(n+jj)
            end do

            CLOSE(unit=9)
            CLOSE(unit=10)

            ! Escribir archivos de datos calculados a partir de modelos invertidos
            WRITE (filename3, '("output_jointGram_gv_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=11,file=filename3,status='unknown') !crear y escribir en archivo
            WRITE(11,*) 'x[m],y[m],z[m],g[mGal]'

            WRITE (filename4, '("output_jointGram_mg_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=12,file=filename4,status='unknown') !crear y escribir en archivo
            WRITE(12,*) 'x[m],y[m],z[m],m[nT]'

            do jj = 1,m,1
                  WRITE(11,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_gv(jj),Am(jj)
                  WRITE(12,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_mg(jj),Am(m+jj)
            end do

            CLOSE(unit=11)
            CLOSE(unit=12)


      end if
end do

CLOSE(unit=16)

! Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
WRITE(0,*) 'Execution time: serial jointGram =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

return
end subroutine jointGram3
!*******************************************************************************


subroutine jointGram4(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                     d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
! Subrutina que calcula inversion conjunta de dos modelos gravimétrico y magnético
! considerando similitud estructural evaluada por gradientes cruzados implementada mediante el enfoque Gramiano.
! Empleando Lapack y Blas como paso intermedio para despues utilizar cuBlas. 

use iso_fortran_env , dp => real64
use auxiliar, only: cov, partialder, laplaciano, coord_voxel
!use model_ini, only: quasi_inf

implicit none

integer, value :: m, n
integer, intent(in) :: num_iters, nx, ny, nz, ciclo
real, intent(in) :: step, err0, strike
real, dimension(10) :: param_reg
!real, intent(in), dimension(nx*ny) :: zobs_gv, zobs_mg
!real, intent(in), dimension(nx*ny) :: d_gv, d_mg, stdevd_gv, stdevd_mg
!real, intent(in), dimension(nx*ny*nz) :: m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
!real(dp), intent(in), dimension(nx*ny,nx*ny*nz) :: A_gv, A_mg
real, allocatable, dimension (:) :: xobs, yobs, zobs_gv, zobs_mg, xcell, ycell, zcell
real, allocatable, dimension(:) :: d_gv, d_mg, stdevd_gv, stdevd_mg, m_gv, m_mg, mapr_gv, mapr_mg, stdevm_gv, stdevm_mg
real(dp), allocatable, dimension(:,:) :: A_gv, A_mg

integer :: ii, jj, kk
real :: conv_gv, conv_mg, rms_gv, rms_mg

real(dp) :: m_n(2*n), Am(2*m), dr(2*m)
real :: Wd_gv(m), Wd_mg(m)

!real, dimension(nx*ny,nx*ny) :: Wd2_gv, Wd2_mg
!real, dimension(nx*ny*nz,nx*ny*nz) :: Wm2_gv, Wm2_mg
real, allocatable, dimension(:,:) :: Wd2_gv, Wd2_mg
real, allocatable, dimension(:,:) :: Wm2_gv, Wm2_mg
!real, dimension(nx*ny*nz,nx*ny*nz) :: bWm2
real, allocatable :: bWm2(:,:)

real(dp), dimension(n,n) :: dx, dy, dz!, opM, opN, opL, opM2, opN2, opL2
real(dp), allocatable, dimension(:,:) :: opM, opN, opL, opM2, opN2, opL2
!real, dimension(nx*ny*nz,nx*ny*nz) :: aL2, gM2, eN2
real, allocatable, dimension(:,:) :: aL2, gM2, eN2
!real(dp), dimension(nx*ny*nz,nx*ny) :: ATWd2_gv, ATWd2_mg
real(dp), allocatable, dimension(:,:) :: ATWd2_gv, ATWd2_mg
real(dp), allocatable, dimension(:,:) :: ATWd2A_gv, ATWd2A_mg
real(dp), dimension(n,n) :: Ia_gv, Ia_mg!, ATWd2A_gv, ATWd2A_mg

real(dp), dimension(n) :: dxm1, dym1, dzm1, dxm2, dym2, dzm2, ter1, ter2, ter3, &
                                 diva, divb, divc, dive, ax, ay, az, bx, by, bz, cx, cy, cz, ex, ey, ez
real(dp), dimension(n) :: dxI1, dyI1, dzI1, dxI2, dyI2, dzI2, ter1J, ter2J

real(dp), dimension(2*n) :: I_new, I, Iconj_new, Iconj, kn_2a, mIg, mJg, Ib
real(dp) :: kn_1, kn_2b, kn_2c, kn
real(dp), dimension(2*n) :: m_new

character(len=42) :: filename00, filename1, filename2, filename3, filename4, filename5

!variables extra para BLAS y LAPACK
real(dp), external :: DDOT
real(dp) :: ALPHA = 1., BETA
!real(dp) :: I2i, Inew2

real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)
PRINT*, 'Inicia inversion conjunta en serie'

!**********************
!reservar memoria para arreglos temporales
!allocate ( d_gv(m),d_mg(m),mapr_gv(n),mapr_mg(n),m_gv(n),m_mg(n),stdevd_gv(m),stdevd_mg(m),stdevm_gv(n),stdevm_mg(n), &
!          A_gv(m,n),A_mg(m,n),xcell(n), ycell(n), zcell(n), xobs(m), yobs(m), zobs_gv(m), zobs_mg(m) )

allocate ( Wd2_gv(m,m),Wd2_mg(m,m),Wm2_gv(n,n),Wm2_mg(n,n),bWm2(n,n), &
opM(n,n),opN(n,n),opL(n,n),opM2(n,n),opN2(n,n),opL2(n,n),aL2(n,n),gM2(n,n),eN2(n,n), &
ATWd2_gv(n,m),ATWd2A_gv(n,n),ATWd2_mg(n,m),ATWd2A_mg(n,n) )

!matrices para derivadas
call partialder(nx,ny,nz,step,'x',dx)
call partialder(nx,ny,nz,step,'y',dy)
call partialder(nx,ny,nz,step,'z',dz)

!Operador Suavidad
call laplaciano(nx,ny,nz,step,opL)
opL2 = MATMUL(transpose(opL),opL)
!call DGEMM('T','N',size(opL,1),size(opL,2),size(opL,1),ALPHA,opL,size(opL,2),opL,size(opL,1),BETA,opL2,size(opL2,1))

!Operador variacion horizontal total
opM = dx*cos(strike) + dy*sin(strike)
opM2 = MATMUL(transpose(opM),opM)
!call DGEMM('T','N',size(opM,1),size(opM,2),size(opM,1),ALPHA,opM,size(opM,2),opM,size(opM,1),BETA,opM2,size(opM2,1))

!Operador variacion vertical
opN = dz
opN2 = MATMUL(transpose(opN),opN)
!call DGEMM('T','N',size(opN,1),size(opN,2),size(opN,1),ALPHA,opN,size(opN,2),opN,size(opN,1),BETA,opN2,size(opN2,1))

!para gv
m_n(1:n) = m_gv
!Obteniendo datos calculados
!Am(1:m) = MATMUL(A_gv,m_gv)
BETA=0.
call DGEMV('N',size(A_gv,1),size(A_gv,2),ALPHA,A_gv,size(A_gv,1),m_n(1:n),1,BETA,Am(1:m),1)

!***********************
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_gv = 1/stdevd_gv
!Calculando matrices de pesos
!subroutine cov(num,type,std_dev,covar)
call cov(m,1,stdevd_gv,Wd2_gv)
call cov(n,1,stdevm_gv,Wm2_gv)
!operadores al cuadrado, pesados por parametros regularizadores
bWm2 = param_reg(3) * Wm2_gv
aL2 = param_reg(1) * opL2
gM2 = param_reg(5) * opM2
eN2 = param_reg(7) * opN2

!**********************
!Para calcular direccion de descenso I(m)
ATWd2_gv = MATMUL(transpose(A_gv),Wd2_gv)
ATWd2A_gv = MATMUL(ATWd2_gv,A_gv)
!call DGEMM('T','N',size(A_gv,1),size(Wd2_gv,2),size(Wd2_gv,1),ALPHA,A_gv,size(A_gv,2),Wd2_gv,size(Wd2_gv,1),BETA,ATWd2_gv,size(ATWd2_gv,1))
!call DGEMM('N','N',size(ATWd2_gv,1),size(A_gv,2),size(A_gv,1),ALPHA,ATWd2_gv,size(ATWd2_gv,1),A_gv,size(A_gv,1),BETA,ATWd2A_gv,size(ATWd2A_gv,1))

Ia_gv = ATWd2A_gv + aL2 + bWm2 + gM2 + eN2            
Ib(1:n) = MATMUL(ATWd2_gv,d_gv) + MATMUL(bWm2,mapr_gv)
!BETA=1.
!call DGEMV('N',size(ATWd2_gv,1),size(ATWd2_gv,2),ALPHA,ATWd2_gv,size(ATWd2_gv,1),d2_gv,1,BETA,Ib(1:n),1)
!call DGEMV('N',size(bWm2,1),size(bWm2,2),ALPHA,bWm2,size(bWm2,1),mapr2_gv,1,BETA,Ib(1:n),1)

!para mg
m_n(n+1:2*n) = m_mg
!Obteniendo datos calculados
!Am(m+1:2*m) = MATMUL(A_mg,m_mg)
BETA=0.
call DGEMV('N',size(A_mg,1),size(A_mg,2),ALPHA,A_mg,size(A_mg,1),m_n(n+1:2*n),1,BETA,Am(m+1:2*m),1)

!***********************
!Recordar W**2 = C-1 = 1/(stdev**2) = 1/varianza (type=1). Pesos W son 1/stdev (type=2).
Wd_mg = 1/stdevd_mg
!Calculando matrices de pesos
!subroutine cov(num,type,std_dev,covar)
call cov(m,1,stdevd_mg,Wd2_mg)
call cov(n,1,stdevm_mg,Wm2_mg)
!operadores al cuadrado, pesados por parametros regularizadores
bWm2 = param_reg(4) * Wm2_mg
aL2 = param_reg(2) * opL2
gM2 = param_reg(6) * opM2
eN2 = param_reg(8) * opN2

!**********************
!Para calcular direccion de descenso I(m)
ATWd2_mg = MATMUL(transpose(A_mg),Wd2_mg)
ATWd2A_mg = MATMUL(ATWd2_mg,A_mg)
!call DGEMM('T','N',size(A_mg,1),size(Wd2_mg,2),size(Wd2_mg,1),ALPHA,A_mg,size(A_mg,2),Wd2_mg,size(Wd2_mg,1),BETA,ATWd2_mg,size(ATWd2_mg,1))
!call DGEMM('N','N',size(ATWd2_mg,1),size(A_mg,2),size(A_mg,1),ALPHA,ATWd2_mg,size(ATWd2_mg,1),A_mg,size(A_mg,1),BETA,ATWd2A_mg,size(ATWd2A_mg,1))

Ia_mg = ATWd2A_mg + aL2 + bWm2 + gM2 + eN2            
Ib(n+1:2*n) = MATMUL(ATWd2_mg,d_mg) + MATMUL(bWm2,mapr_mg)
!BETA=1.
!call DGEMV('N',size(ATWd2_mg,1),size(ATWd2_mg,2),ALPHA,ATWd2_mg,size(ATWd2_mg,1),d2_mg,1,BETA,Ib(n+1:2*n),1)
!call DGEMV('N',size(bWm2,1),size(bWm2,2),ALPHA,bWm2,size(bWm2,1),mapr2_mg,1,BETA,Ib(n+1:2*n),1)


!**********************
!liberar memoria para arreglos temporales
deallocate (m_gv,m_mg,mapr_gv,mapr_mg,stdevd_gv,stdevd_mg,stdevm_gv,stdevm_mg)!d_gv,d_mg,A_gv,A_mg)
deallocate (Wd2_gv,Wd2_mg,Wm2_gv,Wm2_mg,bWm2,opL,opL2,opM,opM2,opN,opN2,aL2,gM2,eN2,ATWd2_gv,ATWd2A_gv,ATWd2_mg,ATWd2A_mg)

!**********************
!Método iterativo Gradiente Conjugado Regularizado RCG.

WRITE (filename00, '("output_misfit_jointGram_c",I1,".txt")' ) ciclo
OPEN(unit=16,file=filename00,status='unknown')
WRITE(16,*) 'iteracion,convergencia%_gv,convergencia%_mg,RMS_gv,RMS_mg'

conv_gv = 1e6 !error inicial (valor arbitrario grande)
conv_mg = 1e6
do ii = 1,num_iters,1
      if ((conv_gv > err0).and.(conv_mg > err0)) then
         
            !Preparativos para calcular I(m).
            !matrices de derivadas
            !dxm1 = MATMUL(dx,m_n(1:n))
            !dym1 = MATMUL(dy,m_n(1:n))
            !dzm1 = MATMUL(dz,m_n(1:n))
            !dxm2 = MATMUL(dx,m_n(n+1:2*n))
            !dym2 = MATMUL(dy,m_n(n+1:2*n))
            !dzm2 = MATMUL(dz,m_n(n+1:2*n))
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),m_n(1:n),1,BETA,dxm1,1)
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),m_n(1:n),1,BETA,dym1,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),m_n(1:n),1,BETA,dzm1,1)
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),m_n(n+1:2*n),1,BETA,dxm2,1)
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),m_n(n+1:2*n),1,BETA,dym2,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),m_n(n+1:2*n),1,BETA,dzm2,1)
            
            !Ig
            ter1 = (dxm1*dxm2)+(dym1*dym2)+(dzm1*dzm2) !producto interno nabla_m1,nabla_m2
            ax = ter1*dxm2 !array operation
            ay = ter1*dym2
            az = ter1*dzm2
            
            ter2 = (dxm2**2) + (dym2**2) + (dzm2**2) !modulo cuadrado de m2
            bx = ter2*dxm1 
            by = ter2*dym1
            bz = ter2*dzm1
            
            cx = ter1*dxm1
            cy = ter1*dym1
            cz = ter1*dzm1

            ter3 = (dxm1**2) + (dym1**2) + (dzm1**2) !modulo cuadrado de m1
            ex = ter3*dxm2
            ey = ter3*dym2
            ez = ter3*dzm2
            
            !diva = MATMUL(dx,ax) + MATMUL(dy,ay) + MATMUL(dz,az)
            !divb = MATMUL(dx,bx) + MATMUL(dy,by) + MATMUL(dz,bz)
            !divc = MATMUL(dx,cx) + MATMUL(dy,cy) + MATMUL(dz,cz)
            !dive = MATMUL(dx,ex) + MATMUL(dy,ey) + MATMUL(dz,ez)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),ax,1,BETA,diva,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),ay,1,BETA,diva,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),az,1,BETA,diva,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),bx,1,BETA,divb,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),by,1,BETA,divb,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),bz,1,BETA,divb,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),cx,1,BETA,divc,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),cy,1,BETA,divc,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),cz,1,BETA,divc,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),ex,1,BETA,dive,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),ey,1,BETA,dive,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),ez,1,BETA,dive,1)        

            mIg(1:n) = param_reg(9) * (diva - divb)
            mIg(n+1:2*n) = param_reg(10) * (divc - dive)


            !Calculando direccion de ascenso I(m) por descenso más pronunciado.
            !I_new(1:n) = MATMUL(Ia_gv,m_n(1:n)) - Ib(1:n) + mIg(1:n)
            !I_new(n+1:2*n) = MATMUL(Ia_mg,m_n(n+1:2*n)) - Ib(n+1:2*n) + mIg(n+1:2*n)
            I_new(1:n) = mIg(1:n) - Ib(1:n)
            I_new(n+1:2*n) = mIg(n+1:2*n) - Ib(n+1:2*n)
            BETA=1.
            call DGEMV('N',size(Ia_gv,1),size(Ia_gv,2),ALPHA,Ia_gv,size(Ia_gv,1),m_n(1:n),1,BETA,I_new(1:n),1)            
            call DGEMV('N',size(Ia_mg,1),size(Ia_mg,2),ALPHA,Ia_mg,size(Ia_mg,1),m_n(n+1:2*n),1,BETA,I_new(n+1:2*n),1)            

            !Reset variables temporales.
            ax=0; ay=0; az=0; bx=0; by=0; bz=0; cx=0; cy=0; cz=0; ex=0; ey=0; ez=0
            diva=0; divb=0; divc=0; dive=0

            !Jg
            !dxI1 = MATMUL(dx,I_new(1:n))
            !dyI1 = MATMUL(dy,I_new(1:n))
            !dzI1 = MATMUL(dz,I_new(1:n))
            !dxI2 = MATMUL(dx,I_new(n+1:2*n))
            !dyI2 = MATMUL(dy,I_new(n+1:2*n))
            !dzI2 = MATMUL(dz,I_new(n+1:2*n))
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),I_new(1:n),1,BETA,dxI1,1)
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),I_new(1:n),1,BETA,dyI1,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),I_new(1:n),1,BETA,dzI1,1)
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),I_new(n+1:2*n),1,BETA,dxI2,1)
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),I_new(n+1:2*n),1,BETA,dyI2,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),I_new(n+1:2*n),1,BETA,dzI2,1)

            ter1J = (dxm2*dxI1)+(dym2*dyI1)+(dzm2*dzI1) !producto interno nabla_m2,nabla_I
            ax = ter1J*dxm2
            ay = ter1J*dym2
            az = ter1J*dzm2
            
            bx = ter2*dxI1
            by = ter2*dyI1
            bz = ter2*dzI1

            ter2J = (dxm1*dxI2)+(dym1*dyI2)+(dzm1*dzI2) !producto interno nabla_m1,nabla_I
            cx = ter2J*dxm1
            cy = ter2J*dym1
            cz = ter2J*dzm1

            ex = ter3*dxI2
            ey = ter3*dyI2
            ez = ter3*dzI2
            
            !divaJ = MATMUL(dx,axJ) + MATMUL(dy,ayJ) + MATMUL(dz,azJ)
            !divbJ = MATMUL(dx,bxJ) + MATMUL(dy,byJ) + MATMUL(dz,bzJ)
            !divcJ = MATMUL(dx,cxJ) + MATMUL(dy,cyJ) + MATMUL(dz,czJ)
            !diveJ = MATMUL(dx,exJ) + MATMUL(dy,eyJ) + MATMUL(dz,ezJ)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),ax,1,BETA,diva,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),ay,1,BETA,diva,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),az,1,BETA,diva,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),bx,1,BETA,divb,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),by,1,BETA,divb,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),bz,1,BETA,divb,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),cx,1,BETA,divc,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),cy,1,BETA,divc,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),cz,1,BETA,divc,1)
            BETA=0.
            call DGEMV('N',size(dx,1),size(dx,2),ALPHA,dx,size(dx,1),ex,1,BETA,dive,1)
            BETA=1.
            call DGEMV('N',size(dy,1),size(dy,2),ALPHA,dy,size(dy,1),ey,1,BETA,dive,1)
            call DGEMV('N',size(dz,1),size(dz,2),ALPHA,dz,size(dz,1),ez,1,BETA,dive,1)

            mJg(1:n) = param_reg(9) * (diva - divb)
            mJg(n+1:2*n) = param_reg(10) * (divc - dive)


            !Direccion de ascenso mediante gradientes conjugados.
            if (ii == 1) then
                  Iconj_new = I_new
            else
                  Iconj_new = I_new + ( (DOT_PRODUCT(I_new,I_new)/DOT_PRODUCT(I,I)) * Iconj )
                  !Inew2 = DDOT(size(I_new),I_new,1,I_new,1)
                  !I2i = DDOT(size(I),I,1,I,1)
                  !Iconj_new = I_new + ( (Inew2/I2i) * Iconj )
            end if

           
            ! Tamaño de paso kn en la dirección I(m).
            kn_1 = DOT_PRODUCT(Iconj_new,I_new)
            !kn_1 = DDOT(size(Iconj_new),Iconj_new,1,I_new,1)

            !kn_2a(1:n) = MATMUL(Ia_gv,Iconj_new(1:n))
            !kn_2a(n+1:2*n) = MATMUL(Ia_mg,Iconj_new(n+1:2*n))
            BETA=0.
            call DGEMV('N',size(Ia_gv,1),size(Ia_gv,2),ALPHA,Ia_gv,size(Ia_gv,1),Iconj_new(1:n),1,BETA,kn_2a(1:n),1)
            call DGEMV('N',size(Ia_mg,1),size(Ia_mg,2),ALPHA,Ia_mg,size(Ia_mg,1),Iconj_new(n+1:2*n),1,BETA,kn_2a(n+1:2*n),1)

            kn_2b = DOT_PRODUCT(Iconj_new,kn_2a)
            !kn_2b = DDOT(size(Iconj_new),Iconj_new,1,kn_2a,1)

            kn_2c = DOT_PRODUCT(Iconj_new,mJg)
            !kn_2c = DDOT(size(Iconj_new),Iconj_new,1,mJg,1)

            kn = kn_1 / (kn_2b + kn_2c +1e-10) !para evitar divisiones por cero

            !Calculo m_new
            m_new = m_n - (kn*Iconj_new) !array operation


            !Condicion de paro
            !Convergencia
            conv_gv = 100 * SQRT( SUM( ((m_new(1:n)-m_n(1:n))**2)/((m_n(1:n)**2)+1e-10) )/n )
            conv_mg = 100 * SQRT( SUM( ((m_new(n+1:2*n)-m_n(n+1:2*n))**2)/((m_n(n+1:2*n)**2)+1e-10) )/n )

            !Actualizando variables para la siguiente iteración
            m_n = m_new
            I = I_new
            Iconj = Iconj_new


            !Obteniendo datos calculados
            !Am(1:m) = MATMUL(A_gv,m_new(1:n))
            !Am(m+1:2*m) = MATMUL(A_mg,m_new(n+1:2*n))
            BETA=0.
            call DGEMV('N',size(A_gv,1),size(A_gv,2),ALPHA,A_gv,size(A_gv,1),m_new(1:n),1,BETA,Am(1:m),1)
            call DGEMV('N',size(A_mg,1),size(A_mg,2),ALPHA,A_mg,size(A_mg,1),m_new(n+1:2*n),1,BETA,Am(m+1:2*m),1)

            !RMS
            dr(1:m) = Am(1:m) - d_gv !array operation
            dr(m+1:2*m) = Am(m+1:2*m) - d_mg
            dr(1:m) = Wd_gv * dr(1:m) !array operation
            dr(m+1:2*m) = Wd_mg * dr(m+1:2*m)
            !dr(1:m) = DOT_PRODUCT((Wd_gv**2),dr(1:m)) !tratando corregir estimacion
            !dr(m+1:2*m) = DOT_PRODUCT((Wd_mg**2),dr(m+1:2*m))

            rms_gv = DOT_PRODUCT(dr(1:m),dr(1:m))
            !rms_gv = DDOT(m,dr(1:m),1,dr(1:m),1)
            rms_gv = SQRT(rms_gv/m)

            rms_mg = DOT_PRODUCT(dr(m+1:2*m),dr(m+1:2*m))
            !rms_mg = DDOT(m,dr(m+1:2*m),1,dr(m+1:2*m),1)
            rms_mg = SQRT(rms_mg/m)

            WRITE(16,'(I3,",",F14.4,",",F14.4,",",F10.4,",",F10.4)') ii,conv_gv,conv_mg,rms_gv,rms_mg

            
            !Escribir resultados
            ! Escribir archivos de modelos invertidos
            WRITE (filename1, '("output_jointGram_gv_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=9,file=filename1,status='unknown') !crear y escribir en archivo
            WRITE(9,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'

            WRITE (filename2, '("output_jointGram_mg_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=10,file=filename2,status='unknown') !crear y escribir en archivo
            WRITE(10,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'

            do jj = 1,n,1
                  WRITE(9,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(jj)
                  WRITE(10,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(n+jj)
            end do

            CLOSE(unit=9)
            CLOSE(unit=10)

            ! Escribir archivos de datos calculados a partir de modelos invertidos
            WRITE (filename3, '("output_jointGram_gv_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=11,file=filename3,status='unknown') !crear y escribir en archivo
            WRITE(11,*) 'x[m],y[m],z[m],g[mGal]'

            WRITE (filename4, '("output_jointGram_mg_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=12,file=filename4,status='unknown') !crear y escribir en archivo
            WRITE(12,*) 'x[m],y[m],z[m],m[nT]'

            do jj = 1,m,1
                  WRITE(11,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_gv(jj),Am(jj)
                  WRITE(12,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_mg(jj),Am(m+jj)
            end do

            CLOSE(unit=11)
            CLOSE(unit=12)

      end if
end do

CLOSE(unit=16)

! Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
WRITE(0,*) 'Execution time: serial jointGram =', (finish-start), 'seg =', (finish-start)/60, 'min'
!close(unit=0)

return
end subroutine jointGram4
!*******************************************************************************


subroutine jointGram5(m,n,nx,ny,nz,step,xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell,d_gv,stdevd_gv,m_gv,stdevm_gv,mapr_gv,A_gv, &
                     d_mg,stdevd_mg,m_mg,stdevm_mg,mapr_mg,A_mg,param_reg,strike,num_iters,err0,ciclo)
! Subrutina que calcula inversion conjunta de dos modelos gravimétrico y magnético
! considerando similitud estructural evaluada por gradientes cruzados implementada mediante el enfoque Gramiano.

use iso_fortran_env , dp => real64
use auxiliar, only: cov, partialder, laplaciano, coord_voxel
!use model_ini, only: quasi_inf

implicit none

integer, value :: m, n 
integer, intent(in) :: num_iters, nx, ny, nz, ciclo
real, intent(in) :: step, err0, strike
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
character(len=42) :: filename00, filename1, filename2, filename3, filename4, filename5
real :: start, finish

!Evaluar eficiencia.
call cpu_time(start)
PRINT*, 'Inicia calculo de modelos en CPU'

!*****************************************
!reservar memoria para arreglos temporales
!allocate ( d_gv(m),d_mg(m),stdevd_gv(m),stdevd_mg(m),m_gv(n),m_mg(n),mapr_gv(n),mapr_mg(n),stdevm_gv(n),stdevm_mg(n), &
!          A_gv(m,n),A_mg(m,n),xcell(n), ycell(n), zcell(n), xobs(m), yobs(m), zobs_gv(m), zobs_mg(m) )
allocate ( dx(n,n),dy(n,n),dz(n,n),opL(n,n),opM(n,n),opL2(n,n),opM2(n,n),opN2(n,n) )

!matrices para derivadas
call partialder(nx,ny,nz,step,'x',dx)
call partialder(nx,ny,nz,step,'y',dy)
call partialder(nx,ny,nz,step,'z',dz)

!Operador Suavidad
call laplaciano(nx,ny,nz,step,opL)
opL2 = MATMUL(transpose(opL),opL)

!Operador variacion horizontal total
opM = dx*cos(strike) + dy*sin(strike)
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


WRITE (filename00, '("output_misfit_jointGram_c",I1,".txt")' ) ciclo
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

            kn = kn_1 / (kn_2b + kn_2c +1e-10) !para evitar divisiones por cero

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
            
            WRITE(16,'(I3,",",F14.4,",",F14.4,",",F10.4,",",F10.4)') ii,conv_gv,conv_mg,rms_gv,rms_mg


            !Escribir resultados
            ! Escribir archivos de modelos invertidos
            WRITE (filename1, '("output_jointGram_gv_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=9,file=filename1,status='unknown') !crear y escribir en archivo
            WRITE(9,*) 'x[m],y[m],z[m],Density_Contrast[g/cm**3]'

            WRITE (filename2, '("output_jointGram_mg_minv_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=10,file=filename2,status='unknown') !crear y escribir en archivo
            WRITE(10,*) 'x[m],y[m],z[m],Magnetization_Contrast[A/m]'

            do jj = 1,n,1
                  WRITE(9,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(jj)
                  WRITE(10,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xcell(jj),ycell(jj),zcell(jj),m_new(n+jj)
            end do

            CLOSE(unit=9)
            CLOSE(unit=10)

            ! Escribir archivos de datos calculados a partir de modelos invertidos
            WRITE (filename3, '("output_jointGram_gv_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=11,file=filename3,status='unknown') !crear y escribir en archivo
            WRITE(11,*) 'x[m],y[m],z[m],g[mGal]'

            WRITE (filename4, '("output_jointGram_mg_anomaly_c",I1,"_iter",I1,".csv")' ) ciclo, ii
            OPEN(unit=12,file=filename4,status='unknown') !crear y escribir en archivo
            WRITE(12,*) 'x[m],y[m],z[m],m[nT]'

            do jj = 1,m,1
                  WRITE(11,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_gv(jj),Am(jj)
                  WRITE(12,'(F10.3,",",F11.3,",",F10.3,",",E11.4)') xobs(jj),yobs(jj),zobs_mg(jj),Am(m+jj)
            end do

            CLOSE(unit=11)
            CLOSE(unit=12)


      end if
end do

CLOSE(unit=16)

! Evaluar eficiencia.
call cpu_time(finish)

OPEN(unit=0,file='output_time.txt',status='old')
WRITE(0,*) 'Execution time: CPU jointGram =', (finish-start), 'seg =', (finish-start)/60, 'min'

deallocate ( dxm1,dym1,dzm1,dxm2,dym2,dzm2,diva,divb,divc,dive, &
             ter1,ter2,ter3,dum6,dum7,dum8,ter1J,ter2J, &
             dxI1,dyI1,dzI1,dxI2,dyI2,dzI2 )
deallocate ( dx,dy,dz,mIg,mJg,I_new)
deallocate ( A_gv,A_mg,d_gv,d_mg, xobs,yobs,zobs_gv,zobs_mg,xcell,ycell,zcell )

!close(unit=0)

return
end subroutine jointGram5
!*******************************************************************************


end module inversion_serial
!*******************************************************************************
! para compilar
!nvfortran -c inversion_serial.f08 -o inversion_serial.o

