!
! ################## Self-interaction term of EWALD sum #######################
! 
SUBROUTINE EWALD_SELF(NS, q, param, sys)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(PM):: param
TYPE(ST):: sys
!
! Internal variables
INTEGER:: i, id
REAL*8 :: PI, q1
!
!
PI = ATAN(1.0)*4.0
param%a = SQRT(PI/(sys%box(1)*sys%box(2)*sys%box(3))**(1./3.))
sys%Uo = 0.0
!
!
DO i = 1, NS%Npt
   id = q(i)%id
   q1 = param%q(id)
   sys%Uo = sys%Uo - q1**2
END DO
sys%Uo = sys%Uo * param%a/SQRT(PI)
RETURN
END SUBROUTINE EWALD_SELF
!
! ##################### LJ potential for water molecule #######################
SUBROUTINE FORCE(NS, q, mc, param, sys)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
!
! Internal variables
INTEGER:: i, j, k, ii, jj
REAL*8:: r2, r2i, r6i, lj, xr(3),  x
!
! Initialize
DO i=1, NS%Npt
   q(i)%ff(:) = 0.0
END DO
sys%Epot = 0.0
!
DO i=1, NS%Nmc-1
   ii = mc(i)%comp(2)
   DO j=i+1, NS%Nmc
      jj = mc(j)%comp(2)
      r2 = 0.0
      DO k=1,3
         x = q(ii)%xx(k) - q(jj)%xx(k)
         xr(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         r2 = r2 + xr(k)**2
      END DO
      IF (r2 < param%rc2) THEN
         r2i= 1./r2
         r6i = param%sig**6*r2i**3
         lj = 48.*r2i*r6i*(r6i - .5)
         DO k = 1,3
            q(ii)%ff(k) = q(ii)%ff(k) + lj*xr(k)*param%eps
            q(jj)%ff(k) = q(jj)%ff(k) - lj*xr(k)*param%eps
         END DO
         sys%Epot = sys%Epot + 4.*r6i*(r6i - 1.)*param%eps - param%ecut
      END IF
   END DO
END DO
!
! long range potential
CALL COULOMB(NS, q, mc, param, sys)
!
RETURN
!
END SUBROUTINE FORCE
!
! ################ Coulomb energy and force using EWALD sum ###################
!
!
SUBROUTINE COULOMB(NS, q, mc, param, sys)
USE DATASTR
IMPLICIT NONE
!
! Charge unit: e = 1.6022E-19 C
! Length unit: A = 1.E-10 m
! Mass unit: amu = 1.6605E-27 kg 
! time unit: 1 = 10fs
! Below C constant is 1/4pi e_0. C times 1/r whose has Angstrom unit will be
! energy in eV
!
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
!
!
REAL*8,  PARAMETER:: C = 14.39983566 ! in eV*A = 1/8.8542E-12C^2/N/m^2/4/pi
INTEGER:: nx, ny, nz, i, j, k, mx, my, mz, id1, id2, Mmax
REAL*8 :: r2, q1, q2, PI, x, xr(3), erfcc, V, rc2, rx, cm, sm, s, m2, hPI, &
     df, r
!
!
Mmax = max_iteration
rc2 = param%rc2
sys%Ur = 0.0
sys%Um = 0.0
PI = ATAN(1.0)*4.0
hPI = SQRT(PI)
V = sys%box(1)*sys%box(2)*sys%box(3)
!
! Direct - real space sum for EWALD energy and force 
!
! Large simulation cell w.r.t cutoff radius
! 
DO i=1, NS%Npt-1
   id1 = q(i)%id
   q1 = param%q(id1)
   DO j=i+1, NS%Npt
      id2 = q(j)%id
      q2 = param%q(id2)
!      rss = RS(Nxid(i),Nxid(j)) <= de Broglie parameter
      r2 = 0.0
      DO k=1,3
         x = q(i)%xx(k) - q(j)%xx(k)
         xr(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         r2 = r2 + xr(k)**2
      END DO
      IF (r2 < rc2) THEN
         r = SQRT(r2)
         sys%Ur = sys%Ur + q1*q2*erfcc(param%a*r)/r 
         df = C*q1*q2*(2.*EXP(-param%a**2*r2)*param%a*r/hPI + &
              erfcc(param%a*r))/r2/r
         DO k=1,3
            q(i)%ff(k) = q(i)%ff(k) + df*xr(k)
            q(j)%ff(k) = q(j)%ff(k) - df*xr(k)
         END DO
      END IF
   END DO
END DO
!
! Reciprocal sum
DO mx = -Mmax,Mmax
   DO my = -Mmax,Mmax
      DO mz = -Mmax,Mmax
         IF ( ((mx ==0) .AND. (my == 0)) .AND. (mz ==0)) GOTO 30
         sm = 0.0
         cm = 0.0
         DO i=1, NS%Npt
            id1 = q(i)%id
            q1 = param%q(id1)
            s  = REAL(mx)*q(i)%xx(1)/sys%box(1) + &
                 REAL(my)*q(i)%xx(2)/sys%box(2) + &
                 REAL(mz)*q(i)%xx(3)/sys%box(3) 
            cm = cm + q1*COS(2.*PI*s)
            sm = sm + q1*SIN(2.*PI*s)
         END DO
         m2 = REAL(mx**2)/sys%box(1)**2 + REAL(my**2)/sys%box(2)**2 + &
              REAL(mz**2)/sys%box(3)**2
         sys%Um = sys%Um + EXP(-m2*PI**2/param%a**2)*(cm**2 + sm**2)/m2
         DO i=1, NS%Npt
            id1 = q(i)%id
            q1 = param%q(id1)
            s  = REAL(mx)*q(i)%xx(1)/sys%box(1) + &
                 REAL(my)*q(i)%xx(2)/sys%box(2) + &
                 REAL(mz)*q(i)%xx(3)/sys%box(3) 
            df = C*2.*q1*EXP(-m2*PI**2/param%a**2)* &
                 (SIN(2.*PI*s)*cm - COS(2.*PI*s)*sm)/m2/V
            q(i)%ff(1) = q(i)%ff(1) + df*REAL(mx)/sys%box(1)
            q(i)%ff(2) = q(i)%ff(2) + df*REAL(my)/sys%box(2)
            q(i)%ff(3) = q(i)%ff(3) + df*REAL(mz)/sys%box(3)
         END DO
30    END DO
   END DO
END DO
sys%Um = sys%Um/PI/V/2.
sys%Epot = sys%Epot + (sys%Ur + sys%Uo + sys%Um)*C
!
RETURN
END SUBROUTINE COULOMB
!
!
FUNCTION erfcc(x)
IMPLICIT NONE
REAL*8:: erfcc, x
REAL*8:: t,z
z = abs(x)
t = 1./(1.+.5*z)
erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
        t*(0.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
        t*(1.48851587+t*(-0.82215223+t*.17087277)))))))))
IF (x.LT.0.) erfcc = 2.-erfcc

RETURN
END FUNCTION
!
! ############# Constant pressure implementation ##########
SUBROUTINE BAROSTAT(NS, q, mc, sys, param)
USE DATASTR
IMPLICIT NONE
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(ST):: sys
TYPE(PM):: param
!
INTEGER:: i, j, k, ii, jj
REAL*8 :: V, P, chi, x, xr(3), ff(3)
V = sys%box(1)*sys%box(2)*sys%box(3)
P = 0.
DO i=1, NS%Nmc-1
   ii = mc(i)%comp(2)
   DO j=i+1, NS%Nmc
      jj = mc(j)%comp(2)
      DO k=1,3
         x = q(ii)%xx(k) - q(jj)%xx(k)
         xr(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         ff(k) = q(ii)%ff(k) - q(jj)%ff(k)
         P = P + xr(k)*ff(k)
      END DO
   END DO
END DO
P = (P + 2.*sys%mv2)/3./V
!GOTO 100
chi = 1. - 0.1*(param%P - P)
chi = chi**(1./3.)
DO i=1, NS%Npt
   DO k=1,3
      q(i)%xx(k) = q(i)%xx(k)*chi
   END DO
END DO
sys%box(:) = sys%box(:)*chi
100 sys%P = P*1.58117E6
!
RETURN
END SUBROUTINE BAROSTAT
         
