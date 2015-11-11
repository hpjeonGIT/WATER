!
! ####### Verlet time integration routine with Thermal Noise #########
! 
! Velocity verlet with RATTLE - stochastic thermostat implemented
!
SUBROUTINE VVerletStoch1(NS, q, mc, param, sys, crit, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Nmax, Ncount, Nbond(3,2)
REAL*8:: x, distance, Errmax, gamma, sigma, s(3), xm(3), mu(3), d2(3), &
     qt(3,3), rr(3,3)
!
! Initialization
sigma = 1.       ! Deviation for Gaussian distribution
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
!
! RATTLE 1st stage
d2(1) = 1.0
d2(2) = 2.6666258
d2(3) = 1.0
xm(1) = param%xm(1)
xm(2) = param%xm(2)
xm(3) = xm(1)
mu(1) = xm(1)*xm(2)/(xm(1)+xm(2))
mu(2) = xm(1)/2.
mu(3) = mu(1)
DO n=1, NS%Nmc
   !
   ! Allocate bond pair, relative distance, and mass   
   id1 = mc(n)%comp(1)
   id2 = mc(n)%comp(2)
   id3 = mc(n)%comp(3)
   Nbond(1,1) = 1
   Nbond(1,2) = 2
   Nbond(2,1) = 1
   Nbond(2,2) = 3
   Nbond(3,1) = 2
   Nbond(3,2) = 3
   DO k=1,3
      qt(1,k) = q(id1)%xv(k)*(1.-0.5*param%alpha*dt/xm(1)) + &
           0.5*dt*q(id1)%ff(k)/xm(1)
      qt(2,k) = q(id2)%xv(k)*(1.-0.5*param%alpha*dt/xm(2)) + &
           0.5*dt*q(id2)%ff(k)/xm(2)
      qt(3,k) = q(id3)%xv(k)*(1.-0.5*param%alpha*dt/xm(1)) + &
           0.5*dt*q(id3)%ff(k)/xm(1)
      x = q(id2)%xx(k) - q(id1)%xx(k)
      rr(1,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H-O bond
      x = q(id3)%xx(k) - q(id1)%xx(k)
      rr(2,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H-O-H' angle
      x = q(id3)%xx(k) - q(id2)%xx(k)
      rr(3,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H'-O bond
   END DO
   !
   ! Iteration loop for gamma
   DO m=1,Nmax
      Ncount = 0
      DO i=1,3
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            s(k) = rr(i,k) + dt*(qt(id2,k) - qt(id1,k))
         END DO
         distance = s(1)**2 + s(2)**2 + s(3)**2 - d2(i)
         gamma = 0.5*mu(i)*distance/dt/ &
              (s(1)*rr(i,1)+s(2)*rr(i,2)+s(3)*rr(i,3))
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1         
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == 3) GOTO 10
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 1st stage"
   STOP
   !
   ! Position update. Intermediate velocity is available now.
10 DO i=1, 3
      DO k=1,3
         x = q(mc(n)%comp(i))%xx(k) + dt*qt(i,k)
         q(mc(n)%comp(i))%xx(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         q(mc(n)%comp(i))%xv(k) = qt(i,k)
      END DO
   END DO
   !
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch1
!
!
SUBROUTINE VVerletStoch2(NS, q, mc, param, sys, crit, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Nmax, Ncount, Nbond(3,2)
REAL*8:: x, distance, Errmax, gamma,  sigma, s(3), xm(3), mu(3), d2(3), &
     qt(3,3), rr(3,3), T0, eta, beta
!
! Initialization
sigma = 1.       ! Deviation for Gaussian distribution
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
beta = SQRT(2.*param%alpha*kb*sys%temp/dt)
!
! RATTLE 1st stage
d2(1) = 1.0
d2(2) = 2.6666258
d2(3) = 1.0
xm(1) = param%xm(1)
xm(2) = param%xm(2)
xm(3) = xm(1)
mu(1) = xm(1)*xm(2)/(xm(1)+xm(2))
mu(2) = xm(1)/2.
mu(3) = mu(1)
DO n=1, NS%Nmc
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, 3
      id1 = mc(n)%comp(i)
      DO k=1,3
         eta = fluct(sigma)
         q(id1)%ff(k) = q(id1)%ff(k) + eta*beta
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   !
   ! Iteration loop for gamma
   Nbond(1,1) = 1
   Nbond(1,2) = 2
   Nbond(2,1) = 1
   Nbond(2,2) = 3
   Nbond(3,1) = 2
   Nbond(3,2) = 3
   DO m=1,Nmax
      Ncount = 0
      DO i=1,3
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            x = q(mc(n)%comp(id2))%xx(k) - q(mc(n)%comp(id1))%xx(k)
            rr(i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond distance
         END DO
         distance = rr(i,1)*(qt(id2,1)-qt(id1,1)) + &
              rr(i,2)*(qt(id2,2)-qt(id1,2)) + rr(i,3)*(qt(id2,3)-qt(id1,3))
         gamma = mu(i)*distance/d2(i)
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == 3) GOTO 20
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 2nd stage"
   STOP
   !
   ! Velocity update (phew! ~~~)
20 DO i=1, 3
      id1 = mc(n)%comp(i)
      DO k=1,3
         q(id1)%xv(k) = qt(i,k)/(1.+0.5*param%alpha*dt/xm(i))
         sys%mv2 = sys%mv2 + xm(i)*q(id1)%xv(k)**2
      END DO
   END DO
   !
END DO
!
RETURN
!
END SUBROUTINE VVerletStoch2
!
! ####### Pure Verlet time integration routine with No Thermostat #########
! 
! Velocity Verlet algorithm with RATTLE - no thermostat implemented
!
SUBROUTINE VVerletNotemp1(NS, q, mc, param, sys, crit, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Nmax, Ncount, Nbond(3,2)
REAL*8:: x, distance, Errmax, gamma,  s(3), xm(3), mu(3), d2(3), qt(3,3),&
     rr(3,3)
!
! Initialization
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
!
! RATTLE 1st stage
d2(1) = 1.0
d2(2) = 2.6666258
d2(3) = 1.0
xm(1) = param%xm(1)
xm(2) = param%xm(2)
xm(3) = xm(1)
mu(1) = xm(1)*xm(2)/(xm(1)+xm(2))
mu(2) = xm(1)/2.
mu(3) = mu(1)
DO n=1, NS%Nmc
   !
   ! Allocate bond pair, relative distance, and mass   
   id1 = mc(n)%comp(1)
   id2 = mc(n)%comp(2)
   id3 = mc(n)%comp(3)
   Nbond(1,1) = 1
   Nbond(1,2) = 2
   Nbond(2,1) = 1
   Nbond(2,2) = 3
   Nbond(3,1) = 2
   Nbond(3,2) = 3
   DO k=1,3
      qt(1,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(1)
      qt(2,k) = q(id2)%xv(k) + 0.5*dt*q(id2)%ff(k)/xm(2)
      qt(3,k) = q(id3)%xv(k) + 0.5*dt*q(id3)%ff(k)/xm(1)
      x = q(id2)%xx(k) - q(id1)%xx(k)
      rr(1,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H-O bond
      x = q(id3)%xx(k) - q(id1)%xx(k)
      rr(2,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H-O-H' angle
      x = q(id3)%xx(k) - q(id2)%xx(k)
      rr(3,k) =  x - sys%box(k)*ANINT(x/sys%box(k)) ! H'-O bond
   END DO
   !
   ! Iteration loop for gamma
   DO m=1,Nmax
      Ncount = 0
      DO i=1,3
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            s(k) = rr(i,k) + dt*(qt(id2,k) - qt(id1,k))
         END DO
         distance = s(1)**2 + s(2)**2 + s(3)**2 - d2(i)
         gamma = 0.5*mu(i)*distance/dt/ &
              (s(1)*rr(i,1)+s(2)*rr(i,2)+s(3)*rr(i,3))
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1         
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == 3) GOTO 10
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 1st stage"
   STOP
   !
   ! Position update. Intermediate velocity is available now.
10 DO i=1, 3
      DO k=1,3
         x = q(mc(n)%comp(i))%xx(k) + dt*qt(i,k)
         q(mc(n)%comp(i))%xx(k) = x - sys%box(k)*ANINT(x/sys%box(k))
         q(mc(n)%comp(i))%xv(k) = qt(i,k)
      END DO
   END DO
   !
END DO
RETURN
END SUBROUTINE VVerletNotemp1
!
SUBROUTINE VVerletNotemp2(NS, q, mc, param, sys, crit, dt)
USE DataStr
IMPLICIT NONE
INTERFACE
   FUNCTION fluct(x)
     IMPLICIT NONE
     REAL*8:: fluct, x, l0, l, r, v1, v2
     REAL*8:: rand1, rand2, ran2
   END FUNCTION fluct
END INTERFACE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(CH):: mc(NS%Nmc)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, m, id1, id2, id3, n, Nmax, Ncount, Nbond(3,2)
REAL*8:: x, distance, Errmax, gamma,  ss(3), xm(3), mu(3), d2(3), qt(3,3), &
     rr(3,3)
!
! 
sys%mv2 = 0.0
Nmax = crit%Nshake
Errmax = crit%Errshake
d2(1) = 1.0
d2(2) = 2.6666258
d2(3) = 1.0
xm(1) = param%xm(1)
xm(2) = param%xm(2)
xm(3) = xm(1)
mu(1) = xm(1)*xm(2)/(xm(1)+xm(2))
mu(2) = xm(1)/2.
mu(3) = mu(1)
DO n=1, NS%Nmc
   !
   ! Allocate bond pair, relative distance, and mass
   DO i=1, 3
      id1 = mc(n)%comp(i)
      DO k=1,3
         qt(i,k) = q(id1)%xv(k) + 0.5*dt*q(id1)%ff(k)/xm(i)
      END DO
   END DO
   !
   ! Iteration loop for gamma
   Nbond(1,1) = 1
   Nbond(1,2) = 2
   Nbond(2,1) = 1
   Nbond(2,2) = 3
   Nbond(3,1) = 2
   Nbond(3,2) = 3
   DO m=1,Nmax
      Ncount = 0
      DO i=1,3
         id1 = Nbond(i,1)
         id2 = Nbond(i,2)
         DO k=1,3
            x = q(mc(n)%comp(id2))%xx(k) - q(mc(n)%comp(id1))%xx(k)
            rr(i,k) = x - sys%box(k)*ANINT(x/sys%box(k)) ! bond distance
         END DO
         distance = rr(i,1)*(qt(id2,1)-qt(id1,1)) + &
              rr(i,2)*(qt(id2,2)-qt(id1,2)) + rr(i,3)*(qt(id2,3)-qt(id1,3))
         gamma = mu(i)*distance/d2(i)
         IF (ABS(gamma) < Errmax ) Ncount = Ncount + 1
         DO k=1,3
            qt(id2,k) = qt(id2,k) - gamma*rr(i,k)/xm(id2)
            qt(id1,k) = qt(id1,k) + gamma*rr(i,k)/xm(id1)
         END DO
      END DO
      IF (Ncount == 3) GOTO 20
   END DO
   PRINT *, "RATTLE is not converged at", n,"th chain for 2nd stage"
   STOP
   !
   ! Velocity update (phew! ~~~)
20 DO i=1, 3
      id1 = mc(n)%comp(i)
      DO k=1,3
         q(id1)%xv(k) = qt(i,k)
         sys%mv2 = sys%mv2 + xm(i)*q(id1)%xv(k)**2
      END DO
   END DO
   !
END DO
RETURN
!
END SUBROUTINE VVerletNotemp2
!
! ############## Random Gaussian(Normal) Distribution Function ################
! 
! For stochastic thermostat, fluctuation dissipation theorem is implemented.
! Basically, random number generator which follows Gaussian distribution is
! needed to implement this thermal noise force.
! Random number is generated using FORTRAN intrinsic fucntion - RANDOM_SEED
! and RANDOM_NUMBER. But during the implementation, it is found that those
! intrinsic functions may not work well under false seeding - to use this
! routine on new machine or new compiler, please make sure that the 
! distribution follows zero-mean with suitable deviation.
!
! This function provides a random number with zero-mean and deviation x 
! along Gaussian distribution.
! <p> = 0.0
! <p**2> = x**2
! 
FUNCTION fluct(x)
  IMPLICIT NONE
  REAL*8:: fluct, x, r, v1, v2
  REAL*8:: rand1, rand2, ran2
  REAL:: ranf
  !
  ! Initialization
  r=1.
  DO WHILE (r.ge.1.)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.*rand1 - 1.
     v2 = 2.*rand2 - 1.
     r = v1*v1+v2*v2
  END DO
  fluct = v1*SQRT(-2.*log(r)/r)*x
  RETURN
END FUNCTION fluct
