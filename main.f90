PROGRAM WATER_SPC
!
! Water molecule modeling using SPC model
!
USE DATASTR
IMPLICIT NONE
!
TYPE(NM):: NS
TYPE(PT), POINTER:: q(:)
TYPE(CH), POINTER:: mc(:)
TYPE(PM):: param
TYPE(ST):: sys
TYPE(TM):: time
TYPE(LM):: crit
REAL*8  :: dt
!
!
! INTERNAL VARIABLES
REAL   :: time0, time1, time2, secnds
INTEGER:: Nloop_max
time1 = secnds(time0)
!
! Open XYZ file for output
OPEN(UNIT=55, file="energy.dat")
WRITE(55,100)
100 FORMAT("# time(fs) instant. temp. Kinetic Eenrgy, Potential Energy, Total energy, Pressure(atm)")
!
! Parsing input data and allocate pointer variables
CALL Init(NS, q, mc, time, param, sys, crit, dt)
!
! Initialization
time%Ndump = INT((time%tdump+dt*.5)/dt)
time%Nrest = INT((time%trest+dt*.5)/dt)
time%Nsamp = INT((time%tsamp+dt*.5)/dt)
time%Nloop = 0
time%tnow = 0.0
Nloop_max = INT((time%tmax+dt*.5)/dt)
CALL Vinit(NS, q, param, sys)
CALL EWALD_SELF(NS, q, param, sys)
!
! Stochastic thermostat loop
IF (param%thermo == 'STOCH ')  THEN
   DO WHILE (time%Nloop < Nloop_max)
      time%tnow = time%tnow + dt
      time%Nloop = time%Nloop + 1
      !
      ! L.2 Time integration with constraint dynamics
      CALL VVerletStoch1(NS, q, mc, param, sys, crit, dt)
      CALL Force(NS, q, mc, param, sys)
      CALL VVerletStoch2(NS, q, mc, param, sys, crit, dt)
      CALL BAROSTAT(NS, q, mc, sys, param)
      IF (MOD(time%Nloop,time%Nrest) == 0) THEN
         CALL RestartBIN(NS, q, time)
      END IF
      IF (MOD(time%Nloop,time%Ndump) == 0) THEN
         Call DUMP(NS, q, time)
      END IF
      IF (MOD(time%Nloop,time%Nsamp) == 0) THEN
         WRITE(55,200) time%tnow*10.18, sys%mv2/REAL(NS%Npt*3-NS%Ncn-6)/kb, &
              sys%mv2, sys%Epot, sys%mv2*0.5+sys%Epot, sys%P
      END IF
END DO
!
! Microcanonical(NVE) ensemble loop
ELSE IF (param%thermo == 'NOTEMP')  THEN
   DO WHILE (time%Nloop < Nloop_max)
      time%tnow = time%tnow + dt
      time%Nloop = time%Nloop + 1
      CALL VVerletNotemp1(NS, q, mc, param, sys, crit, dt)
      CALL Force(NS, q, mc, param, sys)
      CALL VVerletNotemp2(NS, q, mc, param, sys, crit, dt)
      CALL BAROSTAT(NS, q, mc, sys, param)
      !
      IF (MOD(time%Nloop,time%Nrest) == 0) THEN
         CALL RestartBIN(NS, q, time)
      END IF
      IF (MOD(time%Nloop,time%Ndump) == 0) THEN
         Call DUMP(NS, q, time)
      END IF
      IF (MOD(time%Nloop,time%Nsamp) == 0) THEN
         WRITE(55,200) time%tnow*10.18, sys%mv2/REAL(NS%Npt*3-NS%Ncn-6)/kb, &
              sys%mv2, sys%Epot, sys%mv2*0.5+sys%Epot, sys%P
      END IF
   END DO
!
!
ELSE
   PRINT *, "Unknown Thermostat"
END IF
!
PRINT *, sys%box
CLOSE(25)
CLOSE(55)
200 FORMAT(6(ES14.6, 1X))
DEALLOCATE(q, mc)
time2 = secnds(time1)
PRINT *, "Wall time is", time2
!
STOP
!
CONTAINS
!
! ################### Data parsing and memory allocation ####################
! Unit normalization
! Length: 1. means 1Angstrom = 10E-10 m
! Mass: 1. means 1.6605E-27 kg (a.m.u)
! Energy: 1. means 1 eV = 1.6022E-19 J
! Time: 1. means 10.18 fs = 1.018E-14 sec
SUBROUTINE INIT(NS, q, mc, time, param, sys, crit, dt)
USE DATASTR
IMPLICIT NONE
TYPE(PT), POINTER:: q(:)
TYPE(CH), POINTER:: mc(:)
TYPE(TM):: time
TYPE(PM):: param
TYPE(ST):: sys
TYPE(LM):: crit
TYPE(NM):: NS
REAL*8  :: dt
!
! INTERNAL VARIABLES
INTEGER:: i, j, k, length
CHARACTER(len=5):: dummy
REAL*8 :: temp
REAL   :: rand
!
!
![[[[[[[[[[[[[[[[[ "input.dat" parsing and data arrangement ]]]]]]]]]]]]]]]]]
OPEN (UNIT=11, file="control.h2o")
READ(11,*) dummy
!
! Time parameter
! 1.0 = 10.18fs
READ(11,*) dummy
READ(11,*) time%tmax, time%trest, time%tdump, time%tsamp, dt
time%tmax  = time%tmax  / 10.18
time%trest = time%trest / 10.18
time%tdump = time%tdump / 10.18
time%tsamp = time%tsamp / 10.18
dt = dt / 10.18
READ(11,*) dummy
READ(11,*) (sys%box(j), j=1,3)
READ(11,*) dummy
READ(11,*) sys%temp
!
! i.1 Number of all particles
READ(11,*) dummy
READ(11,*) NS%Npt
!
!
READ(11,*) dummy
READ(11,*) param%eps
! kJ/mol -> eV/ea
param%eps = param%eps/96.4852912
!
! 
READ(11,*) dummy
READ(11,*) param%sig
!
! p.3 Square of critical radius for LJ(12-6) potential
READ(11,*) dummy
READ(11,*) param%rc2
param%rc2 = param%rc2**2
!
! p.4 Energy cut-off for LJ(12-6) potential
param%ecut = param%eps*4.*(param%sig**6/param%rc2**3 - 1.)*param%sig**6 &
     /param%rc2**3
!
! p.5 Mass for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%xm(i)
END DO
!
! p.5 Charge for each particle kind
READ(11,*) dummy
DO i=1,Nparam
   READ(11,*) param%q(i)
END DO
!
! p.7 Berendsen thermostat/damping constant
READ(11,*) dummy
READ(11,*) param%thermo
READ(11,*) param%tau, param%alpha
!
! p.8 Iteration limit
READ(11,*) dummy
READ(11,*) crit%Nshake, crit%Errshake
!
READ(11,*) dummy
READ(11,*) param%P
param%P = param%P/1.58117E6
CLOSE(11)
!
!#################### READ particle position and id data ######################
!
OPEN (UNIT=15, file="water.bin", FORM="UNFORMATTED", STATUS = "OLD")
NS%Nmc = NS%Npt / 3
ALLOCATE(q(NS%Npt),mc(NS%Nmc))
!
DO i=1, NS%Nmc
   j = (i-1)*3 + 1
   READ(15) q(j)%xx(1), q(j)%xx(2), q(j)%xx(3)
   q(j)%id = 1
   mc(i)%comp(1) = j
   j = (i-1)*3 + 2
   READ(15) q(j)%xx(1), q(j)%xx(2), q(j)%xx(3)
   q(j)%id = 2
   mc(i)%comp(2) = j
   j = i*3 
   READ(15) q(j)%xx(1), q(j)%xx(2), q(j)%xx(3)
   q(j)%id = 1
   mc(i)%comp(3) = j
END DO
NS%Ncn = NS%Nmc*3
CLOSE(15)
!
END SUBROUTINE INIT
!
!
END PROGRAM WATER_SPC
