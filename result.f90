!
! ############ binary restart file generation ###################
SUBROUTINE RestartBIN(NS, q, time)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(TM):: time
!
!
INTEGER:: i, j, Nfreq
CHARACTER(LEN=10):: NUMBER, FILENAME
DATA NUMBER/'0123456789'/
FILENAME = "REST00.bin"
!
Nfreq = INT(time%Nloop/time%Nrest)
IF(INT(Nfreq/10) < 1 ) THEN
   FILENAME(5:5) = NUMBER(1:1)
   FILENAME(6:6) = NUMBER(Nfreq+1:Nfreq+1)
ELSE
   FILENAME(5:5) = NUMBER(INT(Nfreq/10)+1:INT(Nfreq/10)+1)
   FILENAME(6:6) = NUMBER(Nfreq-INT(Nfreq/10)*10+1:Nfreq-INT(Nfreq/10)*10+1)
END IF
OPEN(UNIT=30,FILE=FILENAME, FORM="UNFORMATTED")
!
DO i=1,NS%Npt
   WRITE(30) (q(i)%xx(j), j=1,3)  
END DO
!
CLOSE(30)
RETURN
END SUBROUTINE RestartBIN
!
! ################### RDF print ####################
SUBROUTINE DUMP(NS, q, time)
USE DATASTR
IMPLICIT NONE
!
! COMMUNICATION VARIABLES
TYPE(NM):: NS
TYPE(PT):: q(NS%Npt)
TYPE(TM):: time
!
!
INTEGER:: i, j, Nfreq
REAL*8 :: Ntau
CHARACTER(LEN=10):: NUMBER, FILENAME
DATA NUMBER/'0123456789'/
FILENAME = "RDF000.dat"
!
Ntau = REAL(time%Ndump)
Nfreq = INT(time%Nloop/time%Ndump)
IF(INT(Nfreq/10) < 1 ) THEN
   FILENAME(5:5) = NUMBER(1:1)
   FILENAME(6:6) = NUMBER(Nfreq+1:Nfreq+1)
ELSE
   FILENAME(5:5) = NUMBER(INT(Nfreq/10)+1:INT(Nfreq/10)+1)
   FILENAME(6:6) = NUMBER(Nfreq-INT(Nfreq/10)*10+1:Nfreq-INT(Nfreq/10)*10+1)
END IF
!
OPEN(UNIT=30,file=FILENAME)
WRITE(30,100)
DO i=1, Ndist+1
   WRITE(30,200) NS%dx*REAL(i), NS%RDF(i)/Ntau
END DO
CLOSE(30)
!
100 FORMAT("# R_ei,      RDF_ei ")
200 FORMAT(2(ES14.6, 1X))
!
RETURN
END SUBROUTINE DUMP
!
! ############ Initial velocity routine ####################
SUBROUTINE Vinit(NS, q, param, sys)
USE DATASTR
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
TYPE(NM)::NS
TYPE(PT)::q(NS%Npt)
TYPE(PM)::param
TYPE(ST)::sys
!
INTEGER:: i, j, k, id, Nr
REAL*8:: xm, lambda, T0, xv(NS%Npt,3)
!
sys%mv2 = 0.0
IF (param%thermo == 'STOCH ') THEN
   Nr = 0
ELSE
   Nr = 4
END IF
DO i=1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   DO j=1, 3
      xv(i,j) = fluct(SQRT(sys%temp))
      sys%mv2 = sys%mv2 + xm*xv(i,j)**2
   END DO
END DO
T0 = sys%mv2/kb/real(3*NS%Npt-NS%Ncn - Nr)
lambda = SQRT(sys%temp/T0)
sys%mv2 = 0.0
DO i = 1, NS%Npt
   id = q(i)%id
   xm = param%xm(id)
   DO j = 1,3
      q(i)%xv(j) = xv(i,j)*lambda
      sys%mv2 = sys%mv2 + xm*q(i)%xv(j)**2
   END DO
END DO
!
RETURN
END SUBROUTINE Vinit
