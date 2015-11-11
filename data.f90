MODULE DATASTR
IMPLICIT NONE
INTEGER, PARAMETER:: Nparam = 2, max_iteration = 10, Ndist = 1000
REAL*8,  PARAMETER:: kb = 8.617343E-5
!
! Particle data type
! xx = current position
! xv = current velocity
! ff = current force
! id = particle type
TYPE PT
   REAL*8:: xx(3), xv(3), ff(3)
   INTEGER:: id
END TYPE PT
!
! Water molecule chain
! comp = component particle
TYPE CH
   INTEGER:: comp(3)
END TYPE CH
!
! Time data type
! Nloop = Number of loops
! Ndump = Number of dumps
! Nsamp = Number of samplings
! tmax = Maximum simulation time
! tdump = xyz file dump time
! tsamp = sampling time
! tnow = Current time
TYPE TM
   INTEGER:: Nloop, Ndump, Nsamp, Nrest
   REAL*8:: tmax, tdump, tsamp, trest, tnow
END TYPE TM
!
! Parameter data type
! eps = Epsilon for LJ(12-6) potential
! sigma = Sigma for LJ(12-6) potential
! rc2 = Square of critical radius
! ecut = Energy cut-off
! xm = Mass of particle
! kt = Torsion stiffness
! tau = Berendsen thermostat parameter
TYPE PM
   REAL*8:: eps, sig, rc2, ecut, xm(Nparam), q(Nparam), tau, alpha, a, P
   CHARACTER*6:: thermo
END TYPE PM
!
! System variable data type
! box = Size of simulation cell
! temp = Temperature of (NVT) simulation
! mv2 = Sum of mv^2(=twice of kinetic energy)
! Epot = Potential energy of system
! Etor = Torsional energy of system
! Ecc  = CHx - CHx interaction energy
TYPE ST
   REAL*8:: box(3), temp, mv2, Epot, Etor, Ecc, Uo, Um, Ur, P
END TYPE ST
!
! Iteration limit data type
! Nshake = Maximum number of SHAKE iterations
! Errshake = Error limit of SHAKE iterations
TYPE LM
   INTEGER:: Nshake
   REAL*8:: Errshake
END TYPE LM
!
! Various number type
! Main objective is transfer specific numbers into print routine
! Npt = Number of particles
! Nmc = Number of molecules
! Ncn = Number of constraints
! Other sampling variables
TYPE NM
   INTEGER:: Npt, Nmc, Ncn
   REAL*8:: temp, mv2, dx, RDF(Ndist+1)
END TYPE NM
END MODULE DATASTR
