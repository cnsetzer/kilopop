MODULE physics_constants

  IMPLICIT NONE

  DOUBLE PRECISION,PARAMETER:: hplanck=  6.62607015d-27 ! Planck constant [erg*s]
  DOUBLE PRECISION,PARAMETER:: kB=       1.38064852d-16 ! Boltzmann constant [erg/K]
  DOUBLE PRECISION,PARAMETER:: parsec=   3.08567758d+18 ! parsec [cm]
  DOUBLE PRECISION,PARAMETER:: clight=   2.99792458d+10 ! speed of light [cm/s]
  DOUBLE PRECISION,PARAMETER:: sigma=    5.67051d-5     ! Stefan-Boltzmann const. [erg/(s*cm2*K4)]
  DOUBLE PRECISION,PARAMETER:: msol=     1.9891d+33     ! solar mass [g]
  DOUBLE PRECISION,PARAMETER:: day_in_s= 8.64d+4        ! one day [s]
  DOUBLE PRECISION,PARAMETER:: Robs=     10.0d0*parsec    ! distance for abs. mags (10 pc,[cm])
  DOUBLE PRECISION,PARAMETER:: Ang=      1.0d-8         ! angstrom [cm]
  DOUBLE PRECISION,PARAMETER:: Pi=       ASIN(1.0d0)*2.0d0    ! pi (3.14159265358979..)

END MODULE physics_constants
