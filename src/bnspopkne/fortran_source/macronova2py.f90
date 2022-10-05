MODULE macronova2py

USE macronova_Pinto_eastman_CNS
USE physics_constants


CONTAINS

SUBROUTINE calculate_luminosity(n, MNE_parameters, func_hrate, func_therm, read_hrate, heating_rates_file, Nt, luminosity)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n, Nt
    !f2py Integer, intent(in):: n, Nt
    DOUBLE PRECISION, INTENT(IN)  :: MNE_parameters(n)
    !f2py DOUBLE PRECISION, intent(in), depend(n) :: MNE_parameters
    DOUBLE PRECISION, INTENT(OUT) :: luminosity(Nt+1,4)
    !f2py DOUBLE PRECISION, intent(out), depend(Nt) :: luminosity
    LOGICAL, INTENT(IN) :: read_hrate, func_hrate, func_therm
    !f2py intent(in) :: read_hrate, func_hrate
    CHARACTER*255, INTENT(IN) :: heating_rates_file
    !f2py intent(in) :: heating_rates_file

    CALL macronova(n, MNE_parameters, func_hrate, func_therm, read_hrate, heating_rates_file, Nt, luminosity)

END SUBROUTINE calculate_luminosity

END MODULE
