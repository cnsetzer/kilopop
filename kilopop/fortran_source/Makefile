F90FLAGS = '-cpp -g -O3 -ffpe-trap=overflow,underflow,invalid -Wall -fbacktrace \
		   -fimplicit-none -fdefault-double-8 -fdefault-real-8 -fopenmp'

LAPACK_DIR = /Users/cnsetzer/software/lib/lapack-3.8.0/
#LAPACK_DIR = /lib64
FPATH = -L$(LAPACK_DIR) -I$(LAPACK_DIR)
FLIBS = -llapack -lblas

all:
	f2py -c $(FPATH) --f90flags=$(F90FLAGS) $(FLIBS) only: calculate_luminosity calculate_heating_rates : module_physics_constants.f90 hratelib.f90 macronova_Pinto_Eastman_CNS.f90 macronova2py.f90 -m m2py

debug:
	f2py -m mac2py -h macronova2py.pyf macronova2py.f90
	f2py -c --debug --build-dir build $(FPATH) --f90flags=$(F90FLAGS) $(FLIBS) only: calculate_luminosity calculate_heating_rates :  macronova2py.pyf module_physics_constants.f90 hratelib.f90 macronova_Pinto_Eastman_CNS.f90 macronova2py.f90

step1:
	f2py -m macronova2py -h macronova2py.pyf macronova2py.f90

step2:
	f2py -c --build-dir build $(FPATH) --f90flags=$(F90FLAGS) $(FLIBS) only: calculate_luminosity calculate_heating_rates :  macronova2py.pyf module_physics_constants.f90 hratelib.f90 macronova_Pinto_Eastman_CNS.f90 macronova2py.f90

clean:
	rm *.pyf
	rm *.so
	rm -rf build/
