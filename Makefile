FC=pgf90
FCFLAGS=-O3
HOMEPATH=$(PWD)




msd_rot_calc: src/fortran/msd_rot_calc.f90 src/fortran/visc_calc.f90 src/fortran/flux_side.f90
	test -f module/path.include || echo "module use $(HOMEPATH)/module" >> ~/.bash_profile
	test -f module/path.include || source ~/.bash_profile
	@echo "prepend_path('PATH', '$(HOMEPATH)/src/fortran/exec')" > module/path.include
	@echo "prepend_path('PATH', '$(HOMEPATH)')" >> module/path.include
	cat module/direct_calc_header module/path.include > module/Dir_Calc_Fluct.lua
	$(FC) $(FCFLAGS) -o src/exec/msd_rot_calc src/fortran/msd_rot_calc.f90
	$(FC) $(FCFLAGS) -o src/exec/visc_calc    src/fortran/visc_calc.f90
	$(FC) $(FCFLAGS) -o src/exec/flux_side    src/fortran/flux_side.f90
	chmod 777 backbone.py



