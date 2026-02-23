#==============================================================
# Makefile for 1D viscous disk evolution code
#==============================================================

#------------------------------
# Compiler and flags
#------------------------------
FC      = gfortran
# Debug: -fcheck=all adds bounds checking (slow, for development only)
#FFLAGS  = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
#         -finit-real=snan -finit-integer=-999999 -llapack -lblas
# Release: no bounds checking, max optimization, OpenMP, native arch
FFLAGS  = -O3 -fopenmp -march=native -llapack -lblas
LDFLAGS = $(FFLAGS)

#------------------------------
# Source files
#------------------------------
SRC = kind_params.f90 \
      constants.f90 \
      mod_global.f90 \
      input_file_mod.f90 \
      units_disk_mod.f90 \
      inflow_source.f90 \
      star_params_mod.f90 \
      mdot_units_mod.f90 \
      bessel_i.f90 \
      initial_mass.f90 \
      analytic_lbp.f90 \
      compare_analytic.f90 \
      compute_LSigma.f90 \
      build_tridiag_coeff.f90 \
      solve_tridiag.f90 \
      diffusion_theta_step.f90 \
      op_params.f90 \
      checkpoint_util_mod.f90 \
      run_control.f90 \
      spline2D_pac.f90 \
      polint2D_pac.f90 \
      opacity_table_mod.f90 \
      wind_mod.f90 \
      relax_mod.f90 \
      lapack_linear_solvers_mod.f90 \
      disk_flux_mod.f90 \
      radiation_params_mod.f90 \
      irradiation_mod.f90 \
      timestep_seed_mod.f90 \
      disk_thermal_mod.f90 \
      thermal_stability_analysis_mod.f90 \
      disk_energy_mod.f90 \
      disk_energy_pde_mod.f90 \
      hot_region_metrics_mod.f90 \
      run_control_mod.f90 \
      evolve_try_mod.f90 \
      setup.f90 \
      state_io_mod.f90 \
      disk_state_update_mod.f90 \
      output_mod.f90 \
      checkpoint_mod.f90 \
      output_run_summary_mod.f90 \
      main.f90

OBJ = $(SRC:.f90=.o)
EXE = vd1d_stage5

#------------------------------
# Rules
#------------------------------
$(EXE): $(OBJ)
	$(FC) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#------------------------------
# S-curve analysis program
#------------------------------
SCURVE_OBJ = $(filter-out main.o, $(OBJ)) scurve_main.o
scurve: $(SCURVE_OBJ)
	$(FC) -o scurve $(SCURVE_OBJ) $(LDFLAGS)

scurve_main.o: scurve_main.f90
	$(FC) $(FFLAGS) -c $<

#------------------------------
# One-cell thermal stability (root find) program
#------------------------------
ONECELL_OBJ = $(filter-out main.o, $(OBJ)) onecell_rootfind_from_diskt.o
onecell_rootfind: $(ONECELL_OBJ)
	$(FC) -o onecell_rootfind $(ONECELL_OBJ) $(LDFLAGS)

onecell_rootfind_from_diskt.o: onecell_rootfind_from_diskt.f90
	$(FC) $(FFLAGS) -c $<

#------------------------------
# Planck opacity test
#------------------------------
TEST_PLANCK_OBJ = $(filter-out main.o, $(OBJ)) test_planck_opacity.o
test_planck: $(TEST_PLANCK_OBJ)
	$(FC) -o test_planck $(TEST_PLANCK_OBJ) $(LDFLAGS)

test_planck_opacity.o: test_planck_opacity.f90
	$(FC) $(FFLAGS) -c $<

#------------------------------
# Debug build (bounds checking, slower)
#------------------------------
FFLAGS_DEBUG = -O2 -g -fopenmp -fcheck=all -fbacktrace -llapack -lblas
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: LDFLAGS = $(FFLAGS_DEBUG)
debug: clean $(EXE)
	@echo "Built with debug flags (-fcheck=all)"

#------------------------------
# Utility targets
#------------------------------
.PHONY: clean show debug test_planck scurve onecell_rootfind

clean:
	rm -f $(OBJ) $(EXE) scurve scurve_main.o onecell_rootfind onecell_rootfind_from_diskt.o test_planck test_planck_opacity.o *.mod

show:
	@echo "Compiler  : $(FC)"
	@echo "Flags     : $(FFLAGS)"
	@echo "Link libs : $(LDFLAGS)"
	@echo "Sources   : $(SRC)"
