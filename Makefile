F90 = gfortran
FFLAGS1 = -O2 -fbounds-check -fbacktrace -pg
FFLAGS2 = -fbounds-check -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,invalid -g -pg
FFLAGS3 = -Wall -O3
FFLAGS4 = -O3

OBJ = obj/
OUTPUT = output/
WAKE = output/wake
WING = output/wing
CP = output/cp
CPFFT = output/cpfft
CEPH_1 = cephes_subset/ellf/
CEPH_2 = cephes_subset/eval/


vlatu:
#	$(F90) $(FFLAGS4) -c $(CEPH_1)ellie.c $(CEPH_1)ellpe.c $(CEPH_1)ellpk.c $(CEPH_1)mtherr.c $(CEPH_1)polevl.c $(CEPH_1)const.c -J $(OBJ)
	$(F90) $(FFLAGS4) -c modules3.f90 functions.f90 pg5.f90 gtmod10.f90 vlatus7.f90 vtxpart5.f90 4dgrid_10.f90 hyper_code.f90 pen_check.f90 influences_5.f90 vlatu_st_15_2.f90 -J $(OBJ)
	$(F90) -o vlatu *.o -I $(OBJ) -llapack -pg

vlatu_harmonic:
	$(F90) $(FFLAGS4) -c modules3.f90 functions_oscillation.f90 pg5.f90 gtmod10.f90 vlatus7.f90 vtxpart5.f90 4dgrid_10.f90 hyper_code.f90 pen_check.f90 influences_5.f90 vlatu_st_15_2.f90 -J $(OBJ)
	$(F90) -o vlatu *.o -I $(OBJ) -llapack -pg

vlatu_steady:
	$(F90) $(FFLAGS4) -c modules3.f90 functions_steady.f90 pg5.f90 gtmod10.f90 vlatus7.f90 vtxpart5.f90 4dgrid_10.f90 hyper_code.f90 pen_check.f90 influences_5.f90 vlatu_st_15_2.f90 -J $(OBJ)
	$(F90) -o vlatu *.o -I $(OBJ) -llapack -pg

test_traj:
	$(F90) $(FFLAGS4) modules3.f90 functions.f90 data_types.f90 test_traj.f90 -o test_traj

test_traj_steady:
	$(F90) $(FFLAGS4) modules3.f90 functions_steady.f90 data_types.f90 test_traj.f90 -o test_traj

test_traj_oscillation:
	$(F90) $(FFLAGS4) modules3.f90 functions_oscillation.f90 data_types.f90 test_traj.f90 -o test_traj

test_traj_interp:
	$(F90) $(FFLAGS4) modules3.f90 functions.f90 test_traj_interp.f90 -o test_traj_interp

clean:
	mv -f vlatu $(OBJ)*.mod *.o *.mod test_traj test_hyper trash
