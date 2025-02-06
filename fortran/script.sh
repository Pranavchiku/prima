set -ex

lfortran --cpp -c ./common/consts.F90 
lfortran --cpp -c ./common/infos.f90 
lfortran --cpp -c ./common/debug.F90 
lfortran --cpp -c ./common/huge.F90 
lfortran --cpp -c ./common/inf.F90 
lfortran --cpp -c ./common/infnan.F90 
lfortran --cpp -c ./common/memory.F90 
lfortran --cpp -c ./common/string.f90 
lfortran --cpp -c ./common/linalg.f90 
lfortran --cpp -c ./common/univar.f90 
lfortran --cpp -c ./common/ratio.f90 
lfortran --cpp -c ./common/redrho.f90 
lfortran --cpp -c ./common/xinbd.f90 
lfortran --cpp -c ./common/history.f90 
lfortran --cpp -c ./common/selectx.f90 
lfortran --cpp -c ./common/checkexit.f90 
lfortran --cpp -c ./common/fprint.f90 
lfortran --cpp -c ./common/message.f90 
lfortran --cpp -c ./common/preproc.f90 
lfortran --cpp -c ./common/pintrf.f90 
lfortran --cpp -c ./common/evaluate.f90 
lfortran --cpp -c ./common/powalg.f90 
lfortran --cpp -c ./common/shiftbase.f90 
lfortran --cpp -c ./cobyla/update.f90 
lfortran --cpp -c ./cobyla/initialize.f90 
lfortran --cpp -c ./cobyla/trustregion.f90 
lfortran --cpp -c ./cobyla/geometry.f90 
lfortran --cpp -c ./cobyla/cobylb.f90 
lfortran --cpp -c ./cobyla/cobyla.f90 
lfortran --cpp -c ./tests/testsuite/param.f90 
lfortran --cpp -c ./tests/testsuite/rand.f90 
lfortran --cpp -c ./tests/testsuite/noise.f90 
lfortran --cpp -c ./tests/testsuite/prob.f90 
lfortran --cpp -c ./tests/testsuite/datetime.f90 
lfortran --cpp -c ./tests/test_lincoa.f90 -I../build/fortran -I./

lfortran --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./
