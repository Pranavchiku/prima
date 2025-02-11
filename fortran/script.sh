set -ex

# take FC from the environment
echo "FC is $FC"

$FC --cpp -c ./common/consts.F90 
$FC --cpp -c ./common/infos.f90 
$FC --cpp -c ./common/debug.F90 
$FC --cpp -c ./common/huge.F90 
$FC --cpp -c ./common/inf.F90 
$FC --cpp -c ./common/infnan.F90 
$FC --cpp -c ./common/memory.F90 
$FC --cpp -c ./common/string.f90 
$FC --cpp -c ./common/linalg.f90 
$FC --cpp -c ./common/univar.f90 
$FC --cpp -c ./common/ratio.f90 
$FC --cpp -c ./common/redrho.f90 
$FC --cpp -c ./common/xinbd.f90 
$FC --cpp -c ./common/history.f90 
$FC --cpp -c ./common/selectx.f90 
$FC --cpp -c ./common/checkexit.f90 
$FC --cpp -c ./common/fprint.f90 
$FC --cpp -c ./common/message.f90 
$FC --cpp -c ./common/preproc.f90 
$FC --cpp -c ./common/pintrf.f90 
$FC --cpp -c ./common/evaluate.f90 
$FC --cpp -c ./common/powalg.f90 
$FC --cpp -c ./common/shiftbase.f90
$FC --cpp -c ./tests/testsuite/param.f90 
$FC --cpp -c ./tests/testsuite/rand.f90 
$FC --cpp -c ./tests/testsuite/noise.f90 
$FC --cpp -c ./tests/testsuite/prob.f90 
$FC --cpp -c ./tests/testsuite/datetime.f90

$FC --cpp -c ./tests/test_bobyqa.f90 -I../build/fortran -I./
$FC --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./

$FC --cpp -c ./tests/test_cobyla.f90 -I../build/fortran -I./
$FC --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./

$FC --cpp -c ./tests/test_lincoa.f90 -I../build/fortran -I./
$FC --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./

$FC --cpp -c ./tests/test_newuoa.f90 -I../build/fortran -I./
$FC --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./

$FC --cpp -c ./tests/test_uobyqa.f90 -I../build/fortran -I./
$FC --cpp -DPRIMA_DEBUGGING=1 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=16 -DPRIMA_REAL_PRECISION=128 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./

