set -ex

# take FC from the environment
echo "FC is $FC"

rm -rf *.o

$FC --cpp -c ./tests/testsuite/param.f90 -I../build/fortran -I./
$FC --cpp -c ./tests/testsuite/rand.f90 -I../build/fortran -I./
$FC --cpp -c ./tests/testsuite/noise.f90 -I../build/fortran -I./
$FC --cpp -c ./tests/testsuite/prob.f90 -I../build/fortran -I./
$FC --cpp -c ./tests/testsuite/datetime.f90 -I../build/fortran -I./

#!/bin/bash

# Source directory (can be changed as needed)
SRC_DIR="../build/fortran/CMakeFiles/primaf.dir/$name"

# Find and copy all *.o files recursively
find "$SRC_DIR" -type f -name "*.o" -exec cp {} . \;

SRC_DIR="../build/fortran/CMakeFiles/primaf.dir/common"

# Find and copy all *.o files recursively
find "$SRC_DIR" -type f -name "*.o" -exec cp {} . \;

echo "All .o files copied to current directory."

$FC --cpp -c ./tests/$test_name -I../build/fortran -I./
$FC --cpp -c -DPRIMA_DEBUGGING=0 -DPRIMA_AGGRESSIVE_OPTIONS=0 -DPRIMA_INTEGER_KIND=4 -DPRIMA_REAL_PRECISION=32 -DPRIMA_QP_AVAILABLE=0 -DPRIMA_TESTDIM="'small'" ./tests/test.F90 -I../build/fortran -I./
$FC *.o -I../build/fortran -I./
