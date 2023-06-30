cd ..
mkdir build_for_grid
rm build_for_grid/*
cd build_for_grid
../../../../../tools/genmake2 -mpi -mods ../code_for_grid -optfile ../../../../../tools/build_options/darwin_amd64_gfortran
make depend
make
cd ../utils

