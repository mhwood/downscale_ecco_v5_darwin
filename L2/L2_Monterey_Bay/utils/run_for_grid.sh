cd ../
rm -r run_for_grid/mnc_*
rm run_for_grid/*
cd run_for_grid
ln -s ../input/* .
ln -s ../namelist_for_grid/* .
ln -s ../build_for_grid/mitgcmuv .
mpirun -np 9 ./mitgcmuv
