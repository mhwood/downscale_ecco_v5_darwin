cd ../
# mkdir run_for_grid
rm run_for_grid/*
rm -r run_for_grid/mnc*
cd run_for_grid
ln -s ../namelist_for_grid/* .
ln -s ../input/* .
ln -s ../build_for_grid/mitgcmuv .
mpirun --oversubscribe -np 12 ./mitgcmuv
