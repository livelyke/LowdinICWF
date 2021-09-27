gfortran -g -O3 -L/lib/x86_64-linux-gnu/  -c bases.f90 modSlaterIndex.f90 params.f90 operators.f90 main.f90
gfortran -g -O3 bases.o main.o params.o modSlaterIndex.o operators.o -o main -lblas -llapack #-L/home/livelyke/bin/SPARSKIT2/ -lskit
