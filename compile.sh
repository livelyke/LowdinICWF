gfortran -g -O3 -Wunused -L/lib/x86_64-linux-gnu/ -c tdEvolve.f90 bases.f90 modSlaterIndex.f90 params.f90 operators.f90 main.f90
gfortran -g -O3 tdEvolve.o bases.o main.o params.o modSlaterIndex.o operators.o -o main -lblas -llapack #-L/home/livelyke/bin/SPARSKIT2/ -lskit
