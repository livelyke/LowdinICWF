#gfortran -g -O3 -Wunused   -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -c tdEvolve.f90 bases.f90 modSlaterIndex.f90 params.f90 operators.f90 main.f90
#gfortran -g -O3  -m64  tdEvolve.o bases.o main.o params.o modSlaterIndex.o operators.o -o main -I"${MKLROOT}/include" -lblas -llapack #-L/home/livelyke/bin/SPARSKIT2/ -lskit
rm *o 
rm *mod
gfortran -c -g -m64 -O3  modSlaterIndex.f90 -o modSlaterIndex.o
gfortran -c -g -m64 -O3  params.f90 -o params.o
gfortran -c -g -m64 -O3  -I"${MKLROOT}/include" bases.f90  -o bases.o
gfortran -c -g -m64 -O3  -I"${MKLROOT}/include" operators.f90  -o operators.o
gfortran -c -g -m64 -O3  -I"${MKLROOT}/include" tdEvolve.f90   -o tdEvolve.o
gfortran -c -g -m64 -O3  -I"${MKLROOT}/include" main.f90   -o main.o
gfortran -g -m64   main.o tdEvolve.o bases.o params.o modSlaterIndex.o operators.o  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -o main #-L${MKLROOT}/lib/intel64  -lblas -llapack

#gfortran -c -g -O3  modSlaterIndex.f90 -o modSlaterIndex.o
#gfortran -c -g -O3  params.f90 -o params.o
#gfortran -c -g -O3   bases.f90  -o bases.o
#gfortran -c -g -O3   operators.f90  -o operators.o
#gfortran -c -g -O3   tdEvolve.f90   -o tdEvolve.o
#gfortran -c -g -O3   main.f90   -o main.o
#gfortran -g  main.o tdEvolve.o bases.o params.o modSlaterIndex.o operators.o /u/klively/bin/OpenBLAS/libopenblas.a  -o main -lpthread #-L${MKLROOT}/lib/intel64  -lblas -llapack


