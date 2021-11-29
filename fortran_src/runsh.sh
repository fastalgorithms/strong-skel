rm -rf helmquadcorr.o
gfortran -O3 -march=native -fPIC -fopenmp -c helmquadcorr.f -o helmquadcorr.o -L/usr/local/lib -lfmm3d -lfmm3dbie
../../mwrap/mwrap -c99complex -list -mex fmm3dbierouts -mb fmm3dbierouts.mw
../../mwrap/mwrap -c99complex -mex -fmm3dbierouts -c fmm3dbierouts.c fmm3dbierouts.mw
/Applications/MATLAB_R2021a.app/bin/mex -v fmm3dbierouts.c helmquadcorr.o -largeArrayDims -DMWF77_UNDERSCORE1 -D_OPENMP -L/usr/local/lib/gcc/11 -output fmm3dbierouts -L/usr/local/lib -lfmm3d -lfmm3dbie -lgomp -lstdc++ -lm -ldl -lgfortran
