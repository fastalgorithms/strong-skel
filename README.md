# strong-skel
A research implementation of quadrature coupled RS-S from Sushnikova et al., "FMM-LU: A fast direct solver for multiscale boundary integral equations in three dimensions", [arXiv](https://arxiv.org/pdf/2201.07325.pdf). This work is based on the previous work of Minden et al., "
A recursive skeletonization factorization based on strong admissibility", [arXiv](https://arxiv.org/pdf/1609.08130.pdf). 


Original contributors:
Victor Minden, Ken L. Ho, Anil Damle, Lexing Ying.

Modifications by:
Daria Sushnikova, Mike O'Neil, Leslie Greengard, and Manas Rachh.

## Installation instructions
- Download the repository and run startup.m
- Download the [FLAM library](https://github.com/klho/FLAM/) and run startup.m in the repository to add the relevant files to path. 
- Download the [FMM3D library](https://github.com/flatironinstitute/FMM3D) and run `make install PREFIX=(FMM3DINSTALLDIR)` in the main directory
- Download the [fmm3dbie library](https://github.com/fastalgoritms/fmm3dbie) and run `make install LIBNAME=fmm3dbie_matlab PREFIX=(FMM3DBIEINSTALLDIR) PREFIX_FMM=(FMM3DINSTALLDIR) BLAS_64=ON` and copy over 
- Go to fortran_src in this repository, and update the compile_script.make file to update FMM3DINSTALLDIR and FMM3DBIEINSTALLDIR and run make -f compile_script.make
- This should generate a mex file. To test a successful installation, run test1.m in the fortran_src folder

## Note on installation
- For mac machines, we recommend setting PREFIX and PREFIX_FMM to blank or to `/usr/local/lib` and then FMM3DINSTALLDIR and FMM3DBIEINSTALLDIR to `/usr/local/lib` as well
- On linux machines, make sure that FMMINSTALLDIR and FMM3DINSTALLDIR are in the `LD_LIBRARY_PATH` environment variable

## Details
Built on the [FLAM library](https://github.com/klho/FLAM/) by Ken L. Ho, this research code implements the factorizations described in "A recursive skeletonization factorization based on strong admissibility" by Minden et al.


The main functions provided are as follows:
- `srskelf.m`: the strong recursive skeletonization factorization (aka "RS-S") for symmetric matrices.  This function uses Cholesky factorizations of diagonal subblocks that can be assumed positive-definite in exact arithmetic.
- `srskelf_asym.m`: RS-S for asymmetric matrices.  This function uses LU factorizations where `srskelf.m` uses Cholesky.
- `srskelf_hybrid.m`: the hybrid recursive skeletonization factorization (aka "RS-WS") for symmetric matrices.  This function is like `srskelf.m`, but additionally interlaces levels of standard ("weak") skeletonization.
- `srskelf_asymhybrid.m`: RS-WS for asymmetric matrices.  This function is like `srskelf_asym.m`, but additionally interlaces levels of standard ("weak") skeletonization.
- `srskelf_asym_new.m`: RS-S for asymmetric matrices, where the proxy points per level are chosen based on the size of the box in wavelengths

Tests are available as MATLAB functions in the directory "test" and may be run with no arguments, assuming that FLAM and the strong-skel directories are in your MATLAB path.  The tests provided are as follows:
- `ie_square.m`: The setup of Example 1 as described in Minden et al.
- `ie_cube.m`: The setup of Example 2 as described in Minden et al.
- `ie_sphere.m`: The setup of Example 3 as described in Minden et al.
-  `ie_fmm3dbie_sqrtscaling.m`: Setup for solving Helmholtz equation with combined field representation and square root scaled matrix entries, where the near quadrature is generated using fmm3dbie. 
