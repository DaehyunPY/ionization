### CMake project file for myprog
cmake_minimum_required (VERSION 2.6.4)
project (myprog)
enable_language (Fortran)





add_executable (a.out main.f90 global.f90 hamiltonian.f90 basis.f90 boundary.f90 mod_zbes.f90)

# library & link directories
include_directories (/usr/include)
include_directories (/usr/local/include)
include_directories ($ENV{HOME}/include)
link_directories (/usr/lib)
link_directories (/usr/local/lib)
link_directories ($ENV{HOME}/lib)

### fgsl library
target_link_libraries (a.out -lfgsl)

### lapack library
# target_link_libraries (a.out -llapack)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	include_directories (-i8)
	include_directories (/opt/intel/mkl/include)
	include_directories (/opt/intel/mkl/include/intel64/ilp64)
	target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_blas95_ilp64.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_lapack95_ilp64.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_intel_ilp64.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_core.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_sequential.a)
	target_link_libraries (a.out -lpthread)
	target_link_libraries (a.out -lm)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	include_directories (-i8)
	include_directories (/opt/intel/mkl/include)
	include_directories (/opt/intel/mkl/include/intel64/ilp64)
	target_link_libraries (a.out /opt/intel/mkl/lib/intel64/libmkl_blas95_ilp64.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/intel64/libmkl_lapack95_ilp64.a)
	target_link_libraries (a.out -WL,--start-group)
	target_link_libraries (a.out /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/intel64/libmkl_core.a)
	target_link_libraries (a.out /opt/intel/mkl/lib/intel64/libmkl_sequential.a)
	target_link_libraries (a.out -Wl,--end-group)
	target_link_libraries (a.out -lpthread)
	target_link_libraries (a.out -lm)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

### my library
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	include_directories (/Users/hyuna917/Dropbox/Library/Fortran/include/)
	target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libconst.a)
	# target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libfft.a)
	target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/liblinear.a)
	target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libplot.a)
	# target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libsde.a)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	include_directories (/users/daehyun/include/)
	target_link_libraries (a.out /users/daehyun/lib/libconst.a)
	# target_link_libraries (a.out /users/daehyun/lib/libfft.a)
	target_link_libraries (a.out /users/daehyun/lib/liblinear.a)
	# target_link_libraries (a.out /users/daehyun/lib/libplot.a)
	# target_link_libraries (a.out /users/daehyun/lib/libsde.a)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")





### make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
			"Choose the type of build, options are: None Debug Release."
			FORCE)
endif (NOT CMAKE_BUILD_TYPE)

### FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	### gfortran
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
	# set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -fimplicit-none -fbounds-check -O -Wuninitialized -ffpe-trap=invalid, zero, overflow -fbacktrace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
	### ifort (untested)
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
	# set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -warn all -check bounds -check uninit -fpe0 -traceback -heap-arrays")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
	### g77
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	message ("No optimized Fortran compiler flags are known, we just try -O2...")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
