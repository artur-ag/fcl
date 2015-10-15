INCLUDE( "FindPackageHandleStandardArgs" )


FIND_PATH( LAPACKE_INCLUDE_DIRS "lapacke.h" )
#PATH_SUFFIXES "lapacke" )

FIND_LIBRARY( LAPACKE_LIBRARY_DIRS
NAMES "liblapacke.so" )
#PATH_SUFFIXES "liblapacke" )

FIND_PACKAGE_HANDLE_STANDARD_ARGS( "LAPACKE" DEFAULT_MSG LAPACKE_INCLUDE_DIRS LAPACKE_LIBRARY_DIRS )
