# Try to find TCLAP
# TCLAP_FOUND
# TCLAP_INCLUDE_DIRS

find_package(PkgConfig)
pkg_check_modules(PC_TCLAP QUIET tclap)

find_path(TCLAP_INCLUDE_DIR tclap/CmdLine.h
	  HINTS ${PC_TCLAP_INCLUDEDIR} ${PC_TCLAP_INCLUDE_DIRS}
	  PATH_SUFFIXES tclap )

set ( TCLAP_INCLUDE_DIRS ${TCLAP_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TCLAP  DEFAULT_MSG
                                  TCLAP_INCLUDE_DIR)

mark_as_advanced(TCLAP_INCLUDE_DIR  )
