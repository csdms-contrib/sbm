cmake_minimum_required (VERSION 2.6)

include (FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
pkg_check_modules (GMODULE2 gmodule-2.0)
include_directories (${GLIB2_INCLUDE_DIRS} ${GMODULE2_INCLUDE_DIRS})
link_directories (${GLIB2_LIBRARY_DIRS} ${GMODULE2_LIBRARY_DIRS})

set (SBM_DATADIR ${CMAKE_INSTALL_PREFIX}/share/sbm)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/sbm_dirs.h.cmake
                ${CMAKE_CURRENT_SOURCE_DIR}/sbm_dirs.h)

find_library (LIBINTL intl HINTS ${CMAKE_LIBINTL_DIR})
if (NOT LIBINTL)
  message (STATUS
           "Looking for libintl... - NOT found.  Try setting CMAKE_LIBINTL_DIR")
endif (NOT LIBINTL)
get_filename_component (LIBINTL_DIR ${LIBINTL} PATH)
message (STATUS "Looking for libintl... - found ${LIBINTL_DIR}")
link_directories (${LIBINTL_DIR})

set (sbm_test_SRCS test_sbm.c)
add_executable (sbm_test ${sbm_test_SRCS})
target_link_libraries (sbm_test m sbm ${GLIB2_LIBRARIES} ${GMODULE2_LIBRARIES})
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbm_test
         DESTINATION bin COMPONENT sbm RENAME sbm-test)

