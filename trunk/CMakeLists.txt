cmake_minimum_required (VERSION 2.6)

add_subdirectory (data)
add_subdirectory (tests)

set (sbm_LIB_SRCS sbm.c)
add_library (sbm ${sbm_LIB_SRCS})
install (TARGETS sbm DESTINATION lib COMPONENT sbm)

set (sbm_SRCS sbm.c)
add_executable (sbm_exe ${sbm_SRCS})
target_link_libraries (sbm_exe m)
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbm_exe
         DESTINATION bin COMPONENT sbm RENAME sbm)

set (sbm_orig_SRCS SBcode_sortmar2005TWO.c)
add_executable (sbm_orig ${sbm_orig_SRCS})
target_link_libraries (sbm_orig m)
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sbm_orig
         DESTINATION bin COMPONENT sbm RENAME sbm-orig)

