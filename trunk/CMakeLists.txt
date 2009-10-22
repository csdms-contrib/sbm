cmake_minimum_required (VERSION 2.6)

set (sbm_orig_SRCS SBcode_sortmar2005TWO.c)
add_executable (sbm_orig ${sbm_orig_SRCS})
target_link_libraries (sbm_orig m)
install (TARGETS sbm_orig DESTINATION bin COMPONENT sbm)

