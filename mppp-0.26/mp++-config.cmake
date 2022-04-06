# Get current dir.
get_filename_component(_MPPP_CONFIG_SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Find the deps. Alter the cmake module path.
set(_MPPP_CONFIG_OLD_MODULE_PATH "${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH "${_MPPP_CONFIG_SELF_DIR}")

# Mandatory dep on GMP.
find_package(mp++_GMP REQUIRED)

# Public optional deps.
if(y)
    find_package(mp++_MPFR REQUIRED)
endif()
set(mp++_WITH_MPFR y)

if(OFF)
    find_package(mp++_MPC REQUIRED)
endif()
set(mp++_WITH_MPC OFF)

if(OFF)
    find_package(Boost 1.60 REQUIRED COMPONENTS serialization)
endif()
set(mp++_WITH_BOOST_S11N OFF)

# Private optional deps.
set(mp++_WITH_ARB OFF)
set(mp++_WITH_QUADMATH y)

# Restore original module path.
set(CMAKE_MODULE_PATH "${_MPPP_CONFIG_OLD_MODULE_PATH}")
unset(_MPPP_CONFIG_OLD_MODULE_PATH)

include(${_MPPP_CONFIG_SELF_DIR}/mp++_export.cmake)

# Clean up.
unset(_MPPP_CONFIG_SELF_DIR)
