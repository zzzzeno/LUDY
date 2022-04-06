#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mp++::mp++" for configuration "Release"
set_property(TARGET mp++::mp++ APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(mp++::mp++ PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmp++.so.12.0"
  IMPORTED_SONAME_RELEASE "libmp++.so.12"
  )

list(APPEND _IMPORT_CHECK_TARGETS mp++::mp++ )
list(APPEND _IMPORT_CHECK_FILES_FOR_mp++::mp++ "${_IMPORT_PREFIX}/lib/libmp++.so.12.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
