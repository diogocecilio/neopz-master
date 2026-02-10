#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "NeoPZ::pz" for configuration "Debug"
set_property(TARGET NeoPZ::pz APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(NeoPZ::pz PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/pz/lib/libpz.so"
  IMPORTED_SONAME_DEBUG "libpz.so"
  )

list(APPEND _cmake_import_check_targets NeoPZ::pz )
list(APPEND _cmake_import_check_files_for_NeoPZ::pz "${_IMPORT_PREFIX}/pz/lib/libpz.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
