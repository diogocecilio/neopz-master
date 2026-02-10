#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "NeoPZ::pz" for configuration "Release"
set_property(TARGET NeoPZ::pz APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(NeoPZ::pz PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/pz/lib/libpz.so"
  IMPORTED_SONAME_RELEASE "libpz.so"
  )

list(APPEND _cmake_import_check_targets NeoPZ::pz )
list(APPEND _cmake_import_check_files_for_NeoPZ::pz "${_IMPORT_PREFIX}/pz/lib/libpz.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
