set(tcmalloc_NAMES
    tcmalloc
    libtcmalloc
    libtcmalloc.so
    tcmalloc4
    libtcmalloc.so.4)

find_library(
  tcmalloc_LIBRARY NO_DEFAULT_PATH
  NAMES ${tcmalloc_NAMES}
  PATHS /lib
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
        /opt/local/lib)

if(tcmalloc_LIBRARY)
  set(tcmalloc_FOUND TRUE)
else()
  set(tcmalloc_FOUND FALSE)
endif()

if(tcmalloc_FOUND)
  message(STATUS "Found tcmalloc: ${tcmalloc_LIBRARY}")
else()
  message(STATUS "Not Found tcmalloc: ${tcmalloc_LIBRARY} with names ${tcmalloc_NAMES}")
endif()

mark_as_advanced(tcmalloc_LIBRARY)

# Output variables generation
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(tcmalloc DEFAULT_MSG tcmalloc_LIBRARY)

if(tcmalloc_FOUND)
  find_package(Threads REQUIRED)
  add_library(tcmalloc::tcmalloc INTERFACE IMPORTED)
  set_property(TARGET tcmalloc::tcmalloc PROPERTY INTERFACE_LINK_LIBRARIES ${tcmalloc_LIBRARY} Threads::Threads)
endif()





add_executable(${PROJECT_NAME}_benchmarks benchmarks.cpp)
target_link_libraries(
  ${PROJECT_NAME}_benchmarks
  PRIVATE ${PROJECT_NAME}
          ${PROJECT_NAME}_bgl
          ${PROJECT_NAME}_test_utils
          tcmalloc::tcmalloc
          GTest::GTest
          GTest::Main)
target_clang_tidy(
  ${PROJECT_NAME}_benchmarks
  ENABLE ${DESCARTES_ENABLE_CLANG_TIDY}
  WARNINGS_AS_ERRORS ${DESCARTES_ENABLE_TESTING}
  CHECKS ${DEFAULT_CLANG_TIDY_CHECKS})

# set(CMAKE_EXE_LINKER_FLAGS -Wl,-no-as-needed -lprofiler -Wl,--as-needed)

install(TARGETS ${PROJECT_NAME}_benchmarks RUNTIME DESTINATION lib/${PROJECT_NAME})
