cmake_minimum_required(VERSION 3.20)
project(netevolve VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- SUNDIALS (CVODE, the default solver) ----------------------------------
# Resolved before the strict warning flags below, so that SUNDIALS' own
# sources (when downloaded) are not compiled with /W4 /WX.
#
#   1. If SUNDIALS_DIR points at a SUNDIALS 5.x install, use it:
#        -DSUNDIALS_DIR=<install>/lib/cmake/sundials
#   2. Otherwise (NETEVOLVE_FETCH_SUNDIALS=ON, the default) download the
#      official SUNDIALS 5.7.0 release, verify its SHA-256, and build just
#      CVODE as static libraries inside the build directory (once per build
#      directory; needs internet access the first time).
#   3. -DNETEVOLVE_WITH_CVODE=OFF builds without SUNDIALS; the default solver
#      is then the in-house `auto` (RK45 + Rosenbrock), which gives different
#      results from CVODE.
option(NETEVOLVE_WITH_CVODE "Build the CVODE solver (SUNDIALS 5.x); it is the default solver" ON)
option(NETEVOLVE_FETCH_SUNDIALS "Download and build SUNDIALS 5.7.0 if no install is found" ON)
if(NETEVOLVE_WITH_CVODE)
  find_package(SUNDIALS 5 CONFIG QUIET)
  if(SUNDIALS_FOUND)
    set(NETEVOLVE_SUNDIALS_SOURCE "installed SUNDIALS ${SUNDIALS_VERSION} from ${SUNDIALS_DIR}")
  elseif(NETEVOLVE_FETCH_SUNDIALS)
    include(FetchContent)
    if(POLICY CMP0135)
      cmake_policy(SET CMP0135 NEW)  # extracted file timestamps = extraction time
    endif()
    # Build only what netevolve needs: CVODE, static, no examples.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
    set(BUILD_CVODE ON CACHE BOOL "" FORCE)
    foreach(pkg ARKODE CVODES IDA IDAS KINSOL)
      set(BUILD_${pkg} OFF CACHE BOOL "" FORCE)
    endforeach()
    set(EXAMPLES_ENABLE_C OFF CACHE BOOL "" FORCE)
    set(EXAMPLES_ENABLE_CXX OFF CACHE BOOL "" FORCE)
    set(EXAMPLES_INSTALL OFF CACHE BOOL "" FORCE)
    set(SUNDIALS_INDEX_SIZE 64 CACHE STRING "" FORCE)
    set(_netevolve_fetch_opts)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
      list(APPEND _netevolve_fetch_opts EXCLUDE_FROM_ALL)  # build only the libraries we link
    endif()
    FetchContent_Declare(sundials
      URL https://github.com/LLNL/sundials/releases/download/v5.7.0/sundials-5.7.0.tar.gz
      URL_HASH SHA256=48da7baa8152ddb22aed1b02d82d1dbb4fbfea22acf67634011aa0303a100a43
      # Short directories: SUNDIALS' nested target directories plus MSBuild's
      # tracking logs can otherwise exceed the 260-character Windows path limit.
      SOURCE_DIR ${CMAKE_BINARY_DIR}/sd-src
      BINARY_DIR ${CMAKE_BINARY_DIR}/sd
      ${_netevolve_fetch_opts})
    if(FETCHCONTENT_SOURCE_DIR_SUNDIALS)
      message(STATUS "Using local SUNDIALS source ${FETCHCONTENT_SOURCE_DIR_SUNDIALS} (no download)")
    else()
      message(STATUS "Fetching SUNDIALS 5.7.0 (first configure of this build directory only)")
    endif()
    FetchContent_MakeAvailable(sundials)
    set(NETEVOLVE_SUNDIALS_SOURCE "SUNDIALS 5.7.0 built from ${sundials_SOURCE_DIR} in ${sundials_BINARY_DIR}")
  else()
    message(FATAL_ERROR
      "SUNDIALS 5.x not found and NETEVOLVE_FETCH_SUNDIALS=OFF. CVODE is netevolve's default solver.\n"
      "Either point CMake at a SUNDIALS 5.x install:\n"
      "  -DSUNDIALS_DIR=<install>/lib/cmake/sundials\n"
      "or allow the download: -DNETEVOLVE_FETCH_SUNDIALS=ON\n"
      "or build without CVODE (default solver becomes the in-house `auto`):\n"
      "  -DNETEVOLVE_WITH_CVODE=OFF")
  endif()
endif()

if(MSVC)
  # Warnings are errors from M2 onwards (spec §0.2).
  # Release already gets /O2 from CMake's default MSVC flags.
  add_compile_options(/W4 /WX /permissive- /EHsc)
endif()

add_compile_definitions(NETEVOLVE_VERSION="${PROJECT_VERSION}")

# --- Core library (everything except main) -----------------------------
add_library(netevolve_core STATIC
  src/network.cpp
  src/model.cpp
  src/integrator.cpp
  src/config.cpp
  src/fitness.cpp
  src/evolve.cpp
)
target_include_directories(netevolve_core PUBLIC src)

if(NETEVOLVE_WITH_CVODE)
  target_link_libraries(netevolve_core PUBLIC
    SUNDIALS::cvode_static SUNDIALS::nvecserial_static
    SUNDIALS::sunmatrixdense_static SUNDIALS::sunlinsoldense_static)
  target_compile_definitions(netevolve_core PUBLIC NETEVOLVE_HAVE_CVODE)
  message(STATUS "CVODE solver enabled: ${NETEVOLVE_SUNDIALS_SOURCE}")
else()
  message(STATUS "CVODE solver disabled (NETEVOLVE_WITH_CVODE=OFF): default solver is the in-house `auto`")
endif()

# --- CLI ---------------------------------------------------------------
add_executable(netevolve src/main.cpp)
target_link_libraries(netevolve PRIVATE netevolve_core)

# --- Tests -------------------------------------------------------------
enable_testing()

add_executable(netevolve_tests
  tests/test_main.cpp
  tests/test_smoke.cpp
  tests/test_network.cpp
  tests/test_model.cpp
  tests/test_integrator.cpp
  tests/test_config.cpp
  tests/test_fitness.cpp
  tests/test_evolve.cpp
  tests/test_prune.cpp
  tests/test_cvode.cpp
  tests/test_timeseries.cpp
)
target_include_directories(netevolve_tests PRIVATE third_party)
target_link_libraries(netevolve_tests PRIVATE netevolve_core)
target_compile_definitions(netevolve_tests PRIVATE
  NETEVOLVE_EXAMPLES_DIR="${CMAKE_SOURCE_DIR}/examples")

add_test(NAME unit_tests COMMAND netevolve_tests)

# `netevolve --version` must exit 0 ...
add_test(NAME cli_version_exit COMMAND netevolve --version)
# ... and print a version string.
add_test(NAME cli_version_output COMMAND netevolve --version)
set_tests_properties(cli_version_output PROPERTIES
  PASS_REGULAR_EXPRESSION "^netevolve [0-9]+\\.[0-9]+\\.[0-9]+")

# --- CLI commands (M4) ---------------------------------------------------
set(EX "${CMAKE_SOURCE_DIR}/examples")
add_test(NAME cli_simulate COMMAND netevolve simulate --model ${EX}/brusselator.ant --t-end 10 --points 11)
set_tests_properties(cli_simulate PROPERTIES PASS_REGULAR_EXPRESSION "^time,S0,S1\n0,1,1\n")
add_test(NAME cli_fitness_brusselator COMMAND netevolve fitness --model ${EX}/brusselator.ant --config ${EX}/default.cfg)
set_tests_properties(cli_fitness_brusselator PROPERTIES PASS_REGULAR_EXPRESSION "fitness (0\\.9[5-9][0-9]*|1)\n")
add_test(NAME cli_unknown_option COMMAND netevolve simulate --model ${EX}/brusselator.ant --bogus 1)
set_tests_properties(cli_unknown_option PROPERTIES WILL_FAIL TRUE)
add_test(NAME cli_bad_model COMMAND netevolve fitness --model ${EX}/default.cfg)
set_tests_properties(cli_bad_model PROPERTIES WILL_FAIL TRUE)

# --- evolve command (M5): tiny run end to end ------------------------------
file(WRITE "${CMAKE_BINARY_DIR}/cli_small.cfg" "population = 10\ngenerations = 2\nseed = 3\n")
add_test(NAME cli_evolve COMMAND netevolve evolve --config ${CMAKE_BINARY_DIR}/cli_small.cfg
                                                  --out ${CMAKE_BINARY_DIR}/cli_evolve_out)
set_tests_properties(cli_evolve PROPERTIES PASS_REGULAR_EXPRESSION "best fitness [0-9.e-]+ .generation [0-9]+")

# --- evolve --numtrials: several independent trials in one run ----------------
add_test(NAME cli_evolve_numtrials COMMAND netevolve evolve --config ${CMAKE_BINARY_DIR}/cli_small.cfg
                                                            --numtrials 2 --out ${CMAKE_BINARY_DIR}/cli_trials_out)
set_tests_properties(cli_evolve_numtrials PROPERTIES PASS_REGULAR_EXPRESSION "trial 2: best fitness .*2 trials in")
add_test(NAME cli_evolve_numtrials_zero COMMAND netevolve evolve --config ${CMAKE_BINARY_DIR}/cli_small.cfg --numtrials 0)
set_tests_properties(cli_evolve_numtrials_zero PROPERTIES WILL_FAIL TRUE)

# --- prune command (M7) ------------------------------------------------------
add_test(NAME cli_prune COMMAND netevolve prune --model ${EX}/brusselator.ant
                                                --out ${CMAKE_BINARY_DIR}/cli_pruned.ant)
set_tests_properties(cli_prune PROPERTIES PASS_REGULAR_EXPRESSION "pruned 0 of 4 reactions.*long-run check .to 500.: passed")

# --- Acceptance test (M7): label "slow", skipped unless NETEVOLVE_RUN_SLOW=1 --
add_executable(netevolve_acceptance tests/acceptance.cpp)
target_link_libraries(netevolve_acceptance PRIVATE netevolve_core)
# std::getenv is fine here; MSVC's C4996 would otherwise be an error under /WX.
target_compile_definitions(netevolve_acceptance PRIVATE _CRT_SECURE_NO_WARNINGS)
add_test(NAME acceptance_seeds_1_to_5
         COMMAND netevolve_acceptance ${EX}/default.cfg ${CMAKE_BINARY_DIR}/acceptance_out)
set_tests_properties(acceptance_seeds_1_to_5 PROPERTIES LABELS slow SKIP_RETURN_CODE 77 TIMEOUT 3600)
