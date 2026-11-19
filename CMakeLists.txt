cmake_minimum_required(VERSION 3.20)
project(netevolve VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
