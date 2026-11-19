# Progress

## Current milestone

**M0 (skeleton and toolchain) is complete.** It's waiting for the human's approval before M1 starts.

## M0: Skeleton and toolchain

Done:
- Created the directory layout from SPEC §0.3. Module source files will be added in the milestones that need them.
- `CMakeLists.txt` builds `netevolve.exe` and `netevolve_tests.exe` with `/W4 /permissive- /EHsc`, using C++20. Release gets `/O2` from CMake's defaults.
- Vendored doctest v2.4.12 at `third_party/doctest/doctest.h`, downloaded from the official GitHub tag.
- `netevolve --version` prints `netevolve 0.1.0`. The version comes from `project(VERSION)` in CMake.
- Added `CLAUDE.md`, `.gitignore` (which excludes `build/` and `out/`), and a git repo.

Tests (Release, from an empty `build/`): **3/3 passed**
- `unit_tests`: the doctest smoke case
- `cli_version_exit`: `--version` exits with 0
- `cli_version_output`: the output matches `^netevolve X.Y.Z`

## Decisions and notes on ambiguity

- The spec is copied verbatim to `docs/SPEC.md`. The original stays at `../docs/netevolve-spec.md`.
- The `--version` check is two CTest tests: one checks the exit code, the other matches the output. A regex match by itself would ignore the exit code.

## Open issues

- No GitHub remote yet. Creating one needs the human's go-ahead.
