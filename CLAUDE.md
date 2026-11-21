# NetEvolve

An evolutionary generator of mass-action reaction networks, written in C++20. The full specification is in [docs/SPEC.md](docs/SPEC.md); read it before starting work. Current status is in [PROGRESS.md](PROGRESS.md).

## Rules (summary of SPEC §0.1)

- Work one milestone at a time, in order. Don't start M(n+1) until every M(n) test passes in a **Release** build.
- Write tests first or alongside the code. When reporting, give pass/fail counts.
- Never weaken, skip or delete a test to make it pass. If a test looks wrong, stop and explain why to the human.
- Keep `PROGRESS.md` current: the milestone, what was done, test status, and open issues.
- Commit at the end of each milestone as `M<n>: <title>`. Smaller commits are fine.
- Use only the C++20 standard library plus vendored doctest, and SUNDIALS 5.x (CVODE), which the PI approved on 2026-11-20 as the default solver. Ask before adding anything else.
- Determinism: draw everything from one seeded `std::mt19937_64`. The same seed and config must give bit-identical results.
- Keep it simple: plain structs and free functions, nothing the spec doesn't ask for.
- When the spec is ambiguous, pick the simplest reading and note it in `PROGRESS.md`.
- Treat warnings as errors from M2 onwards.

## Build and test (verified on this machine)

`cmake` and `ctest` aren't on PATH, so use the full path. From Git Bash, in the repo root:

```bash
"/c/Program Files/CMake/bin/cmake.exe" -S . -B build -G "Visual Studio 17 2022" -A x64
"/c/Program Files/CMake/bin/cmake.exe" --build build --config Release -- -m
"/c/Program Files/CMake/bin/ctest.exe" --test-dir build -C Release --output-on-failure
```

From PowerShell, use `& "C:\Program Files\CMake\bin\cmake.exe" ...` with the same arguments.

- Compiler: MSVC from Visual Studio 2022 Community (v143).
- The Visual Studio generator is multi-config, so always pass `--config` / `-C`.
- Outputs go to `build/Release/netevolve.exe` and `build/Release/netevolve_tests.exe`.
- The M7 acceptance test has the CTest label `slow`. It's skipped (exit code 77) unless `NETEVOLVE_RUN_SLOW=1` is set: `NETEVOLVE_RUN_SLOW=1 "/c/Program Files/CMake/bin/ctest.exe" --test-dir build -C Release -L slow -V`. It takes about 15 s on 8 threads.

## Layout

The CLI is in `src/`, the doctest tests in `tests/`, example models and configs in `examples/`, and doctest v2.4.12 in `third_party/doctest/doctest.h`.

## Workspace note

This repo sits in `D:\zyrr\Research\UW\RoadRunnerProject2026\netevolve`, next to an unrelated libroadrunner workspace (`../roadrunner`, `../build`, etc.). Don't touch anything outside `netevolve/`.
