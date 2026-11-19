# NetEvolve — Specification for an Evolutionary Network Generator (C++)

**Audience:** a Claude Code session running in agentic mode, and the researcher supervising it.
**Goal:** a command-line tool that evolves mass-action biochemical reaction networks toward a target behaviour, with sustained oscillation as the first target. The algorithm follows the approach of Deckard & Sauro (2004, *ChemBioChem* 5:1423) and Paladugu et al. (2006, *IEE Proc. Syst. Biol.* 153:223).

---

## 0. Instructions to the Claude session (read first)

You are building this project agentically: you write the code, run the builds, run the tests, and fix failures yourself. The human supervises and approves milestones.

### 0.1 Working rules

1. **Work one milestone at a time, in order.** Do not start milestone *n+1* until every test for milestone *n* passes in a Release build.
2. **Tests first or alongside.** Each milestone lists required tests. Implement them, run them, and paste a short summary of the results (pass/fail counts) when you report.
3. **Never mark a milestone complete while tests fail.** Do not weaken, skip or delete a test to make it pass. If you believe a test is wrong, stop and explain why to the human.
4. **Keep `PROGRESS.md` up to date** at the repo root. Record the current milestone, what was done, test status, and open issues. A future session should be able to resume from it.
5. **Create `CLAUDE.md`** at the repo root in Milestone 0. It should hold a short summary of these rules, the build and test commands that actually worked on this machine, and a pointer to this spec (`docs/SPEC.md`).
6. **Use git.** Initialise a repository in M0 and commit at the end of each milestone with message `M<n>: <title>`. Commit smaller steps as well when useful.
7. **Dependencies:** standard C++20 library only, plus the test framework below. Do not add other libraries (Boost, Eigen, SUNDIALS, JSON libraries, etc.) without asking the human first.
8. **Determinism:** every stochastic component draws from a single seeded RNG (`std::mt19937_64`). The same seed and config must give bit-identical results on the same build.
9. **Keep it simple.** Prefer plain structs and free functions over class hierarchies. Do not add plugin systems, abstract factories, or configurability that the spec doesn't ask for.
10. If something in this spec is ambiguous, choose the simplest reasonable reading, note it in `PROGRESS.md`, and carry on. Only stop to ask when the choice would be expensive to undo.

### 0.2 Environment

- **OS:** Windows 10/11, x64.
- **Compiler:** MSVC from Visual Studio 2022 (toolset v143), using C++20 (`/std:c++20`).
- **Build system:** CMake (≥ 3.20) with the `Visual Studio 17 2022` generator.
- **Shell:** the session's shell may be PowerShell or Git Bash. Check which one you have before writing commands.

**Locating CMake.** CMake ships with VS 2022 but is often not on `PATH`. If `cmake --version` fails, find it with `vswhere`:

```powershell
$vs = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -property installationPath
$cmake = "$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
& $cmake --version
```

The `ctest.exe` executable is in the same folder. Record the working paths in `CLAUDE.md`.

**Standard commands** (run from the repo root):

```
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Use Debug builds when you need to debug, but every milestone is accepted on Release.

**Test framework:** [doctest](https://github.com/doctest/doctest), as a single vendored header at `third_party/doctest/doctest.h`. If you cannot download it, ask the human to place the file there. Register test executables with CTest via `add_test`.

**Compiler flags:** use `/W4 /permissive- /EHsc`, plus `/O2` in Release. Treat warnings as errors only from M2 onwards.

### 0.3 Repository layout

```
netevolve/
  CMakeLists.txt
  CLAUDE.md
  PROGRESS.md
  docs/SPEC.md            <- this file
  src/
    main.cpp              CLI entry point
    network.hpp/.cpp      network representation, random generation, Antimony I/O
    model.hpp/.cpp        ODE right-hand side and Jacobian from a network
    integrator.hpp/.cpp   ODE solvers
    fitness.hpp/.cpp      objective functions
    evolve.hpp/.cpp       mutation, selection, generation loop
    config.hpp/.cpp       configuration parsing
    rng.hpp               seeded RNG wrapper
  tests/
    test_*.cpp
  examples/
    brusselator.ant
    damped.ant
    default.cfg
  third_party/doctest/doctest.h
```

---

## 1. Problem description

### 1.1 What is being evolved

A **network** consists of:

- `F` **floating species** `S0 … S(F-1)`. These are dynamic, with initial concentrations given.
- `B` **boundary species** `X0 … X(B-1)`. These are held at fixed concentrations and act as sources and sinks. Open systems are required: a closed mass-action network relaxes to equilibrium and cannot sustain oscillations.
- `R` **reactions**. Each is irreversible mass action with rate constant `k > 0`.

Allowed reaction **types** (reactants → products), where any participant may be floating or boundary:

| Type | Form | Rate |
|---|---|---|
| UniUni | A → P | k·A |
| UniBi | A → P + Q | k·A |
| BiUni | A + B → P | k·A·B |
| BiBi | A + B → P + Q | k·A·B |

Same-species stoichiometry is allowed (e.g. `A + A → P`, rate `k·A²`; `A → A + A`, which is autocatalytic). The ODE engine (M2) must support general integer stoichiometry so that test models such as the Brusselator (which contains a trimolecular step) can be loaded, even though mutation only generates the four types above.

**Validity constraints** on any network produced by generation or mutation:
- No reaction has identical reactant and product multisets.
- Every reaction involves at least one floating species.
- `R_min ≤ R ≤ R_max` (configurable).
- No reaction duplicates another reaction (same reactants and products). Duplicates should be merged by summing their rate constants, or rejected.

### 1.2 The dynamics

For floating species concentrations **s**, the dynamics are d**s**/dt = **N v**(**s**). Here **N** is the F×R stoichiometry matrix restricted to floating species, and **v** is the vector of mass-action rates (boundary species enter **v** as constants). The analytic Jacobian **J** = **N** ∂**v**/∂**s** is also required for the stiff solver.

### 1.3 The evolutionary algorithm

Outline, following the early papers:

```
initialise population of P random valid networks
repeat for G generations (or until fitness >= target):
    evaluate fitness of every network
    sort by fitness (higher is better)
    keep the top E% unchanged (elitism)
    fill the remaining P - E slots:
        pick a parent from the top fraction T% (uniformly)
        child = copy(parent)
        apply one mutation (see below)
        if child invalid: retry mutation (bounded attempts) or recopy
    log generation statistics
```

**Mutation operators.** Each child receives exactly one mutation, chosen with configurable probabilities:
- **Rate-constant change** (most frequent, default p = 0.7): pick a reaction and multiply its `k` by `exp(σ·N(0,1))`, then clamp to `[k_min, k_max]`.
- **Add reaction** (default p = 0.15): random type, random participants, and log-uniform `k`.
- **Delete reaction** (default p = 0.10): only if `R > R_min`.
- **Rewire** (default p = 0.05): replace one participant of one reaction with another species.

Initial concentrations of floating species are fixed (default 1.0 for each, configurable) and are not evolved in v1. Boundary concentrations are fixed (default 1.0).

### 1.4 Objective: peak-and-trough oscillation fitness

Simulate from t = 0 to `t_end` and sample at `n_points` uniform times. Discard the first `transient_fraction` (default 0.5) of the samples. For each floating species, take the remaining trajectory `y` and compute:

1. **Extrema detection.** A sample is a peak if `y[i-1] < y[i] ≥ y[i+1]`; troughs are defined symmetrically. Keep only alternating peak/trough pairs whose height difference exceeds `min_rel_amplitude × mean(y)` (default 0.01). This filters out numerical ripple.
2. **Count** `n` = the number of complete peak-trough pairs.
3. **Damping.** Compare the amplitude of the first pair, `a_first`, with the last pair, `a_last`: `d = clamp((a_first - a_last)/a_first, 0, 1)`.
4. **Species score:** `f = min(n, n_target)/n_target × (1 - d)`, with `n_target` defaulting to 5.

The **network fitness** is the maximum species score over all floating species, so it lies in [0, 1].

**Fitness is 0** if any of the following happen: the integrator fails or exceeds its step budget, any concentration becomes negative (below `-1e-8`), any concentration exceeds `max_conc` (default 1e6), or the result is non-finite.

Design the fitness code so that other objectives (e.g. the steady-state "square root" target from the 2004 paper) can be added later as a separate function. **Do not** build a general objective framework now.

---

## 2. Configuration and CLI

A plain `key = value` config file, with `#` for comments. No JSON library.

```
seed = 42
population = 200
generations = 500
elite_fraction = 0.1
parent_fraction = 0.3
n_floating = 4
n_boundary = 2
r_min = 3
r_max = 12
r_initial = 6
k_min = 0.01
k_max = 100
mut_sigma = 0.5
p_rate = 0.7
p_add = 0.15
p_delete = 0.10
p_rewire = 0.05
t_end = 100
n_points = 1000
transient_fraction = 0.5
n_target = 5
min_rel_amplitude = 0.01
max_conc = 1e6
target_fitness = 0.99
output_dir = out
threads = 0          # 0 = use all cores
```

Command-line interface:

```
netevolve evolve   --config FILE [--seed N] [--out DIR]
netevolve simulate --model FILE.ant [--t-end T] [--points N] [--csv OUT.csv]
netevolve fitness  --model FILE.ant [--config FILE]
netevolve prune    --model FILE.ant [--config FILE] [--out FILE.ant]
netevolve --version
```

Command-line flags override values from the config file. Unknown keys or flags are errors and are reported clearly.

**Outputs of `evolve`** (in `output_dir`):
- `log.csv`: `generation,best,mean,median,best_reactions` per generation.
- `best_genNNNN.ant`: the best network, written whenever best fitness improves.
- `best_final.ant` and `best_final.csv` (time course of the best network).
- `run.cfg`: the fully resolved config, including the seed actually used.

---

## 3. Antimony I/O

Networks are exported, and loaded, in a restricted subset of Antimony. This makes results directly loadable in Tellurium for independent checking. The exported form is:

```
// netevolve: generation 213, fitness 0.97
const X0, X1;
J0: X0 -> S1; k0*X0;
J1: S1 + S2 -> S0; k1*S1*S2;
J2: S0 -> S2 + S2; k2*S0;
J3: S2 -> X1; k3*S2;
k0 = 1.23; k1 = 0.45; k2 = 3.1; k3 = 0.8;
S0 = 1; S1 = 1; S2 = 1;
X0 = 1; X1 = 1;
```

The **parser** only needs to accept this subset:
- `const` declarations;
- reactions written as `name: reactants -> products; rate;`, where integer coefficients like `2 S1` are allowed, and the rate is a product of a parameter and species names;
- numeric assignments;
- `//` comments.

It may ignore the written rate law and rebuild mass action from the reactants, but it must check that the two agree. A mismatch is a parse error. Anything outside the subset is a clear error, not silently ignored.

---

## 4. Milestones and tests

Each milestone lists its deliverables and the tests that must pass. The tolerances given are requirements.

### M0 — Skeleton and toolchain
**Deliverables:** directory layout; `CMakeLists.txt` building `netevolve.exe` and a `netevolve_tests` executable; vendored doctest; `CLAUDE.md`; `PROGRESS.md`; git repo with `.gitignore` (excluding `build/`, `out/`).
**Tests:**
- A trivial doctest case passes under `ctest`.
- `netevolve --version` prints a version string and exits with code 0 (add this as a CTest test).
- A clean configure, build and test from an empty `build/` succeeds.

### M1 — Network representation, RNG, Antimony I/O
**Deliverables:** `Network` struct; validity checker; random network generator; Antimony writer and parser.
**Tests:**
- Random generation with a fixed seed yields identical networks across two calls.
- 10,000 random networks all pass the validity checker.
- Round trip: write → parse → write produces identical text for 1,000 random networks.
- Parsing `examples/brusselator.ant` gives the expected stoichiometry and constants.
- Malformed input (unknown species, missing semicolon, rate law not matching the reactants) produces an error with a line number.

### M2 — ODE model and Jacobian
**Deliverables:** construction of **N**, the rate function `v(s)`, the RHS `ds/dt`, and the analytic Jacobian.
**Tests:**
- Hand-computed RHS for a 3-reaction network matches to 1e-12.
- The analytic Jacobian matches a central finite-difference Jacobian to a relative error of 1e-6 for 100 random networks at random positive states.
- Boundary species do not appear as state variables and have zero derivative.
- Mass conservation: for a closed network `A → B → A`, the RHS sums to 0 (within 1e-12).

### M3 — Integrators
**Deliverables:**
- An adaptive explicit solver: Dormand–Prince RK45 with embedded error control, `rtol`/`atol`, and dense or interpolated output at the requested sample times.
- A stiff solver: an L-stable Rosenbrock method, e.g. ROS3P or Rodas-type, using the analytic Jacobian, with an LU solve written in-house (systems are small, F ≤ ~20).
- A step budget and a wall-clock guard per simulation.

Default policy: try RK45, and switch to the stiff solver if it signals stiffness or exceeds the step budget.

**Tests:**
- `dy/dt = -k y` matches the analytic solution to a relative error of 1e-6 at t = 10, for both solvers.
- Lotka–Volterra (expressed as mass action, with boundary prey food): the conserved quantity drifts by less than 1e-5 relative over t = 50.
- The Robertson stiff problem is solved by the Rosenbrock solver to t = 1e4 within the step budget, with `y1+y2+y3 = 1` to 1e-6.
- The Brusselator (A = 1, B = 3; unstable focus, so a limit cycle) settles to a period within 2% of a reference value. Obtain the reference independently (e.g. a very tight-tolerance run) and record it in the test with a comment saying how it was obtained.
- A deliberately explosive network (e.g. `S → 2S` with large k) triggers the `max_conc` guard rather than hanging.

### M4 — Fitness function
**Deliverables:** extrema detection and the peak-and-trough fitness described in §1.4, plus the `simulate` and `fitness` CLI commands.
**Tests on synthetic signals** (no ODEs):
- `sin(t)` sampled over 10 periods → `n ≥ n_target`, `d ≈ 0`, fitness ≈ 1.
- `e^{-0.3t} sin(t)` → fitness substantially below 1 (assert < 0.5 with default settings).
- A constant signal and a monotone ramp → fitness 0.
- `sin(t)` plus noise at 0.1% of the amplitude → extrema are not double-counted.

**Tests on models:**
- `examples/brusselator.ant` → fitness ≥ 0.95.
- `examples/damped.ant`, a network converging to a stable focus → fitness < 0.5.
- A network that blows up → fitness exactly 0.

### M5 — Mutation and the evolutionary loop
**Deliverables:** the four mutation operators, the selection and elitism loop, logging, and the `evolve` command (single-threaded).
**Tests:**
- Every mutation operator, applied 10,000 times to random networks, always yields valid networks.
- Rate mutation keeps `k` in `[k_min, k_max]`.
- Elitism: with a fixed seed, best fitness is non-decreasing across generations.
- Determinism: two runs with the same seed produce byte-identical `log.csv`.
- A toy fitness (e.g. "reward networks with exactly 7 reactions") is solved within 50 generations. This tests the loop independently of the ODE machinery.

### M6 — Parallel evaluation and performance
**Deliverables:** parallel fitness evaluation using `std::thread` or OpenMP (`/openmp` is available in MSVC). Results must not depend on the thread count. To ensure this, all random draws happen in the serial part; fitness evaluation is pure.
**Tests:**
- `log.csv` is byte-identical for `threads = 1` and `threads = 8`.
- Report evaluations per second for the default config in `PROGRESS.md`. There is no pass threshold; this is for information only.

### M7 — End-to-end oscillator evolution and validation
**Deliverables:**
- The `prune` command. It removes each reaction in turn and accepts the removal if fitness stays ≥ 0.95 × the original, repeating until no further removal is possible. This mirrors the observation in the early papers that evolved networks often contain dispensable reactions.
- A long-run check: re-simulate the best network to `5 × t_end` and confirm that the amplitude in the last window is within 10% of the window at `t_end`. This distinguishes sustained oscillation from slow damping.

**Acceptance test** (registered in CTest with the label `slow`, and excluded from the default run):
- With the default config and seeds 1–5, at least 3 of the 5 runs reach fitness ≥ 0.95 within 500 generations, and the pruned result passes the long-run check.
- The exported `.ant` loads in Tellurium. This step is manual and done by the human: `te.loada(open('best_final.ant').read()).simulate(0, 500, 5000)` should show sustained oscillation.

If the acceptance test fails, do not tune the tests. Report the fitness curves and propose changes to the parameters (population, mutation rates, `t_end`, reaction limits) to the human.

---

## 5. Out of scope for v1

Evolving initial concentrations; reversible or non-mass-action rate laws; SBML export; GUI; objectives other than oscillation; crossover. These are possible later extensions. Do not scaffold for them.

## 6. Definition of done

- M0–M7 complete, with all non-slow tests passing in Release.
- The acceptance test has been run and its result recorded in `PROGRESS.md`.
- `README.md` covers how to build, how to run an evolution, the output files, and how to inspect results in Tellurium.
