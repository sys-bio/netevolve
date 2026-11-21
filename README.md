# NetEvolve

NetEvolve is a command-line tool that evolves mass-action biochemical reaction networks toward a target behaviour. The first target is sustained oscillation. The algorithm follows Deckard & Sauro (2004, *ChemBioChem* 5:1423) and Paladugu et al. (2006, *IEE Proc. Syst. Biol.* 153:223).

It's written in C++20. The default ODE solver is CVODE from [SUNDIALS](https://computing.llnl.gov/projects/sundials) 5.x, linked statically. Everything else uses only the standard library, and the in-house solvers let you build without SUNDIALS (see Building). Tests use [doctest](https://github.com/doctest/doctest), which is vendored. The full specification is in [docs/SPEC.md](docs/SPEC.md).

## Status

| Milestone | Contents | State |
|---|---|---|
| M0 | Skeleton, toolchain | done |
| M1 | Network representation, RNG, Antimony I/O | done |
| M2 | ODE model and analytic Jacobian | done |
| M3 | Integrators (RK45, RODAS3) | done |
| M4 | Oscillation fitness; `simulate` and `fitness` commands | done |
| M5 | Mutation and evolutionary loop; `evolve` command | done |
| M6 | Parallel evaluation | done |
| M7 | `prune`, long-run check, acceptance test | done |

Details and test results for each milestone are in [PROGRESS.md](PROGRESS.md).

## Building

You need Windows x64, Visual Studio 2022 (MSVC v143), CMake 3.20 or later (3.28+ recommended), and internet access the first time you configure. From the repo root:

```bash
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The executable is `build/Release/netevolve.exe`. If `cmake` isn't on your PATH, use the copy installed with CMake (`C:\Program Files\CMake\bin\cmake.exe`) or the one bundled with Visual Studio.

### SUNDIALS (CVODE)

NetEvolve's default ODE solver is CVODE from SUNDIALS 5.x. You don't need to install it: CMake looks for SUNDIALS in this order.

1. **An existing install.** If you already have SUNDIALS 5.x, for example from the libRoadRunner dependencies, point CMake at it. It must be static libraries built with the release DLL runtime (`/MD`):
   ```bash
   cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DSUNDIALS_DIR=<install>/lib/cmake/sundials
   ```
2. **Automatic download (the default).** Otherwise, the first configure downloads the official [SUNDIALS 5.7.0 release](https://github.com/LLNL/sundials/releases/tag/v5.7.0) from GitHub and checks its SHA-256 hash. It then builds only CVODE, as static libraries, inside the build folder (`build/sd-src`, `build/sd`). This adds about 15 s to the first configure and about 30 s to the first build. Later builds reuse it, and nothing is installed outside the build folder. The result is byte-for-byte identical to using a separately installed SUNDIALS 5.7.0. To forbid the download, pass `-DNETEVOLVE_FETCH_SUNDIALS=OFF`; configuration then stops with instructions if no install is found.
3. **No SUNDIALS.** `-DNETEVOLVE_WITH_CVODE=OFF` builds without CVODE. The default solver is then the in-house `auto` (RK45 + Rosenbrock), which is about 2.5× slower and less robust and gives different results.

`netevolve --version` shows which solvers a binary has. CMake prints which SUNDIALS it used, for example `CVODE solver enabled: SUNDIALS 5.7.0 downloaded and built in .../build/sd`.

Notes:
- **Keep the build path short.** SUNDIALS' build has deeply nested folders, and Visual Studio's build logs can exceed Windows' 260-character path limit if the repo sits deep in a folder tree. If you get `FileTracker : error FTK1011 ... could not create the new file tracking log file`, move the repo to a shorter path, or use a short build folder such as `-B C:/nb`.
- **Old build folders.** A build folder configured before CVODE became the default (or with `-DNETEVOLVE_WITH_CVODE=OFF`) keeps that setting cached. Pass `-DNETEVOLVE_WITH_CVODE=ON` once when you reconfigure it, or start a new build folder.
- **Offline builds.** To build offline, download `sundials-5.7.0.tar.gz` once and pass `-DFETCHCONTENT_SOURCE_DIR_SUNDIALS=<extracted sundials-5.7.0 folder>`, or install SUNDIALS and use `SUNDIALS_DIR`.

## Usage

```
netevolve evolve   --config FILE [--seed N] [--out DIR] [--numtrials N]
netevolve simulate --model FILE.ant [--t-end T] [--points N] [--csv OUT.csv]
netevolve fitness  --model FILE.ant [--config FILE]
netevolve prune    --model FILE.ant [--config FILE] [--out FILE.ant]
netevolve --version
```


### Run an evolution

```bash
netevolve evolve --config examples/default.cfg --seed 1 --out runs/seed1
```

Flags override the config file. Progress is printed once per generation. The output directory then contains:

| File | Contents |
|---|---|
| `log.csv` | `generation,best,mean,median,best_reactions`, one row per generation (generation 0 is the random initial population) |
| `best_genNNNN.ant` | the best network, written whenever the best fitness improves |
| `best_final.ant` | the best network at the end of the run |
| `best_final.csv` | time course of the best network (`t_end`, `n_points` from the config) |
| `run.cfg` | the fully resolved config, including the seed used |

A run stops after `generations` generations, or earlier once the best fitness reaches `target_fitness`. The same seed and config give byte-identical outputs, whatever the thread count.

Fitness evaluation runs on `threads` threads (0 means all cores). On the development machine, an 8-thread laptop CPU, the in-house `auto` solver evaluates about 1,450 networks per second on one thread and about 9,000 per second on 8, so a full 500-generation run takes roughly 10–15 seconds. CVODE, the default, is about 2.5× faster per evaluation (measured on a different machine).

### Run several independent trials

```bash
netevolve evolve --config examples/default.cfg --seed 1 --numtrials 50 --out runs/batch
```

`--numtrials N` runs N independent evolutions in one command. The seed is set once: it seeds a master `mt19937_64`, and trial i runs with that generator's i-th output as its seed. Trial i writes the usual outputs to `runs/batch/trial<i>/`, and its `run.cfg` records its own seed, so any trial can be rerun alone with `--seed`. `runs/batch/summary.csv` has one row per trial: `trial,seed,best_fitness,best_generation,reactions,evaluations,long_run`. Without `--numtrials`, `evolve` runs a single trial exactly as before.

### Prune a network

Evolved networks often contain reactions they don't need. `prune` tries removing each reaction in turn and keeps a removal if the fitness stays at least 0.95 × the original. It repeats this until no more reactions can be removed, and never goes below `r_min`.

```bash
netevolve prune --model runs/seed1/best_final.ant --config examples/default.cfg --out runs/seed1/pruned.ant
```
```
pruned 2 of 8 reactions: fitness 1 -> 1
long-run check (to 500): passed, S0 amplitude 138.195 -> 144.057 (4.24171%)
```

Without `--out`, the pruned model goes to stdout and the summary to stderr.

**Long-run check.** This distinguishes sustained oscillation from slow damping, and both `prune` and `evolve` report it. The network is re-simulated to 5 × `t_end`. For its best-scoring species, the mean peak-to-trough amplitude in the last window must be within 10% of the amplitude in the window ending at `t_end`. Both windows have the fitness window's length.

### Simulate a model

Prints the time course as CSV (`time,S0,S1,...`), or writes it to a file with `--csv`. The defaults are `t_end = 100` and 1000 samples.

```bash
netevolve simulate --model examples/brusselator.ant --t-end 50 --points 501 --csv bruss.csv
```

### Score a model

Simulates the model with the settings from the config file (or the defaults) and prints each species' score and the network fitness:

```bash
netevolve fitness --model examples/brusselator.ant --config examples/default.cfg
```
```
status: ok
S0: pairs 7, damping 0.000376191, score 0.999624
S1: pairs 6, damping 0.000955597, score 0.999044
re-check (rtol 1e-10, atol 1e-14): status ok, score 0.999621
fitness 0.999621
```

The per-species scores come from the default-tolerance simulation. The `re-check` line appears when that score reaches `recheck_threshold`. The final `fitness` is the value used for selection and pruning.

### Evolve toward a time series

Instead of oscillation, the target can be a time series read from a text file. Set `objective = timeseries` and `target_file` in the config ([examples/timeseries.cfg](examples/timeseries.cfg)):

```bash
netevolve evolve --config examples/timeseries.cfg --seed 1 --numtrials 10 --out runs/ts
```

The target file has two columns, time and concentration, separated by spaces, tabs or a comma. `#` starts a comment, and one text header line is allowed. Times must be strictly increasing. Alternating peak and trough values at the times they should occur work well. [examples/target_peaks.txt](examples/target_peaks.txt) has the troughs and peaks of the Brusselator's S0. A relative `target_file` is relative to the config file. It is stored as an absolute path in each run's `run.cfg`.

Each network is simulated from t = 0, with every species starting at its initial value (1 for evolved networks), and sampled at exactly the target times. For each species, the error is the RMS difference from the targets divided by the mean |target| (NRMSE), and the score is 1 / (1 + NRMSE). That gives 1 for a perfect match, about 0.91 for a 10% error and 0.5 for a 100% error. The network's fitness is its best-matching species. The target doesn't say which species, so evolution may use any of them. The tight-tolerance re-check, `prune` and `--numtrials` work as for oscillation. The long-run check is specific to oscillation and is reported as n/a. `netevolve fitness --config examples/timeseries.cfg --model M.ant` prints each species' NRMSE.

Only the values at the target times are scored. A network can pass through a "peak" value without having a maximum there, so dense targets (for example including points between the peaks and troughs) constrain the shape more tightly.

Exit codes: 0 for success, 1 for a runtime error (for example a bad model file, or a simulation that failed in `simulate`), and 2 for command-line misuse.

## How it works

**Networks.** A network has `F` floating species `S0..S(F-1)`, `B` boundary species `X0..X(B-1)` held at fixed concentrations, and irreversible mass-action reactions. Evolution generates four reaction types: A → P, A → P + Q, A + B → P and A + B → P + Q. The ODE engine handles any integer stoichiometry.

**Simulation.** The tool integrates ds/dt = N·v(s). The solver is set by the `solver` config key:
- `cvode` (the default): SUNDIALS CVODE, variable-order BDF with the analytic Jacobian and a dense linear solver. It runs one internal step at a time, so the guards below are checked after every step. Samples are interpolated with CVODE's dense output.
- `auto`: the in-house adaptive Dormand–Prince RK45. If RK45 detects stiffness or runs out of steps, it switches to RODAS3, an L-stable Rosenbrock method that uses the analytic Jacobian.
- `rk45` or `rosenbrock` alone.

In a 50-trial comparison (`--seed 1 --numtrials 50`), CVODE ran in 110 s against 258 s for `auto`, about 22,000 against 8,000 evaluations per second, with the same yield (15 against 16 trials passing the long-run check). RK45 alone was 20–85× slower than `auto`. `max_steps` (default 50,000) is the step budget per simulation. A simulation fails if it exceeds its step budget or wall-clock limit, if a concentration rises above `max_conc` or falls below −1e-8, or if a value becomes non-finite.

**Evolution.** Each generation is evaluated and sorted by fitness. The top `elite_fraction` are kept unchanged. The rest of the next generation are children of parents drawn uniformly from the top `parent_fraction`. Each child gets exactly one mutation:
- a rate-constant change, k × exp(σ·N(0,1)), clamped to [k_min, k_max];
- adding a random reaction;
- deleting a reaction;
- rewiring one participant to another species.

The chances follow `p_rate : p_add : p_delete : p_rewire`. If a mutation would make the network invalid, it's retried up to 10 times; after that the child is a copy of its parent. Every random draw comes from one `mt19937_64` seeded with `seed`, and all of them happen in the serial part of the loop. Only the pure fitness evaluations run in parallel.

**Oscillation fitness.** This follows SPEC §1.4:
1. Discard the first `transient_fraction` of the samples.
2. For each species, find alternating peaks and troughs that differ by more than `max(min_rel_amplitude × mean|y|, 1e-6)`. The absolute floor of 1e-6 stops integrator noise in a species that has decayed to about 0 from counting as oscillation.
3. Count the complete pairs `n`, and compute a damping factor `d` from the first and last pair amplitudes.
4. The species score is `min(n, n_target)/n_target × (1 − d)`. The network fitness is the best species score, in [0, 1]. A failed simulation scores 0.
5. **Tight-tolerance re-check.** If the score is at least `recheck_threshold` (default 0.9), the network is simulated again with `recheck_rtol` / `recheck_atol` (defaults 1e-10 / 1e-14), and its fitness is the lower of the two scores. This applies to selection, `prune`, `fitness` and the acceptance test. Without it, evolution finds networks whose oscillation exists only because of integrator error at the default tolerance (rtol 1e-6). In a 50-trial test, 2 of 37 high-scoring winners were such artifacts and 3 more had inflated scores. CVODE (libRoadRunner) failed to integrate most of them. [examples/integrator_artifact.ant](examples/integrator_artifact.ant) is one: with `solver = auto` it scores 1 at the default tolerance and 0 at the tight one, and CVODE scores it 0 at both. The re-check stays on with CVODE. Set `recheck_threshold` above 1 to turn the re-check off.
6. **Concentration floor.** If a species that oscillates in the scored window (at least one peak-trough pair) drops below `osc_conc_floor` (default 1e-12) there, the network's fitness is 0. This covers the whole network, not just that species, because other species can oscillate cleanly *because of* the one crashing. Such dips are physically meaningless (far below one molecule) and numerically unreliable. [examples/deep_dip.ant](examples/deep_dip.ant), an evolved network whose species fall to about 1e-22 every cycle, scored 1.0 with CVODE before the floor. Solvers disagreed on it, and RK45 and LSODA failed. Species that just decay toward 0 without oscillating don't trigger the floor. Set `osc_conc_floor = 0` to turn it off.

## Configuration

A config file uses plain `key = value` lines, with `#` starting a comment. [examples/default.cfg](examples/default.cfg) lists every key with its default value. Unknown keys and malformed values are errors, reported with their line number.

## Model files (Antimony subset)

Networks are read and written in a restricted subset of [Antimony](https://tellurium.readthedocs.io/en/latest/antimony.html), so they load straight into Tellurium:

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

The parser accepts:
- species named `S<n>` (floating) or `X<n>` (boundary, which must be declared `const`);
- coefficients such as `2 S0`;
- `//` comments.

Each rate law must be one parameter multiplied by the reactant species, matching mass action. A mismatch, or any syntax outside the subset, is an error with a line number.

Examples:
- [examples/brusselator.ant](examples/brusselator.ant): the Brusselator with A = 1, B = 3. It has a limit cycle, and its fitness is about 1.
- [examples/damped.ant](examples/damped.ant): the same network with B = 1.5. It has a stable focus, and its fitness is 0.

## Tests

```bash
ctest --test-dir build -C Release --output-on-failure
```

This runs everything except the acceptance test, which has the label `slow` and is reported as skipped. The acceptance test runs five full evolutions (seeds 1–5, default config, about 15 s on 8 threads). It requires at least 3 of them to reach fitness ≥ 0.95, with the pruned result passing the long-run check. Run it with:

```bash
NETEVOLVE_RUN_SLOW=1 ctest --test-dir build -C Release -L slow -V
```

It keeps each run's output, plus `pruned.ant`, in `build/acceptance_out/seedN/`. (In PowerShell, set the variable with `$env:NETEVOLVE_RUN_SLOW = "1"` first.)

## Checking results in Tellurium

```python
import tellurium as te
r = te.loada(open('best_final.ant').read())
r.simulate(0, 500, 5000)
r.plot()
```

Evolved oscillators are often spiky relaxation oscillators, with species dipping to about 1e-6 on every cycle. libRoadRunner's CVODE can stop on these with `CV_TOO_MUCH_WORK`. If that happens, tighten the tolerances before simulating:

```python
r.integrator.absolute_tolerance = 1e-14
r.integrator.relative_tolerance = 1e-10
```

This fixed seed 5 of the acceptance run. Some networks behave the other way round: they run with the defaults and stall with tight tolerances. NetEvolve's own solver and SciPy's Radau/LSODA agreed on all of them (see PROGRESS.md, M7).

## Repository layout

```
src/        main.cpp, network, model, integrator, fitness, evolve, config, rng.hpp
tests/      doctest unit tests (test_*.cpp) and the acceptance test (acceptance.cpp)
examples/   example models and the default config
docs/       SPEC.md
third_party/doctest/doctest.h
```
