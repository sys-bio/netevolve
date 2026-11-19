# NetEvolve

NetEvolve is a command-line tool that evolves mass-action biochemical reaction networks toward a target behaviour. The first target is sustained oscillation. The algorithm follows Deckard & Sauro (2004, *ChemBioChem* 5:1423) and Paladugu et al. (2006, *IEE Proc. Syst. Biol.* 153:223).

It's written in C++20 and uses only the standard library. Tests use [doctest](https://github.com/doctest/doctest), which is vendored. The full specification is in [docs/SPEC.md](docs/SPEC.md).

## Status

| Milestone | Contents | State |
|---|---|---|
| M0 | Skeleton, toolchain | done |
| M1 | Network representation, RNG, Antimony I/O | done |
| M2 | ODE model and analytic Jacobian | done |
| M3 | Integrators (RK45, RODAS3) | done |
| M4 | Oscillation fitness; `simulate` and `fitness` commands | done |
| M5 | Mutation and evolutionary loop; `evolve` command | done |
| M6 | Parallel evaluation | not started |
| M7 | `prune`, long-run check, acceptance test | not started |

Details and test results for each milestone are in [PROGRESS.md](PROGRESS.md).

## Building

You need Windows x64, Visual Studio 2022 (MSVC v143) and CMake 3.20 or later. From the repo root:

```bash
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

If `cmake` isn't on your PATH, use the copy installed with CMake (`C:\Program Files\CMake\bin\cmake.exe`) or the one bundled with Visual Studio. The executable is `build/Release/netevolve.exe`.

## Usage

```
netevolve evolve   --config FILE [--seed N] [--out DIR]
netevolve simulate --model FILE.ant [--t-end T] [--points N] [--csv OUT.csv]
netevolve fitness  --model FILE.ant [--config FILE]
netevolve prune    --model FILE.ant [--config FILE] [--out FILE.ant]
netevolve --version
```

`prune` arrives in M7.

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

A run stops after `generations` generations, or earlier once the best fitness reaches `target_fitness`. The same seed and config give byte-identical outputs.

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
fitness 0.999624
```

Exit codes: 0 for success, 1 for a runtime error (for example a bad model file, or a simulation that failed in `simulate`), and 2 for command-line misuse.

## How it works

**Networks.** A network has `F` floating species `S0..S(F-1)`, `B` boundary species `X0..X(B-1)` held at fixed concentrations, and irreversible mass-action reactions. Evolution generates four reaction types: A → P, A → P + Q, A + B → P and A + B → P + Q. The ODE engine handles any integer stoichiometry.

**Simulation.** The tool integrates ds/dt = N·v(s) with adaptive Dormand–Prince RK45. If RK45 detects stiffness or runs out of steps, it switches to RODAS3, an L-stable Rosenbrock method that uses the analytic Jacobian. A simulation fails if it exceeds its step budget or wall-clock limit, if a concentration rises above `max_conc` or falls below −1e-8, or if a value becomes non-finite.

**Evolution.** Each generation is evaluated and sorted by fitness. The top `elite_fraction` are kept unchanged. The rest of the next generation are children of parents drawn uniformly from the top `parent_fraction`. Each child gets exactly one mutation:
- a rate-constant change, k × exp(σ·N(0,1)), clamped to [k_min, k_max];
- adding a random reaction;
- deleting a reaction;
- rewiring one participant to another species.

The chances follow `p_rate : p_add : p_delete : p_rewire`. If a mutation would make the network invalid, it's retried up to 10 times; after that the child is a copy of its parent. Every random draw comes from one `mt19937_64` seeded with `seed`.

**Oscillation fitness.** This follows SPEC §1.4:
1. Discard the first `transient_fraction` of the samples.
2. For each species, find alternating peaks and troughs that differ by more than `max(min_rel_amplitude × mean|y|, 1e-6)`. The absolute floor of 1e-6 stops integrator noise in a species that has decayed to about 0 from counting as oscillation.
3. Count the complete pairs `n`, and compute a damping factor `d` from the first and last pair amplitudes.
4. The species score is `min(n, n_target)/n_target × (1 − d)`. The network fitness is the best species score, in [0, 1]. A failed simulation scores 0.

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

## Checking results in Tellurium

```python
import tellurium as te
r = te.loada(open('examples/brusselator.ant').read())
r.simulate(0, 500, 5000)
r.plot()
```

## Repository layout

```
src/        main.cpp, network, model, integrator, fitness, evolve, config, rng.hpp
tests/      doctest unit tests (test_*.cpp)
examples/   example models and the default config
docs/       SPEC.md
third_party/doctest/doctest.h
```
