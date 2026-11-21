# NetEvolve Web: Specification for a Rust/WebAssembly Port

**Audience:** a Claude Code session working agentically, and the researcher supervising it.

**Goal:** a browser application, hosted on GitHub Pages, that evolves mass-action biochemical reaction networks toward a target behaviour. It covers two targets: sustained oscillation, and matching a user-supplied time series.

**Reference implementation:** the C++ NetEvolve in this repository, branch `cvode-timeseries`. Part A describes the algorithm as that implementation actually behaves, including details the original `docs/SPEC.md` left open or that have changed since. **Where this document and the C++ code disagree, the C++ code is authoritative. Report the disagreement; don't silently pick one.**

Background: Deckard & Sauro (2004, *ChemBioChem* 5:1423); Paladugu et al. (2006, *IEE Proc. Syst. Biol.* 153:223).

---

## 0. Instructions to the Claude session (read first)

### 0.1 Working rules

1. **Work one milestone at a time, in order** (§C). Don't start milestone *n+1* until every test for milestone *n* passes in a release build (`cargo test --release`).
2. **Write tests first or alongside the code.** Report pass/fail counts.
3. **Never weaken, skip or delete a test to make it pass.** If a test looks wrong, stop and explain why.
4. **Keep `PROGRESS.md` current:** the milestone, what was done, test status and open issues. Create `CLAUDE.md` in the first milestone. It should hold these rules in short form, the build/test/serve commands that actually worked, and a pointer to this file.
5. **Commit at the end of each milestone** as `R<n>: <title>`.
6. **Dependencies:** keep them few, and record every crate with its reason in `PROGRESS.md`. Ask before adding anything heavy (a UI framework, a linear-algebra library beyond what the ODE solver needs).
7. **Determinism:** the same seed and config give bit-identical results within this implementation, whatever the number of workers. Bit-identity with the C++ build is *not* required (see §A.11).
8. **Keep it simple:** plain structs and free functions. Build nothing this spec doesn't ask for.
9. **Ambiguities:** pick the simplest reading, note it in `PROGRESS.md`, and continue. Ask only when a choice would be expensive to undo.

### 0.2 Environment and project layout

A Cargo workspace with three crates:

```
netevolve-web/
  Cargo.toml              workspace
  CLAUDE.md  PROGRESS.md  README.md
  docs/SPEC_RUST_WEB.md   <- this file (copied)
  core/                   library: everything in Part A; no browser dependencies
  cli/                    native command-line tool over `core` (for testing and
                          comparison with the C++ build; mirrors its commands)
  web/                    wasm-bindgen bindings, Web Worker glue, static site
  fixtures/               copied from the C++ repo's examples/ (see §C.0)
```

- **`core`** must compile both natively and for `wasm32-unknown-unknown`. It must not use threads, the file system or the clock directly. Wall-clock limits go through an injected time source (§A.5).
- **The site** is static files (HTML, JS, the `.wasm` file) deployable to GitHub Pages by a GitHub Actions workflow. No server code.

---

# Part A: The algorithm (language-neutral)

## A.1 Networks

A **network** has:
- `F` floating species `S0 … S(F-1)`, with initial concentrations (default 1.0 each);
- `B` boundary species `X0 … X(B-1)`, fixed concentrations (default 1.0);
- `R` irreversible mass-action reactions, each with a rate constant `k > 0`.

**Species index:** a single integer. `0..F-1` are floating, and `F..F+B-1` are boundary.

**Reaction:** `reactants` and `products`, each a list of `(species, coeff ≥ 1)` terms in **canonical form**: sorted by species index, one term per species. The rate is `v = k · Π reactants s_i^coeff`; boundary species enter as constants. General integer stoichiometry must be supported, because the Brusselator fixture has a trimolecular step. Generation and mutation only produce four shapes: A→P, A→P+Q, A+B→P and A+B→P+Q, where any participant may repeat (e.g. A+A→P, A→A+A).

**Validity**, checked after every generation or mutation:
- **Size:** `F ≥ 1`, and `r_min ≤ R ≤ r_max`.
- **Terms and values:** terms are canonical with valid indices, every `k` is finite and > 0, and concentrations are finite and ≥ 0.
- **Each reaction:** reactants ≠ products, and at least one participant is floating.
- **No duplicates:** no two reactions have identical reactant *and* product lists.

## A.2 Random numbers

There is one generator per run: **MT19937-64** (`std::mt19937_64`), seeded with the 64-bit `seed`. The helpers must be implemented exactly like this, not with a library's distributions, whose algorithms vary:

| Helper | Definition |
|---|---|
| `next()` | next 64-bit engine output |
| `uniform01()` | `(next() >> 11) * 2^-53`, in [0, 1) |
| `uniform(lo, hi)` | `lo + (hi - lo) * uniform01()` |
| `uniform_int(lo, hi)` (inclusive) | `span = hi - lo + 1`; `limit = U64_MAX - (U64_MAX % span)`; repeat `x = next()` while `x >= limit`; return `lo + x % span` |
| `normal()` | repeat `u1 = uniform01()` while `u1 <= 0`; `u2 = uniform01()`; return `sqrt(-2 ln u1) · cos(2π u2)` |
| `log_uniform(lo, hi)` | `exp(uniform(ln lo, ln hi))` |

The `rand_mt` crate's `Mt64` is expected to match `std::mt19937_64`. Verify this against two values, both checked with MSVC's `std::mt19937_64`:
- the 10000th output for the default seed 5489 is `9981545732273789042` (the C++ standard's own check);
- the first output for seed 42 is `13930160852258120406`.

**Order matters.** Draws happen in exactly the order given in A.3 and A.7. Fitness evaluation never draws random numbers.

## A.3 Random generation

`random_reaction(F, B, k_min, k_max)`, with the draws in this order:
1. `type = uniform_int(0, 3)`: 0 UniUni, 1 UniBi, 2 BiUni, 3 BiBi.
2. `n_react` = 2 for types 2 and 3, else 1. `n_prod` = 2 for types 1 and 3, else 1.
3. Draw `n_react` reactant species, then `n_prod` product species, each `uniform_int(0, F+B-1)`.
4. Make both sides canonical.
5. `k = log_uniform(k_min, k_max)`.

The result may be invalid.

`random_network`: start with `init = 1.0`×F and `boundary = 1.0`×B. For each of the `r_initial` reactions, draw `random_reaction` up to 1000 times until it's *acceptable*: reactants ≠ products, it involves a floating species, and it doesn't duplicate an existing reaction. If no draw succeeds, that's an error.

## A.4 ODE model

`ds/dt = N · v(s)` over floating species only, where `N[i][j]` = (products − reactants) coefficient of floating species i in reaction j.

The **analytic Jacobian** is `J = N · ∂v/∂s`. For `v_j = rate_factor_j · Π s_i^c_i`, the derivative is `∂v_j/∂s_l = rate_factor_j · c_l · s_l^(c_l−1) · Π_{i≠l} s_i^c_i`. Here `rate_factor_j` is `k_j` times the product of its boundary reactant concentrations raised to their coefficients.

## A.5 Simulation

`simulate(model, y0, sample_times, solver_options)`. Sample times are strictly increasing and start at 0. The state at t = 0 is recorded first. Two entry points:
- `simulate(t_end, n_points)`: uniform samples `t_i = t_end·i/(n_points−1)`, with the last one exactly `t_end`.
- `simulate_at(times)`: arbitrary times ≥ 0. If `times[0] > 0`, t = 0 is integrated from but not returned.

**Solver options and defaults:**
- `rtol = 1e-6`, `atol = 1e-9` (normal evaluation);
- `max_steps = 50000` accepted steps per simulation;
- `max_wall_seconds = 10`;
- `max_conc = 1e6`, `min_conc = −1e-8`.

**Guards:** checked after **every accepted step** and at every sample, not only at samples. Each failure has a status:
- **MaxConc:** a value is > `max_conc`.
- **Negative:** a value is < `min_conc`.
- **NonFinite:** a value is NaN or infinite.
- **StepBudget:** more than `max_steps` steps.
- **WallClock:** more than `max_wall_seconds`. Check the clock every 64 step attempts, through an injected time source (`performance.now()` in a worker, `Instant` natively).
- **StepTooSmall:** the solver can't make progress (step-size underflow, repeated convergence failures).

A failed simulation returns the samples up to the failure and its status.

**Solver requirement:** you need a **stiff** solver. This is the main technical risk of the port. Read §B.1 before choosing.

**Background (C++ findings, PROGRESS.md):**
- **RK45 alone isn't viable:** stiff candidates made it 20–85× slower.
- **The in-house Rosenbrock (RODAS3) is unreliable** at tight tolerances (negative concentrations, step-budget failures on valid networks).
- **CVODE (variable-order BDF with the analytic Jacobian)** became the default: about 2.5× faster than RK45 + Rosenbrock, and robust.

## A.6 Objectives

The fitness lies in [0, 1], and higher is better. The `objective` config key selects `oscillation` (default) or `timeseries`.

### A.6.1 Oscillation fitness

Simulate to `t_end` (default 100) at `n_points` (default 1000) uniform samples. For each floating species with trajectory `y`, compute the following.

1. **Window:** `start = floor(transient_fraction · len(y))` (default 0.5). `w = y[start..]`.
2. **Threshold:** `max(min_rel_amplitude · mean(|w|), 1e-6)`, with `min_rel_amplitude` default 0.01. The absolute floor 1e-6 stops solver noise in a decayed species from counting.
3. **Extrema** (`find_extrema(w, threshold)`), exactly as follows:
   - **Raw candidates:** for i in 1..len−2, a *peak* if `w[i−1] < w[i] >= w[i+1]`, a *trough* if `w[i−1] > w[i] <= w[i+1]`.
   - **Hysteresis:** `pending = first candidate`. For each next candidate `c`:
     - if `c` is the same kind as `pending`, replace `pending` when `c` is more extreme (higher peak or lower trough);
     - otherwise, if `|c − pending| > threshold`, commit `pending` to the output and set `pending = c`;
     - otherwise ignore `c`.
   - **Edges:**
     - Append `pending` only if the signal *after* it moves away by more than the threshold: some later `w[j]` with `pending − w[j] > thr` for a peak, `w[j] − pending > thr` for a trough.
     - Then drop the first committed extremum unless the signal *before* it moves away by more than the threshold, by the same test.
4. **Pairs:** `n = floor(len(extrema)/2)`. If n = 0, the score is 0.
5. **Damping:** `a_first = |e0 − e1|`, `a_last = |e(2n−2) − e(2n−1)|`, `d = clamp((a_first − a_last)/a_first, 0, 1)`. A growing amplitude gives d = 0, with no penalty (known weakness).
6. **Species score:** `min(n, n_target)/n_target · (1 − d)`, with `n_target` default 5.

**Network fitness** = the maximum species score, except in these cases:
- **Failed simulation:** the fitness is 0.
- **Concentration floor** (`osc_conc_floor`, default 1e-12, 0 = off): if *any* species with n ≥ 1 has `min(w) < osc_conc_floor`, the network's fitness is 0. This rule covers the whole network, because other species can oscillate cleanly *because of* the one crashing. Species that decay without oscillating don't trigger it.

### A.6.2 Time-series fitness

**Target file:** text lines with two numbers, time and concentration, separated by spaces, tabs or a comma.
- `#` starts a comment, and blank lines are skipped.
- At most one text header line is allowed, before the data. A line is a header if its first field doesn't parse as a number.
- Times must be ≥ 0, strictly increasing and finite, and the last time must be > 0.
- Values must be finite.
- Anything else is an error with a line number.

**Score:** simulate with `simulate_at(target times)`. For each floating species:

```
NRMSE = sqrt(mean_i (sim_i − target_i)^2) / mean_i |target_i|
```

If the mean |target| is 0, the denominator is 1 instead. The species score is `1/(1 + NRMSE)`. **Network fitness** = the best species score; a failed simulation scores 0. No extremum is required at the target times, only the values are scored.

### A.6.3 Tight-tolerance re-check (both objectives)

If a network's score is ≥ `recheck_threshold` (default 0.9), simulate it again with `rtol = recheck_rtol` (1e-10) and `atol = recheck_atol` (1e-14). The re-check uses the same solver method, `max_steps = 5,000,000` and `max_wall_seconds = 30`. **The fitness is the lower of the two scores.** A threshold > 1 turns the re-check off.

This exists because evolution exploits solver error. With the in-house solvers, about 5% of high-scoring C++ winners were integrator artifacts.

## A.7 Evolution

**Population:**
- `P = population` (200).
- `n_elite = round(elite_fraction·P)` (0.1). It's at least 1 if `elite_fraction > 0`, and at most `P−1`.
- `n_parents = clamp(round(parent_fraction·P), 1, P)` (0.3).

**Loop:**

```
rng = MT19937-64(seed)
pop = [random_network(rng) for _ in 0..P]            # generation 0
best = none
for gen = 0, 1, 2, ...:
    evaluate every individual not yet evaluated      # pure; may be parallel
    stable-sort pop by fitness, descending
    log stats(gen); if pop[0].fitness > best.fitness (strictly): best = pop[0] (record gen)
    if best.fitness >= target_fitness or gen >= generations: stop
    next = pop[0..n_elite]                           # elites keep their fitness
    while len(next) < P:
        parent = pop[uniform_int(0, n_parents−1)]
        next.push(make_child(parent))                # unevaluated
    pop = next
```

**`make_child`:** up to 10 attempts. Each attempt copies the parent, applies a mutation chosen by `choose_mutation`, and accepts the child if the mutation succeeded and the child is valid. If all 10 attempts fail, the child is a copy of the parent and is re-evaluated.

**`choose_mutation`:**
- `u = uniform01()·(p_rate + p_add + p_delete + p_rewire)`.
- Subtract the probabilities in that order. The first that takes `u` below 0 wins.
- If rounding leaves `u ≥ 0`, take the last kind with p > 0.
- Defaults: 0.7, 0.15, 0.10, 0.05.

**Mutations,** with draws in this order:

| Kind | Draws and rule | Fails if |
|---|---|---|
| Rate | `j = uniform_int(0, R−1)`; `k_j = clamp(k_j · exp(mut_sigma · normal()), k_min, k_max)` (`mut_sigma` 0.5) | R = 0 |
| Add | `random_reaction(...)`; append | R ≥ r_max, or not acceptable |
| Delete | `j = uniform_int(0, R−1)`; remove it | R ≤ r_min |
| Rewire | `j = uniform_int(0, R−1)`; expand reaction j into one slot per molecule, reactant slots first then product slots, in canonical order; `slot = uniform_int(0, slots−1)`; `r = uniform_int(0, F+B−2)`; if `r >= old species` then `r += 1` (a different species, uniformly); rebuild canonical sides, keep k | F+B < 2, or the result isn't acceptable against the *other* reactions |

**Generation stats** (`log.csv`: `generation,best,mean,median,best_reactions`):
- `median`: middle element, or the mean of the two middle elements, of the sorted fitnesses.
- `best_reactions`: reaction count of `pop[0]`.

**Multiple trials** (`numtrials = N`): a *master* MT19937-64 seeded once with `seed`. Trial i (1-based) runs the loop above with `seed = i-th master output`, so any trial can be rerun alone from its recorded seed.

## A.8 Long-run check (oscillation only)

This separates sustained oscillation from slow damping.
1. **Species:** take the species with the best score in the normal fitness simulation.
2. **Long simulation:** simulate to `5·t_end` with the same sample spacing: `intervals = round(5·(n_points−1))`, samples = intervals + 1. Scale the step budget and wall-clock limit by 5.
3. **Windows:** each is `round((1 − transient_fraction)·(n_points−1)) + 1` samples long. The first ends at sample `n_points−1` (t_end); the last ends at the final sample.
4. **Amplitude:** in each window, the mean `|e_i − e_(i−1)|` over consecutive extrema from `find_extrema`, with the same threshold rule as A.6.1. It's 0 if there are fewer than 2 extrema.
5. **Pass:** both amplitudes > 0, and `|a_last − a_ref|/a_ref ≤ 0.10`.

## A.9 Prune

1. Compute the fitness of the original network; `required = 0.95 ×` that fitness.
2. Repeat full passes until a pass removes nothing. In each pass, for i = 0, 1, …: try removing reaction i. Keep the removal if `fitness(candidate) >= required`, and don't advance i, because the next reaction moves into slot i. Otherwise advance i.
3. Never go below `r_min` reactions.

Prune uses the objective's full fitness, including the re-check.

## A.10 Antimony subset (import/export)

The writer produces exactly this form (the `const` line is omitted if B = 0):

```
// <optional header>
const X0, X1;
J0: X0 -> S1; k0*X0;
J1: S1 + S1 -> S0 + X1; k1*S1*S1;
k0 = 1.23; k1 = 0.45;
S0 = 1; S1 = 1;
X0 = 1; X1 = 1;
```

- **Sides:** each side lists one entry per molecule (`S1 + S1`, not `2 S1`), in canonical order. Generation and mutation never produce an empty side. If the parser accepts one, the writer writes it as empty, e.g. `J0: S0 -> ; k0*S0;`. Check the C++ parser's behaviour and match it.
- **Rate law:** `k<i>` followed by `*<species>` for each reactant molecule.
- **Numbers:** the shortest representation that round-trips. Note that Rust's `{}` never uses exponents (it prints `0.00000000000001`). Either match C++ `to_chars` (e.g. `1e-14`) or document the difference. The parser must accept both.

**The parser** accepts exactly:
- `//` comments;
- `const` lists;
- reactions `name: reactants -> products; rate;`, where coefficients like `2 S1` are allowed;
- numeric assignments.

It rebuilds mass action from the reactants and **checks that the written rate law matches**. A mismatch is an error. Unknown species, a missing `;`, duplicate reaction names and anything outside the subset are errors with a line number.

Exported files must load unchanged in Tellurium: `te.loada(text)`.

## A.11 Configuration

Configs are `key = value` lines; `#` starts a comment. An unknown key or a malformed value is an error with a line number. Every key with its default:

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
recheck_threshold = 0.9
recheck_rtol = 1e-10
recheck_atol = 1e-14
solver = <the port's stiff solver>     # C++: cvode | auto | rk45 | rosenbrock
max_steps = 50000
osc_conc_floor = 1e-12
objective = oscillation                # or timeseries
target_file =                          # web: an uploaded file instead
threads = 0                            # web: number of workers, 0 = hardwareConcurrency
```

**Validation:**
- **Population and generations:** `population ≥ 2`; `generations ≥ 0`.
- **Fractions:** `0 ≤ elite_fraction < 1`; `0 < parent_fraction ≤ 1`.
- **Network shape:** `n_floating ≥ 1`; `n_boundary ≥ 0`; `1 ≤ r_min ≤ r_initial ≤ r_max`.
- **Rate constants and mutation:** `0 < k_min < k_max`, finite; `mut_sigma ≥ 0`; all p ≥ 0 with a sum > 0.
- **Simulation:** `t_end > 0`; `n_points ≥ 3`; `0 ≤ transient_fraction < 1`; `n_target ≥ 1`.
- **Thresholds and limits:** `min_rel_amplitude ≥ 0`; `max_conc > 0`; `recheck_threshold ≥ 0`; recheck rtol and atol > 0; `max_steps ≥ 1`; `osc_conc_floor ≥ 0`.
- **Objective:** `timeseries` requires a target.

**Determinism:** within the port, a given seed and config must give bit-identical results whatever the worker count. Across languages, expect *statistically* equivalent behaviour, not identical runs. Even with identical RNG streams, `exp`, `ln` and `cos` may differ in the last bit between C++ and Rust/wasm math libraries, and different ODE solvers give different trajectories. §C.8 defines how to compare against C++.

---

# Part B: The Rust/WebAssembly application

## B.1 ODE solver in WebAssembly (decide in R3; the main risk)

SUNDIALS/CVODE is C. Linking it into a `wasm32-unknown-unknown` Rust build isn't practical, so the port needs a **pure-Rust stiff solver**. Options, in order of preference:

1. **A maintained pure-Rust BDF or SDIRK crate.** `diffsol` is the first candidate to evaluate. In a spike, confirm that it:
   - compiles to `wasm32-unknown-unknown`;
   - accepts a user-supplied dense Jacobian;
   - supports stepping one step at a time, or a stop time plus dense output, so the per-step guards (A.5) can be applied;
   - passes the M3-equivalent tests below.

   Record the result in `PROGRESS.md`.
2. **Port the C++ in-house solvers** (`src/integrator.cpp`: DOPRI5 + RODAS3, with an in-house LU solver). This is known to work, but known to be weaker. If you choose it, keep the re-check and the floor. Also investigate the RODAS3 tight-tolerance failures (the open issue in PROGRESS.md), for example by adding positivity handling, and consider a variable-order BDF later.
3. **SUNDIALS through Emscripten.** This works in principle but needs a different wasm target and toolchain. Avoid it unless 1 and 2 fail.

Whatever you choose:
- **Guards:** the per-step guards and statuses of A.5 must be applied.
- **Interpolation:** samples come from the solver's dense output, not from extra steps.
- **Validation:** the solver must pass §C.3 and the §C.8 comparison.

## B.2 Parallelism in the browser

GitHub Pages can't set the `Cross-Origin-Opener-Policy` and `Cross-Origin-Embedder-Policy` headers. So `SharedArrayBuffer` and wasm threads aren't available without workarounds (such as a `coi-serviceworker` shim). **Use a pool of Web Workers with message passing instead:**
- **Main-thread role:** the main thread, or one coordinator worker, runs the serial part of A.7: random draws, sorting and mutation.
- **Evaluation:** each generation, it sends the unevaluated networks to `W` evaluation workers in index-tagged batches. Each worker loads the same wasm module and returns `(index, fitness)`.
- **Order:** results are written by index, so the run is independent of `W` and of completion order. This is the browser version of the C++ M6 rule.
- **Responsiveness:** the UI thread must never run an evaluation, so the page stays responsive. Provide **Stop**, which finishes the current generation and keeps the results.

Per-network data sent to workers should be compact: a typed array or a small serialised struct, not Antimony text.

## B.3 User interface

A single page, no login, no server. **Inputs:**
- **Parameters:** a form for all config keys (A.11), with defaults, inline validation messages, and "load/save config" as a `.cfg` text file compatible with the C++ tool.
- **Objective:** oscillation or time series. For time series, upload or paste the two-column target, and show it as a plot.
- **Run settings:** seed, number of trials and number of workers.

**While running:**
- a progress line per generation (generation, best, mean, reactions, evaluations per second);
- a live plot of best and mean fitness against generation;
- for multiple trials, a table that fills in as trials finish (A.7 trials and the `summary.csv` columns).

**Results:**
- **Best network:** its Antimony text, with copy and download buttons.
- **Time course:** a plot of every species, with the target points overlaid for the time-series objective.
- **Fitness breakdown:** per-species scores (pairs and damping, or NRMSE), and the re-check result.
- **Oscillation checks:** the long-run check result, and **Prune** with before/after reaction counts and fitness.
- **Downloads:** `log.csv`, `best_final.ant`, `best_final.csv`, `run.cfg` and, for trials, `summary.csv`. These use the same formats as the C++ tool, so results can be compared and loaded in Tellurium. A zip of everything is optional.

**Also:**
- A **"Load model"** panel: paste or upload Antimony, then simulate and plot it, and score it (the C++ `simulate` and `fitness` commands).
- **Example presets:** the Brusselator model, the target-peaks time series, and the default config.

Keep the UI plain: HTML, CSS and a small amount of JS or TypeScript. Use a lightweight plotting library such as uPlot. Ask before adopting a framework.

## B.4 Deployment

- **Build:** `wasm-pack` (or `wasm-bindgen` CLI) plus a minimal static bundler, or none.
- **Publishing:** a GitHub Actions workflow builds the release wasm and publishes the static site to GitHub Pages on pushes to the main branch.
- **No secrets:** none are needed.
- **README:** state the URL, how to run locally (`python -m http.server` in the built folder), and the browser support (current Chrome, Firefox, Safari, Edge).

---

# Part C: Milestones and tests

## C.0 Fixtures

Copy these from the C++ repo's `examples/` folder into `fixtures/`:
- **`brusselator.ant`:** a genuine limit cycle (A = 1, B = 3).
- **`damped.ant`:** a stable focus (B = 1.5).
- **`integrator_artifact.ant`:** scores 1.0 with loose in-house solvers and 0 when solved accurately.
- **`deep_dip.ant`:** species dip to about 1e-22 each cycle; must score 0 under the floor.
- **`target_peaks.txt`:** Brusselator S0 troughs and peaks.
- **`timeseries.cfg`** and **`default.cfg`**.

## C.1 R0: Skeleton
- **Deliverables:** the workspace, the three crates, CI running `cargo test` and a wasm build, and `CLAUDE.md` and `PROGRESS.md`.
- **Tests:** a trivial test passes; the `core` crate builds for `wasm32-unknown-unknown`; the CLI prints `--version`.

## C.2 R1: Networks, RNG, Antimony, config
- **RNG:** MT19937-64 matches the standard check (A.2). A fixed seed gives fixed, recorded first values for every helper.
- **Generation:** 10,000 random networks are all valid, and a fixed seed gives identical networks.
- **Antimony:**
  - Write → parse → write is identical for 1,000 random networks.
  - Every fixture parses.
  - Errors with line numbers for: an unknown species, a missing `;`, a rate law that doesn't match.
- **Config:** a round trip; errors with line numbers; `default.cfg` equals the built-in defaults; validation rejects each rule in A.11.

## C.3 R2–R3: Model and solver
- **Model:**
  - The RHS of a hand-worked 3-reaction network matches to 1e-12.
  - The analytic Jacobian matches central differences to a relative error of 1e-6 on 100 random networks.
  - Boundary species aren't state variables.
  - The closed network A→B→A sums to 0.
- **Solver** (the §B.1 choice):
  - `dy/dt = −0.5y`: relative error < 1e-6 at all samples to t = 10, with rtol 1e-9 to 1e-10.
  - Robertson to t = 1e4: `y1+y2+y3 = 1` to 1e-6. At t = 1e4, y = (1.073004285e-1, 4.800166973e-7, 8.926990914e-1), each to a relative error of 1e-3 (SciPy reference).
  - Brusselator period over t ∈ [100, 200], from upward crossings of S0 = 1: within 2% of **7.1569199212**.
  - Lotka–Volterra conserved quantity: relative drift < 1e-5 over t = 50, with rtol 1e-10.
  - Guards: an explosive `S0 → S0 + S0` with k = 50 gives MaxConc in < 1 s; a step budget of 50 gives StepBudget with exactly 50 steps; `min_conc = 0.5` on the decay gives Negative.
  - `simulate_at` matches the exact decay at arbitrary times, with and without t = 0.
  - Determinism: identical results on repeated runs.

## C.4 R4: Fitness
- **Synthetic signals:**
  - sin(t) over 10 periods → n ≥ 5, d < 1e-3, fitness > 0.99.
  - `e^{−0.3t} sin t` → < 0.5.
  - A constant and a ramp → 0.
  - sin plus 0.1% noise → no double counting.
- **Floor:** a clean oscillator plus an oscillating species dipping to 1e-15 → 0. Plus a decaying species instead → no effect.
- **Fixtures, with default settings:**
  - Brusselator ≥ 0.95, and its re-check agrees to within 0.01.
  - damped < 0.5.
  - deep_dip = 0.
  - integrator_artifact < 0.5.
  - A blow-up network = 0.
- **Time series:**
  - Parsing (formats, comments, header, CRLF) and every error case.
  - The exact NRMSE formula on the decay with target = 2× exact.
  - Brusselator vs `target_peaks.txt` > 0.999 (S0 is the best match; C++ gives 0.99997).
  - damped vs `target_peaks.txt` < 0.7 (C++: 0.524).

## C.5 R5: Evolution
- Every mutation kind, applied 10,000 times, gives valid networks. Rate mutation stays in [k_min, k_max]. The chosen mutation kinds match their probabilities to within 1% over 100,000 draws.
- Elitism: best fitness never decreases.
- A toy fitness ("exactly 7 reactions") is solved within 50 generations.
- Determinism: the same seed gives a byte-identical `log.csv`, and a different seed gives a different best network.
- Trials: trial seeds come from the master generator; prefixes are stable; rerunning a trial alone reproduces it; `summary.csv` has one row per trial.
- Long-run check: the Brusselator passes; a slowly damped Brusselator (B = 1.999, near the fixed point) fails.
- Prune: it removes exactly the padding reactions from a padded Brusselator, removes nothing from the minimal one, and respects r_min.

## C.6 R6: Workers
- **Same results with any number of workers:** a run gives byte-identical `log.csv` with 1, 2 and 4 workers. Test this in a headless browser (e.g. `wasm-pack test --headless`), or natively against the same batching logic.
- **Stop:** Stop ends the run cleanly after the current generation.
- **Speed:** report evaluations per second against the C++ build's (about 22,000 per second with CVODE on 8 threads, on the PI's machine). There is no pass threshold.

## C.7 R7: User interface and deployment
- Every B.3 feature works in Chrome and Firefox. Downloaded files load in the C++ tool (`netevolve fitness --model`) and in Tellurium.
- The GitHub Pages deployment works from a clean push.

## C.8 R8: Validation against the C++ reference (acceptance)

Both use the default config and master seed 1.
- **Yield:** with 50 trials, the share of trials passing the long-run check should be statistically consistent with C++. C++ with CVODE gave 15/50 (13/50 with the floor); with the in-house solvers, 16/50. Use a two-proportion test at the 5% level and report it. Also report the share reaching fitness ≥ 0.95 (C++: 33/50).
- **Independent check of winners:** every winner with fitness ≥ 0.95 is re-scored by the native C++ build (`netevolve fitness --model best_final.ant`). At least 90% must also score ≥ 0.95 there. List the disagreements.
- **Time series:** 10 trials on `timeseries.cfg`. C++ gave best scores of 0.98, 0.91 and 0.88, with the rest at about 0.55 (a broad local optimum). The port should behave similarly; report the distribution.

If R8 fails, don't adjust the tests. Report the evidence and propose causes, for example solver differences or a deviation from A.2/A.3 in the RNG or mutation order.

---

## D. Known limitations carried over from the C++ version

- **Growing amplitudes aren't penalised** (d is clamped at 0). Only `max_conc` and the long-run check catch outward spirals.
- **Only the first and last pairs set the damping.** Irregular oscillation in between isn't penalised.
- **A short window:** `t_end = 100` with `n_target = 5` requires a period < 10. Slow oscillators, or ones still settling at t_end, can be misjudged.
- **Time series:** only the values at the target times are scored, the target species isn't named, and initial conditions aren't evolved.
- **Yield:** about 30% of trials give a sustained oscillator. In about 70% of time-series trials the best species gets stuck drifting slowly through the middle of the target range instead of reaching the peaks and troughs. More generations, a larger population or a higher `r_max` probably help; this is untested.

## E. Out of scope

Crossover; evolving initial concentrations; reversible or non-mass-action kinetics; SBML; user accounts or server storage; mobile-specific layouts beyond being usable.
