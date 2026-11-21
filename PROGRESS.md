# Progress

## Current milestone

**M7 (end-to-end oscillator evolution and validation) is complete. All milestones M0–M7 are done.** One step remains for the human: loading the results into Tellurium manually (SPEC M7). See "Manual Tellurium step" below.

Remote: https://github.com/zyrrron/netevolve (private), branch `main`.

## After M7: multiple trials and the tight-tolerance re-check (2026-11-20)

Both were added at the PI's request and go beyond SPEC.

**`evolve --numtrials N`** (`run_trials_command`, `trial_seeds` in `src/evolve.cpp`):
- Runs N independent trials in one command.
- `--seed` seeds one master `mt19937_64` once, and trial i uses its i-th output as its seed.
- Trial i writes to `<out>/trial<i>/`, and its `run.cfg` holds its own seed. Rerunning a trial alone with that seed gives byte-identical output (tested).
- `<out>/summary.csv` has one row per trial.
- Without the flag, `evolve` is unchanged: seed 1 still stops at generation 72 after 13,160 evaluations.

**Why the re-check was needed:**
- **Test batch:** a 50-trial batch, `--seed 1 --numtrials 50` (default config, before the re-check), gave 37/50 trials with fitness ≥ 0.95 and 17/50 passing the long-run check.
- **Re-scoring at tight tolerance:** the 37 winners were re-scored at rtol 1e-10 / atol 1e-14 with three solver setups: Auto, Rosenbrock only and RK45 only.
  - 28 were robust.
  - 2 were integrator artifacts. Trial 1 dropped from 1.0 to 0.33. Trial 18 dropped from 1.0 to 0 and went to negative concentrations; it had also **passed** the long-run check.
  - 3 had inflated scores, 1.0 dropping to 0.8.
  - The remaining 4 differed only in the Rosenbrock-only solve or in amplitude.
- **CVODE (Tellurium 2.2.11.1, libRoadRunner 2.9.1):** it failed to integrate 8 of the 9 flagged networks, at default and/or tight settings, mostly with `CV_CONV_FAILURE` (h down to 1e-156). It reproduced all 4 robust controls, matching NetEvolve's amplitudes to 4–5 significant figures. Conclusion: the flagged models are numerically bad, and the default tolerance let evolution select them.

**Tight-tolerance re-check** (`oscillation_fitness` in `src/fitness.cpp`; config keys `recheck_threshold = 0.9`, `recheck_rtol = 1e-10`, `recheck_atol = 1e-14`):
- **Rule:** a network scoring ≥ the threshold is re-simulated at the tight tolerances, and its fitness is the lower of the two scores. A failed re-check simulation scores 0. A threshold > 1 turns the re-check off.
- **Where it applies:** evolution, `prune` and the acceptance test. The `fitness` command prints the re-check line.
- **Fixture:** `examples/integrator_artifact.ant` (trial 18 above) scores 1 at the default tolerance and 0 at the tight one.
- **Tests:** 3 new doctest cases plus config validation checks. No existing test was changed.

**Results with the re-check:**
- **Acceptance test:** **4 of 5 passed**, up from 3 of 5. Seed 2 now evolves a different network, which passes (213.6 → 213.6).
- **Same 50-trial batch:**
  - 37/50 trials had fitness ≥ 0.95, and 16/50 passed the long-run check.
  - Run time was 257.5 s, against 240.6 s before (+7%).
  - All 37 winners keep their fitness at tight tolerance with both Auto and RK45.
  - The yield before the re-check, counting only winners that also passed the tolerance check, was 14/50.

**CVODE is now the default solver** (PI decision, 2026-11-20):
- **Why:** RK45 alone is not viable. In a 50-trial batch with `solver = rk45` and `max_steps = 500000`, trials ran 20–85× slower than `auto` (projected 2.5–3 h against 4 min). The batch was stopped after 3 trials, which gave the same results as `auto`. Stiff candidates make RK45 take hundreds of thousands of steps.
- **Implementation** (`run_cvode` in `src/integrator.cpp`):
  - SUNDIALS 5.7.0, CVODE BDF with the analytic Jacobian, `SUNDenseMatrix` and `SUNLinSol_Dense`.
  - `CV_ONE_STEP` with a stop time at t_end, so the max_conc, min_conc and non-finite guards, the step budget and the wall-clock limit apply after every step, as for the in-house solvers.
  - Samples come from `CVodeGetDky`. A non-finite RHS returns a recoverable error.
  - CVODE errors map to `StepBudget` (`CV_TOO_MUCH_WORK`), `NonFinite` (RHS failures) or `StepTooSmall` (others).
- **Build:**
  - `NETEVOLVE_WITH_CVODE` is ON by default. If SUNDIALS isn't found, configuration stops with instructions; `-DNETEVOLVE_WITH_CVODE=OFF` builds without it.
  - Config key `solver = cvode | auto | rk45 | rosenbrock`. The default is `cvode`, or `auto` in builds without CVODE, so `default.cfg` leaves it commented. Also `max_steps`.
  - `--version` lists the solvers available.
  - On the PI's machine: `-DSUNDIALS_DIR=D:/buildroadrunner/libroadrunner-deps/install-Release/lib/cmake/sundials`. That install is static and `/MD`, matching our Release build.
- **Tests:**
  - `tests/test_cvode.cpp`: decay against the exact solution (7e-9), Robertson against the SciPy reference, Brusselator against RK45 (4e-6), guards, determinism, and the artifact scoring 0.
  - Without CVODE, it checks that requesting CVODE fails cleanly.
  - Test changes: (1) `integrator_artifact.ant` test now pins `solver = auto`, because the artifact fools the in-house solvers and CVODE rejects it even at the default tolerance; (2) my own `solver` validation test used "cvode" as an unknown name and now uses "lsoda"; (3) the default-solver check expects CVODE when available.
  - Results: 79/79 doctest cases with CVODE and 74/74 without; CTest 12/12 in both.
- **Acceptance with CVODE:** 4 of 5 passed. Seed 5 fails the long-run check (0.43 → 0, a damped transient).
- **Same 50-trial batch with CVODE:**
  - 110 s against 258 s for `auto`; 2.46 M evaluations, about 22,400 per second against about 8,000.
  - 33/50 trials had fitness ≥ 0.95, and 15/50 passed the long-run check (16 for `auto`; no real difference at n = 50).
  - 32 of the 33 winners were confirmed by the in-house tight solvers. The exception, trial 19, scores 1.0 with CVODE and with Rosenbrock, 0.92 with SciPy Radau, and fails with RK45 (negative) and with LSODA. Its S0 and S1 drop to about 1e-22 every cycle: numerically pathological and physically meaningless. It fails the long-run check anyway.
  - Proposed, not done: an optional concentration floor (e.g. reject oscillators whose species dip below 1e-12).

**Concentration floor** (`osc_conc_floor = 1e-12`, PI request; `oscillation_report` in `src/fitness.cpp`):
- **Rule:** if a species with at least one peak-trough pair in the scored window drops below the floor there, the network's fitness is 0.
- **Why the whole network:** in deep_dip.ant, S3 oscillates cleanly *because* S0 and S1 crash, so a per-species rule would not have rejected it.
- **What doesn't trigger it:** species that decay without oscillating.
- **Fixtures and tests:**
  - `examples/deep_dip.ant` is CVODE trial 19. Its minimum is about 1e-22 (Radau); with CVODE at the default tolerance, S0 reaches -3.2e-10.
  - Tests: synthetic signals, the config key, the Brusselator unaffected, deep_dip rejected with CVODE (and passing with the floor off), and integrator_artifact.ant rejected by the floor alone.
  - The earlier re-check test now sets `osc_conc_floor = 0` to isolate the re-check, because the floor catches that network first.
- **Acceptance** (CVODE + floor): **3 of 5** passed, down from 4 of 5. Seeds 2 and 5 stall at fitness 0.2 and 0.9.
- **50-trial batch** (CVODE + floor, master seed 1): 109 s. 33/50 trials had fitness ≥ 0.95 (same as without the floor), and **13/50** passed the long-run check, against 15/50 without the floor. At n = 50 the binomial SD is about 3.2, so this isn't a clear loss, but the floor may block some evolutionary paths that pass through deep-dip networks. Worth watching with more trials.

**Time-series objective** (PI request; `parse_target`, `evaluate_target`, `network_fitness` in `src/fitness.cpp`, and `simulate_at` in `src/integrator.cpp`):
- **Config:** `objective = oscillation | timeseries` and `target_file`.
  - The target file has two columns, time and concentration, with whitespace or comma separators, `#` comments and one optional header line. Times must be strictly increasing.
  - A relative `target_file` resolves against the config file's directory and is stored absolute, so `run.cfg` works from the output directory.
- **Score:** simulate from t = 0 and sample at exactly the target times (`simulate_at`). Per species, NRMSE = RMS error / mean |target|, and score = 1 / (1 + NRMSE). The network fitness is the best species.
  - The re-check applies, and so do `prune` and `--numtrials`.
  - The long-run check is n/a for this objective. The trial summary shows "n/a", and the batch line reports how many trials reached target_fitness.
- **Ambiguities resolved simply:**
  - The file names no species, so the best match is used, as in the oscillation fitness.
  - Only the values at the target times are scored; no extremum is required there.
  - Initial conditions are those of the network (all 1 for evolved networks).
- **Example:** `examples/target_peaks.txt` holds the Brusselator S0 troughs and peaks on t ∈ [0, 40], and `examples/timeseries.cfg` uses it. The Brusselator scores 0.99997 against it (S0), and damped.ant scores 0.52.
- **Tests:** `tests/test_timeseries.cpp` covers file formats and errors, `simulate_at` for all in-house solvers, the exact NRMSE formula, the Brusselator target, config validation and relative paths, and evolve plus `--numtrials` end to end.
- **10 trials, master seed 1, default population and generations:** 60 s, and **0/10 reached 0.99**.
  - Best: trial 6 at 0.980 (NRMSE 0.020). Its S3 matches every target to within about 0.07, a genuine oscillator found from scratch.
  - Trials 2 and 9 reached 0.91 and 0.88.
  - The other 7 are stuck at about 0.54–0.57 (NRMSE ≈ 0.83). That is the score of a slowly drifting line near the target's mean, a broad local optimum.
  - Every trial was still improving at generation ~500 and was at or near r_max = 12 reactions. More generations, a larger population or a higher r_max would probably help.
- **Test counts:** 90/90 doctest cases with CVODE and 84/84 without; CTest 12/12 in both.

**SUNDIALS is downloaded automatically** (PI request; `CMakeLists.txt`):
- **Lookup order:** `SUNDIALS_DIR` / `find_package`; otherwise FetchContent of the official v5.7.0 release tarball with a pinned SHA-256 (`NETEVOLVE_FETCH_SUNDIALS=ON` by default); otherwise a configure error with instructions.
- **What gets built:** only CVODE, static, with no examples. SUNDIALS is resolved *before* `add_compile_options(/W4 /WX ...)`, so its sources aren't built with our warnings-as-errors.
- **Short folders:** source and build go in `<build>/sd-src` and `<build>/sd`, because the default `_deps/sundials-build` plus MSBuild's .tlog paths exceeded 260 characters (FTK1011) in a deep build directory.
- **Measured:** a fresh build directory takes 15 s to configure (including the download) and 34 s to build. All 12 CTest tests and all 90 doctest cases pass.
- **Identical results:** seed 3 gives byte-identical `log.csv`, `best_final.ant` and `best_final.csv` with the downloaded SUNDIALS and with the PI's installed SUNDIALS 5.7.0.
- **Offline:** `-DFETCHCONTENT_SOURCE_DIR_SUNDIALS=<extracted folder>` uses a local source and downloads nothing (verified).
- `.gitignore` now ignores `build_*/` as well as `build/`.

**Open issue: RODAS3 at tight tolerances** (it affects only `solver = auto` / `rosenbrock` now).
- In the re-check batch, Rosenbrock alone scores 0 on 5 winners (trials 16, 17, 21, 22, 39) where Auto and RK45 agree on about 1.0. It gives negative concentrations on 16 and 39 and exceeds the step budget on 21.
- On trial 21 (first batch), CVODE agreed with Auto/RK45. So this looks like a weakness of the Rosenbrock solver (no positivity handling?), not of the networks.
- Auto switches to Rosenbrock on stiffness, so this could wrongly reject stiff genuine oscillators. It needs investigating.

## M7: End-to-end oscillator evolution and validation

Done:
- **`prune`** (`src/evolve.cpp`, CLI `prune --model F [--config F] [--out F]`):
  - Tries removing each reaction in order and keeps a removal if the fitness stays ≥ 0.95 × the *original* network's fitness. It repeats full passes until nothing more can be removed, and never goes below `r_min`.
  - Writes the pruned model with a header comment and reports how many reactions were removed, the fitness before and after, and the long-run check.
- **Long-run check** (`long_run_check` in `src/fitness.cpp`):
  - Re-simulates to 5 × t_end with the same sample spacing.
  - Measures the amplitude of the best-scoring species as the mean |peak − trough| of the filtered extrema (`window_amplitude`, same threshold as the fitness). It uses two windows of length (1 − transient_fraction) × t_end: one ending at t_end, one ending at 5 × t_end.
  - Passes if both amplitudes are > 0 and differ by at most 10%.
  - The step budget and wall-clock limit scale with the simulated time (× 5).
  - `evolve` and `prune` both print the result.
- **Acceptance test** (`tests/acceptance.cpp`, executable `netevolve_acceptance`, CTest `acceptance_seeds_1_to_5`, label `slow`):
  - Skipped (exit code 77, CTest "Skipped") unless `NETEVOLVE_RUN_SLOW=1`, so it isn't part of the default run.
  - Runs the default config for seeds 1–5, prunes each best network, runs the long-run check, writes `build/acceptance_out/seedN/` (including `pruned.ant`) and prints a table.

Tests (Release): **9/9 CTest tests passed in the default run, with the acceptance test reported as Skipped. Doctest: 67/67 cases and 1,575,605 assertions passed.**
- **Prune:**
  - Removes exactly the two padding reactions (negligible k) from a 6-reaction Brusselator, leaving the 4 core reactions in order.
  - Removes nothing from the minimal Brusselator.
  - Respects `r_min` (1, 3 and 4 tested).
  - With a toy fitness, the 0.95 threshold is relative to the original network's fitness: k values 50, 46, 3, 1 → removes the 3 and the 1, keeping 0.96.
- **Long-run check:**
  - Brusselator passes (−0.04%).
  - **A slowly damped Brusselator** (B = 1.999, just past the Hopf point, starting near the fixed point) **has fitness 0.972 but fails the long-run check (−24.8%)**. This is exactly the slow damping the check exists to catch.
  - The stable focus and a blow-up network both fail.
- **`window_amplitude`:** measures 4.0 for 5 + 2 sin.
- **CLI (CTest):** `prune` on the Brusselator removes 0 of 4 reactions and passes the long-run check.

### Acceptance test result (recorded 2026-11-18)

`NETEVOLVE_RUN_SLOW=1 ctest -C Release -L slow`, default config, 8 threads, 12 s in total: **PASSED, 3 of 5 runs** (3 needed).

| seed | best fitness | at generation | evaluations | pruned fitness | reactions | long-run amplitude (best species) | result |
|---|---|---|---|---|---|---|---|
| 1 | 1.0000 | 72 | 13,160 | 1.0000 | 7 → 5 | 1.574 → 0 (−100%) | FAIL (long-run) |
| 2 | 1.0000 | 112 | 20,360 | 1.0000 | 7 → 6 | 0.4819 → 0.3811 (−20.9%) | FAIL (long-run) |
| 3 | 1.0000 | 58 | 10,640 | 1.0000 | 7 → 7 | 491.7 → 484.7 (−1.4%) | PASS |
| 4 | 1.0000 | 93 | 16,940 | 1.0000 | 9 → 9 | 104.4 → 104.4 (−0.0%) | PASS |
| 5 | 0.9943 | 30 | 5,600 | 0.9974 | 11 → 7 | 0.1541 → 0.1540 (−0.1%) | PASS |

Every run reached fitness ≥ 0.95, within 30–112 generations; runs stop at `target_fitness = 0.99`. What each long-run result means:
- **Seed 1, correctly rejected.** S0 and S2 grow without bound (S2 is about 297 at t = 500 and still rising). The growing S2 damps the oscillating species through S1 + S2 → S2, so the oscillation is a transient, not sustained.
- **Seed 2, a real limit cycle rejected by the check as specified.** At t = 100 the trajectory is still settling onto the cycle from above: S1's peak goes 1.3 → 0.96 → 0.945, then stays constant from t ≈ 200 to 500. The check compares the t = 50–100 window with the last window, which gives −21%. This is a false negative caused by t_end being short relative to the convergence time. Neither the check nor the threshold was changed.
- **Seeds 3, 4 and 5, sustained oscillators, confirmed independently:**
  - **SciPy** `solve_ivp`, Radau and LSODA (rtol 1e-10, atol 1e-13), on equations built directly from `pruned.ant`, gives the same value ranges as NetEvolve in every window from t = 50 to 500, and constant amplitude from t ≈ 200.
  - **Tellurium 2.2.13.1** (libRoadRunner, CVODE) with the spec's call `te.loada(...).simulate(0, 500, 5000)` shows sustained oscillation for seeds 3 and 4, with the same amplitude at t = 450–500 as at t = 50–100. For seed 5, CVODE stops with `CV_TOO_MUCH_WORK` at t ≈ 242 with default settings, and still stops even with `maximum_num_steps = 1e6`. With `absolute_tolerance = 1e-14` and `relative_tolerance = 1e-10` it completes and shows sustained oscillation. Those same tight settings make CVODE stall on seed 3 at t = 16.6 (step size about 1e-15), while seed 3 runs fine with the defaults. CVODE is fragile on these spiky relaxation oscillators, whose species dip to about 1e-6 on every cycle. NetEvolve and SciPy handle all of them.

Defect found and fixed during the acceptance run: the long-run check first used the normal per-simulation step budget (50,000 steps) for a simulation five times as long. Seeds 3 and 4 hit it (seed 3 at t = 358) and failed with an amplitude of 0, although their oscillation is genuine. The budget and wall-clock limit now scale with the simulated duration. With that fix the result went from 1 of 5 to 3 of 5. No test or threshold was changed. The acceptance table now prints the simulation status when a long run fails, so this can't hide behind "0 → 0" again.

**The margin is thin: exactly 3 of 5.** If you want more headroom without changing the test, the options are parameter changes, as SPEC M7 suggests:
- a longer `t_end`, e.g. 200, so the fitness window starts after the approach to the limit cycle (this would likely turn seed 2 into a pass);
- a `target_fitness` closer to 1, so runs don't stop at the first network scoring 0.99;
- adding the long-run check to the selection itself, instead of only checking afterwards.

None of these is applied.

### Manual Tellurium step (for the human, SPEC M7)

```python
import tellurium as te
r = te.loada(open(r'build/acceptance_out/seed3/best_final.ant').read())   # also seed4
r.simulate(0, 500, 5000); r.plot()
```

For seed 5, set `r.integrator.absolute_tolerance = 1e-14` and `r.integrator.relative_tolerance = 1e-10` before `simulate`, or CVODE stops with `CV_TOO_MUCH_WORK`.

## M6: Parallel evaluation and performance

Done:
- `evaluate_pending()` in `src/evolve.cpp` evaluates the individuals with no fitness yet on `resolve_threads(cfg.threads)` `std::thread`s (0 means `hardware_concurrency`).
  - Workers take indices from an atomic counter, and each writes only its own individual's fitness. The fitness function is pure: each call builds its own `Model` and `SimResult`.
  - Every random draw (initial population, parent choice, mutation) stays in the serial part of the loop, so results can't depend on the thread count or on scheduling.
  - An exception thrown in a worker is captured and rethrown in the caller.
  - Threads are created per generation. That overhead is negligible next to about 180 evaluations of roughly 0.7 ms each.
- `EvolutionResult::n_evaluations` counts fitness evaluations; cached elites aren't re-evaluated. `evolve` prints the number of evaluations, the elapsed time and the evaluations per second.

Tests (Release): **8/8 CTest tests passed. Doctest: 59/59 cases and 1,575,573 assertions passed.** The M6 cases passed on 3 repeated runs.
- `log.csv`, `best_final.ant` and `best_final.csv` are byte-identical for `threads = 1`, `8` and `0` (all cores). The test uses a real oscillation run (seed 27, population 100, 4 generations) that reaches a non-zero best fitness, so the logs aren't trivially all zeros.
- `evolve()` with a varying toy fitness gives identical logs and best networks for 1, 2, 3, 8 and 17 threads. The evaluation count is exact: 120 + 30 × 108.
- An exception in a worker thread reaches the caller.

**Evaluations per second** (information only; default config, 40 generations, 7,400 evaluations, Release, 8 hardware threads):

| threads | evals/s | speed-up |
|---|---|---|
| 1 | 1,455 | 1.0× |
| 2 | 2,729 | 1.9× |
| 4 | 5,049 | 3.5× |
| 8 | 9,247 | 6.4× |

**Preview of a full default run** (seed 42, 8 threads), for information ahead of M7. The best fitness stayed at 0 for generations 0–135, then climbed to **1.0 by generation 140**, and the run stopped at the target. It took 25,400 evaluations in 3.2 s. The best network is an 8-reaction relaxation oscillator; for example, S1 swings between 0.18 and 256. Re-simulated to t = 500, the amplitude in [450, 500] matches [50, 100]. So the oscillation is genuine, not an artifact, and the default search can find oscillators despite the flat start.

## M5: Mutation and the evolutionary loop

Done (`src/evolve.hpp/.cpp`, CLI `evolve`):
- **Mutation operators:**
  - Rate: k × exp(σ·N(0,1)), clamped to [k_min, k_max].
  - Add: a random reaction with log-uniform k, refused at r_max.
  - Delete: refused at r_min.
  - Rewire: one molecule slot moves to a different species, keeping the reaction type.
  - `choose_mutation` picks the operator with probabilities p_rate : p_add : p_delete : p_rewire.
  - `make_child` makes up to 10 attempts, choosing the kind afresh each time, then falls back to copying the parent.
- **Generation loop (`evolve`):**
  - Generation 0 is the random initial population, followed by `generations` more, stopping early at `target_fitness`.
  - Each generation: evaluate, stable-sort (best first), log, record any improvement, then build the next generation from the elites plus children of parents drawn uniformly from the top `parent_fraction`.
  - Elites keep their cached fitness, which is safe because evaluation is pure.
  - Fitness is passed in as a function, so the toy test runs the loop without the ODE machinery.
- **`run_evolve_command`:** writes `run.cfg`, `log.csv` (flushed every generation), `best_genNNNN.ant` on every improvement, `best_final.ant` and `best_final.csv`. CLI: `evolve --config FILE [--seed N] [--out DIR]`.

Tests (Release): **8/8 CTest tests passed. Doctest: 56/56 cases and 1,575,544 assertions passed.**
- Each of the four operators, applied 10,000 times to random networks, always gives a valid network with k in range. A refused mutation leaves the network unchanged. Applied rates: Rate 100%, Add 76%, Delete 90%, Rewire 85%.
- `make_child` gives valid children 10,000 times out of 10,000, and more than 99% of them differ from their parent.
- Operator details: Add and Delete change R by ±1, and Rewire keeps every reaction's type. Add is refused at r_max and Delete at r_min.
- Rate mutation with σ = 3 for 10,000 steps keeps k in [k_min, k_max], reaching both bounds.
- Mutation kinds are chosen in the configured proportions (200,000 draws).
- **Toy fitness** (exactly 7 reactions, starting from 3): solved at generation 5, well within the limit of 50.
- **Elitism:** best fitness never decreases, both with the oscillation fitness and with a bumpy toy fitness over 60 generations.
- **Determinism:**
  - Two `evolve` runs with the same seed give byte-identical `log.csv`, `best_final.ant` and `best_final.csv`. A different seed gives a different best network.
  - With a fitness that varies, log rows are identical for the same seed and differ for another seed.
- **Output files:** `log.csv` has the header plus one row per generation. `run.cfg` re-reads to the same config. `best_final.ant` loads back as the best network. `best_gen0000.ant` and the best generation's file both exist.
- **CLI (CTest):** a tiny `evolve` run finishes and reports the best fitness.

**Two defects found while checking M5 results.** Both made random networks score a perfect 1.0 without oscillating:
1. **Integrator noise near zero.** A species with no production decays to about 0 and then jitters at the 1e-9 level. A threshold relative to mean|y| (~1e-11) counted that as a sustained oscillation. Fix: an absolute floor on the peak-trough difference, `min_abs_amplitude = 1e-6`, applied as `max(min_rel_amplitude × mean|y|, 1e-6)`. Regression test in `test_fitness.cpp`.
2. **RODAS3 dense output overshoot.** For a very stiff, quasi-steady species (relaxation rate ~4e4), cubic Hermite interpolation between large Rosenbrock steps overshot by about 1e-6, creating fake oscillations. Fix: RODAS3 no longer interpolates; it shortens steps to land exactly on sample times, and the controller keeps its preferred step size. Regression test in `test_integrator.cpp`: S3 stays within 1% of its quasi-steady value and decreases monotonically. Every M3 test still passes. RODAS3 step counts went up: Robertson 422 → 1,326, Brusselator (t = 200, 20,001 samples) 11,610 → 23,635. Both are within the budget of 50,000.

**Fitness landscape, a risk for M7.** With both fixes in, random networks rarely score above 0. Among 20,000 random networks (default shape), the mean fitness is 4e-5 to 1e-4, and the best are 0.25 and 0.34: genuine damped oscillators. A 30-generation run with population 200 (seed 1) stayed at 0 throughout. The spec's fitness gives no gradient until some species shows at least one large peak-trough pair after t = 50. M7 will tell whether the default 200 × 500 search finds oscillators. If not, the spec says to report the curves and propose parameter changes, not to tune the tests.

**Speed:** about 1,260 fitness evaluations per second on one thread (20,000 networks in 15.8 s).

## M4: Fitness function

Done:
- `src/fitness.hpp/.cpp`:
  - `find_extrema` detects peaks and troughs with an amplitude filter.
  - `signal_score` scores one species as in SPEC §1.4, after discarding the transient.
  - `oscillation_report`, `evaluate_oscillation` and `oscillation_fitness` compute the network fitness: the best species score, or 0 if the simulation fails.
  - `fitness_params(Config)` maps config keys to fitness settings.
- `src/config.hpp/.cpp`: the `key = value` config parser from SPEC §2, with all keys. It handles errors with line numbers, `validate_config`, and `write_config` (for `run.cfg` in M5). `examples/default.cfg` is the spec's config.
- CLI: `simulate` (CSV to stdout or `--csv`) and `fitness` (per-species breakdown plus the fitness value). `parse_flags` rejects unknown options, missing values and repeated options. Exit codes are 0 for success, 1 for runtime errors and 2 for usage errors.
- `examples/damped.ant`: the Brusselator with B = 1.5, a stable focus with eigenvalues −0.25 ± 0.968i.
- `README.md`: building, usage, the fitness definition, the Antimony subset and Tellurium.

Tests (Release): **7/7 CTest tests passed. Doctest: 44/44 cases and 777,521 assertions passed.**
- **Synthetic signals:**
  - sin(t) over 10 periods: n ≥ 5, d < 1e-3, score > 0.99.
  - e^{−0.3t} sin(t): score 0.003, below 0.5.
  - A constant signal and increasing and decreasing ramps all score 0.
  - sin(t) with ±0.1% uniform noise: the same pair count as the clean signal across 20 trials.
  - The same at 20,000 samples, where noise creates more than 1,100 raw extrema: the filter keeps exactly the 10 real ones.
- **Models:**
  - `brusselator.ant` scores 0.9996 (S0: 7 pairs, d = 4e-4).
  - `damped.ant` scores 0.
  - A blow-up network stops with `MaxConc` and scores exactly 0.
- **Extrema unit tests:** alternation, merging of ripples, flat tops, and ripple at the window edges.
- **Config:** `default.cfg` matches the defaults. Also covered: comments, round trip, and errors for an unknown key, a bad value or a missing `=`, each with its line number.
- **CLI (CTest):** `simulate` prints the expected CSV header and first row; `fitness` reports ≥ 0.95 for the Brusselator; an unknown option and a bad model file both fail.

Defect found and fixed during M4: `find_extrema` committed the last pending extremum unconditionally. With noise, a ripple at the end of the window, where the signal was still rising, counted as a peak. The dense-noise test caught it. An edge extremum now counts only if the signal really moves past the threshold beyond it (after it, for the last one; before it, for the first).

## M3: Integrators

Done (`src/integrator.hpp/.cpp`):
- **RK45:** Dormand–Prince 5(4), as in Hairer's DOPRI5. It uses FSAL, a weighted RMS error norm with `rtol`/`atol`, the 4th-order DOPRI5 continuous extension for samples between steps, and the DOPRI5 stiffness test (15 hits of h·|λ| > 3.25).
- **Stiff solver:** RODAS3 (Sandu et al. 1997, KPP formulation). It's a 4-stage Rosenbrock method of order 3, with an embedded order-2 error estimate, and it's L-stable and stiffly accurate. It uses the analytic Jacobian and an in-house LU with partial pivoting (`lu_factor` / `lu_solve`), and cubic Hermite interpolation for samples.
- **Guards:** a step budget (`max_steps`, per solver), a wall-clock limit (`max_wall_seconds`, checked every 64 attempted steps), `max_conc`, `min_conc` (−1e-8), non-finite values, and step-size underflow. They are checked at every accepted step and every sample.
- **Auto policy:** run RK45 first. If it reports stiffness, the step budget or step underflow, carry on from the last accepted point with RODAS3. Samples already taken are kept.
- **Sample times:** `simulate()` samples at t_i = t_end·i/(n_points−1), including both ends. On failure it returns the samples up to that point together with a `SimStatus`.

Tests (Release): **3/3 CTest tests passed. Doctest: 27/27 cases and 777,386 assertions passed.**
- **LU:** 200 random systems with n ≤ 20 have residuals below 1e-9, and singular matrices are detected.
- **dy/dt = −ky at t = 10:** relative error 1.1e-9 (RK45) and 1.8e-9 (RODAS3), against a limit of 1e-6. Interior samples are just as accurate.
- **Lotka–Volterra over t = 50:** the conserved quantity drifts by at most 1.5e-10 relative (limit 1e-5), for both RK45 and Auto.
- **Robertson to t = 1e4 with RODAS3:** 422 accepted steps (budget 50,000). y1+y2+y3 stays within 1.3e-15 of 1 (limit 1e-6). The final state matches the SciPy reference to within 2.4e-6 relative. The Auto policy switches to RODAS3 on this problem.
- **Brusselator period:** 7.15692 with every solver, against the reference 7.1569199212, a relative error of about 2e-7 (limit 2%).
- **Explosive `S0 → 2 S0` (k = 50):** every method stops with `MaxConc` in under a second.
- **Guards:** the step budget, wall clock and concentration floor each stop a simulation with the right status.
- **Determinism:** two identical runs give bit-identical samples.

Bug found and fixed during M3: `lu_solve` applied the row swaps interleaved with forward elimination, but `lu_factor` swaps whole rows, L multipliers included. The new LU test caught it with residuals around 7. The solver tests hadn't, because 2–3-species systems rarely pivot in a way that exposes it.

## M2: ODE model and Jacobian

Done:
- `src/model.hpp/.cpp`:
  - `make_model(Network)` builds the F×R stoichiometry matrix **N**, restricted to floating species. It also folds each reaction's `k` and its boundary-species factors into one constant, `rate_factor`.
  - `rates`, `rhs` (ds/dt = N v) and `jacobian` (J = N ∂v/∂s, analytic, F×F row-major).
  - Every function supports general integer stoichiometry, which the Brusselator's trimolecular step needs.
- Warnings are now errors (`/WX`), as the spec requires from M2.

Tests (Release): **3/3 CTest tests passed. Doctest: 18/18 cases and 776,937 assertions passed.**
- The RHS of a 3-reaction network matches values worked out by hand to within 1e-12. Rates and N are checked too.
- The Brusselator's RHS and Jacobian match the textbook formulas to within 1e-12.
- The analytic Jacobian matches a central finite-difference Jacobian for 100 random networks at random positive states. The worst relative error was 1.5e-9 (required: < 1e-6). A separate case checks cubic rate laws.
- Boundary species aren't state variables: the state, N's rows and the Jacobian are all F-sized, X enters only as a constant, and boundary concentrations never change.
- Mass conservation: for closed A → B → A, the RHS sums to 0 within 1e-12 at 1,000 random states.
- As a sanity check, dropping the stoichiometric coefficient from the Jacobian on purpose made 3 Jacobian tests fail. The change was reverted.

## M1: Network representation, RNG, Antimony I/O

Done:
- `src/rng.hpp`: `ne::Rng` wraps `std::mt19937_64`. The uniform, integer, normal and log-uniform helpers are hand-written instead of `<random>` distributions, whose output differs between standard libraries. Results are therefore reproducible across compilers, not just on one build.
- `src/network.hpp/.cpp`:
  - The `Network`, `Reaction` and `Term` structs.
  - `check_validity`, `random_reaction` and `random_network`.
  - `write_antimony`, `parse_antimony` and `load_antimony_file`.
- `examples/brusselator.ant`: A = 1, B = 3.
- `netevolve_core` static library, linked by both executables.

Tests (Release): **3/3 CTest tests passed. Doctest: 11/11 cases and 775,806 assertions passed.**
- RNG: the same seed gives the same stream, and helpers stay in range.
- Random generation is deterministic for a fixed seed.
- 10,000 random networks pass the validity checker. F, B and R vary from network to network.
- The validity checker rejects each kind of invalid network.
- Round trip (write → parse → write) gives identical text for 1,000 random networks, and the parsed network is bit-identical.
- The writer reproduces the example in SPEC §3 exactly.
- `examples/brusselator.ant` parses to the expected stoichiometry, constants and concentrations.
- Malformed input gives an error with the correct line number. Cases covered: an unknown species in an equation or a rate law, a missing `;` in three positions, a rate law that doesn't match its reactants (two variants), two parameters in a rate law, a boundary species not declared `const`, a parameter with no value, unsupported syntax, and a duplicate reaction.

## Decisions and notes on ambiguity

### M0
- The spec is copied verbatim to `docs/SPEC.md`. The original stays at `../docs/netevolve-spec.md`.
- The `--version` check is two CTest tests: one checks the exit code, the other matches the output. A regex match by itself would ignore the exit code.

### M1
- **Species indexing:** a single index addresses every species. `0..F-1` are floating species `S*` and `F..F+B-1` are boundary species `X*`. Reaction sides are stored in canonical form, sorted by species with one `Term {species, coeff}` per species, so comparing two reactions is an exact `==`.
- **Antimony names:**
  - Species must be named `S<n>` (floating) or `X<n>` (boundary), with no leading zeros. Any other name on a reaction side is an "unknown species" error.
  - `X*` must be declared `const`, and only `X*` may be.
  - F and B are taken from the highest index used. Every species from 0 to that index must have an initial value, so gaps are an error.
- **Rate law:** exactly one parameter, of any name that isn't a species, multiplied by species names. The species must match the reactant multiset exactly. Two reactions may share a parameter. A parameter that is assigned but never used is an error, to catch typos.
- **Duplicates:** the parser rejects them rather than merging them, and random generation redraws until there are none.
- **Writer:**
  - Coefficients are written as repeated species (`S2 + S2`), as in the SPEC §3 example.
  - Rate constants are written in the shortest form that reads back as exactly the same double (`std::to_chars`), which makes the round trip bit-exact.
  - The parser also accepts the `2 S1` form.
- **Not supported:** negative numbers, empty reaction sides and unnamed reactions are parse errors.
- **Validity:** at least one floating species is required, and `k` must be finite and > 0. `check_validity` doesn't check `k` against `[k_min, k_max]`, because a loaded model may use any positive k. The generator and mutation keep `k` within range.
- **Brusselator boundaries:** the sink products D and E share one boundary species (`X2`). Its concentration doesn't affect the dynamics.

### M2
- **Jacobian error measure:** the finite-difference comparison uses the matrix-norm relative error, max|J_analytic − J_fd| / max|J_analytic|. Entry-wise relative error isn't meaningful for entries that are zero or cancel to almost zero. The FD step is h = 1e-5·max(1, |s|).
- **Boundary-derivative test:** "boundary species have zero derivative" is tested structurally. Boundary species have no state entry, no row in N and no Jacobian column, and the model never writes their concentrations.

### M3
- **Stiff method:** RODAS3 rather than ROS3P, because ROS3P is A-stable but not L-stable, and the spec asks for L-stable. RODAS3's coefficients are simple rationals, and it's stiffly accurate.
- **Brusselator reference:** T = 7.1569199212, from SciPy `solve_ivp` with DOP853 and with Radau (rtol 1e-12, atol 1e-14; they agree to 1e-10). It's the mean spacing of upward crossings of x = 1 over t ∈ [200, 400], starting from (1, 1). The script lives outside the repo; the method is recorded in a comment in the test.
- **Robertson reference** at t = 1e4: (0.1073004285, 4.800166973e-7, 0.8926990914), from SciPy Radau and LSODA (rtol 1e-12). This check goes beyond what the spec asks.
- **Step budget:** it applies per solver, so in Auto mode RK45 and RODAS3 each get `max_steps`. RODAS3 continues from where RK45 stopped instead of restarting.
- **Wall-clock guard and determinism:** results depend on the wall clock only when the guard actually fires, which a normal simulation shouldn't. The default limit is 10 s per simulation. If this ever matters for byte-identical logs, raise the limit.
- **Dense output:** for RODAS3 it's a cubic Hermite (3rd order, the same as the method). Both interpolants preserve linear conservation laws.

### M4
- **Amplitude threshold:** `min_rel_amplitude × mean(|y|)` instead of `× mean(y)`. For concentrations, which are never negative, the two are identical. For the zero-mean `sin(t)` test signals that the spec asks for, `mean(y) ≈ 0` would make the filter useless.
- **Filtering:** "keep only alternating pairs whose height difference exceeds the threshold" is implemented as hysteresis over the raw extrema. Same-kind candidates merge into the more extreme one, and an opposite-kind candidate confirms the pending extremum only if it differs by more than the threshold.
- **Pairing:** pairs are consecutive extrema (e0, e1), (e2, e3), …, starting from the first confirmed extremum, whatever its kind. So `n = floor(#extrema / 2)`, and a pair's amplitude is |e_2k − e_2k+1|.
- **Noise in the test:** "noise at 0.1% of the amplitude" is read as uniform noise in ±0.001 for a unit-amplitude sine.
- **Solver settings for fitness:** rtol 1e-6 and atol 1e-9 (defaults), max_conc from the config. The config doesn't expose solver tolerances, because the spec doesn't list them.
- **`simulate` on failure:** it writes the samples computed so far, prints a warning and exits with 1.
- **`evolve` and `prune`:** they print "not implemented yet" and exit with 1 until M5 and M7.

### M5
- **Generation count:** generation 0 is the random initial population, and `generations = G` gives up to G+1 evaluated generations, so G+1 rows in `log.csv`. "Solved within 50 generations" is checked as generation index ≤ 50.
- **Elite count:** `round(elite_fraction × P)`, at least 1 if elite_fraction > 0, and at most P−1. Parent pool: `round(parent_fraction × P)`, clamped to [1, P].
- **Child validity:** "retry mutation (bounded attempts) or recopy" means up to 10 attempts, each with a freshly chosen kind, then an unchanged copy of the parent.
- **Rewire:** chooses a molecule slot uniformly over all reactant and product molecules of a uniformly chosen reaction, and replaces it with a uniformly chosen *different* species.
- **Absolute amplitude floor (1e-6):** a deliberate addition to SPEC §1.4. It's a constant in `OscillationParams`, not a config key, because SPEC §2 fixes the key list. It's cheap to change or remove.
- **RODAS3 sampling:** steps land exactly on every sample time. This caps the RODAS3 step at the sample spacing (0.1 by default), so a stiff simulation takes at least `n_points` steps.

### M6
- **Threads:** the pool is created per generation instead of being kept alive, because that's simpler and the overhead is negligible.
- **Wall-clock guard:** this is the only possible source of thread-count dependence. If heavy CPU load ever made a simulation exceed `max_wall_seconds` (10 s, against a typical 0.7 ms), its fitness would become 0. It hasn't happened in any run so far.

### M7
- **"Fitness stays ≥ 0.95 × the original":** the reference is the input network's fitness, fixed for the whole prune, not the fitness after earlier removals.
- **Prune floor:** prune respects `r_min` from the config (default 3), so the pruned network is still valid by SPEC §1.1.
- **Long-run amplitude:** the mean |peak − trough| over the filtered extrema in each window, for the species with the best fitness score. Plain max − min would depend on how the samples happen to fall on sharp spikes.
- **Long-run windows:** both have the fitness window's length, (1 − transient_fraction) × t_end = 50, and end at t_end and at 5 × t_end.
- **Acceptance pass rule:** a run passes if it reaches fitness ≥ 0.95 *and* its pruned network passes the long-run check. At least 3 of 5 such runs are needed.
- **Excluding the acceptance test from the default run:** the test exits with code 77 (`SKIP_RETURN_CODE`) unless `NETEVOLVE_RUN_SLOW=1` is set, so a plain `ctest -C Release` reports it as skipped, not failed.

## Open issues

- **Acceptance margin:** the test passes 3 of 5, exactly the threshold. Seed 2 is a genuine limit cycle that the spec's long-run check rejects because it's still settling at t_end = 100. Parameter options are listed in the M7 notes; none is applied.
- **Tellurium/CVODE:** fragile on some evolved spiky oscillators; see M7 for the settings that work.

