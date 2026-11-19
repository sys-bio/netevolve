# Progress

## Current milestone

**M3 (integrators) is complete.** It's waiting for the human's approval before M4 starts.

Remote: https://github.com/zyrrron/netevolve (private), branch `main`.

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

## Open issues

- None.
