# Progress

## Current milestone

**M1 (network representation, RNG, Antimony I/O) is complete.** It's waiting for the human's approval before M2 starts.

Remote: https://github.com/zyrrron/netevolve (private), branch `main`.

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

## Open issues

- None.
