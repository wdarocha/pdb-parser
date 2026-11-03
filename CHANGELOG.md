# Changelog

## 0.4.0

Parser outputs are unchanged. In the remap stage every file is renamed; the contents of
`ddgp` and `backbone` outputs are unchanged, and `redmdgp` `I_*` and `T_*` change (below).

### Breaking changes

- Remapped suffixes: `*_ddgp_order_<N>.dat` (was `*_ddgp_hc_order_<N>.dat`) and
  `*_redmdgp_hc_order_<N>.dat` (was `*_redmdgp_order_<N>.dat`). `hc` (hand-crafted) now
  marks only the `redmdgp` patterns.

### New

- `order_type: dmdgp` requests the fixed DMDGP order `N CA C` of `backbone` outputs (NMR or
  X-ray), without `order_id`. It cannot be combined with `ddgp` or `redmdgp`, and
  `backbone_plus_hydrogens` outputs reject it. `ddgp` or `redmdgp` on `backbone` outputs
  still write that order, with a warning suggesting `order_type: dmdgp`.

### Changes in `redmdgp` outputs

- Every copy of an atom has four predecessors, one of them a `[0, 0]` with the previous
  occurrence: when that occurrence is `i-3`, `i-4` is added at `d(i-3, i-4)`; otherwise the
  zero is an extra row (before, only a copy placed by an interval had it). `I_*` grows by
  up to about 2%.
- `T_*` signs: a torsion of 0 has sign `+1` (in every `T_*`; `0` means an unknown sign), and
  the copy of `C(k-1)` placed from `C(k) N(k) CA(k)` (orders 1 and 2) gets the sign of
  `-phi(k)` with the `phi` interval of `A_*`, instead of sign `0`.

### Documentation

- A step-by-step `order_type: dmdgp` example (`docs/remap_stage.md`, README);
  `data/instance_remap.cfg` lists the remappable outputs per method.
- `redmdgp` pattern legend: N-terminal hydrogens labelled `h1 h2 h3` (placement order, was
  `n1 n2 n3`), and a marked copy does not branch the search.
- `data/params.cfg`: `max_distance` also applies to `precise`, and `interval_experimental`
  uses neither the epsilons nor `max_distance`.

## 0.3.0

Except for the fixes below, the interval formulas, van der Waals radii, planar
templates, angle-to-distance conversion and remapping orders are unchanged. For
the same NMR input, atom set and seed file, the values match 0.2.0 up to the new
12-decimal precision; torsion-derived values now come from 12-decimal angles.

### Breaking changes

- The remapping stage was renamed from *reorder* to *remap*: CLI
  `instance-remap`, config `data/instance_remap.cfg`, module
  `pdb_parser.instance_remap`, package `pdb_parser.remapping`, output directory
  `data/outputs/<pdb_id>/remapped/` and suffix `_ddgp_hc_order_<N>.dat`. The
  `run_pipeline.sh` modes are `remap` and `remap-all`, and every mode except
  `parser` requires `output_mode: parser_outputs`.
- Output files use aligned, space-separated columns and 12 decimals for every
  distance and angle:
  - parser `X_*`: `%5d %-4s %4d %3s %9.3f %9.3f %9.3f` (0.2.0 could glue a
    coordinate ≤ -100 Å to the previous one);
  - `A_*`: `%4d %3s` followed by six `%17.12f` (previously `%9.4f`);
  - `I_*`: `%5d %5d %6d %6d %16.12f %16.12f %4.4s %4.4s %s %s` (previously
    `%20.16f`);
  - `T_*`: `%5d %5d %5d %5d %2d %16.12f %16.12f`;
  - remapped `X_*`: `%9.3f %9.3f %9.3f`.

  Files written by earlier releases remain readable, except a 0.2.0 `X_*` with a
  glued coordinate, which must be regenerated.
- The configuration is validated as a whole before any structure is processed,
  and an error exits with status `2`. Rejected: unknown or repeated keys,
  non-finite numbers, non-positive distances (even when the selected mode does not
  use them), a negative `torsion_angle_width`, `percentage_backbone_torsion_angles`
  outside 0--100 (it was clamped), and NOE bounds that do not satisfy
  `vdwr_hh < noe_strong <= noe_medium <= noe_weak`.
- `distance_constraints: precise` forces `vdw_constraints: no` and
  `tightening_upper_bounds: no`.
- Input rules, which can change the atom set and therefore every output:
  - residues numbered ≤ 0 are dropped with a warning;
  - residue numbers must be consecutive, and consecutive residues must be bonded
    (C(i)--N(i+1) ≤ 2.0 Å), so chain breaks hidden by the numbering are rejected;
  - for each residue number, the residue without insertion code is kept
    (otherwise the first code); the others are discarded with a warning;
  - alternate locations are resolved deterministically (blank conformer,
    otherwise the alphabetically first);
  - atoms that repeat the coordinates of an earlier atom are dropped with a
    warning;
  - residues without a covalent template and ambiguous chain ids are rejected.
- `backbone_plus_hydrogens` and `backbone_plus_neighbors` keep the atoms that are
  present (and accept `HN`) instead of failing when H, HA, O or CB is missing. In
  prolines, `backbone_plus_hydrogens` keeps the delta hydrogen nearer the `HA` of
  the previous residue (0.2.0 always kept `HD2`).
- The CLIs keep processing after a failed structure but exit with status `1`; a
  structure skipped by the remap stage also counts. Failures are reported as
  `[SKIP] <pdb_id>: ...` and warnings go to stderr.
- The `tmp/` directory of each structure is removed when it finishes; use
  `--keep-tmp` to keep it.
- The example `data/params.cfg` uses `interval_experimental`,
  `torsion_angle_width: 60.0` and `tightening_upper_bounds: yes`.

### New

- X-ray structures: `method: auto | nmr | xray`, where `auto` reads `EXPDTA`. For
  X-ray, the long-range restraints are CA--CA distances, `backbone_plus_hydrogens`
  becomes `backbone_plus_neighbors`, `interval_experimental` becomes
  `interval_centered`, and a one-model structure uses model 1 (each with a warning).
- `output_mode: filtered_atoms | parser_outputs`, and a
  `metadata_<pdb_id>_model<m>_chain<c>.json` file read by the remap stage.
- Remapping of `backbone` outputs (`*_dmdgp_order.dat`; the order is ignored).
- `order_type: ddgp | redmdgp`: `redmdgp` writes each `order_id` as an iDMDGP
  re-order (Lavor et al., Discrete Appl. Math. 256 (2019) 91-104, Def. 3), placing
  every vertex from the three previous ones and repeating atoms when needed
  (`*_redmdgp_order_<N>.dat`).
- Several remap orders in one call: `order_id: 1-10` or `1,3,9`, and
  `order_type: ddgp, redmdgp`. The parser files are read once per structure, and a
  failed combination does not stop the others. `run_pipeline.sh --mode remap-all`
  makes a single call.
- `tightening_upper_bounds: yes`: every upper bound of `I_*`, including the `999`
  van der Waals bounds, becomes its shortest path over the finite bounds.
- `percentage_hh_distance_constraints` (NMR): keeps a reproducible random
  percentage of the H--H restraints. The H--H edges the remap needs for
  discretization are always kept.
- Outputs are written transactionally, invalid cached PDB files are downloaded
  again, and downloads retry transient network failures.
- Dependencies are bounded by major version (`docs/dependencies.md`), with a
  `[dev]` extra for the tests; CI runs on Python 3.10 to 3.13.

### Fixes

- The covalent templates lacked aromatic ring bonds (PHE, TYR, HIS, TRP) and bonded
  the TRP indole hydrogen to the wrong atom; `full_chain` with tightening failed.
- Lower bounds follow one rule in every mode: the base value (`vdwr_hh` or the sum
  of van der Waals radii) when it is below the reference distance `d`, otherwise
  the first of 90%, 80%, ..., 10% of it below `d`, otherwise `0.9 d`. In 0.2.0 an
  NOE lower bound could exceed `d`, so the reference structure violated its own
  instance.
- The remap stage failed for a first residue with three N-terminal hydrogens but
  without `H2` or `H3`, for cis prolines in orders 7 to 10, and for atoms with
  repeated coordinates. A missing discretization interval is now reported as a
  data issue naming the atoms.
- Numbers are read back exactly (`float_precision="round_trip"`).

### Performance

The parser and the remap stage are several times faster, with the same outputs:
the merge, validators, van der Waals bounds and NMR and X-ray restraints are
vectorized (the X-ray restraint time grew with the cube of the chain length), the
merge keeps the tables in memory and rounds each value once, and MDAnalysis and
SciPy are imported only when needed.

### Compatibility notes

- Seed files from 0.2.0 remain valid and reproduce the same random streams for the
  same atom set. Seed values outside `1`--`2^63 - 1` are rejected.
