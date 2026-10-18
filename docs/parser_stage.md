# Parser Stage Guide

The parser stage is the public `pdb-parser` CLI. It consumes a PDB-id list plus
`data/params.cfg`. With `output_mode: parser_outputs`, it writes parser-stage
`X_*`, `A_*`, and `I_*` files; with `output_mode: filtered_atoms`, it writes
only the persistent filtered-coordinate TSV and stops.

Current public `atom_selection` values:

- `full_chain`
- `backbone`
- `backbone_plus_hydrogens`
- `backbone_plus_neighbors`

Public `method` values are `auto`, `nmr`, and `xray`. `auto` accepts only NMR
and X-ray diffraction according to `EXPDTA`. For X-ray,
`backbone_plus_hydrogens` becomes `backbone_plus_neighbors` with a warning;
the other selections remain unchanged.

The parser/remapping method policy is narrower than the parser selection API:
`backbone` may be remapped for NMR or X-ray, while
`backbone_plus_hydrogens` may be remapped only for NMR. `full_chain` and
`backbone_plus_neighbors` remain parser-only layouts.

`filtered_atoms` is not an input for the remap stage; remapping requires the
`X_*`, `A_*`, and `I_*` files produced by `parser_outputs`.

## Input contract

Each execution selects exactly one model and one chain; `full_chain` means all
retained atoms from that chain, never atoms from multiple chains. Residues
numbered `<= 0` are dropped with a warning; the remaining ones must have
strictly consecutive numbers (the first number need not be 1). The positional
`X_*`, `A_*`, and `I_*` formats cannot encode insertion codes, so for each
residue number the blank-code residue is kept (otherwise the alphabetically
first code) and the others are discarded with a warning. Alternate locations
use a deterministic policy: blank is preferred, otherwise the lexicographically
first conformer is retained for each residue/atom. After these rules,
consecutively numbered residues must be bonded: a C(i)--N(i+1) distance above
2.0 Å is a chain break and the structure is rejected (this includes inserted
residues whose removal leaves a gap that the numbering hides). Protein residues must have
a covalent template in `RESIDUE_BONDS`; unsupported residue names are rejected.

Configuration keys are unique. Numeric values must be finite; distance
cutoffs, interval widths, and VDW radii must be greater than zero, even when the
selected mode does not use them; `torsion_angle_width` must be zero or greater. For NMR
experimental intervals, the bounds must satisfy
`vdwr_hh < noe_strong <= noe_medium <= noe_weak`.

The same one-chain, consecutive-residue, insertion-code, alternate-location,
and residue-template rules are applied before parser artifacts are written.

The CLI processes the full PDB-id list even if individual structures fail. It
returns `0` when all jobs complete, `1` when at least one structure is skipped
or fails, and `2` for a configuration error before jobs start.
Unexpected ordinary exceptions are also reported for the affected PDB and do
not stop the remaining jobs. Interrupts and explicit process exits still
abort immediately.

## Commands

Installed entry point:

```bash
pdb-parser data/pdb_ids.txt data/params.cfg data/pdb data/seeds data/outputs
```

Direct module execution:

```bash
python -m pdb_parser.pdb_parser data/pdb_ids.txt data/params.cfg data/pdb data/seeds data/outputs
```

The intermediate files of each structure (`filtered_atoms.dat`,
`distance_constraints_*.dat` and `angular_constraints_5.dat`) are written to
`data/outputs/<pdb_id>/tmp/` and removed when the structure finishes. Pass
`--keep-tmp` to keep them:

```bash
pdb-parser --keep-tmp data/pdb_ids.txt data/params.cfg data/pdb data/seeds data/outputs
```

When a structure fails, its `tmp/` directory is always removed and its previous
outputs are restored.

## Positional arguments

| argument | description |
|----------|-------------|
| `data/pdb_ids.txt` | text file containing the list of PDB identifiers to process |
| `data/params.cfg` | configuration file controlling the parser behavior |
| `data/pdb` | directory where PDB structures will be stored |
| `data/seeds` | directory where per-PDB seed files will be read and written |
| `data/outputs` | directory where the generated constraint files will be written |

## Main outputs

For each processed `pdb_id`, files are written to:

```text
data/outputs/<pdb_id>/
```

Main parser outputs:

- `X_<pdb_id>_model<model_number>_chain<chain_id>.dat`
  Structure file with one atom per line.
  Columns: `atom_id atom_name resid resname x y z`
  Format: `%5d %-4s %4d %3s %9.3f %9.3f %9.3f`
- `A_<pdb_id>_model<model_number>_chain<chain_id>.dat`
  Angular constraint file with one residue per line.
  Columns: `resid resname omega_center omega_radius phi_center phi_radius psi_center psi_radius`
  Format: `%4d %3s %17.12f %17.12f %17.12f %17.12f %17.12f %17.12f`
- `I_<pdb_id>_model<model_number>_chain<chain_id>.dat`
  Distance-constraint file.
  Columns: `atom_id_i atom_id_j resid_i resid_j d_l d_u atom_name_i atom_name_j resname_i resname_j`
  Format: `%5d %5d %6d %6d %16.12f %16.12f %4s %4s %s %s`

## Seed files

Per-PDB seed files are stored in:

```text
data/seeds/<pdb_id>/seed_<pdb_id>.dat
```

The parser creates or reuses the file only when the current execution
actually needs stochastic streams. The file stores only the seed fields that
can affect the current execution:

- `seed_interval_centered_distances` only when `distance_constraints = interval_centered`
- `seed_talos_angle_centers` only when the effective torsion interval width is greater than `0`
- `seed_phi_psi_mask` only when the `phi/psi` selection is genuinely random
- `seed_hh_distance_subset` only for NMR with `0 < percentage_hh_distance_constraints < 100`

If the current execution is fully deterministic, no seed file is created.
