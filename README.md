# pdb-parser

Lightweight Python toolkit for parsing **Protein Data Bank (PDB)**
Nuclear Magnetic Resonance (NMR) and X-ray structures and generating **distance and torsion-angle constraints**
tailored for **Distance Geometry Problem (DGP)** workflows.

The package extracts structural information from PDB files and produces
constraint sets that can be used in algorithms for protein structure
determination.

Typical applications include:

-   generation of **distance-constraint lists**
-   generation of **backbone torsion angles**
-   extraction of **3D coordinates**

------------------------------------------------------------------------

## 🔬 Project structure

    pdb-parser/
    │
    ├── data/                  # configuration and input lists
    │   ├── params.cfg
    │   ├── instance_remap.cfg
    │   └── pdb_ids.txt
    │
    ├── src/pdb_parser/
    │   ├── pdb_parser.py      # CLI for parsing and writing output files
    │   ├── instance_remap.py # CLI for remapping previously written outputs
    │   ├── geometry/          # distance geometry utilities
    │   ├── io/                # PDB parsing and filtering routines
    │   ├── pipeline/          # main parsing pipeline
    │   ├── remapping/         # remap order patterns and dataframe reindexing
    │   └── utils/             # auxiliary utilities
    │
    ├── pyproject.toml
    ├── README.md
    └── .gitignore

------------------------------------------------------------------------

## Installation

This project requires Python 3.10 or newer. Clone the repository:

    git clone https://github.com/wdarocha/pdb-parser.git
    cd pdb-parser

Create and activate a virtual environment if desired, then install the package in editable mode:

    python -m pip install -e .

This installs the runtime dependencies declared in `pyproject.toml`: `MDAnalysis`,
`numpy`, `pandas`, `scipy`, and `requests`.

Dependency compatibility is intentionally bounded by major version: NumPy
1.23--2.x, pandas 1.5--2.x, SciPy 1.9--1.x, MDAnalysis 2.4--2.x, and requests
2.28--2.x. The `pyproject.toml` file is the single packaging metadata source.
The supported Python versions, rationale for each bound, and the minimum
verification matrix are recorded in [`docs/dependencies.md`](docs/dependencies.md).

To run the test suite and the lint step used by CI, install the development extras:

    python -m pip install -e .[dev]

This exposes the command-line entry points:

    pdb-parser --help
    instance-remap --help

If you prefer not to install the package, run the modules directly from the source tree with `PYTHONPATH=src`:

    PYTHONPATH=src python -m pdb_parser.pdb_parser --help
    PYTHONPATH=src python -m pdb_parser.instance_remap --help

------------------------------------------------------------------------

## Quick check

From the repository root, the following commands should display the available command-line options:

    PYTHONPATH=src python -m pdb_parser.pdb_parser --help
    PYTHONPATH=src python -m pdb_parser.instance_remap --help
    ./run_pipeline.sh --help

The full local verification (lint, tests, and the 70% coverage gate used by CI) is:

    make check

or, step by step:

    python -m ruff check .
    python -m pytest -q

------------------------------------------------------------------------

## 📂 Input configuration

The parser stage is controlled by `data/params.cfg`.

The instance remapping stage is controlled by `data/instance_remap.cfg`.

Both files use the same `key: value` syntax and support `#` comments.

Unknown keys are rejected by the current public configuration interface.
Repeated keys are rejected as ambiguous. Numeric parameters must be finite;
distance cutoffs, interval widths, and VDW radii must be strictly positive,
even when the selected mode does not use them; `torsion_angle_width` must be
zero or greater.

Below is a complete description of the parameters currently supported by each stage.

------------------------------------------------------------------------

## Public interface in this release

The current public interface is intentionally narrower than the internal codebase.

Publicly supported through configuration files:

- parser CLI: `pdb-parser`
- remapping CLI: `instance-remap`
- parser config file: `data/params.cfg`
- remapping config file: `data/instance_remap.cfg`
- experimental methods: `auto`, `nmr`, `xray`
- output modes: `filtered_atoms`, `parser_outputs`
- atom selections accepted by public parser config:
  - `full_chain`
  - `backbone`
  - `backbone_plus_hydrogens`
  - `backbone_plus_neighbors`
- distance-constraint modes accepted by public parser config:
  - `precise`
  - `interval_centered`
  - `interval_experimental`
- remapping `order_id` values accepted by public config:
  - `1..10`, alone or as a list of values and ranges (`1-10`, `1,3,9`)
- remapping `order_type` values: `ddgp` (default), `redmdgp`, or both (`ddgp, redmdgp`),
  for `backbone_plus_hydrogens`; one remap call writes every listed combination
- `order_type: dmdgp` for `backbone` outputs: the fixed DMDGP order `N CA C`, without `order_id`

Internal code paths may contain additional options for future releases, but they are not part of the supported configuration surface of this release.

### Input contract

Each parser execution selects exactly one model and one chain. `full_chain` means
all atoms from that selected chain; it never combines chains. Residues numbered
`<= 0` are dropped with a warning; the remaining residue numbers may start at
any positive value, but must be strictly consecutive. The `X_*`, `A_*`, and
`I_*` formats do not encode insertion codes, so for each residue number the
residue with a blank insertion code is kept (otherwise the alphabetically first
code) and the other insertion-coded residues are discarded with a warning. For alternate locations, a blank conformer is preferred; when no
blank conformer exists, the lexicographically first conformer is selected for
each residue/atom identity. After these rules, consecutively numbered residues
must also be bonded: a C(i)--N(i+1) distance above 2.0 Å is a chain break and the
structure is rejected (for example, residues missing from the model, or
insertion-coded residues whose removal leaves two unbonded residues).

Parser outputs require the standard amino-acid residue templates present in
`RESIDUE_BONDS` (ALA, GLY, PRO, SER, THR, VAL, LEU, ILE, MET, CYS, ASN, GLN,
ASP, GLU, HIS, PHE, TYR, TRP, LYS, and ARG). A residue outside this set is
rejected rather than producing incomplete covalent topology.

------------------------------------------------------------------------

## Experimental-method policy

`auto` reads `EXPDTA` and supports only NMR and X-ray diffraction. NMR keeps
the existing H--H long-range policy (or N/C proxies in backbone mode). X-ray
uses backbone angular distances for adjacent residues and CA--CA long-range
distances for non-adjacent residues. An X-ray request for
`interval_experimental` is converted to `interval_centered` with a warning.

For X-ray coordinates, `backbone_plus_hydrogens` is converted to
`backbone_plus_neighbors`, because hydrogens are normally absent.

The `noe_strong`, `noe_medium`, `noe_weak`, `vdwr_hh`, and
`percentage_hh_distance_constraints` parameters are NMR-only. For X-ray, `interval_centered` uses `epsilon_long` and
`max_distance`; its CA--CA lower bound is internal (`2 × vdw_c`).

------------------------------------------------------------------------

## Compatibility policy

The following behaviors are treated as compatibility-sensitive and change only in a documented release:

- the scientific logic used to build the generated instances
- the accepted public values already documented for this release
- the naming pattern of parser outputs `X_*`, `A_*`, `I_*`
- the naming pattern of remapped outputs `X_*`, `I_*`, `T_*`
- the current text-file structure of generated instance files
- the current directory layout under `data/outputs/<pdb_id>/`
- the meaning of `order_id` values `1..10`

See [`CHANGELOG.md`](CHANGELOG.md) for changes between releases.

------------------------------------------------------------------------

## Future and experimental surface

The repository may contain implementation hooks that support future releases, experiments, or internal workflows. Those code paths should be interpreted with the following rule:

- implemented does not automatically mean publicly supported
- undocumented input values are not part of the current release contract
- internal options may change between releases until they are explicitly documented in this `README`

For practical use, treat this document and the public config validation rules as the release contract. Everything else should be considered internal or provisional until promoted here.

------------------------------------------------------------------------

## Architecture overview

The project is organized as a small pipeline with two public stages.

1. Parser stage
   - downloads or reads one PDB structure
   - validates the selected experimental method, model and chain
   - filters the public atom set
   - generates `X_*`, `A_*`, and `I_*`

2. Remapping stage
   - reads parser outputs for one structure
   - applies one or more fixed `order_id` patterns
   - writes remapped `X_*`, `I_*`, and `T_*`

Main module responsibilities:

- `src/pdb_parser/pdb_parser.py`
  - public CLI for the parser stage
- `src/pdb_parser/instance_remap.py`
  - public CLI for the remapping stage
- `src/pdb_parser/pipeline/`
  - stage orchestration and public pipeline flow
- `src/pdb_parser/io/`
  - file conversion, filtering, PDB access, and constraint writers
- `src/pdb_parser/remapping/`
  - order patterns, clique construction, and dataframe reindexing
- `src/pdb_parser/geometry/`
  - geometric helper functions used by the pipeline
- `src/pdb_parser/utils/`
  - config readers, filesystem helpers, and seed management

------------------------------------------------------------------------

## Operational Guides

The detailed operational material now lives in dedicated short guides inside `docs/`:

- parser stage: [`docs/parser_stage.md`](docs/parser_stage.md)
- remap stage: [`docs/remap_stage.md`](docs/remap_stage.md)
- pipeline wrapper: [`docs/pipeline_wrapper.md`](docs/pipeline_wrapper.md)

Use those guides for command examples, positional arguments, output layout, and seed behavior.

At a high level:

- `pdb-parser` reads `data/pdb_ids.txt` and `data/params.cfg`, downloads or reuses PDB files, and writes either a filtered-coordinate file (`output_mode: filtered_atoms`) or parser-stage `X_*`, `A_*`, and `I_*` files (`output_mode: parser_outputs`) under `data/outputs/<pdb_id>/`.
- `instance-remap` reads parser-stage outputs plus `data/instance_remap.cfg`, applies one or more fixed `order_id` values (as DDGP orders or, with `order_type: redmdgp`, as iDMDGP re-orders) for `backbone_plus_hydrogens`, or the fixed DMDGP order `N CA C` for `backbone` (`order_type: dmdgp`), and writes remapped `X_*`, `I_*`, and `T_*` files under `data/outputs/<pdb_id>/remapped/`.
- `run_pipeline.sh` orchestrates parser-only, remap-only, full, and `remap-all` executions from the source tree.

The command-line tools keep processing the remaining PDB identifiers after a
per-structure failure, but return exit status `1` when one or more identifiers
fail. Configuration errors return status `2` before processing begins.

The remap stage supports parser outputs generated with:

- `atom_selection: backbone`
- `atom_selection: backbone_plus_hydrogens`

`backbone_plus_neighbors` is a public parser selection and is also the X-ray
effective selection when `backbone_plus_hydrogens` is requested. It is not a
remapping input. `full_chain` and `filtered_atoms` likewise cannot be remapped.

Remapping is supported for `backbone` with NMR or X-ray outputs and for
`backbone_plus_hydrogens` with NMR outputs only. The X-ray fallback is
therefore parser-only.

For `backbone`, the remapped `T_*` uses the simplified chain rule based on the
three previous vertices and parser angles from `A_*`. In this layout,
`instance_remap.cfg` still selects `model_number` and `chain_id`, the order is
requested with `order_type: dmdgp` (no `order_id` needed), and the remapped outputs
always use the suffix `*_dmdgp_order.dat` (step-by-step example in
[`docs/remap_stage.md`](docs/remap_stage.md#example-the-dmdgp-order-of-backbone-outputs)). For
`backbone_plus_hydrogens`, the existing DDGP/iDDGP clique construction remains
in place.

The generated files are compatibility-sensitive contracts:

- parser-stage `X_*` keeps atom ids, atom names, residue ids, residue names, and coordinates together;
- parser-stage `A_*` keeps one row per residue with `omega/phi/psi` centers and radii;
- parser-stage `I_*` keeps merged distance constraints with atom and residue metadata;
- remapped `X_*` remains coordinate-only in remapped atom order;
- remapped `I_*` remains the same logical constraint set, reindexed to the remapped order;
- remapped `T_*` remains the seven-column clique matrix consumed by the downstream DDGP workflow.

------------------------------------------------------------------------

## Reference And Theory

The public operational surface is intentionally summarized here.

Short theoretical summary:

- the parser builds backbone-focused `X_*`, `A_*`, and `I_*` files from one NMR or X-ray model/chain;
- the remapping stage converts those files into DDGP/iDDGP-oriented `X_*`, `I_*`, and `T_*` outputs;
- the current public parser atom selections are `full_chain`, `backbone`,
  `backbone_plus_hydrogens`, and `backbone_plus_neighbors`;
- the current public remapping choices are the fixed `order_id=1..10`, with `order_type: ddgp` (default) or `redmdgp`.

Detailed reference material was moved to:

- [docs/remap_reference.md](docs/remap_reference.md)

That document keeps:

- parser parameter reference and examples;
- distance/torsion modeling notes;
- remapping configuration details;
- internal-residue order tables used by the public interface.

------------------------------------------------------------------------

## 📖 Citation

If this code is useful in your research, please cite the repository (also available via the **“Cite this repository”** button on GitHub thanks to the included `CITATION.cff`) and cite the preprint associated with the ordering used in your experiments.

The remapping code supports fixed `order_id` values from `1` to `10`. At the
moment, the preprint references documented in this repository are:

- `order_id=1`: `A hybrid combinatorial-continuous strategy for solving molecular distance geometry problems` (arXiv:2510.19970).
- `order_id=9`: `An Angle-Based Algorithmic Framework for the Interval Discretizable Distance Geometry Problem` (arXiv:2508.09143).

If your work uses one of these two orderings, cite the corresponding preprint. If your work uses both, cite both preprints.

### Preprint for `order_id=1`

**Leonardo D. Secchin, Wagner da Rocha, Mariana da Rosa, Leo Liberti, Carlile Lavor.**  
*A hybrid combinatorial-continuous strategy for solving molecular distance geometry problems.*  
**arXiv:2510.19970**, 2025.  
https://arxiv.org/abs/2510.19970

### BibTeX for `order_id=1`
```bibtex
@misc{secchin2025hybrid,
  title        = {A hybrid combinatorial-continuous strategy for solving molecular distance geometry problems},
  author       = {Leonardo D. Secchin and Wagner da Rocha and Mariana da Rosa and Leo Liberti and Carlile Lavor},
  year         = {2025},
  eprint       = {2510.19970},
  archivePrefix= {arXiv},
  primaryClass = {math.OC},
  url          = {https://arxiv.org/abs/2510.19970}
}
```

### Preprint for `order_id=9`

**Wagner da Rocha, Carlile Lavor, Leo Liberti, Leticia de Melo Costa, Leonardo D. Secchin, Therese E. Malliavin.**
*An Angle-Based Algorithmic Framework for the Interval Discretizable Distance Geometry Problem.*  
**arXiv:2508.09143**, 2025.  
https://arxiv.org/abs/2508.09143

### BibTeX for `order_id=9`
```bibtex
@misc{darocha2025anglebased,
  title        = {An Angle-Based Algorithmic Framework for the Interval Discretizable Distance Geometry Problem},
  author       = {Wagner da Rocha and Carlile Lavor and Leo Liberti and Leticia de Melo Costa and Leonardo D. Secchin and Therese E. Malliavin},
  year         = {2025},
  eprint       = {2508.09143},
  archivePrefix= {arXiv},
  primaryClass = {q-bio.BM},
  url          = {https://arxiv.org/abs/2508.09143}
}
```

------------------------------------------------------------------------

## 📜 License

This repository is licensed under the [MIT License](./LICENSE).
© 2025–2026 Wagner da Rocha

------------------------------------------------------------------------

## 👤 Author

Developed and maintained by [Wagner da Rocha](https://github.com/wdarocha).

------------------------------------------------------------------------

## 🙏 Acknowledgments

Special thanks to [Leonardo D. Secchin](https://github.com/leonardosecchin) for the valuable support provided during the development of the code.
