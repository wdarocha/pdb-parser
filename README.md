# pdb-parser

Lightweight Python toolkit for parsing **Protein Data Bank (PDB)**
structures and generating **distance and torsion-angle constraints**
tailored for **Distance Geometry Problem (DGP)** workflows.

The package extracts structural information from PDB files and produces
constraint sets that can be used in algorithms for protein structure
determination, including **Discretizable Distance Geometry Problem (DDGP
/ _i_DDGP)** frameworks.

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
    │   └── pdb_ids.txt
    │
    ├── src/pdb_parser/
    │   ├── geometry/          # distance geometry utilities
    │   ├── io/                # PDB parsing and filtering routines
    │   ├── pipeline/          # main parsing pipeline
    │   ├── reordering/        # vertex ordering strategies
    │   └── utils/             # auxiliary utilities
    │
    ├── pyproject.toml
    ├── README.md
    ├── CITATION.cff
    └── .gitignore

------------------------------------------------------------------------

# 📂 Input configuration

The behavior of the parser is controlled by a configuration file located
at

    data/params.cfg

This file defines how the PDB structures are processed and how the
geometric constraints used in the pipeline are generated.

Below is a complete description of all parameters currently supported.

------------------------------------------------------------------------

# Running the parser

The parser is executed as

    python3 src/pdb_parser/pdb_parser.py data/pdb_ids.txt data/params.cfg data/pdb data/outputs

The arguments are:

| argument | description |
|----------|-------------|
| `data/pdb_ids.txt` | text file containing the list of PDB identifiers to process |
| `data/params.cfg` | configuration file controlling the parser behavior |
| `data/pdb` | directory where PDB structures will be stored |
| `data/outputs` | directory where the generated constraint files will be written |

------------------------------------------------------------------------

# PDB list

The file

    data/pdb_ids.txt

contains a list of PDB identifiers, one per line.

Example

    1TOS
    1UAO

Each identifier corresponds to a structure that will be retrieved and
processed by the pipeline.

------------------------------------------------------------------------

# Configuration parameters (`params.cfg`)

## Model selection

Selects which model from the PDB structure will be used.

    model_number: 1

Model numbering follows the PDB convention and is **1-based**.

------------------------------------------------------------------------

## Chain selection

Specifies which chain in the selected model will be processed.

    chain_id: A

Only atoms belonging to this chain are considered.

------------------------------------------------------------------------

## Atom selection strategy

Defines which atoms are extracted from the PDB structure.

    atom_selection: backbone_plus_hydrogens

Currently supported option:

| option | description |
|-------|-------------|
| `backbone_plus_hydrogens` | backbone atoms (N, CA, C) plus hydrogens directly bonded to them |

------------------------------------------------------------------------

# Distance constraint model

Defines how **NMR-derived distance constraints** are generated.

    distance_constraints: interval_centered

Note that **covalent distances, planar constraints, and peptide-group
distances are always treated as precise**.

------------------------------------------------------------------------

## Synthetic distance intervals

When using

    distance_constraints: interval_centered

distance intervals are generated around the reference distance extracted
from the PDB.

The reference distance, $d_{ij}$, is perturbed as

$$
d_{ij}^* \sim \mathcal{N}\left(d_{ij},\left(\frac{\varepsilon_{ij}}{8}\right)^2\right)
$$

and the resulting interval is

$$
\mathcal{D}_{ij} =
\left[
\max\left(d_{ij}^* - \frac{\varepsilon_{ij}}{2},\ v_{\mathrm{dw}}\right),
\
\min\left(d_{ij}^* + \frac{\varepsilon_{ij}}{2},\ 5 \ \mathrm{Å}\right)
\right]
$$

where the interval width satisfies $\varepsilon_{ij} = 8\sigma$, corresponding to $\pm4\sigma$ around the mean.

------------------------------------------------------------------------

## Interval parameters

The interval width depends on whether the atoms belong to nearby
residues.

    epsilon_short: 1.0
    epsilon_long: 2.0
    min_distance: 2.4
    max_distance: 5.0

| parameter | meaning |
|-----------|---------|
| `epsilon_short` | interval width for atoms in the same or adjacent residues |
| `epsilon_long` | interval width for atoms in non-adjacent residues |
| `min_distance` | minimum allowed lower bound for distance intervals |
| `max_distance` | maximum allowed upper bound for distance intervals |

Suggested values are

    epsilon_short = 1.0 Å
    epsilon_long  = 2.0 Å
    max_distance  = 5.0 Å

------------------------------------------------------------------------

## van der Waals constraints

Lower-bound constraints based on van der Waals radii can optionally be
included.

    vdw_constraints: yes

Options:

| value | meaning |
|------|---------|
| `yes` | include van der Waals lower-bound constraints |
| `no`  | ignore van der Waals constraints |

------------------------------------------------------------------------

# Torsion-angle intervals

Torsion angles are derived from the PDB structure and converted into
intervals.

Given a reference torsion angle $\tau_{i}$, taken from the PDB structure, a perturbed value is
sampled as

$$
\tau_i^* \sim \mathcal{N}\left(\tau_i,\left(\frac{\Delta\tau_i}{8}\right)^2\right)
$$

and the resulting interval is

$$
\left[
\tau_i^* - \frac{\Delta\tau_i}{2},
\
\tau_i^* + \frac{\Delta\tau_i}{2}
\right],
$$

where `torsion_angle_width` defines the total interval width. This corresponds to $\Delta\tau = 8\sigma$

------------------------------------------------------------------------

## Backbone torsion selection

The percentage of backbone torsion angles ($\phi/\psi$) that will be included as
interval constraints is controlled by

    percentage_backbone_torsion_angles: 100.0

| value | behavior |
|------|----------|
| `100` | all backbone torsion angles are used |
| `<100` | a random subset of torsion angles is used |

Angles that are **not selected** receive the default range

$$
\left[
\tau_i^* - \frac{\Delta\tau_i}{2},
\
\tau_i^* + \frac{\Delta\tau_i}{2}
\right],
$$
------------------------------------------------------------------------

## ⚙️ Installation instructions

Clone the repository:

    git clone https://github.com/wdarocha/pdb-parser.git
    cd pdb-parser

Install in editable mode:

    pip install -e .

------------------------------------------------------------------------

## ▶️ Usage

After configuring the parameters in `data/params.cfg` and listing the
PDB identifiers in `data/pdb_ids.txt`, run the parsing pipeline.

Example:

``` 
python3 src/pdb_parser/pdb_parser.py data/pdb_ids.txt data/params.cfg data/pdb data/outputs

```

------------------------------------------------------------------------

## 📖 Citation

If this code is useful in your research, please cite the preprint below (also available via the **“Cite this repository”** button on GitHub thanks to the included `CITATION.cff`):

**W. da Rocha, C. Lavor, L. Liberti, L. de Melo Costa, L. D. Secchin, T. E. Malliavin.**  
*An Angle-Based Algorithmic Framework for the Interval Discretizable Distance Geometry Problem.*  
**arXiv:2508.09143**, 2025.  
https://arxiv.org/abs/2508.09143

### BibTeX
```bibtex
@misc{darocha2025,
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
© 2025 Wagner Alan Aparecido da Rocha

------------------------------------------------------------------------

## 👤 Author

Developed and maintained by [Wagner Alan Aparecido da Rocha](https://github.com/wdarocha).

------------------------------------------------------------------------

## 🙏 Acknowledgments

Special thanks to [Leonardo D.Secchin](https://github.com/leonardosecchin) for the valuable support provided during the development of the code.
