# pdb-parser

Lightweight Python toolkit for parsing **Protein Data Bank (PDB)**
structures and generating **distance and torsion-angle constraints**
tailored for **Distance Geometry Problem (DGP)** workflows.

The package extracts structural information from PDB files and produces
constraint sets that can be used in algorithms for protein structure
determination, including **Discretizable Distance Geometry Problem (DDGP
/ iDDGP)** frameworks.

Typical applications include:

-   generation of **distance-constraint lists**
-   extraction of **backbone torsion angles**
-   preparation of datasets for **branch-and-prune distance geometry
    algorithms**
-   preprocessing of **NMR-derived structural information**

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

## 📂 Input configuration

The parser behavior is controlled by a configuration file located in

    data/params.cfg

This file specifies:

-   which **model** from the PDB file will be used
-   which **chain** will be processed
-   the **atom selection strategy**
-   the **distance-constraint model**
-   parameters for synthetic or experimental **NMR distance intervals**
-   torsion-angle interval generation

Example:

    model_number: 1
    chain_id: A
    atom_selection: backbone_plus_hydrogens
    distance_constraints: interval_centered

------------------------------------------------------------------------

## PDB list

The file

    data/pdb_ids.txt

contains the list of PDB identifiers to process.

Example:

    1TOS
    1UAO

Each entry corresponds to a structure that will be retrieved and
processed by the pipeline.

------------------------------------------------------------------------

## Distance constraint models

The parser supports multiple strategies for constructing distance
constraints.

## Precise distances

Distances extracted directly from the PDB structure.

    distance_constraints: precise

## Synthetic intervals

Synthetic intervals centered around the reference distance.

$$
\mathcal{D}_{ij} =
\left[
\max\left\{d_{ij}^* - \frac{\varepsilon_{ij}}{2},\ v_{\mathrm{dw}}\right\},
\;
\min\left\{d_{ij}^* + \frac{\varepsilon_{ij}}{2},\ d_{\max}\right\}
\right]
$$

where $\varepsilon_{ij}$ corresponds to the parameter `interval_width`.

## Experimental NOE intervals

Distance bounds derived from NOESY peak intensity classes:

  NOE class   upper bound
  ----------- -------------
  strong      2.5 Å
  medium      3.5 Å
  weak        5.0 Å

------------------------------------------------------------------------

## Torsion-angle intervals

Torsion angles are derived from the PDB structure and converted into
intervals.

Given a reference torsion angle $\tau_i$, a perturbed value is sampled as

$$
\tau_i^* \sim \mathcal{N}\left(\tau_i,\left(\frac{\omega}{8}\right)^2\right)
$$

and the resulting interval is

$$
\left[
\tau_i^* - \frac{\Delta\tau_i}{2},
\;
\tau_i^* + \frac{\Delta\tau_i}{2}
\right],
$$

where $\Delta\tau_i$ corresponds to the parameter `torsion_angle_width`.

------------------------------------------------------------------------

## ⚙️ Build Instructions

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

``` python
from pdb_parser.pipeline import parser_pipeline

parser_pipeline.run()
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
