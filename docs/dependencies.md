# Dependency compatibility policy

This project supports Python 3.10 or newer; CI runs the full check on Python 3.10,
3.11, 3.12 and 3.13. Runtime dependencies are bounded
to prevent an unreviewed major release from changing numerical behavior or
the file-reading APIs used by the parser:

| Dependency | Supported range | Role |
| --- | --- | --- |
| NumPy | `>=1.23,<3` | arrays, random generators, numerical operations |
| pandas | `>=1.5,<3` | tabular X/I/A processing |
| SciPy | `>=1.9,<2` | sparse shortest-path tightening |
| MDAnalysis | `>=2.4,<3` | PDB model, chain, and coordinate access |
| requests | `>=2.28,<3` | PDB downloads |

The ranges are declared in `pyproject.toml`, the single packaging metadata
source used by the build frontend. A dependency
upgrade requires running the full test suite, the coverage gate (currently
70%), shell syntax checks, and the parser/remapper integration tests. Changes
to numerical behavior or generated X/I/A/T formats require an explicit
compatibility review before widening a range.

The lower bounds represent the oldest supported API baseline. The upper
bounds are exclusive major-version bounds; patch and minor updates within a
range are accepted after the normal verification checks. No dependency is
vendored, and production installations should use a lock or constraints file
owned by the deploying application when reproducibility beyond this policy is
required.
