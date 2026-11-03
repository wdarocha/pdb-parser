# Pipeline Wrapper Guide

The repository provides `run_pipeline.sh` as a source-tree wrapper around the
public parser and remap CLIs.

## Default command

Without `--mode`, the script runs the parser and then the remap stage (`--mode full`):

```bash
./run_pipeline.sh data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
```

## Positional arguments

The script always accepts the same six positional arguments:

| argument | description |
|----------|-------------|
| `data/pdb_ids.txt` | text file containing the list of PDB identifiers |
| `data/params.cfg` | configuration file for the parser stage |
| `data/instance_remap.cfg` | configuration file for the remap stage |
| `data/pdb` | directory where downloaded PDB structures will be stored |
| `data/seeds` | directory where per-PDB seed files will be read and written |
| `data/outputs` | directory where parser and remapped outputs are written |

## Modes

```bash
./run_pipeline.sh --mode full data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
./run_pipeline.sh --mode parser data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
./run_pipeline.sh --mode remap data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
./run_pipeline.sh --mode remap-all data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
```

- `full`: runs parser and then remap with the configured `order_id`
- `parser`: runs only the parser stage
- `remap`: runs only the remap stage with the configured `order_id`
- `remap-all`: runs the remap stage repeatedly for all order ids selected by `--orders`

The underlying CLIs continue processing after an individual PDB failure, but
return `1` if any PDB fails. The wrapper therefore also stops with a nonzero
status for a failed stage; configuration errors return `2`.
Unexpected ordinary exceptions are reported as per-PDB processing failures;
interrupts and explicit process exits still abort immediately.

This is the `run_cli_jobs` contract: `FileNotFoundError` and `ValueError` are
reported as data issues, `RuntimeError` and other `OSError` subclasses as
processing errors, and every other ordinary `Exception` as an unexpected
processing error. Each failure is associated with its PDB identifier and the
remaining identifiers are still attempted. The function returns the number
of failed jobs, which the public CLIs map to exit status `1`.

## `--orders`

`remap-all` defaults to `1-10`. You can pass fixed values, ranges, or both:

```bash
./run_pipeline.sh --mode remap-all --orders 1,3,9 data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
./run_pipeline.sh --mode remap-all --orders 1-4,9 data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
```

The script writes the expanded list (for example `order_id: 1,2,3,4,9`) in a
temporary copy of `instance_remap.cfg` and runs a single remap call, which reads
the parser files once per structure and writes every order. With `--log-dir`,
this call logs to `remap_all.log`.

`--orders` has no effect for `order_type: dmdgp` (`backbone` outputs): the list is
checked, and the single DMDGP order is written once. Use `--mode remap` there.

## Intermediate files

`--keep-tmp` is passed to the parser stage, which then keeps the intermediate
files in `<out_dir>/<pdb_id>/tmp/` instead of removing them:

```bash
./run_pipeline.sh --keep-tmp data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
```

## Logging

To write per-stage logs instead of streaming command output directly:

```bash
./run_pipeline.sh --mode remap-all --log-dir logs data/pdb_ids.txt data/params.cfg data/instance_remap.cfg data/pdb data/seeds data/outputs
```
