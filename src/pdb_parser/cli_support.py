"""Shared output helpers for the public CLI entry points."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Iterable, Sequence
from typing import Protocol, TypeVar

from pdb_parser.messages import print_error

CliContext = TypeVar("CliContext")
CliArgs = argparse.Namespace
CliJob = Callable[[str], None]


class CliParser(Protocol):
	"""Minimal parser interface required by the shared CLI runner."""

	def parse_args(self, args: Sequence[str] | None = None) -> CliArgs: ...


CliParserFactory = Callable[[], CliParser]


class ContextLoader(Protocol[CliContext]):
	"""Callable loader used by public CLI wrappers."""

	def __call__(self, args: CliArgs) -> CliContext: ...


class JobRunnerFactory(Protocol[CliContext]):
	"""Callable factory that binds one loaded context to a per-PDB job runner."""

	def __call__(self, context: CliContext) -> CliJob: ...


def print_config_error(message: str) -> None:
	"""Print a standardized configuration error message to stderr."""
	print_error(f"Configuration error: {message}")


def print_skip(pdb_id: str, message: str) -> None:
	"""Print a standardized skip message for one PDB identifier."""
	print(f"[SKIP] {pdb_id}: {message}")


def print_job_banner(pdb_id: str) -> None:
	"""Print a visual banner for one PDB job."""
	print("------------------------------------------------------")
	print(f"------------------------ {pdb_id} ------------------------")
	print("------------------------------------------------------")


def run_cli_jobs(
	pdb_ids: Iterable[str],
	run_job: CliJob,
) -> int:
	"""Run every job, reporting ordinary exceptions and continuing the batch.

	Expected data and processing errors keep their stable labels. Other
	``Exception`` instances are reported as unexpected processing errors for
	the affected PDB. ``KeyboardInterrupt`` and ``SystemExit`` still abort.
	"""
	failures = 0
	for pdb_id in pdb_ids:
		print_job_banner(pdb_id)
		try:
			run_job(pdb_id)
		except (FileNotFoundError, ValueError) as exc:
			failures += 1
			print_skip(pdb_id, f"Data issue: {exc}")
		except (RuntimeError, OSError) as exc:
			failures += 1
			print_skip(pdb_id, f"Processing error: {exc}")
		except Exception as exc:
			failures += 1
			print_skip(pdb_id, f"Unexpected processing error: {exc}")
	return failures


def run_cli_command(
	argv: list[str] | None,
	build_parser: CliParserFactory,
	load_context: ContextLoader[CliContext],
	build_job_runner: JobRunnerFactory[CliContext],
) -> int:
	"""Run one public CLI using shared parse/load/dispatch behavior."""
	arg_parser = build_parser()
	args = arg_parser.parse_args(argv)

	try:
		context = load_context(args)
	except (FileNotFoundError, ValueError, OSError) as exc:
		print_config_error(str(exc))
		return 2

	failures = run_cli_jobs(context["pdb_ids"], build_job_runner(context))
	return 1 if failures else 0
