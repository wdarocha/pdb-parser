from __future__ import annotations

from pdb_parser.cli_support import (
	print_config_error,
	print_error,
	print_job_banner,
	print_skip,
	run_cli_jobs,
)


def test_print_error_writes_to_stderr(capsys) -> None:
	print_error("bad input")
	captured = capsys.readouterr()
	assert captured.out == ""
	assert captured.err == "[ERROR] bad input\n"


def test_print_config_error_writes_expected_stderr(capsys) -> None:
	print_config_error("missing order_id")
	captured = capsys.readouterr()
	assert captured.out == ""
	assert captured.err == "[ERROR] Configuration error: missing order_id\n"


def test_print_skip_writes_expected_stdout(capsys) -> None:
	print_skip("1ABC", "missing file")
	captured = capsys.readouterr()
	assert captured.out == "[SKIP] 1ABC: missing file\n"
	assert captured.err == ""


def test_print_job_banner_writes_three_lines(capsys) -> None:
	print_job_banner("1ABC")
	captured = capsys.readouterr()
	assert captured.out.splitlines() == [
		"------------------------------------------------------",
		"------------------------ 1ABC ------------------------",
		"------------------------------------------------------",
	]
	assert captured.err == ""


def test_run_cli_jobs_runs_all_jobs_and_prints_banners(capsys) -> None:
	processed: list[str] = []

	def run_job(pdb_id: str) -> None:
		processed.append(pdb_id)

	assert run_cli_jobs(["1ABC", "2XYZ"], run_job) == 0

	captured = capsys.readouterr()
	assert processed == ["1ABC", "2XYZ"]
	assert "[SKIP]" not in captured.out
	assert "1ABC" in captured.out
	assert "2XYZ" in captured.out


def test_run_cli_jobs_converts_exceptions_to_skip_messages(capsys) -> None:
	def run_job(pdb_id: str) -> None:
		if pdb_id == "2XYZ":
			raise ValueError("missing file")

	assert run_cli_jobs(["1ABC", "2XYZ"], run_job) == 1

	captured = capsys.readouterr()
	assert "[SKIP] 2XYZ: Data issue: missing file" in captured.out
	assert captured.err == ""


def test_run_cli_jobs_reports_unexpected_exceptions_and_continues(capsys) -> None:
	processed: list[str] = []

	def run_job(pdb_id: str) -> None:
		processed.append(pdb_id)
		if pdb_id == "2XYZ":
			raise KeyError("unexpected failure")

	assert run_cli_jobs(["1ABC", "2XYZ", "3LMN"], run_job) == 1

	captured = capsys.readouterr()
	assert processed == ["1ABC", "2XYZ", "3LMN"]
	assert "[SKIP] 2XYZ: Unexpected processing error: 'unexpected failure'" in captured.out
	assert captured.err == ""


def test_cli_support_uses_typed_parser_contracts() -> None:
	from pdb_parser.cli_support import CliArgs, CliParser, CliParserFactory

	assert CliArgs.__name__ == "Namespace"
	assert CliParser is not None
	assert CliParserFactory is not None
