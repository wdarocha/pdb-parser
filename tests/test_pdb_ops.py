from __future__ import annotations

from pathlib import Path

import pytest
import requests

from pdb_parser.io.pdb_ops import (
	detect_experimental_method,
	download_pdb,
	extract_model_chain,
	is_nmr_structure,
	list_chains_for_model,
)


@pytest.mark.parametrize(("content", "expected"), [
	("EXPDTA NMR\nEXPDTA NMR\n", "nmr"),
	("EXPDTA X-RAY DIFFRACTION\nEXPDTA X-RAY DIFFRACTION\n", "xray"),
])
def test_detect_experimental_method_is_deterministic_for_repeated_records(
	tmp_path: Path, content: str, expected: str
) -> None:
	path = tmp_path / "method.pdb"
	path.write_text(content, encoding="utf-8")
	assert detect_experimental_method(path) == expected


def test_detect_experimental_method_rejects_combined_supported_methods(tmp_path: Path) -> None:
	path = tmp_path / "method.pdb"
	path.write_text("EXPDTA NMR\nEXPDTA X-RAY DIFFRACTION\n", encoding="utf-8")
	with pytest.raises(ValueError, match="Ambiguous EXPDTA"):
		detect_experimental_method(path)


def test_detect_experimental_method_rejects_unsupported_method(tmp_path: Path) -> None:
	path = tmp_path / "method.pdb"
	path.write_text("EXPDTA ELECTRON MICROSCOPY\n", encoding="utf-8")
	with pytest.raises(ValueError, match="Unsupported or missing EXPDTA"):
		detect_experimental_method(path)


def test_download_pdb_rejects_invalid_identifier(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="Invalid PDB id"):
		download_pdb("BAD-ID", tmp_path)


def test_download_pdb_retries_a_network_failure_and_reports_it(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	calls, sleeps = [], []

	def fake_get(*args, **kwargs):
		calls.append(args)
		raise requests.RequestException("network down")

	monkeypatch.setattr("pdb_parser.io.pdb_ops.requests.get", fake_get)
	monkeypatch.setattr("pdb_parser.io.pdb_ops.time.sleep", sleeps.append)

	# A network failure is a processing error (OSError), not a data issue.
	with pytest.raises(ConnectionError, match="Failed to download PDB '1ABC' from RCSB after 3 attempts"):
		download_pdb("1ABC", tmp_path)
	assert len(calls) == 3
	assert sleeps == [2.0, 5.0]


def test_download_pdb_recovers_from_a_transient_timeout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
	class Response:
		text = "HEADER    TEST\nEXPDTA    X-RAY DIFFRACTION\nATOM      1  N   ALA A   1       0.000   0.000   0.000\n"

		def raise_for_status(self) -> None:
			return None

	answers = [requests.Timeout("read timed out"), Response()]

	def fake_get(*args, **kwargs):
		answer = answers.pop(0)
		if isinstance(answer, Exception):
			raise answer
		return answer

	monkeypatch.setattr("pdb_parser.io.pdb_ops.requests.get", fake_get)
	monkeypatch.setattr("pdb_parser.io.pdb_ops.time.sleep", lambda seconds: None)

	assert download_pdb("1ABC", tmp_path).read_text(encoding="utf-8") == Response.text


def test_download_pdb_does_not_retry_an_unknown_id(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
	calls = []

	class NotFound:
		status_code = 404

		def raise_for_status(self) -> None:
			raise requests.HTTPError("404 Client Error: Not Found", response=self)

	def fake_get(*args, **kwargs):
		calls.append(args)
		return NotFound()

	monkeypatch.setattr("pdb_parser.io.pdb_ops.requests.get", fake_get)

	with pytest.raises(ValueError, match="Failed to download PDB '9ZZZ' from RCSB: 404"):
		download_pdb("9ZZZ", tmp_path)
	assert len(calls) == 1


def test_download_pdb_rejects_invalid_response_without_replacing_cache(
	tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
	cache = tmp_path / "1ABC.pdb"
	cache.write_text("old-valid-cache", encoding="utf-8")

	class Response:
		text = "HTML error page"

		def raise_for_status(self) -> None:
			return None

	monkeypatch.setattr("pdb_parser.io.pdb_ops.requests.get", lambda *args, **kwargs: Response())
	with pytest.raises(ValueError, match="invalid or incomplete"):
		download_pdb("1ABC", tmp_path)
	assert cache.read_text(encoding="utf-8") == "old-valid-cache"


def test_is_nmr_structure_returns_false_without_expdta_nmr(tmp_path: Path) -> None:
	pdb_file = tmp_path / "sample.pdb"
	pdb_file.write_text("HEADER TEST\nEXPDTA X-RAY DIFFRACTION\n", encoding="utf-8")

	assert is_nmr_structure(pdb_file) is False


def test_list_chains_for_model_rejects_out_of_range_model(tmp_path: Path) -> None:
	pdb_file = tmp_path / "sample.pdb"
	pdb_file.write_text(
		(
			"MODEL        1\n"
			"ATOM      1  N   ALA A   1      11.104  13.207   9.447  1.00 20.00           N\n"
			"ATOM      2  CA  ALA A   1      12.560  13.201   9.699  1.00 20.00           C\n"
			"ATOM      3  C   ALA A   1      13.113  14.609   9.981  1.00 20.00           C\n"
			"TER\n"
			"ENDMDL\n"
			"END\n"
		),
		encoding="utf-8",
	)

	with pytest.raises(ValueError, match="model_number out of range"):
		list_chains_for_model(pdb_file, 2)


def test_extract_model_chain_reports_missing_chain_with_available_ids(tmp_path: Path) -> None:
	pdb_file = tmp_path / "sample.pdb"
	pdb_file.write_text(
		(
			"MODEL        1\n"
			"ATOM      1  N   ALA A   1      11.104  13.207   9.447  1.00 20.00           N\n"
			"ATOM      2  CA  ALA A   1      12.560  13.201   9.699  1.00 20.00           C\n"
			"ATOM      3  C   ALA A   1      13.113  14.609   9.981  1.00 20.00           C\n"
			"TER\n"
			"ENDMDL\n"
			"END\n"
		),
		encoding="utf-8",
	)

	with pytest.raises(ValueError, match="Available chain identifiers: A"):
		extract_model_chain(pdb_file, 1, "B")
