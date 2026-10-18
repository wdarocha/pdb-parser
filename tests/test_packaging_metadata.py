from __future__ import annotations

import runpy
import subprocess
import sys
import zipfile
from email.parser import Parser
from pathlib import Path

from packaging.requirements import Requirement

from pdb_parser import __version__

try:
	import tomllib
except ModuleNotFoundError:  # Python 3.10
	import tomli as tomllib


def _requirement_key(raw_requirement: str) -> str:
	requirement = Requirement(raw_requirement)
	marker = f"; {requirement.marker}" if requirement.marker else ""
	return requirement.name.lower() + str(requirement.specifier) + marker


def test_manifest_excludes_non_package_content_from_source_distribution() -> None:
	manifest_text = Path("MANIFEST.in").read_text(encoding="utf-8")

	assert "prune tests" in manifest_text
	assert "prune data" in manifest_text
	assert "recursive-include src *.py" in manifest_text


def test_gitignore_covers_local_python_artifacts_and_runtime_outputs() -> None:
	gitignore_text = Path(".gitignore").read_text(encoding="utf-8")

	assert "__pycache__/" in gitignore_text
	assert "*.egg-info/" in gitignore_text
	assert "data/outputs/" in gitignore_text
	assert "data/seeds/" in gitignore_text
	assert "data/pdb/" in gitignore_text


def test_makefile_exposes_local_verification_targets() -> None:
	makefile_text = Path("Makefile").read_text(encoding="utf-8")

	assert ".PHONY: install-dev test coverage lint check clean" in makefile_text
	assert "\ntest:\n" in makefile_text
	assert "\nlint:\n" in makefile_text
	assert "\ncheck: lint test\n" in makefile_text
	assert "\nclean:\n" in makefile_text


def test_wheel_publishes_declared_entry_points_and_dependencies(tmp_path: Path) -> None:
	repo_root = Path(__file__).resolve().parents[1]
	wheel_dir = tmp_path / "wheel"
	subprocess.run(
		[
			sys.executable,
			"-m",
			"pip",
			"wheel",
			"--no-deps",
			"--no-build-isolation",
			"--wheel-dir",
			str(wheel_dir),
			str(repo_root),
		],
		check=True,
		cwd=repo_root,
	)
	wheel_path = next(wheel_dir.glob("pdb_parser-*.whl"))
	with zipfile.ZipFile(wheel_path) as wheel:
		metadata_name = next(name for name in wheel.namelist() if name.endswith("/METADATA"))
		entry_points_name = next(
			name for name in wheel.namelist() if name.endswith("/entry_points.txt")
		)
		metadata = Parser().parsestr(wheel.read(metadata_name).decode("utf-8"))
		entry_points = wheel.read(entry_points_name).decode("utf-8")

	project = tomllib.loads((repo_root / "pyproject.toml").read_text(encoding="utf-8"))["project"]
	assert metadata["Name"] == project["name"]
	assert metadata["Version"] == __version__
	assert metadata["Summary"] == project["description"]
	assert metadata["Requires-Python"] == project["requires-python"]
	published_dependencies = {
		_requirement_key(dependency)
		for dependency in metadata.get_all("Requires-Dist", [])
	}
	declared_dependency_values = list(project["dependencies"])
	for extra, dependencies in project.get("optional-dependencies", {}).items():
		declared_dependency_values.extend(
			f'{dependency}; extra == "{extra}"' for dependency in dependencies
		)
	declared_dependencies = {
		_requirement_key(dependency)
		for dependency in declared_dependency_values
	}
	assert published_dependencies == declared_dependencies
	declared_entry_points = {
		f"{name} = {target}" for name, target in project["scripts"].items()
	}
	published_entry_points = {
		line.strip()
		for line in entry_points.splitlines()
		if line.strip() and not line.startswith("[")
	}
	assert published_entry_points == declared_entry_points


def test_git_does_not_track_generated_python_artifacts() -> None:
	repo_root = Path(__file__).resolve().parents[1]
	result = subprocess.run(
		[
			"git",
			"ls-files",
			"src/pdb_parser/**/__pycache__/*",
			"tests/__pycache__/*",
			"src/pdb_parser.egg-info/*",
		],
		check=True,
		cwd=repo_root,
		capture_output=True,
		text=True,
	)

	assert result.stdout.strip() == ""


def test_operational_docs_are_present_and_linked_from_readme() -> None:
	readme_text = Path("README.md").read_text(encoding="utf-8")

	for relative_path in (
		"docs/parser_stage.md",
		"docs/remap_stage.md",
		"docs/pipeline_wrapper.md",
	):
		assert Path(relative_path).exists()
		assert relative_path in readme_text


def test_legacy_setup_adapter_uses_the_same_project_metadata() -> None:
	repo_root = Path(__file__).resolve().parents[1]
	project = tomllib.loads((repo_root / "pyproject.toml").read_text(encoding="utf-8"))
	setup_module = runpy.run_path(str(repo_root / "setup.py"), run_name="_metadata_test")
	legacy = setup_module["_legacy_kwargs"]()
	canonical = project["project"]
	assert legacy["name"] == canonical["name"]
	assert legacy["description"] == canonical["description"]
	assert legacy["python_requires"] == canonical["requires-python"]
	assert legacy["install_requires"] == canonical["dependencies"]
	assert legacy["extras_require"] == canonical["optional-dependencies"]
	assert legacy["entry_points"]["console_scripts"] == [
		f"{name}={target}" for name, target in canonical["scripts"].items()
	]
