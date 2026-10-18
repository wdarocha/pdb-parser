"""Compatibility adapter for legacy ``python setup.py`` invocations.

All project metadata is declared in ``pyproject.toml``. This adapter only
translates that canonical PEP 621 data for older setuptools versions that do
not read the ``[project]`` table when build isolation is disabled.
"""

from pathlib import Path

from setuptools import find_packages, setup

try:
	import tomllib
except ModuleNotFoundError:  # Python 3.10
	import tomli as tomllib


def _legacy_kwargs() -> dict[str, object]:
	"""Read canonical project metadata for a legacy setuptools invocation.

	``pyproject.toml`` is the sole source for distribution metadata, runtime
	dependencies, extras, and console scripts. The version remains dynamic and
	is read from the attribute named by its canonical setuptools configuration.
	"""
	root = Path(__file__).parent
	project_file = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
	project = project_file["project"]
	version_attr = project_file["tool"]["setuptools"]["dynamic"]["version"]["attr"]
	return {
		"name": project["name"],
		"version": _read_package_attribute(root, version_attr),
		"description": project["description"],
		"long_description": (root / project["readme"]).read_text(encoding="utf-8"),
		"long_description_content_type": "text/markdown",
		"python_requires": project["requires-python"],
		"install_requires": project["dependencies"],
		"extras_require": project.get("optional-dependencies", {}),
		"entry_points": {"console_scripts": [
			f"{name}={target}" for name, target in project["scripts"].items()
		]},
		"package_dir": {"": "src"},
		"packages": find_packages(where="src", include=["pdb_parser", "pdb_parser.*"]),
		"include_package_data": False,
	}


def _read_package_attribute(root: Path, dotted_attribute: str) -> str:
	"""Read a simple package attribute named by setuptools' dynamic config."""
	module_name, separator, attribute = dotted_attribute.rpartition(".")
	if not separator or not module_name or not attribute:
		raise RuntimeError(f"Invalid dynamic package attribute: {dotted_attribute!r}")
	module_path = root / "src" / Path(*module_name.split("."))
	if module_path.is_dir():
		module_path /= "__init__.py"
	else:
		module_path = module_path.with_suffix(".py")
	for line in module_path.read_text(encoding="utf-8").splitlines():
		if line.startswith(f"{attribute} = "):
			return line.split("=", 1)[1].strip().strip('"')
	raise RuntimeError(f"Unable to find {dotted_attribute}")


def _setuptools_reads_project_table() -> bool:
	"""Return whether this setuptools reads ``[project]`` from ``pyproject.toml`` (61+)."""
	import setuptools

	major = setuptools.__version__.split(".", 1)[0]
	return major.isdigit() and int(major) >= 61


if __name__ == "__main__":
	if _setuptools_reads_project_table():
		# Passing the same metadata again only makes setuptools warn that pyproject.toml
		# overrides it (an error under PYTHONWARNINGS=error).
		setup()
	else:
		setup(**_legacy_kwargs())
