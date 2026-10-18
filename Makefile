PYTHON ?= python3

.PHONY: install-dev test coverage lint check clean

install-dev:
	$(PYTHON) -m pip install -e .[dev]

test:
	$(PYTHON) -m pytest -q

coverage:
	$(PYTHON) -m pytest --cov=pdb_parser --cov-report=term-missing -q

lint:
	$(PYTHON) -m ruff check .

check: lint test
	$(PYTHON) -m pytest --cov=pdb_parser --cov-report=term-missing --cov-fail-under=70 -q

clean:
	find . -type d -name '__pycache__' -prune -exec rm -rf {} +
	find . -type d -name '*.egg-info' -prune -exec rm -rf {} +
	for generated_dir in data/outputs data/seeds data/pdb; do \
		if test -d "$$generated_dir"; then find "$$generated_dir" -mindepth 1 ! -name '.gitkeep' -exec rm -rf {} +; fi; \
	done
	rm -rf .pytest_cache .ruff_cache htmlcov .coverage build dist
