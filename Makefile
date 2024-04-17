DIRS_PYTHON := src tests docs

.PHONY: help
help:
	@echo "Usage: make <target>"
	@echo
	@echo "Targets:"
	@echo "  help            This help (default target)"
	@echo "  deps            Install all dependencies"
	@echo "  ci-docs-deps    Install all dependencies for docs in CI"
	@echo "  format          Format source code"
	@echo "  lint            Run lint checks"
	@echo "  example_run     Run example"
	@echo "  run						 Run the application"
	@echo "  test            Run tests"
	@echo "  ci-test         Run tests in CI"
	@echo "  ci              Install dependencies, run lints and tests"
	@echo "  docs            Generate the documentation"
	@echo "  ci-docs		 Generate the documentation in CI"

.PHONY: deps
deps:
	pipenv install --dev

.PHONY: ci-docs-deps
ci-docs-deps:
	python -m pip install --upgrade --no-cache-dir pip setuptools
	python -m pip install --upgrade --no-cache-dir sphinx readthedocs-sphinx-ext
	python -m pip install --no-cache-dir -r docs/requirements.txt

.PHONY: format
format:	\
	format-isort \
	format-black

.PHONY: format-isort
format-isort:
	pipenv run isort --profile=black $(DIRS_PYTHON)

.PHONY: format-black
format-black:
	pipenv run black --line-length 100 $(DIRS_PYTHON)

.PHONY: lint
lint: \
	lint-isort \
	lint-black \
	lint-flake8 \
	lint-mypy

.PHONY: lint-isort
lint-isort:
	pipenv run isort --profile=black --check-only --diff $(DIRS_PYTHON)

.PHONY: lint-black
lint-black:
	pipenv run black --check --line-length 100 --diff $(DIRS_PYTHON)

.PHONY: lint-flake8
flake8:
	pipenv run flake8 --max-line-length 100 $(DIRS_PYTHON)

.PHONY: lint-mypy
lint-mypy:
	pipenv run mypy --check-untyped-defs $(DIRS_PYTHON)

#pipenv run python -m src.main 4-113568536-G-GA --genome_release hg19
.PHONY: example_run
example_run:
	pipenv run python -m src.cli "4-113568536-G-GA" --genome-release hg19

.PHONY: run
run:
ifdef GR
	pipenv run python -m src.cli "$(VAR)" --genome-release $(GR)
else
	pipenv run python -m src.cli "$(VAR)"
endif

.PHONY: test
test:
	pipenv run pytest tests/

.PHONY: ci-test
ci-test:
	pipenv run pytest \
		--cov-report term-missing \
		--cov-report lcov \
		--cov=src \
		tests/

.PHONY: ci
ci: \
	deps \
	lint \
	ci-test

.PHONY: docs
docs:
	PYTHONPATH=$(PWD) pipenv run -- make -C docs clean html

.PHONY: ci-docs
ci-docs:
	make -C docs clean html
