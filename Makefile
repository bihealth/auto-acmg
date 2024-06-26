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
	@echo "  run			 Run the application"
	@echo "  test-remote     Run remote tests"
	@echo "  test            Run tests"
	@echo "  test-all        Run all tests"
	@echo "  ci-unit-test    Run unit tests in CI"
	@echo "  ci-e2e-test     Run end-to-end tests in CI"
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
	format-black \
	format-isort

.PHONY: format-isort
format-isort:
	pipenv run isort --profile=black --line-length 100 $(DIRS_PYTHON)

.PHONY: format-black
format-black:
	pipenv run black --line-length 100 $(DIRS_PYTHON)

.PHONY: lint
lint: \
	lint-black \
	lint-isort \
	lint-flake8 \
	lint-mypy

.PHONY: lint-isort
lint-isort:
	pipenv run isort --profile=black --line-length 100 --check-only --diff $(DIRS_PYTHON)

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

.PHONY: test-remote
test-remote:
	pipenv run pytest \
		-m "remote" \
		tests/

.PHONY: test
test:
	pipenv run pytest \
		-m "not remote" \
		tests/

.PHONY: test-all
test-all:
	pipenv run pytest \
		tests/

.PHONY: ci-unit-test
ci-unit-test:
	pipenv run pytest \
		-m "not remote" \
		--cov-report term-missing \
		--cov-report lcov \
		--cov=src \
		tests/

.PHONY: ci-e2e-test
ci-e2e-test:
	pipenv run pytest \
		-m "remote" \
		--capture=no \
		tests/

.PHONY: docs
docs:
	PYTHONPATH=$(PWD) pipenv run -- make -C docs clean html

.PHONY: ci-docs
ci-docs:
	make -C docs clean html 2>&1 | tee sphinx-output.log | grep -q "ERROR" && exit 1 || exit 0
