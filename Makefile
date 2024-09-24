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
	@echo "  serve		     Run the API application"
	@echo "  bench           Run the benchmark"
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
	format-ruff \
	format-isort

.PHONY: format-ruff
format-ruff:
	pipenv run ruff format --line-length 100 $(DIRS_PYTHON)

.PHONY: format-isort
format-isort:
	pipenv run isort --profile=black --line-length 100 $(DIRS_PYTHON)

.PHONY: lint
lint: \
	lint-ruff \
	lint-isort \
	lint-mypy

.PHONY: lint-ruff
lint-ruff:
	pipenv run ruff check --line-length 100 $(DIRS_PYTHON)

.PHONY: lint-isort
lint-isort:
	pipenv run isort --profile=black --line-length 100 --check-only --diff $(DIRS_PYTHON)

.PHONY: lint-mypy
lint-mypy:
	pipenv run mypy --check-untyped-defs $(DIRS_PYTHON)

.PHONY: serve
serve:
	pipenv run uvicorn src.main:app --host 0.0.0.0 --port 8080 --reload --workers 8

.PHONY: bench
bench:
	pipenv run python -m src.bench.comparison_v4

.PHONY: test-remote
test-remote:
	pipenv run pytest \
		-m "remote" \
		--record-mode=new_episodes \
		tests/

.PHONY: test
test:
	pipenv run pytest \
		-m "not remote" \
		--record-mode=new_episodes \
		tests/

.PHONY: test-all
test-all:
	pipenv run pytest \
		--record-mode=new_episodes \
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
		--record-mode=new_episodes \
		--capture=no \
		tests/

.PHONY: docs
docs:
	PYTHONPATH=$(PWD) pipenv run -- make -C docs clean html

.PHONY: ci-docs
ci-docs:
	make -C docs clean html 2>&1 | tee sphinx-output.log | grep -q "ERROR" && exit 1 || exit 0

.PHONY: jupyternotebook
jupyternotebook:
	pipenv run \
		jupyter notebook

.PHONY: jupyterlab
jupyterlab:
	cp src/bench/results_analysis.ipynb tmp.ipynb && \
	PYTHON=. pipenv run \
		jupyter lab \
			--ip=0.0.0.0 --allow-root --NotebookApp.custom_display_url=http://127.0.0.1:8888 \
			tmp.ipynb

.PHONY: dump-openapi
dump-openapi:
	pipenv run python -c 'from src.utils import dump_openapi_yaml; dump_openapi_yaml("openapi.yaml")'

.PHONY: ci-dump-openapi
ci-dump-openapi:
	pipenv run python -c 'from src.utils import dump_openapi_yaml; dump_openapi_yaml("tmp/openapi.yaml")'
