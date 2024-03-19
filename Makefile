DIRS_PYTHON := src tests

.PHONY: help
help:
	@echo "Usage: make <target>"
	@echo
	@echo "Targets:"
	@echo "  help    This help (default target)"
	@echo "  deps    Install all dependencies"
	@echo "  format  Format source code"
	@echo "  lint    Run lint checks"
	@echo "  example_run    Run example"
	@echo "  test    Run tests"
	@echo "  ci      Install dependencies, run lints and tests"
	@echo "  docs    Generate the documentation"
	@echo "  mksuperuser Create a superuser"
	@echo "  serve   Run the (development) server"
	@echo "  jupyterlab Run jupyterlab"
	@echo "  celery  Run celery"
	@echo "  migrate Create alembic versions and upgrade"

.PHONY: deps
deps:
	pipenv install --dev

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
	MYPYPATH=$(PWD)/stubs pipenv run mypy $(DIRS_PYTHON)

.PHONY: example_run
example_run:
	pipenv run python -m src.main 4-113568536-G-GA --genome_release hg19

.PHONY: test
test:
	pipenv run pytest tests/

.PHONY: test-ci
test-ci:
	pipenv run pytest \
		--cov-report term-missing \
		--cov-report lcov \
		--cov=src \
		tests/

.PHONY: ci
ci: \
	deps \
	lint \
	test-ci

# .PHONY: docs
# docs:
# 	PYTHONPATH=$(PWD) pipenv run -- make -C ../docs clean html

# .PHONY: mksuperuser
# mksuperuser:
# 	PYTHONPATH=. pipenv run python app/backend_pre_start.py
# 	PYTHONPATH=. pipenv run python app/initial_data.py

# .PHONY: serve
# serve:
# 	pipenv run uvicorn app.main:app --host 0.0.0.0 --port 8080 --reload --workers 8

# .PHONY: celery
# celery:
# 	PYTHONPATH=. pipenv run \
# 		watchmedo auto-restart --directory=./ --pattern=*.py --recursive -- \
# 		celery -A app.worker worker --loglevel=debug --beat -Q main-queue

# .PHONY: jupyterlab
# jupyterlab:
# 	cp utils/minimal.ipynb tmp.ipynb && \
# 	PYTHON=. pipenv run \
# 		jupyter lab \
# 			--ip=0.0.0.0 --allow-root --NotebookApp.custom_display_url=http://127.0.0.1:8888 \
# 			tmp.ipynb
