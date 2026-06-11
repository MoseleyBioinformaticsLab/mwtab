# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**mwtab** is a Python library for reading, writing, and validating files in the mwTab format used by the [Metabolomics Workbench](https://www.metabolomicsworkbench.org/) to archive Mass Spectrometry (MS) and NMR experimental data. It also provides a CLI (`mwtab`) and REST API integration for fetching studies.

## Commands

```bash
# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run a single test file
pytest tests/test_mwtab.py

# Run a specific test
pytest tests/test_mwtab.py::test_function_name

# Run tests with coverage
pytest --cov=mwtab

# Build docs
cd docs && make html

# Run the CLI
python -m mwtab --help
mwtab convert <path> <output>
mwtab validate <path>
mwtab download <study_id>
mwtab extract <path>
```

## Architecture

Source lives under `src/mwtab/`. The library is organized around a central data model with distinct layers for I/O, validation, and extraction.

### Data Flow

```
Input sources (file/URL/archive/REST)
  → fileio.py          — opens streams from local paths, URLs, .zip/.tar.gz
  → tokenizer.py       — low-level key-value pair generator from mwTab text
  → mwtab.py           — builds MWTabFile (nested dict-like object)
  → validator.py       — validates structure, consistency, metabolite data
  → converter.py       — serializes to/from JSON
  → mwextract.py       — exports metadata/metabolites as CSV or JSON
```

### Core Modules

| Module | Purpose |
|--------|---------|
| `mwtab.py` | `MWTabFile` — the central data structure. Behaves like a dict keyed on section names (e.g. `SUBJECT`, `METABOLITES`). Provides `study_id`, `analysis_id`, `header` properties. Handles reading/writing via `read()`/`write()`. |
| `validator.py` | Schema-based and consistency validation. Checks required fields, value formats, and cross-section consistency. Entry point: `MWTabFileValidator`. |
| `mwschema.py` | JSON schema definitions for MS and NMR file variants. Defines allowed sections, required keys, and NA value patterns. |
| `metadata_column_matching.py` | Regex-based matching for standard column names (e.g., `kegg_id`, `pubchem_id`) and value validation. Used by the validator. |
| `mwrest.py` | REST API wrapper for querying Metabolomics Workbench. `GenericMWURL` builds query URLs; results feed into `fileio`. |
| `fileio.py` | Unified file opener for local files, directories, URLs, and zip/tar archives. Returns generators of `MWTabFile` objects. |
| `tokenizer.py` | Parses the line-based mwTab text format into `(key, value)` tuples. Handles section start/end markers. |
| `converter.py` | Converts `MWTabFile` ↔ JSON. Handles single files and batch conversion across archives. |
| `mwextract.py` | Extracts metadata and metabolite tables from `MWTabFile` objects and writes CSV/JSON output. |
| `duplicates_dict.py` | Custom dict that preserves duplicate keys using a `{{{_N_}}}` suffix naming scheme. Used when parsing sections that allow repeated keys. |
| `cli.py` | docopt-based CLI. Commands map to `converter`, `validator`, `mwrest`, and `mwextract` functions. |

### MWTabFile Structure

An `MWTabFile` is a `dict`-like object. Top-level keys are section names from the mwTab format (e.g., `"SUBJECT"`, `"COLLECTION"`, `"MS_METABOLITE_DATA"`). Each section is itself a dict of key-value pairs parsed from the file. The `#SUBJECT_SAMPLE_FACTORS` block and data tables are stored as lists of dicts under special keys.

### Adding Validation Rules

New validation logic belongs in `validator.py`. Schema constraints (allowed values, required keys) belong in `mwschema.py`. Column name/value patterns for metabolite tables belong in `metadata_column_matching.py`.

### Extending the CLI

CLI commands are defined in `cli.py` using docopt. The docstring at the top of the file is the docopt spec. Each sub-command dispatches to the appropriate library module.

## Testing

- Tests live in `tests/`; fixtures and helpers in `tests/fixtures.py`.
- Example mwTab and JSON files used as test data are in `tests/example_data/`.
- CI runs pytest across Python 3.10–3.14 on Linux, Windows, and macOS.
- Coverage is tracked via CodeCov.

## Versioning

Version is managed by `setuptools_scm` from git tags. Do not manually edit `src/mwtab/_version.py` — it is auto-generated.
