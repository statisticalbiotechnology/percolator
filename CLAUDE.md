# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Percolator is a C++17 command-line application for postprocessing shotgun proteomics data using semi-supervised machine learning (SVM-based). It produces PSM (Peptide Spectrum Match) scores and FDR estimates. Key executables: `percolator`, `qvality`, and format converters (e.g., `sqt2pin`).

## Build Commands

Dependencies (Ubuntu): `libboost-filesystem-dev`, `libboost-system-dev`, `libboost-thread-dev`. Eigen and GoogleTest are fetched automatically via CMake FetchContent.

```bash
# Standard build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DGOOGLE_TEST=1 ..
make -j$(nproc)

# Quick build script (wraps cmake+make for all platforms)
./quickbuild.sh -s <source_dir> -b <build_dir> -r <release_dir>

# Full build script with installation of dependencies
admin/builders/xx_build.sh -s <parent dir to source_dir> -b <build_dir> -r <release_dir>
# Example: ./admin/builders/ubuntu64_build.sh -s .. -b ../build -r ../release

```

## Testing

```bash
# Run all tests (from build directory)
make test
# or
ctest -V

# Run a single unit test binary directly
./src/UnitTest_Percolator   # adjust binary name as needed
```

Unit tests live in `tests/unit_tests/percolator/` (GoogleTest, files named `UnitTest_*.cpp`).  
System/integration tests live in `data/system_tests/percolator/` (Python scripts).

## Architecture

### Core Processing Pipeline (`src/`)

The main flow: `main.cpp` → `Caller` → `CrossValidation` → model training → scoring → FDR estimation.

- **`Caller`** (`Caller.cpp/h`) — Top-level orchestration: parses options, loads data, drives cross-validation, writes output.
- **`CrossValidation`** (`CrossValidation.cpp/h`) — k-fold cross-validation loop for SVM training.
- **`Scores` / `ScoreHolder`** — PSM score storage and manipulation.
- **`SetHandler` / `DataSet` / `PSMDescription`** — Data loading and feature management. Input is PIN (tab-delimited Peptide-spectrum match Interchange Format).
- **`svm.cpp`** — Integrated libsvm for SVM training (large, ~80KB, treat as third-party).

### Regression / Calibration

Recent work has focused on PEP (Posterior Error Probability) calibration models:
- **`ISplineTRRRegressor`** (`ISplineTRRRegressor.h`) — Constrained isotonic spline regression using PAVA (Pool Adjacent Violators Algorithm). Uses Eigen for linear algebra.
- **`IsotonicPEP`** (`IsotonicPEP.cpp/h`) — Isotonic regression for PEP estimation.
- **`MonotoneRegressor`** — Earlier monotone regression model.
- **`LogisticRegression`** — Logistic regression model.
- **`PosteriorEstimator`** — Orchestrates PEP/q-value computation after scoring.

### Protein-Level Inference

`src/picked_protein/` — Protein inference module (separate subdirectory with its own CMakeLists.txt).  
`ProteinProbEstimator` — Protein-level FDR estimation.

### Normalizers

`NoNormalizer`, `StdvNormalizer`, `UniNormalizer` — Feature normalization strategies selected at runtime.

### Converters

`src/converters/` — Format conversion tools (SQT, MS2, etc. → PIN format).

## C++ Style Guide Summary

Follow `docs/StyleGuide.md`. Key rules:
- **Naming**: `camelCase` variables/functions, `PascalCase` classes, `camelCase_` member variables (trailing underscore), `ALL_CAPS` constants, `lowercase` namespaces. Struct members do NOT get the trailing underscore.
- **Formatting**: 2-space indentation, 80-char line limit, opening brace on same line.
- **Memory**: Prefer `std::unique_ptr`/`std::shared_ptr`; avoid raw `new`/`delete`.
- **Error handling**: Use exceptions, not return codes.
- **Functions**: Keep under ~40 lines; pass objects by `const&`.
- **Include guards**: Use traditional `#ifndef` guards (not `#pragma once`).
- **File names**: PascalCase (e.g., `MyClass.cpp`).

Note: existing files have not all been reformatted yet — consistency within a file matters more than enforcing the guide on untouched code.

## Release Versioning

Version is set in `CommonCMake.txt` as MAJOR.MINOR.PATCH (e.g., `v3.08.01`). Release tags use format `rel-3-08-01`. Current version: `v3.08.01`.
