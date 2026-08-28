# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.0] - 2026-08-28

### Added

- Platform-specific compiler flags.
- Expanded documentation.

### Fixed

- Fixed C++ __restrict keyword to avoid MSVC compiler issues on Windows.

### Changed

- SciPy is now required for a basic install.
- Removed matplotlib from full install dependencies since GERBLS does not call it directly.
- Removed the SavGol filter wrapper from the API to avoid endorsing a specific detrender.

## [0.7.1] - 2025-08-06

### Added

- Significantly improved documentation.

## [0.7.0] - 2025-07-21

### Added

- Automatic period-dependent data downsampling.
- Numerous minor improvements.

## [0.6.2] - 2025-05-28

### Added

- Added an option to manipulate maximum tested transit duration.

## [0.6.1] - 2025-05-23

### Added

- gerbls.bls_run() now produces more useful output.

## [0.6.0] - 2025-05-22

### Added

- Initial public release of GERBLS enabling basic usage.