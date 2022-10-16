# tpxo-tide-prediction - Changelog

## v0.2.3 (2022-10-17)
Release for Zenodo DOI. 
No futher code changes.

## v0.2.2 (2022-06-07)

Minor bugfix.

### Added

- new parameter file *params_track_between_grid_nodes.txt*

### Changed

- `utils.py`
   - `subset_region()`:
      - fix indexing when both minimum and maximum longitude are in between two grid nodes of loaded constituent netCDF files

## v0.2.1 (2022-06-02)

Minor bugfixes.

### Changed

- `predict_tide.py`
   - `read_parameter_file()`
      - fix error after parsing single time input from command line
- `utils.py`
   - `subset_region()`:
      - selecting only single position

## v0.2.0 (2022-04-24)

Fixed several bugs occuring for files crossing _dateline_, *equator*, _prime meridian_, or *north pole*.
Changed interpolation method from _cubic_ to _linear_ (according to [xarray.Dataset.interp]([xarray.Dataset.interp](https://docs.xarray.dev/en/stable/generated/xarray.Dataset.interp.html))).

### Added

- new parameter example files (covering special cases, e.g. *dateline* crossings)
- added test for `longitude_to_180` and `test_longitude_to_360` functions

### Changed

- `predict_tide.py`
   - changed default interpolation method to *linear*
   - raise `MemoryError` if interpolated tide array is too huge to allocate
   - use **only** associated timestamp for tide predition if `mode = 'track'` (more efficient!)
- `utils.py`
   - `subset_region()`:
      - reworked indexing method
      - changed default number of nodes to pad subset region from 10 to **3**
   - `read_h_netCDFs()`:
      - changed subset selection from _labeled_ (`.sel`) to _integer-based_ (`.isel`)
      - changed default number of nodes to pad subset region from 10 to **3**

## v0.1.1 (2022-02-27)

Patch to fix missing consideration of tracks/regions extending over the prime meridian (zero longitude) resulting in similar error messages: `index XXXXX is out of bounds for axis 0 with size YYYYY`

### Added

- [CHANGELOG.md](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/CHANGELOG.md)

### Changed

- updated [`read_h_netCDFs()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L590) to handle new output from [`subset region()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L460)
- changed `offset` parameter default value from 10 to **3**

### Fixed

- [`subset region()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L460): now considering prime meridian for longitude subset

## v0.1.0 (2022-02-11)

First published version.
