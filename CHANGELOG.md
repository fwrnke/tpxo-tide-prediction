# tpxo-tide-prediction | Change Log
## [0.1.1] - 2022-02-27

Patch to fix missing consideration of tracks/regions extending over the prime meridian (zero longitude) resulting in similar error messages: `index XXXXX is out of bounds for axis 0 with size YYYYY`

### Added
- [CHANGELOG.md](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/CHANGELOG.md)

### Changed

- updated [`read_h_netCDFs()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L590) to handle new output from [`subset region()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L460)
- changed `offset` parameter default value from 10 to **3**

### Fixed

- [`subset region()`](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/tpxo_tide_prediction/utils.py#L460): now considering prime meridian for longitude subset

## [0.1] - 2022-02-11

First published version.
