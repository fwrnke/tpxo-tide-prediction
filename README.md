[![DOI](https://zenodo.org/badge/423613235.svg)](https://zenodo.org/badge/latestdoi/423613235)

# tpxo-tide-prediction

This package allows the user to predict tidal elevations computed based on the [TPXO9-atlas](https://www.tpxo.net/global/tpxo9-atlas) models (**netCDF** format) provided by [Gary Egbert & Svetlana Erofeeva from the Oregon State University](https://www.tpxo.net/home) (on request for academic purposes).

## Disclaimer

This package was build to compute tidal corrections for *user-specified timestamps and geographical positions*. It uses the  latest 1/30 degree resolution fully global solution [TPXO9-atlas-v4/v5](https://www.tpxo.net/global/tpxo9-atlas) tidal model.

**NOTE: The `TPXO9-atlas-v*` netCDF files have to be downloaded _separately_ from the OSU [webpage](https://www.tpxo.net/global) (free registration for academic usage)!**


The code is based on the [original OTPS software package](https://www.tpxo.net/otps) written in Fortran-90 and developed at the [OSU](https://www.tpxo.net/home) with further inspirations from:

- [pyTMD](https://pytmd.readthedocs.io/en/latest/) (Python-based tidal prediction software) by T. C. Sutterley
- [OtisTidalPrediction](https://github.com/jklymak/OtisTidalPrediction) by jklymak
- [tidal-prediction-python](https://gitlab.com/jblarsen/tidal-prediction-python) by Jesper Baasch-Larsen

The original OTPS license can be found in the file [LICENSE_OSU.pdf](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/LICENSE_OSU.pdf).

## Current capabilities

**NOTE**: This package only supports the prediction of tidal elevations (**h**) and **<u>not</u>** tidal transport (**u**, **v**).

At the moment, both the `TPXO9-atlas-v4` and `-v5` tidal models are supported. The following tidal constituents are available:

- M2, S2, K1, O1, N2, P1, K2, Q1, 2N2, K2, M4, MF, MM, MN4, MS4, 
- S1 (**only v5**) 

Behind the scenes, `tpxo-tide-prediction` is based on the great `xarray` package that enables easy I/O of netCDF files. To minimize computational time and resources, the global constituent files are clipped to the requested coordinate region plus a buffer of **3** nodes (_default_) in order to avoid edge effects during the following *linear* interpolation. The selected constituents are sampled at the coordinates location of the interpolated grid. 

**NOTE**: The chosen *linear* interpolation and used indexing method result in predicted tides that are slightly different compared to the output of the OSU Tidal Prediction Software ([OTPS](https://www.tpxo.net/otps)) using a *bi-linear* interpolation method!

## Installation

Generally, I recommend using a `conda` environment with `pip` installed. To install this package directly from GitHub, run the following line:

```python
pip install git+https://github.com/fwrnke/tpxo-tide-prediction.git
```

You can also download the repository as a ZIP archive. To install `tpxo-tide-prediction` locally using `pip` after unzipping the downloaded archive:

```bash
>>> cd ./tpxo-tide-prediction  # root directory of unzipped package
>>> pip install [-e] .         # -e: setuptools "develop mode"
```

## Usage

### Command line interface (CLI)

`tpxo-tide-prediction` was designed to be mainly used from the command line and provides an interface that can be used by running:

```python
predict_tide {model_dir} {model_dir} [optional parameters]
```

This script requires two positional arguments:

- `model_dir`: Input directory of the model files (e.g. `./data`).
- `params_file`: Input parameter file with columns [LAT, LON, UTC TIME (YYYY-MM-DDThh:mm:ss)]. Could be provided as [LAT, LON] if `--time` is set or as [UTC TIME] if `--lat` and `--lon` are set.

Optionally, the following parameter can be used:

- `--help`, `-h`: show help
- `--constituents`: Available tidal constituents supported by TPXO9 atlas model.
   - **default**: m2, s2, n2, k2, k1, o1, p1, q1
- `--correct_minor`: Correct for minor tidal constituents.
- `--lat`: Constant latitude for each timestep in parameter file. Expecting only TIME column.
- `--lon`: Constant longitude for each timestep in parameter file. Expecting only TIME column.
- `--time`, `-t`: One of the following options (expecting only LAT and LON columns):
   - constant time for every coordinate position (YYYY-MM-DDThh:mm:ss)
   - start and end time [START END] 
   - start and end time with stepsize [START END STEP] (as _X_(D)ay, _X_(M)onth, _X_(Y)ear, _X_(h)ours, _X_(m)inutes, or _X_(s)econds, with _X_: integer number of units, e.g. "5s")
- `--output_file`: Output file path for predicted tides.
- `--mode`, `-m`: Output either _time x position matrix_ (**full**) or only at matching _time, lat, lon_ locations (**track**).

------

### Python script

Additionally, it is possible to just import the essential functions `read_parameter_file`, `tide_predict`, and `write_tides` when using these methods in a custom script:

#### Import package

```python
from tpxo_tide_prediction import (
    read_parameter_file,
    tide_predict,
    write_tides
    )
```

#### Read input parameter from file

```python
# (A) read inputs from parameter file
lat, lon, times = read_parameter_file('path/to/parameter/file.txt')

# (B) read only `time` from parameter file and provide fixed location
lat = -36.446349
lon = 175.166068
lat, lon, times = read_parameter_file('path/to/parameter/file', lat, lon)
```

#### Compute tidal elevation

```python
# lat, lon:  int, float, np.ndarray
# times:     str, np.ndarray(dtype='datetime64[s]')

# (A) default configuration
tide = tide_predict('path/to/TPXO9-atlas-model', lat, lon, times)

# (B) custom configuration (less constituents, no correction of minor constituents)
tide = tide_predict('path/to/TPXO9-atlas-model', lat, lon, times,
                    constituents=['m2','s2','n2','k2'], correct_minor=False)
```

#### Write computed tides to formatted output file

```python
# (A) custom output file
# full output (tide at every time at every location)
write_tides('path/to/output/file.tide', tide, mode='full')
# or `track` output (lat-lon-time triplets)
write_tides('path/to/output/file.tide', tide, mode='track')

# (B) create output file from parameter file
basepath, filename = os.path.split('path/to/parameter/file.txt')
basename, suffix = os.path.splitext(filename)
write_tides(basepath, basename, tide, mode='full')
```

## Example parameter files

The folder `./examples` contains different input parameter files for the use of `tpxo-tide-prediction` via the CLI (`predict_tide`). 

- [params_almost-constant-lat_almost-constant-lon_time.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_almost-constant-lat_almost-constant-lon_time.txt): slightly varying geographic location and timesteps
- [params_constant-lat_constant-lon_time.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_constant-lat_constant-lon_time.txt): single location and timesteps
- [params_lat_lon.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_lat_lon.txt): different locations **<u>but</u>** no timesteps (provide times using `--time` parameter) 
- [params_time.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_time.txt): only timesteps (provide coordinates using `--lat` and `--lon` parameters) 
- [params_lat_lon_constant-time.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_lat_lon_constant-time.txt): different locations with constant time
- [params_lat_lon_time_track.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_lat_lon_time_track.txt): different locations with different timesteps (**track**, negative latitudes, positive longitudes)
- [params_lat_lon_time_track-2.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_lat_lon_time_track-2.txt): different locations with different timesteps (**track**, negative latitudes and longitudes)
- [params_lat_lon_time_OTPSnc.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_lat_lon_time_OTPSnc.txt): various lat, lon, time combinations adapted from input file shipped with [OTPSnc](https://www.tpxo.net/otps)
- [params_tracks_dateline.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_tracks_dateline.txt): track (lat, lon, time) crossing the dateline (180°E/180°W)
- [params_tracks_equator.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_tracks_equator.txt): track (lat, lon, time) crossing the equator
- [params_tracks_north_pole.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_tracks_north_pole.txt): track (lat, lon, time) running close to the north pole
- [params_tracks_zero_meridian.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_tracks_zero_meridian.txt):  track (lat, lon, time) crossing the prime meridian (0° longitude)
- [params_tracks_between_grid_nodes.txt](https://github.com/fwrnke/tpxo-tide-prediction/blob/main/examples/params_tracks_between_grid_nodes.txt): track (lat, lon, time) with all longitude values located in between two grid nodes of TPXO netCDF files
