"""
Compute tidal corrections for user-specified timestamps and geographical positions.
Using TPXO9-atlas-v4/v5 as tidal model with 1/30 degree resolution fully global solution,
obtained by combining 1/6 degree base global solution TPXO9.v1 and
thirty 1/30 degree resolution local solutions for all coastal areas,
including Arctic and Antarctic.
The model is available from: https://www.tpxo.net/global/tpxo9-atlas

This code is based on the original OTPS software package written in Fortran-90
and developed by Gary Egbert & Svetlana Erofeeva from the Oregon State University.

Further inspirations are taken from:
    pyTMD  - Python-based tidal prediction software
        by T. C. Sutterley, https://pytmd.readthedocs.io/en/latest/
    OtisTidalPrediction
        by jklymak, https://github.com/jklymak/OtisTidalPrediction
    tidal-prediction-python
        by Jesper Baasch-Larsen, https://gitlab.com/jblarsen/tidal-prediction-python

@author: fwrnke
@mail:   fwrnke@mailbox.org
@date:   2021-10-26

"""
import os
import glob
import argparse
import numpy as np

from .tidal_constituents import (
    CONST_ID,
    OMEGA_d,
    PHASE_mkB
    )
from .utils import (
    read_h_netCDFs,
    nodal,
    infer_minor
    )

__all__ = [
    'tide_predict',
    'read_parameter_file',
    'write_tides',
    ]

#%%
def define_input_args():
    parser = argparse.ArgumentParser(
        description='Compute tidal corrections for user-specified timestamps and geographical positions')
    parser.add_argument('model_dir', type=str, help='Input directory of tidal model files.')
    parser.add_argument('params_file', type=str,
                        help='Input file with columns [LAT,LON,UTC TIME (YYYY-MM-DDThh:mm:ss)].\
                              Could be provided as [LAT,LON] if "--time" is set or \
                              as [UTC TIME] if "--lat" and "--lon" are set.')
    parser.add_argument('--constituents', '-c', nargs='+',
                        choices=['m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1',
                                 'm4', 'mf', '2n2', 'mm', 'mn4', 'ms4'],
                        default=['m2','s2','n2','k2','k1','o1','p1','q1'],
                        help='Available tidal constituents supported by TPXO9 atlas model.')
    parser.add_argument('--correct_minor', action='store_true',
                        help='Correct for minor tidal constituents.')
    parser.add_argument('--lat', type=float,
                        help='Constant latitude for each timestep in parameter file. Expecting only TIME column.')
    parser.add_argument('--lon', type=float,
                        help='Constant longitude for each timestep in parameter file. Expecting only TIME column.')
    parser.add_argument('--time', '-t', nargs='*', #type=str,
                        help='Constant time for every coordinate position (YYYY-MM-DDThh:mm:ss)\
                              OR start and end time [START END] \
                              OR start and end time with stepsize [START END STEP] \
                              (as X(D)ay, X(M)onth, X(Y)ear, X(h)ours, X(m)inutes, or X(s)econds\
                               with X: integer number of units, e.g. "5s").\
                              Expecting only LAT and LON columns.')
    parser.add_argument('--mode', '-m', type=str, choices=['full', 'track'], default='full',
                        help='Output either time x position matrix (full) or only at matching "time, lat, lon" locations (track).')
    parser.add_argument('--export_params', action='store_true',
                        help='Export input parameter (lat, lon, time) together with tides.')
    parser.add_argument('--output_file', '-o', type=str,
                        help='Output file path for predicted tides.')
    return parser


def tide_predict(model_dir, lat, lon, times,
                 constituents:list=['m2','s2','n2','k2','k1','o1','p1','q1'],
                 correct_minor=True,
                 method='linear',
                 mode='full'):
    """
    Predict tides for given timestamps and geographical positions.

    Parameters
    ----------
    model_dir : str
        Path to model directory.
    lat : np.ndarray
        Array of latitudes.
    lon : np.ndarray
        Array of longitudes.
    times : np.ndarray
        Array of timestamps.
    constituents : list, optional
        Tidal constituents to use for prediction.
    correct_minor : bool, optional
        Correct tide for minor constituents (default: True).
    method : str
        Interpolation method. Choose from 'nearest' or 'linear'.
    mode : str
        Either `full` tide matrix or only tides along `track`.

    Returns
    -------
    h_pred : np.ndarray
        Predicted tides with shape (timestamp x positions).

    """
    # sanity checks
    if not os.path.isdir(model_dir):
        raise IOError(f'>{model_dir}< is not a valid file directory!')
    
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    
    # get number of positions
    npts = lon.size

    #--- prepare times
    # calc modified julian day (MJD) since 1992-01-01
    if isinstance(times, np.ndarray):
        time_init = times[0]
    elif np.issubdtype(times.dtype, np.datetime64):
        time_init = times
    else:
        raise ValueError('Type of input times not known. Please use np.datetime64 format!')
    mjd = (time_init.astype('datetime64[D]') - np.datetime64('1992-01-01', 'D')).astype('float') + 48622.0
    # calc timesteps in seconds since 1992-01-01
    timesteps = (times - np.datetime64('1992-01-01','s')).astype('float')
    # get number of timesteps
    ntimes = times.size

    #--- prepare input tidal constituents
    constituents = [c.lower() for c in constituents if c.lower() in CONST_ID]
    ncon = len(constituents)

    # get index of each constituent for selecting constants
    const_idx = []
    for constid in constituents:
        try:
            const_idx.append(CONST_ID.index(constid.lower()))
        except ValueError as err:
            print('[WARNING]    ', err)

    #--- read input tidal constituents
    # read OTIS tidal elevation files
    paths_h = glob.glob(model_dir + '/h_*.nc')
    h = read_h_netCDFs(paths_h, lat, lon, constituents, method=method)

    # build (complex) harmonic constant array from netCDFs
    hc = np.empty((ncon, npts), dtype='complex')
    for k, c in enumerate(h.attrs['constituents']):
        hc[k,:] = (h[f'{c}_hRe'].data + h[f'{c}_hIm'].data * 1.0j) / 1000   # convert mm to m

    # compute nodal corrections
    pf, pu = nodal(mjd, constituents)

    #--- predict tides
    try:
        # init output array
        h_pred = np.zeros((1, npts)) if mode == 'track' else np.zeros((ntimes, npts))
    except MemoryError as error:
        raise MemoryError(str(error) + '. Please use mode="track" or split input dataset.')
    
    # loop over each constituent for every (coordinate) postition
    for nn in range(ncon):
        # get index of constituent to get corresponding `OMEGA_d` and `PHASE_mkB` value
        j = const_idx[nn]
        # coordinate position loop for single constituent
        for ii in range(npts):
            mask = ii if mode == 'track' else None
            # sum over all tidal constituents
            h_pred[:,ii] = h_pred[:,ii] + pf[nn] * hc.real[nn,ii] * \
                np.cos(OMEGA_d[j] * timesteps[mask].squeeze() + PHASE_mkB[j] + pu[nn]) - \
 	            pf[nn] * hc.imag[nn,ii] * np.sin(OMEGA_d[j] * timesteps[mask].squeeze() + PHASE_mkB[j] + pu[nn])
                 
    # infer corrections for minor constituents
    if correct_minor:
        dh = infer_minor(hc, constituents, mjd, timesteps)
        h_pred += dh

    return h_pred.squeeze()


def read_parameter_file(params_file, lat=None, lon=None, times=None):
    """
    Read custom parameter file comprising of the following columns:
        (1)  UTC timestamps, LAT, LON
        (2)  LAT, LON           (constant TIME set from command line)
        (3)  UTC timestamps     (constant LAT, LON set from command line)

    Parameters
    ----------
    params_file : str
        Path to parameter file.
    lat : float, None
        Constant latitude for tide predictions or None.
    lon : float, None
        Constant longitude for tide predictions or None.
    times : str, None
        Constant time for tide predictions or None.

    Returns
    -------
    lat : np.ndarray
        Array of latitude(s) for tidal predictions.
    lon : np.ndarray
        Array of longitudes(s) for tidal predictions.
    times : np.ndarray
        Array of time(s) for tidal predictions.

    """
    with open(params_file, 'r') as f:
        lat_list = []
        lon_list = []
        times_list = []
        for line in f:
            line = line.rstrip().strip().split(',')
            ncols = len(line)

            if ncols == 1:
                times_list.append(np.datetime64(line[0]))
            elif ncols == 2:
                lat_list.append(float(line[0]))
                lon_list.append(float(line[1]))
            elif ncols == 3:
                lat_list.append(float(line[0]))
                lon_list.append(float(line[1]))
                times_list.append(np.datetime64(line[2]))

    # command line parameter supersedes parameter file input
    if ((lat is not None) and (lon is not None)) or \
        ((len(lat_list) == 0) and (len(lat_list) == 0) and \
         (lat is not None) and (lon is not None)):
        lat = np.array([lat], dtype='float')
        lon = np.array([lon], dtype='float')
    elif len(lat_list) > 1 and len(lon_list) > 1:
        lat = np.asarray(lat_list)
        lon = np.asarray(lon_list)
    else:
        raise ValueError('No LAT and LON found in parameter files and no coordinates provided. Check your inputs!')

    # check if single postion was provided
    lat_uniq = np.unique(lat)
    lon_uniq = np.unique(lon)
    if lat_uniq.size == 1 and lon_uniq.size == 1:
        lat = lat_uniq
        lon = lon_uniq
    assert lat.size == lon.size, 'Input `lat` and `lon` arrays must have identical sizes and shapes!'

    # command line parameter supersedes parameter file input
    if (times is not None) or (len(times_list) == 0 and times is not None):
        # check if times is iterable and has more than 1 element
        if isinstance(times, (list, tuple, np.ndarray)) and len(times) > 1:
            if len(times) == 2:
                step = np.timedelta64(1, 'D') # set default stepsize
            else:
                number, unit = times[-1][:-1], times[-1][-1]
                step = np.timedelta64(number, unit) # set custom stepsize
            # set time range
            times = np.arange(np.datetime64(times[0]), np.datetime64(times[1]) + step , step)
        else:
            # single time --> constant time for every data point
            times = np.array([np.datetime64(t) for t in times])
    elif len(times_list) > 1:
        times = np.asarray(times_list)
    else:
        raise ValueError('No times found in parameter files and no constant time provided. Check your inputs!')

    # filter for unique time values (e.g. when provided constant time in parameter file)
    times_uniq = np.unique(times)
    if times_uniq.size == 1:
        times = times_uniq

    return lat, lon, times


def write_tides(tide,
                mode,
                output_file=None,
                basepath=None,
                basename=None,
                lat=None,
                lon=None,
                times=None,
                export_params=False):
    """
    Write predicted tides to text file.

    Parameters
    ----------
    tide : np.ndarray
        Predicted tides to output.
    mode : str
        Output mode, eithter 'full' (complete array) or 'track' (only along line).
    output_file : str, None
        Path of output file if provided.
        Either `output_file` or `basepath` and `basename` must be provided.
    basepath : str, None
        Directory of parameter file.
        Either `output_file` or `basepath` and `basename` must be provided.
    basename : str, None
        Parameter filename without suffix.
        Either `output_file` or `basepath` and `basename` must be provided.
    lat : np.ndarray, None
        Latitude used for tide predictions.
    lon : np.ndarray, None
        Longitude used for tide predictions.
    times : np.ndarray, None
        Timestamps used for tide predictions.
    export_params : bool, False
        If True, export LAT, LON, TIMESTAMP ("track") or TIMESTAMP ("full") with tide.

    """
    prefix = np.datetime_as_string(np.datetime64('today')).replace(':','')
    
    if tide.ndim == 1:
        nrows, ncols = tide.size, 1
    elif tide.ndim == 2:
        nrows, ncols = tide.shape
    
    if mode == 'full':
        output_file = output_file if output_file is not None else os.path.join(basepath, f'{prefix}_{basename}_tides.txt')
        header = ','.join([f'pos{i+1}_m' for i in range(ncols)])
        
        if export_params and times is not None and lat.size > 1:
            header = 'timestamp,' + header
            fmt = '%s' + ',%.6f' * ncols
            data = np.empty((nrows, ncols + 1), dtype='object')
            data[:,0]  = np.datetime_as_string(times)
            data[:,1:] = np.atleast_2d(tide).T if tide.ndim == 1 else tide
        elif export_params and times is not None and lat.size == 1:
            header = 'lat,lon,timestamp,tide_m'
            fmt = '%.6f,%.6f,%s,%.6f'
            data = np.empty((nrows, 4), dtype='object')
            data[:,0] = lat
            data[:,1] = lon
            data[:,2] = np.datetime_as_string(times)
            data[:,3] = tide.squeeze()
        elif not export_params:
            # header = 'Predicted tides for input parameter file provided as timestamp x position (rows x cols)'
            fmt = '%.6f'
            data = tide
        np.savetxt(output_file, data, fmt=fmt, delimiter=',', newline='\n', header=header, comments='')
    elif mode == 'track':
        output_file = output_file if output_file is not None else os.path.join(basepath, f'{prefix}_{basename}_tides-along-track.txt')

        if export_params and all([var is not None for var in [lat, lon, times]]):
            header = 'lat,lon,timestamp,tide_m'
            fmt = '%.6f,%.6f,%s,%.6f'
            data = np.empty((nrows, 4), dtype='object')
            data[:,0] = lat
            data[:,1] = lon
            data[:,2] = np.datetime_as_string(times)
            data[:,3] = tide#[idx_track]
        elif not export_params:
            header = 'tide_m'
            fmt = '%.6f'
            data = tide#[idx_track]
        np.savetxt(output_file, data, fmt=fmt, delimiter=',', newline='\n', header=header, comments='')


def main(input_args=None):
    parser = define_input_args()
    args = parser.parse_args(input_args)

    basepath, filename = os.path.split(args.params_file)
    basename, suffix = os.path.splitext(filename)

    # read inputs from parameter file
    lat, lon, times = read_parameter_file(args.params_file, args.lat, args.lon, args.time)

    if args.mode == 'track' and not (lat.size == lon.size == times.size):
        raise ValueError(
            'For mode `track` the input parameter `lat`, `lon` and `times` must have the same size!'
            )

    #  predict tides at given times and positions
    tide = tide_predict(args.model_dir, lat, lon, times, args.constituents,
                        correct_minor=args.correct_minor, method='linear', mode=args.mode)

    # save predicted tides to output file
    write_tides(tide, args.mode, args.output_file, basepath, basename, lat, lon, times, args.export_params)

#%%
if __name__ == '__main__':
    main()
