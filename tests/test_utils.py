import os
import glob

import numpy as np
import pytest

from tpxo_tide_prediction.utils import (
    astrol,
    nodal,
    infer_minor,
    longitude_to_360,
    longitude_to_180,
    subset_region,
    # read_grd_netCDF,
    read_h_netCDFs,
    # open_u_nc
    )

@pytest.mark.parametrize('mjd, expected_astrol', 
                         [(59824, (232.30027963747852, 161.1341070731487, 285.71910461747075, 46.61372421736121)),
                          (47780, (295.7810745175011, 169.99730323314907, 23.974989297470643, 324.38893009736114)),
                          (47960, (147.53244091750094, 347.413828033149, 44.02762469747063, 314.8572514973612))
                          ])
def test_astrol(mjd, expected_astrol):
    # python -c "from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes; print(calc_astrol_longitudes(47780))"
    assert astrol(mjd) == expected_astrol


@pytest.mark.parametrize('mjd, constituents, expected_nodal_corrections',
                         [(59824, ['m2', 's2', 'k1', 'o1'],
                           (np.array([[0.97470489],
                                      [1.        ],
                                      [1.08530079],
                                      [1.13778514]]),
                            np.array([[-0.02728905],
                                      [ 0.        ],
                                      [-0.10156566],
                                      [ 0.11658498]]))
                           ),
                          (47780, ['m2', 's2', 'k1', 'o1', 'n2', 'p1', 'k2', 'q1', '2n2', 'mu2', 'nu2', 'l2'],
                           (np.array([[0.97006719],
                                   [1.        ],
                                   [1.0967259 ],
                                   [1.15652401],
                                   [0.97006719],
                                   [1.        ],
                                   [1.26019397],
                                   [1.15802746],
                                   [0.97006719],
                                   [0.97006719],
                                   [0.97006719],
                                   [0.79374876]], dtype='float64'),
                            np.array([[ 0.02188946],
                                   [ 0.        ],
                                   [ 0.08008796],
                                   [-0.09161668],
                                   [ 0.02188946],
                                   [ 0.        ],
                                   [ 0.16874858],
                                   [-0.09510551],
                                   [ 0.02188946],
                                   [ 0.02188946],
                                   [ 0.02188946],
                                   [ 0.12353511]], dtype='float64'))
                           ),
                          ]
                         )
def test_nodal(mjd, constituents, expected_nodal_corrections, allclose):
    # python -c "from pyTMD.load_nodal_corrections import load_nodal_corrections; print(load_nodal_corrections(59824, ['m2', 's2', 'k1', 'o1']))"
    assert allclose(nodal(mjd, constituents), expected_nodal_corrections)


def test_infer_minor(allclose):
    paths_h = glob.glob(os.path.join(os.path.dirname(__file__), '../data/v4/h_*.nc'))
    if paths_h == []:
        p = os.path.join(os.path.dirname(__file__), '../data/v4/')
        raise FileNotFoundError(f'Could not find model netCDF files in path: {p}')
    lat = -36.5
    lon = -175
    constituents = ['m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1']
    h = read_h_netCDFs(paths_h, lat, lon, constituents)
    
    
    time_start = np.datetime64('2019-12-08','D')
    time_end   = time_start + np.timedelta64(10, 'D') #np.datetime64('1989-11-21','D')
    dtime = np.timedelta64(1, 'h')
    times = np.arange(time_start, time_end, dtime)
    mjd = (time_start.astype('datetime64[D]') - np.datetime64('1992-01-01', 'D')).astype('float') + 48622.0
    timesteps = (times - np.datetime64('1992-01-01','s')).astype('float')
    
    # build (complex) harmonic constant array from netCDFs
    ncon = len(constituents)
    npts = 1 if isinstance(lat, (int, float)) else len(lat)
    hc = np.empty((ncon, npts), dtype='complex')
    for k, c in enumerate(h.attrs['constituents']):
        hc[k,:] = (h[f'{c}_hRe'].data + h[f'{c}_hIm'].data * 1.0j) / 1000   # convert mm to m
    
    # infer corrections for minor constituents
    dh = infer_minor(hc, constituents, mjd, timesteps)
    
    assert allclose(dh, np.atleast_2d(np.repeat([-0.00280321], len(times))).T)


@pytest.mark.parametrize('lon, expected_lon', 
                         [(    0,      0),
                          (-180 ,   -180),
                          (180.5, -179.5),
                          (276.2214, -83.7786)
                          ])
def test_longitude_to_180(lon, expected_lon):
    lon_180 = longitude_to_180(lon)
    assert np.allclose(lon_180, expected_lon)
    
    
@pytest.mark.parametrize('lon, expected_lon', 
                         [(    0,      0),
                          (-180 ,    180),
                          (-179.5, 180.5),
                          (276.2214, 276.2214)
                          ])
def test_longitude_to_360(lon, expected_lon):
    lon_360 = longitude_to_360(lon)
    assert np.allclose(lon_360, expected_lon)