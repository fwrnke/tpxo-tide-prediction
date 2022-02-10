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
                                   [0.79374876]]),
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
                                   [ 0.12353511]]))
                           ),
                          ]
                         )
def test_nodal(mjd, constituents, expected_nodal_corrections, allclose):
    # python -c "from pyTMD.load_nodal_corrections import load_nodal_corrections; print(load_nodal_corrections(59824, ['m2', 's2', 'k1', 'o1']))"
    assert allclose(nodal(mjd, constituents), expected_nodal_corrections)

# def test_infer_minor():
#     times = np.arange(np.datetime64('1989-11-09','D'), np.datetime64('1989-11-21','D'))
#     timesteps = (times - np.datetime64('1992-01-01','s')).astype('float')
#     ncon = 5
#     npts = 20
#     hc = np.empty((ncon, npts), dtype='complex')
    
    