"""
Utility functions used for tidal prediction based on TPXO9-atlas models.

"""

import math
import bisect

import numpy as np
import xarray as xr

from .tidal_constituents import CONST_ID


def astrol(mjd):
    """
    Computes the basic astronomical mean longitudes  s, h, p, N.
    Note N is not N', i.e. N is decreasing with time.
    These formulae are for the period 1990 - 2010, and were derived
    by David Cartwright (personal comm., Nov. 1990).
    mjd is UTC in decimal MJD.
    
    All longitudes returned in degrees.
    R. D. Ray    Dec. 1990
    
    Non-vectorized version.
    
    """
    circle = 360.0
    
    T = mjd - 51544.4993
    
    # mean longitude of moon
    s = 218.3164 + 13.17639648 * T
    
    # mean longitude of sun
    h = 280.4661 +  0.98564736 * T
    
    # mean longitude of lunar perigee
    p =  83.3535 +  0.11140353 * T
    
    # mean longitude of ascending lunar node
    N = 125.0445 -  0.05295377 * T
    
    s = np.mod(s, circle)
    h = np.mod(h, circle)
    p = np.mod(p, circle)
    N = np.mod(N, circle)

    return s, h, p, N


def nodal(mjd:float, constit:list):
    """
    Calculates the nodal corrections for tidal constituents.

    Parameters
    ----------
    mjd : float
        Modified Julian Day (MJD) since 1992-01-01.
    constit : list
        List of constituents to use.

    Returns
    -------
    pf, pu : np.ndarray
        Nodal corrections for the constituents.

    """
    # index_labels = ['m2', 's2', 'k1', 'o1', 'n2', 'p1', 'k2', 'q1', '2n2', 'mu2', 
    #                 'nu2', 'l2', 't2', 'j1', 'm1', 'oo1', 'rho1', 'mf', 'mm', 'ssa', 
    #                 'm4', 'ms4', 'mn4', 'mk3', 's6', '2sm2']
    # index = [29, 34, 18, 11, 26, 16, 36,  9, 24, 25, 
    #          27, 32, 33, 22, 13, 23, 10,  4,  2,  1,
    #          44, 45, 43, 49, 41, 50, 39]
    
    pp = 282.94     # solar perigee at epoch 2000
    rad = math.pi/180
    
    hour = (mjd - int(mjd)) * 24.0
    t1 = 15.0 * hour
    t2 = 30.0 * hour
    
    # get the basic astronomical mean longitudes
    s, h, p, omega = astrol(mjd)
    
    # list of all constituents available for this function
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1', 
              'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1', 
              'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 
              'm2a', 'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 
              'mns2', '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 
              's5', 'm6', 's6', 's7', 's8']
    
    sinn = np.sin(omega*rad)
    cosn = np.cos(omega*rad)
    sin2n = np.sin(2*omega*rad)
    cos2n = np.cos(2*omega*rad)
    sin3n = np.sin(3*omega*rad)
    
    # # arg not needed!
    # arg = np.empty((53), dtype=np.float64)
    # arg[0]  = h - pp                    # Sa
    # arg[1]  = 2*h                       # Ssa
    # arg[2]  = s - p                     # Mm
    # arg[3]  = 2*s - 2*h                 # MSf
    # arg[4]  = 2*s                       # Mf
    # arg[5]  = 3*s - p                   # Mt
    # arg[6]  = t1 - 5*s + 3*h + p - 90   # alpha1
    # arg[7]  = t1 - 4*s + h + 2*p - 90   # 2Q1
    # arg[8]  = t1 - 4*s + 3*h - 90       # sigma1
    # arg[9]  = t1 - 3*s + h + p - 90     # q1
    # arg[10] = t1 - 3*s + 3*h - p - 90   # rho1
    # arg[11] = t1 - 2*s + h - 90         # o1
    # arg[12] = t1 - 2*s + 3*h + 90       # tau1
    # arg[13] = t1 - s + h + 90           # M1
    # arg[14] = t1 - s + 3*h - p + 90     # chi1
    # arg[15] = t1 - 2*h + pp - 90        # pi1
    # arg[16] = t1 - h - 90               # p1
    # arg[17] = t1 + 90                   # s1
    # arg[18] = t1 + h + 90               # k1
    # arg[19] = t1 + 2*h - pp + 90        # psi1
    # arg[20] = t1 + 3*h + 90             # phi1
    # arg[21] = t1 + s - h + p + 90       # theta1
    # arg[22] = t1 + s + h - p + 90       # J1
    # arg[23] = t1 + 2*s + h + 90         # OO1
    # arg[24] = t2 - 4*s + 2*h + 2*p      # 2N2
    # arg[25] = t2 - 4*s + 4*h            # mu2
    # arg[26] = t2 - 3*s + 2*h + p        # n2
    # arg[27] = t2 - 3*s + 4*h - p        # nu2
    # arg[28] = t2 - 2*s + h + pp         # M2a
    # arg[29] = t2 - 2*s + 2*h            # M2
    # arg[30] = t2 - 2*s + 3*h - pp       # M2b
    # arg[31] = t2 - s + p + 180          # lambda2
    # arg[32] = t2 - s + 2*h - p + 180    # L2
    # arg[33] = t2 - h + pp               # t2
    # arg[34] = t2                        # S2
    # arg[35] = t2 + h - pp + 180         # R2
    # arg[36] = t2 + 2*h                  # K2
    # arg[37] = t2 + s + 2*h - pp         # eta2
    # arg[38] = t2 - 5*s + 4.0*h + p      # MNS2
    # arg[39] = t2 + 2*s - 2*h            # 2SM2
    # arg[40] = 1.5*arg[29]               # M3
    # arg[41] = arg[18] + arg[29]         # MK3
    # arg[42] = 3*t1                      # S3
    # arg[43] = arg[26] + arg[29]         # MN4
    # arg[44] = 2*arg[29]                 # M4
    # arg[45] = arg[29] + arg[34]         # MS4
    # arg[46] = arg[29] + arg[36]         # MK4
    # arg[47] = 4*t1                      # S4
    # arg[48] = 5*t1                      # S5
    # arg[49] = 3*arg[29]                 # M6
    # arg[50] = 3*t2                      # S6
    # arg[51] = 7.0*t1                    # S7
    # arg[52] = 4*t2                      # S8    
    
    f = np.empty((53), dtype=np.float64)
    f[0]  = 1                                     # Sa
    f[1]  = 1                                     # Ssa
    f[2]  = 1 - 0.130*cosn                        # Mm
    f[3]  = 1                                     # MSf
    f[4]  = 1.043 + 0.414*cosn                    # Mf
    f[5]  = np.sqrt((1+.203*cosn+.040*cos2n)**2 + (.203*sinn+.040*sin2n)**2)  # Mt
    f[6]  = 1                                     # alpha1
    f[7]  = np.sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  # 2Q1
    f[8]  = f[7]                                  # sigma1
    f[9]  = f[7]                                  # q1
    f[10] = f[7]                                  # rho1
    f[11] = np.sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 + (0.189*sinn-0.0058*sin2n)**2)    # O1
    f[12] = 1                                     # tau1
    tmp1  = 1.36*np.cos(p*rad)+.267*np.cos((p-omega)*rad)  # Ray's
    tmp2  = 0.64*np.sin(p*rad)+.135*np.sin((p-omega)*rad)
    f[13] = np.sqrt(tmp1**2 + tmp2**2)                 # M1
    f[14] = np.sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  # chi1
    f[15] = 1                                     # pi1
    f[16] = 1                                     # P1
    f[17] = 1                                     # S1
    f[18] = np.sqrt((1.+.1158*cosn-.0029*cos2n)**2 + (.1554*sinn-.0029*sin2n)**2)       # K1
    f[19] = 1                                     # psi1
    f[20] = 1                                     # phi1
    f[21] = 1                                     # theta1
    f[22] = np.sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  # J1
    f[23] = np.sqrt((1.0+0.640*cosn+0.134*cos2n)**2 + (0.640*sinn+0.134*sin2n)**2 )      # OO1
    f[24] = np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2)    # 2N2
    f[25] = f[24]                                 # mu2
    f[26] = f[24]                                 # N2
    f[27] = f[24]                                 # nu2
    f[28] = 1                                     # M2a
    f[29] = f[24]                                 # M2
    f[30] = 1                                     # M2b
    f[31] = 1                                     # lambda2
    temp1 = 1.-0.25*np.cos(2*p*rad) - 0.11*np.cos((2*p-omega)*rad) - 0.04*cosn
    temp2 = 0.25*np.sin(2*p) + 0.11*np.sin((2*p-omega)*rad) + 0.04*sinn
    f[32] = np.sqrt(temp1**2 + temp2**2)          # L2
    f[33] = 1                                     # t2
    f[34] = 1                                     # S2
    f[35] = 1                                     # R2
    f[36] = np.sqrt((1.+.2852*cosn+.0324*cos2n)**2 + (.3108*sinn+.0324*sin2n)**2)  # K2
    f[37] = np.sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  # eta2
    f[38] = f[29]**2                              # MNS2
    f[39] = f[29]                                 # 2SM2
    f[40] = 1   # wrong                           # M3
    f[41] = f[18]*f[29]                           # MK3
    f[42] = 1                                     # S3
    f[43] = f[29]**2                              # MN4
    f[44] = f[43]                                 # M4
    f[45] = f[43]                                 # MS4
    f[46] = f[29]*f[36]                           # MK4
    f[47] = 1                                     # S4
    f[48] = 1                                     # S5
    f[49] = f[29]**3                              # M6
    f[50] = 1                                     # S6
    f[51] = 1                                     # S7
    f[52] = 1                                     # S8
    
    u = np.empty((53), dtype=np.float64)
    u[ 0] = 0                                    # Sa
    u[ 1] = 0                                    # Ssa
    u[ 2] = 0                                    # Mm
    u[ 3] = 0                                    # MSf
    u[ 4] = -23.7*sinn + 2.7*sin2n - 0.4*sin3n   # Mf
    u[ 5] = np.arctan(-(.203*sinn+.040*sin2n)/ (1+.203*cosn+.040*cos2n))/rad   # Mt
    u[ 6] = 0                                    # alpha1
    u[ 7] = np.arctan(.189*sinn/(1.+.189*cosn))/rad      # 2Q1
    u[ 8] = u[7]                                 # sigma1
    u[ 9] = u[7]                                 # q1
    u[10] = u[7]                                 # rho1
    u[11] = 10.8*sinn - 1.3*sin2n + 0.2*sin3n    # O1
    u[12] = 0                                    # tau1
    u[13] = np.arctan2(tmp2,tmp1)/rad            # M1
    u[14] = np.arctan(-.221*sinn/(1.+.221*cosn))/rad     # chi1
    u[15] = 0                                    # pi1
    u[16] = 0                                    # P1
    u[17] = 0                                    # S1
    u[18] = np.arctan((-.1554*sinn+.0029*sin2n)/ (1.+.1158*cosn-.0029*cos2n))/rad   # K1
    u[19] = 0                                    # psi1
    u[20] = 0                                    # phi1
    u[21] = 0                                    # theta1
    u[22] = np.arctan(-.227*sinn/(1.+.169*cosn))/rad     # J1
    u[23] = np.arctan(-(.640*sinn+.134*sin2n)/ (1.+.640*cosn+.134*cos2n))/rad  # OO1
    u[24] = np.arctan((-.03731*sinn+.00052*sin2n)/ (1.-.03731*cosn+.00052*cos2n))/rad  # 2N2
    u[25] = u[24]                                # mu2
    u[26] = u[24]                                # N2
    u[27] = u[24]                                # nu2
    u[28] = 0                                    # M2a
    u[29] = u[24]                                # M2
    u[30] = 0                                    # M2b
    u[31] = 0                                    # lambda2
    u[32] = np.arctan(-temp2/temp1)/rad          # L2
    u[33] = 0                                    # t2
    u[34] = 0                                    # S2
    u[35] = 0                                    # R2
    u[36] = np.arctan(-(.3108*sinn+.0324*sin2n)/ (1.+.2852*cosn+.0324*cos2n))/rad  # K2
    u[37] = np.arctan(-.436*sinn/(1.+.436*cosn))/rad     # eta2
    u[38] = u[29]*2                              # MNS2
    u[39] = u[29]                                # 2SM2
    u[40] = 1.5*u[29]                            # M3
    u[41] = u[29] + u[18]                        # MK3
    u[42] = 0                                    # S3
    u[43] = u[29]*2                              # MN4
    u[44] = u[43]                                # M4
    u[45] = u[29]                                # MS4
    u[46] = u[29]+u[36]                          # MK4
    u[47] = 0                                    # S4
    u[48] = 0                                    # S5
    u[49] = u[29]*3                              # M6
    u[50] = 0                                    # S6
    u[51] = 0                                    # S7
    u[52] = 0                                    # S8
    
    # filter input constituents based on list of available ones
    constit = [c for c in constit if c.lower() in cindex]
    # get number of constituents to include
    nconstit = len(constit)
    # init output arrays
    pu = np.zeros((nconstit,1))
    pf = np.ones((nconstit,1))
    
    # add nodal corrections for tidal constituents from input list
    for i, cons in enumerate(constit):
        if cons.lower() in CONST_ID:
            ii = cindex.index(cons)
            #print(f'[INFO]   nodal():  add < {cons.lower()} >')
            pf[i,:] = f[ii]
            pu[i,:] = u[ii] * rad
        else:
            print(f'[WARNING]   nodal(): < {cons.lower()} > not part of primary tidal constituents!')
    
    return pf, pu
    

def infer_minor(z, constituents:list, mjd:float, timesteps):
    """
    Calculate the tidal corrections for minor constituents inferred using
    major constituents.

    Parameters
    ----------
    z : np.ndarray
        Complex harmonic constituents (constituents x points).
    constituents : list
        List of tidal constituents.
    mjd : float
        Modified Julian Day (MJD) since 1992-01-01.
    timesteps : np.ndarray
        Array of timesteps in seconds since 1992-01-01.

    Returns
    -------
    dh : np.ndarray
        Tidal height from minor constituents.

    """
    rad = math.pi/180
    PP = 282.8
    cid8 = ['q1','o1','p1','k1','n2','m2','s2','k2']
    ncid8 = len(cid8)
    
    # number of constituents & number of points (locations)
    nc, npts = z.shape
    # number of timesteps
    nt = len(np.atleast_1d(timesteps))
    # number of data points to calculate
    n = nt if ((npts == 1) & (nt > 1)) else npts
    # init output array
    dh = np.zeros((n))
    
    # re-order constituents to correspond to cid8
    ncon = len(constituents)
    z8 = np.zeros((ncid8, npts), dtype='complex')
    ni = 0
    for i in range(ncid8):
        for j in range(ncon):
            if constituents[j].lower() == cid8[i]:
                z8[i] = z[j]
                if i not in [2, 7]:
                    ni += 1

    if ni < 6:
        raise ValueError('Not enough constituents for inference!')
    
    # list of minor constituents
    minor = ['2q1','sigma1','rho1','m12','m11','chi1','pi1','phi1','theta1',
             'j1','oo1','2n2','mu2','nu2','lambda2','l2','l2','t2']
    # only add minor constituents that are not on the list of major values
    minor_indices = [i for i,m in enumerate(minor) if m not in constituents]
    
    # relationship between major and minor constituent amplitude and phase
    zmin = np.empty((18, n), dtype='complex')
    zmin[0]  = 0.263 * z8[0] - 0.0252 * z8[1]    # 2Q1
    zmin[1]  = 0.297 * z8[0] - 0.0264 * z8[1]    # sigma1
    zmin[2]  = 0.164 * z8[0] + 0.0048 * z8[1]    # rho1 +
    zmin[3]  = 0.0140 * z8[1] + 0.0101 * z8[3]   # M1
    zmin[4]  = 0.0389 * z8[1] + 0.0282 * z8[3]   # M1
    zmin[5]  = 0.0064 * z8[1] + 0.0060 * z8[3]   # chi1
    zmin[6]  = 0.0030 * z8[1] + 0.0171 * z8[3]   # pi1
    zmin[7]  = -0.0015 * z8[1] + 0.0152 * z8[3]  # phi1
    zmin[8]  = -0.0065 * z8[1] + 0.0155 * z8[3]  # theta1
    zmin[9]  = -0.0389 * z8[1] + 0.0836 * z8[3]  # J1 +
    zmin[10] = -0.0431 * z8[1] + 0.0613 * z8[3]  # OO1 +
    zmin[11] = 0.264 * z8[4] - 0.0253 * z8[5]    # 2N2 +
    zmin[12] = 0.298 * z8[4] - 0.0264 * z8[5]    # mu2 +
    zmin[13] = 0.165 * z8[4] + 0.00487 * z8[5]   # nu2 +
    zmin[14] = 0.0040 * z8[5] + 0.0074 * z8[6]   # lambda2
    zmin[15] = 0.0131 * z8[5] + 0.0326 * z8[6]   # L2 +
    zmin[16] = 0.0033 * z8[5] + 0.0082 * z8[6]   # L2 +
    zmin[17] = 0.0585 * z8[6]                    # t2 +
    
    hour = (mjd - int(mjd)) * 24.0
    t1 = 15.0 * hour
    t2 = 30.0 * hour
    
    # get the basic astronomical mean longitudes
    S, H, P, omega = astrol(mjd)
    
    # determine equilibrium tidal arguments
    arg = np.empty((18, n), dtype=np.float64)
    arg[0,:] = t1 - 4. * S + H + 2. * P - 90.  # 2Q1
    arg[1,:] = t1 - 4. * S + 3. * H - 90.  # sigma1
    arg[2,:] = t1 - 3. * S + 3. * H - P - 90.  # rho1
    arg[3,:] = t1 - S + H - P + 90.  # M1
    arg[4,:] = t1 - S + H + P + 90.  # M1
    arg[5,:] = t1 - S + 3. * H - P + 90.  # chi1
    arg[6,:] = t1 - 2. * H + PP - 90.  # pi1
    arg[7,:] = t1 + 3. * H + 90.  # phi1
    arg[8,:] = t1 + S - H + P + 90.  # theta1
    arg[9,:] = t1 + S + H - P + 90.  # J1
    arg[10,:] = t1 + 2. * S + H + 90.  # OO1
    arg[11,:] = t2 - 4. * S + 2. * H + 2. * P  # 2N2
    arg[12,:] = t2 - 4. * S + 4. * H  # mu2
    arg[13,:] = t2 - 3. * S + 4. * H - P  # nu2
    arg[14,:] = t2 - S + P + 180.  # lambda2
    arg[15,:] = t2 - S + 2. * H - P + 180.  # L2
    arg[16,:] = t2 - S + 2. * H + P  # L2
    arg[17,:] = t2 - H + PP  # t2

    # determine nodal corrections f and u
    sinn = np.sin(omega * rad)
    cosn = np.cos(omega * rad)
    sin2n = np.sin(2. * omega * rad)
    cos2n = np.cos(2. * omega * rad)
    
    # init f
    f = np.ones((18, n), dtype=np.float64)
    f[0,:] = np.sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 + (0.189*sinn - 0.0058*sin2n)**2)
    f[1,:] = f[0]
    f[2,:] = f[0]
    f[3,:] = np.sqrt((1.0 + 0.185 * cosn)**2 + (0.185 * sinn)**2)
    f[4,:] = np.sqrt((1.0 + 0.201 * cosn)**2 + (0.201 * sinn)**2)
    f[5,:] = np.sqrt((1.0 + 0.221 * cosn)**2 + (0.221 * sinn)**2)
    f[9,:] = np.sqrt((1.0 + 0.198 * cosn)**2 + (0.198 * sinn)**2)
    f[10,:] = np.sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 + (0.640*sinn + 0.134*sin2n)**2 )
    f[11,:] = np.sqrt((1.0 - 0.0373 * cosn)**2 + (0.0373 * sinn)**2)
    f[12,:] = f[11]
    f[13,:] = f[11]
    f[15,:] = f[11]
    f[16,:] = np.sqrt((1.0 + 0.441 * cosn)**2 + (0.441 * sinn)**2)
    
    # init u
    u = np.zeros((18, n), dtype=np.float64)
    u[0,:] = np.arctan2(0.189*sinn - 0.0058*sin2n, 1.0 + 0.189*cosn - 0.0058*sin2n)/rad
    u[1,:] = u[0]
    u[2,:] = u[0]
    u[3,:] = np.arctan2(0.185 * sinn, 1.0 + 0.185 * cosn) / rad
    u[4,:] = np.arctan2(-0.201 * sinn, 1.0 + 0.201 * cosn) / rad
    u[5,:] = np.arctan2(-0.221 * sinn, 1.0 + 0.221 * cosn) / rad
    u[9,:] = np.arctan2(-0.198 * sinn, 1.0 + 0.198 * cosn) / rad
    u[10,:] = np.arctan2(-0.640*sinn - 0.134*sin2n, 1.0 + 0.640*cosn + 0.134*cos2n)/rad
    u[11,:] = np.arctan2(-0.0373 * sinn, 1.0 - 0.0373 * cosn) / rad
    u[12,:] = u[11]
    u[13,:] = u[11]
    u[15,:] = u[11]
    u[16,:] = np.arctan2(-0.441 * sinn, 1.0 + 0.441 * cosn) / rad
    
    # compute sum of minor tidal constituents
    for k in minor_indices:
        th = (arg[k,:] + u[k,:]) * rad 
        dh += zmin.real[k,:] * f[k,:] * np.cos(th) - zmin.imag[k,:] * f[k,:] * np.sin(th)
    
    # reshape array to fit input "z" and allow for numpy's broadcasting when adding
    # case: timeseries but single position --> shape: (n, 1)
    # case: timeseries AND multiple positions --> shape: (1, 20)
    dh = dh.reshape((-1,z.shape[-1]))
    
    return dh

def longitude_to_360(lon):
    """ Convert longitude (DD) to 0 to 360 degree range. """
    return lon % 360

def longitude_to_180(lon):
    """ Convert longitude (DD) to -180 to 180 degree range. """
    return ((lon + 180) % 360) - 180


def subset_region(ds, lat, lon, coords_lat='lat_z', coords_lon='lon_z', offset=10):
    """
    Return coordinate slices that match existing coordinate locations in xr.Dataset.
    Slices are determined using input coordinate arrays and offset (# of nodes).

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset (netCDF).
    lat : np.ndarray
        Array of latitudes where to sample netCDF (1d).
    lon : np.ndarray
        Array of longitudes where to sample netCDF (1d.
    coords_lat : str, optional
        Name of netCDF coordinate corresponding to latidude. The default is 'lat_z'.
    coords_lon : str, optional
        Name of netCDF coordinate corresponding to longitude. The default is 'lon_z'.
    offset : int, optional
        Offset value (# of nodes in input dataset) used to extend coordinate boundaries. 
        The default is 10.

    Returns
    -------
    slices
        Latitude and longitude coordinate slices.

    """
    # get min/max from coords
    lat_min, lat_max = lat.min(), lat.max()
    lon_min, lon_max = lon.min(), lon.max()
    
    lat_min_pad = float(ds[coords_lat][bisect.bisect(ds[coords_lat], lat_min) - offset].values)
    lon_min_pad = float(ds[coords_lon][bisect.bisect(ds[coords_lon], lon_min) - offset].values)
    
    lat_max_pad = float(ds[coords_lat][bisect.bisect(ds[coords_lat], lat_max) + offset].values)
    lon_max_pad = float(ds[coords_lon][bisect.bisect(ds[coords_lon], lon_max) + offset].values)

    return slice(lat_min_pad, lat_max_pad), slice(lon_min_pad, lon_max_pad)


def read_grd_netCDF(path, lat, lon):  #TODO: update!
    """
    Read bathymetry at provided coordinates (lat, lon) from grid netCDF.

    Parameters
    ----------
    path : str
        File path of grid netCDF.
    lat : np.ndarray
        Array of latitudes where to sample netCDF (1d).
    lon : np.ndarray
        Array of longitudes where to sample netCDF (1d).

    Returns
    -------
    None.

    """
    # convert longitude from -180/180 to 0/360 range
    lon = longitude_to_360(lon)
    
    # remove duplicate coordinates (for interpolation)
    lat_uniq, lat_idx, lat_inv = np.unique(
        lat, return_index=True, return_inverse=True)
    lon_uniq, lon_idx, lon_inv = np.unique(
        lon, return_index=True, return_inverse=True)
    
    # create DataArrays for indexing
    lat_da = xr.DataArray(lat, dims='pos')
    lon_da = xr.DataArray(lon, dims='pos')
    
    # open netCDF
    grd = xr.open_dataset(path)
    
    # split into different datasets
    grd_z = grd[[v for v in grd.keys() if 'z' in v]]
    grd_u = grd[[v for v in grd.keys() if 'u' in v]]
    grd_v = grd[[v for v in grd.keys() if 'v' in v]]
    
    
    # get list of variables containing coordinates
    grd_z = grd_z.swap_dims({'nx':'lon_z', 'ny':'lat_z'})
    grd_u = grd_u.swap_dims({'nx':'lon_u', 'ny':'lat_u'})
    grd_v = grd_v.swap_dims({'nx':'lon_v', 'ny':'lat_v'})
    
    lat_slice, lon_slice = subset_region(grd_z, lat, lon, offset=10)
    
    # sample grid at coordinate locations using cubic interpolation
    # using individual slices for h, u, v because their coordinate grids differ!
    lat_slice, lon_slice = subset_region(grd_z, lat, lon, offset=10)
    # (1) select subset of full grid based on slices
    # (2) interpolate values within subset
    # (3) sample interpolated grid at coordinate positions
    hz_interp = grd_z.sel(lon_z=lon_slice, lat_z=lat_slice).interp(lon_z=lon, lat_z=lat, method='cubic').sel(lon_z=lon_da, lat_z=lat_da)
    
    lat_slice, lon_slice = subset_region(grd_u, lat, lon, offset=10, coords_lat='lat_u', coords_lon='lon_u')
    hu_interp = grd_u.sel(lon_u=lon_slice, lat_u=lat_slice).interp(lon_u=lon, lat_u=lat, method='cubic').sel(lon_u=lon_da, lat_u=lat_da)
    
    lat_slice, lon_slice = subset_region(grd_v, lat, lon, offset=10, coords_lat='lat_v', coords_lon='lon_v')
    hv_interp = grd_v.sel(lon_v=lon_slice, lat_v=lat_slice).interp(lon_v=lon, lat_v=lat, method='cubic').sel(lon_v=lon_da, lat_v=lat_da)
    
    return hz_interp, hu_interp, hv_interp
    

def read_h_netCDFs(paths, lat, lon, constituents:list, offset:int=10, method='cubic'):
    """
    Read and sample tidal elevation constituent files at given coordinates.
    A subset of the input data is selected using coordinate extent increased by
    provided offset (as nodes in input netCDF files).
    This subset is subsequently interpolated (cubic) 
    and sampled at the coordinate locations.

    Parameters
    ----------
    paths : list
        List of input file paths.
    lat : np.ndarray
        Array of latitudes where to sample netCDF (1d).
    lon : np.ndarray
        Array of longitudes where to sample netCDF (1d).
    constituents : list
        List of tidal constituents.
    offset : int
        Offset value (nodes in input netCDFs) used to extend coordinate boundaries
        before extracting subset and interpolating.
    method : str
        Interpolation method. Choose from 'nearest', 'linear', or 'cubic'.

    Returns
    -------
    h_interp : xr.Dataset
        Dataset of tidal constituents as DataArrays sampled from input netCDFs
        at each lat/lon position.

    """
    # convert longitude from -180/180 to 0/360 range
    lon = longitude_to_360(lon)
    
    # remove duplicate coordinates (for interpolation)
    lat_uniq, lat_idx, lat_inv = np.unique(
        lat, return_index=True, return_inverse=True)
    lon_uniq, lon_idx, lon_inv = np.unique(
        lon, return_index=True, return_inverse=True)

    # create DataArrays for indexing
    lat_da = xr.DataArray(lat, dims='pos')
    lon_da = xr.DataArray(lon, dims='pos')
    
    # downcase  constituent names
    constituents = [c.lower() for c in constituents]
    
    # loop over file paths and load only specified ones
    ds_cons = []
    cons_loaded = []
    for con in constituents:
        check = [con in p for p in paths]
        if any(check):
            path = paths[check.index(True)]
            print(f'[INFO]    Loading elevation data for < {con} >')
            #print(path)
            ds = xr.open_dataset(path)
            ds = ds.drop_vars('con').swap_dims({'nx':'lon_z', 'ny':'lat_z'})
            var_names = list(ds.data_vars) #list(ds.var.keys())
            ds = ds.rename_vars(dict(zip(var_names, [f'{con}_{v}' for v in var_names])))
            
            ds_cons.append(ds)
            cons_loaded.append(con)
    
    # merge individual constituents into single dataset
    h = xr.merge(ds_cons)
    h.attrs['constituents'] = cons_loaded
      
    # sample constituents at coordinate locations using cubic interpolation 
    lat_slice, lon_slice = subset_region(h, lat, lon, offset=offset)
    # h_not_interp = h.sel(lon_z=lon_slice, lat_z=lat_slice)
    h_interp = h.sel(
        lon_z=lon_slice, lat_z=lat_slice
        ).interp(
            lon_z=lon_uniq, lat_z=lat_uniq, method=method
            ).sel(lon_z=lon_da, lat_z=lat_da)
    
    if len(constituents) != len(cons_loaded):
        print('[WARNING]    Could not find constituents: ', 
              list(set(constituents) - set(cons_loaded)))
    
    return h_interp


def open_u_nc(paths, lat, lon, constituents:list):
    #TODO: smart way to separate u and v from each netCDF...
    pass
