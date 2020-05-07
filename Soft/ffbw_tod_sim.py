import os

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from ctypes import *

from util.calibration_model import Mission_Model_Parameters, FocalPlane
from util.h5_util import H5IO
from util.bs2px_func import bs2px, bs2px_bm, bs2px4dadb, bs2px_bm_v2, src2bmcoord
from util.bmsys import beam_map_func
from util.beam_characterization import beam_mapmaking

class ScanData(H5IO):
    def __init__(self):
        self.t = None
        self.az = None
        self.el = None

class TOD(H5IO):
    def __init__(self):
        pass
        # self.bs_altaz_az = None
        # self.bs_altaz_el = None
        # self.bs_altaz_dk = 0.

class BM(object):
    def __init__(self):
        pass

def ground_pickup(tod, chidx, pltopt=False):
    """
    Add the ground emission, it depend on the az and el.
    first we can use some polynomail function on az direction;
    in the elevation level, we assume the 1/sin(el) function.
    :param tod: tod object
    :param chidx: the index of the channel
    :param pltopt: plot option
    :return: modify the tod object
    """
    p = np.poly1d([1, 1, 1])
    da = 1./np.sin(tod.pxa_altaz_az[chidx])*p(tod.pxa_altaz_az[chidx])
    db = 1. / np.sin(tod.pxb_altaz_az[chidx]) * p(tod.pxb_altaz_az[chidx])
    tod.da[chidx] += da/np.max(da)*0.5
    tod.db[chidx] += db/np.max(db)*0.5
    if pltopt:
        plt.plot(tod.da[0])
        plt.show()

def add_noise(tod):
    """
    add white noise to the tod, the noise level is about 300uk*sqrt(s)
    as we know that, the source emission is about 260K, which is very larger than noise.
    :param tod:
    :return:
    """
    pass

# @time_consuming
def scanfile_gen(mmp, checkopt=False):
    """
    generate the scan file which contains the t, az, el offset for each scanset.
    it includes left elnod, constant elevation scan, right elnod.
    we get the az, el using the integration from the acc value; or we can just load the existing scan file.
    :param mmp: mission model parameters.
    :param checkopt: whether we want to show the az and el data to check the result is right.
    :return: scandata object including t, az, el attributes.
    """
    scanfile_name = '../Data/scanfile.hdf5'
    if os.path.exists(scanfile_name):
        scandata = ScanData()
        scandata.loadh5(filename=scanfile_name, opt=['t','az','el'])
    else:

        # construct the ces
        # ===========================================================================================
        # construct the ces start
        tr_ces_start = mmp.tacc_ces+mmp.tc1_ces+mmp.tacc_ces
        t_ces_start = np.arange(0, tr_ces_start+1./mmp.f_sample, 1./mmp.f_sample)
        ind1 = np.logical_and(t_ces_start >= 0, t_ces_start < mmp.tacc_ces)
        ind2 = np.logical_and(t_ces_start >= mmp.tacc_ces, t_ces_start < mmp.tacc_ces + mmp.tc1_ces)
        ind3 = np.logical_and(t_ces_start >= mmp.tacc_ces + mmp.tc1_ces, t_ces_start < mmp.tacc_ces + mmp.tc1_ces + mmp.tacc_ces)
        acc_ces_start = np.zeros_like(t_ces_start)
        acc_ces_start[ind1] = mmp.acc_ces_v
        acc_ces_start[ind2] = 0.
        acc_ces_start[ind3] = -mmp.acc_ces_v
        v_ces_start = [np.sum(acc_ces_start[0:i])/mmp.f_sample for i in range(len(acc_ces_start))]
        az_ces_start = np.array([np.sum(v_ces_start[0:i]) / mmp.f_sample for i in range(len(v_ces_start))])
        # acc_ces_start = np.insert(acc_ces_start, 0, 0)
        # v_ces_start = np.cumsum(acc_ces_start[0:-1]/mmp.f_sample)
        # az_ces_start = np.cumsum(v_ces_start/mmp.f_sample)
        # plt.plot(t_ces_start, az_ces_start)
        # plt.show()
        # construct the ces end
        az_ces_end = az_ces_start[-1::-1]

        # construct one full scan
        tr_ces_fs = mmp.tacc_ces+2*mmp.tc2_ces+2*mmp.tacc_ces+2*mmp.tc2_ces+mmp.tacc_ces
        t_ces_fs = np.arange(0, tr_ces_fs, 1./mmp.f_sample)
        ind1 = np.logical_and(t_ces_fs>=0, t_ces_fs<mmp.tacc_ces)
        ind2 = np.logical_and(t_ces_fs >= mmp.tacc_ces, t_ces_fs < mmp.tacc_ces+2*mmp.tc2_ces)
        ind3 = np.logical_and(t_ces_fs >= mmp.tacc_ces+2*mmp.tc2_ces, t_ces_fs < mmp.tacc_ces + 2 * mmp.tc2_ces+2*mmp.tacc_ces)
        ind4 = np.logical_and(t_ces_fs >= mmp.tacc_ces + 2 * mmp.tc2_ces+2*mmp.tacc_ces, t_ces_fs < mmp.tacc_ces + 2 * mmp.tc2_ces + 2 * mmp.tacc_ces+2*mmp.tc2_ces)
        ind5 = np.logical_and(t_ces_fs >= mmp.tacc_ces + 2 * mmp.tc2_ces + 2 * mmp.tacc_ces+2*mmp.tc2_ces,
                              t_ces_fs <= mmp.tacc_ces + 2 * mmp.tc2_ces + 2 * mmp.tacc_ces + 2 * mmp.tc2_ces+mmp.tacc_ces)
        acc_ces_fs = np.zeros_like(t_ces_fs)
        acc_ces_fs[ind1] = -mmp.acc_ces_v
        acc_ces_fs[ind2] = 0.
        acc_ces_fs[ind3] = mmp.acc_ces_v
        acc_ces_fs[ind4] = 0.
        acc_ces_fs[ind5] = -mmp.acc_ces_v
        # acc_ces_fs = np.insert(acc_ces_fs, 0, 0)
        # print(len(acc_ces_fs))
        # plt.plot(acc_ces_fs)
        # plt.show()
        v_ces_fs = [np.sum(acc_ces_fs[0:i])/mmp.f_sample for i in range(len(acc_ces_fs))]
        az_ces_fs = np.array([az_ces_start[-1]+np.sum(v_ces_fs[0:i])/mmp.f_sample for i in range(len(v_ces_fs))])
        # acc_ces_fs = np.insert(acc_ces_fs, 0, 0)
        # v_ces_fs = np.cumsum(acc_ces_fs[:-1]/mmp.f_sample)
        # az_ces_fs = az_ces_start[-1] + np.cumsum(v_ces_fs / mmp.f_sample)
        # print(len(v_ces_fs))
        # plt.plot(v_ces_fs1)
        # plt.plot(v_ces_fs)
        # plt.show()
        # plt.plot(az_ces_fs)
        # plt.plot(az_ces_fs1)
        # plt.show()
        az_ces_fss = np.tile(az_ces_fs[:], mmp.n_fullscans)
        # plt.plot(t_ces_fs, az_ces_fs)
        # az_ces = np.concatenate((az_ces_start, az_ces_fss[1:], az_ces_end[1:]))
        az_ces = az_ces_fss
        el_ces = np.zeros_like(az_ces)
        t_ces = np.arange(len(az_ces))/mmp.f_sample
        # plt.plot(t_ces, az_ces)
        # ===================================================================================================

        # combine all the data
        # ===================================================================================================
        zeros = np.zeros(int(3*mmp.f_sample))
        # convert it to radians
        az = az_ces*np.pi/180
        t = np.arange(len(az))/mmp.f_sample
        # ===================================================================================================
        scandata = ScanData()
        scandata.t = t
        scandata.az = az
        scandata.az = az+np.random.normal(0, mmp.p_rms, size=len(az))
        scandata.el = 0+np.random.normal(0, mmp.p_rms, size=len(az))
        # save the t, az, el into hdf5 file
        # scandata.saveh5(filename=scanfile_name, opt=['t', 'az', 'el'])
    if checkopt:
        # plt.plot(t, az)
        # print(az[0], az[-1])
        # plt.show()
        plt.plot(scandata.t, scandata.el)
        plt.xlim([0, 100])
        plt.show()
    return scandata

def top_simulation(tod, mmp, dk):
    """
    time ordered pointinng simulation
    :param tod: tod object
    :param mmp: mission model parameters object
    :param dk: deck angle
    :return: modified tod object
    """
    scandata = scanfile_gen(mmp)
    az = mmp.az_ref + scandata.az
    tod.bs_altaz_az = np.tile(az, mmp.n_steps)
    el = mmp.el_ref+np.arange(mmp.n_steps)*mmp.el_step - 5.*np.pi/180
    tod.bs_altaz_el = np.repeat(el, len(az))
    tod.bs_altaz_dk = np.ones_like(tod.bs_altaz_az)*dk*np.pi/180
    # tod.bs_altaz_dk = np.random.normal(0, 1, len(tod.bs_altaz_el))
    # tod.bs_altaz_dk[0] = 0.
    if 0:
        # this is used for plot the TOP
        plt.plot(tod.bs_altaz_az*180/np.pi)
        plt.xlim([0,10000])
        plt.xlabel('Sample')
        plt.ylabel('Azimuth')
        plt.show()
        plt.plot(tod.bs_altaz_el * 180 / np.pi)
        plt.xlabel('Sample')
        plt.ylabel('Elevation')
        plt.show()
    return


def tod_sim_old(tod, fp, mmp, dk):
    """
    tod simulation, this version will not be used anymore.
    :param tod: tod object
    :param fp: focal plane structure object
    :param mmp: mission model parameters object
    :param dk: deck angle value
    :return: modified tod object.
    """
    # initialize the detector readout with 0
    tod.da = np.zeros_like(tod.px_altaz_az, dtype=np.float64)
    tod.db = np.zeros_like(tod.px_altaz_az, dtype=np.float64)
    # source modulation scale factor, 20Hz from 1 to 0.1
    scale_fac0 = np.ones_like(tod.da[0])
    f = 5
    idxtmp = np.arange(f, len(scale_fac0), f)
    # idxtmp = (np.array(len(scale_fac)))
    # scale_fac0[idxtmp] = 0.1
    mmp.nside = 2048

    # find the source pixel idx
    src_vec = hp.ang2vec(np.pi/2-mmp.el_ref, mmp.az_ref)
    src_pix_idx = hp.query_disc(mmp.nside, src_vec, 1. / 60 * np.pi / 180)
    # print(len(src_pix_idx))
    # mt = np.zeros(hp.nside2npix(mmp.nside))
    # mt[src_pix_idx] = 1.
    # hp.mollview(mt)
    # plt.show()
    # map_src = np.ones_like(src_pix_idx)
    map_src = np.ones(len(src_pix_idx), dtype=np.float64)
    for ch_idx in range(fp.lenfp):
    # for ch_idx in range(10):
        print('detector index', ch_idx)
        # calculate the beam map
        bm = beam_map_func(None, fp, mmp, ch_idx)
        # find the index which are close to the point source
        # idx = np.logical_and(np.abs(tod.px_altaz_az[ch_idx]-mmp.az_ref)<8*mmp.fwhm, np.abs(tod.px_altaz_el[ch_idx]-mmp.el_ref)<8*mmp.el_ref)
        idx = np.where(np.logical_and(np.abs(tod.px_altaz_az[ch_idx] - mmp.az_ref) < 9 * mmp.fwhm,
                             np.abs(tod.px_altaz_el[ch_idx] - mmp.el_ref) < 9 * mmp.fwhm))[0]
        print(len(idx))
        tod_tmp = TOD()
        tod_tmp.bs_altaz_az = tod.bs_altaz_az[idx]
        tod_tmp.bs_altaz_el = tod.bs_altaz_el[idx]
        tod_tmp.bs_altaz_dk = tod.bs_altaz_dk[idx]
        scale_fac = np.ascontiguousarray(scale_fac0[idx], dtype=np.float64)
        # tod_tmp.bs_altaz_dk = np.zeros_like(tod_tmp.bs_altaz_az)
        # bs2px_bm(tod_tmp, bm, pltopt=False)
        from time import time
        print(time())
        # bs2px_bm_v2(tod_tmp, bm)
        # the first method project the beam to sky
        bs2px_bm(tod_tmp, bm)
        # or we can project the sky to beam coordinate, which will be much faster
        print(time())
        # convert the beam tod (bm_az, bm_el) into pixel tod array
        pix_idx_array = hp.ang2pix(mmp.nside, np.pi/2-tod_tmp.bm_el, tod_tmp.bm_az, lonlat=False)
        bm.length = len(tod_tmp.bm_el)
        da = np.zeros(len(idx), dtype=np.float64)
        lib = cdll.LoadLibrary('./clib/bm_convolve.so')
        func = lib.beam_convolve
        pix_idx_array = np.ascontiguousarray(pix_idx_array, dtype=np.int64)
        src_pix_idx = np.ascontiguousarray(src_pix_idx, dtype=np.int64)
        map_src = np.ascontiguousarray(map_src, dtype=np.float64)
        # src modulation
        bm.Ba = np.ascontiguousarray(bm.Ba, dtype=np.float64)
        func(c_void_p(pix_idx_array.ctypes.data), c_void_p(src_pix_idx.ctypes.data), c_void_p(bm.Ba.ctypes.data), \
             c_void_p(map_src.ctypes.data), c_void_p(scale_fac.ctypes.data), c_void_p(da.ctypes.data), \
             c_int(len(idx)), c_int(bm.length), c_int(len(map_src)))

        tod.da[ch_idx, idx] = da
        tod.db[ch_idx, idx] = da
        if 0:
            # this can be use to check wheter our simulation is right or not
            plt.plot(tod.da[ch_idx])
            # plt.xlim([0, 300000])
            # plt.plot(tod.px_altaz_az[ch_idx], tod.da[ch_idx])
            # plt.plot(tod.px_altaz_el[ch_idx]*180/np.pi, tod.da[ch_idx])
            plt.show()
    ground_pickup(tod, pltopt=True)
    filename = mmp.toddir+'FF_Calibration_dk{}_{}.h5'.format(dk, mmp.srcopt)
    tod.saveh5(filename, opt=['da', 'db', 'bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk'])

def tod_sim_faster(tod, fp, mmp, dk4filename):
    """
    tod simulation
    :param tod: tod object
    :param fp: focal plane structure object
    :param mmp: mission model parameters object
    :param dk4filename: deck angle value which will be used in the filename
    :return: modified tod object.
    """
    # Initialize the detector readout as 0
    tod.da = np.zeros_like(tod.px_altaz_az, dtype=np.float64)
    tod.db = np.zeros_like(tod.px_altaz_az, dtype=np.float64)
    # set the source emission property
    # source modulation scale factor,
    scale_fac0 = np.ones_like(tod.da[0])*0.1
    if mmp.srcopt!='const':
        # here, the srcopt is the chopping frequency
        f = mmp.srcopt
        step = int(mmp.f_sample/f)
        idx0 = np.arange(len(scale_fac0))
        idx1 = idx0%step
        idxtmp = idx1<(int(step/2))
        scale_fac0[idxtmp] = 1
        # plt.plot(scale_fac0[0:200])
        # plt.show()

    # find the source pixel idx, if the source only 1 pixel, here we will not use it
    # here we use the source position as az_ref and el_ref
    mmp.nside = 2048
    src_vec = hp.ang2vec(np.pi / 2 - mmp.el_ref, mmp.az_ref)
    src_pix_idx = hp.query_disc(mmp.nside, src_vec, 1. / 60 * np.pi / 180)
    map_src = np.ones(len(src_pix_idx), dtype=np.float64)

    for ch_idx in range(fp.lenfp):
        print('detector index', ch_idx)
        # this is temperary, in the real case, for each detector we need a beam file
        bmfilename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_0.0000.txt'
        bmfilename = None
        # calculate the beam map
        bm = beam_map_func(bmfilename, fp, mmp, ch_idx, sopt=True)
        # find the index which are close to the point source
        az = tod.px_altaz_az[ch_idx]
        el = tod.px_altaz_el[ch_idx]
        dk = tod.px_dk[ch_idx]
        # this idx is for choose the az and el which is close to the beam center, here we choose a rectangular,
        # or we can use circular.
        r_critical = mmp.n_sigma/2*mmp.sigma
        idx = (((az-mmp.az_ref)**2+(el-mmp.el_ref)**2)<r_critical**2)
        az1 = az[idx]
        el1 = el[idx]
        dk1 = dk[idx]
        # we can project the sky to beam coordinate, which will be much faster
        # convert the src al el to beam coordinate
        theta = (fp.pxa_theta[ch_idx] + fp.pxb_theta[ch_idx]) / 2
        az_bmcoord, el_bmcoord = src2bmcoord(mmp, az1, el1, dk1, theta)
        # find the idx for the source
        az_idx = ((az_bmcoord-mmp.xmin)/mmp.xstep).astype(np.int)
        el_idx = ((el_bmcoord - mmp.ymin) / mmp.ystep).astype(np.int)
        print(np.max(el_idx), np.max(az_idx))
        tod.da[ch_idx,idx] = bm.Ba[el_idx, az_idx]*scale_fac0[idx]
        tod.db[ch_idx, idx] = bm.Bb[el_idx, az_idx] * scale_fac0[idx]
        ground_pickup(tod, ch_idx)
        if 0:
            # this can be use to check wheter our simulation is right or not
            plt.plot(tod.da[ch_idx])
            # plt.xlim([0, 300000])
            # plt.plot(tod.px_altaz_az[ch_idx], tod.da[ch_idx])
            # plt.plot(tod.px_altaz_el[ch_idx]*180/np.pi, tod.da[ch_idx])
            plt.show()

    filename = mmp.toddir + 'FF_Calibration_dk{}_{}_faster.h5'.format(dk4filename, mmp.srcopt)
    tod.saveh5(filename, opt=['da', 'db', 'bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk'])


if __name__=='__main__':
    # load the mission model parameters object
    mmp = Mission_Model_Parameters(src_center_opt=False)
    # load the focal plane object
    fp = FocalPlane()
    # set the dk angles
    dks = np.array([0, 45, 90, 135, 180])
    for dk in dks:
        # for each deck angle, we simulate the tod, then save to a file
        print('dk={}'.format(dk))
        # initialize the tod object
        tod = TOD()
        # simulate the TOP
        top_simulation(tod, mmp, dk)
        # calculate the trajectory of each focal plane pixel, which is very time consuming
        # bs2px(tod, fp, pltopt=False)
        bs2px4dadb(tod, fp)
        # simulate the TOD for each detector
        # tod_sim(tod, fp, mmp, dk)
        tod_sim_faster(tod, fp, mmp, dk)
