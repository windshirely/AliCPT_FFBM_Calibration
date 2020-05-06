import os

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from ctypes import *
from scipy.optimize import curve_fit

from util.calibration_model import Mission_Model_Parameters, FocalPlane
from util.h5_util import H5IO
from util.bs2px_func import bs2px, bs2px_with_dk, bs2px4dadb
from util.bmsys import bm_func, bmplot, beam_map_func
from util.beam_characterization import beam_mapmaking, BM, beam_mapmaking_on_detector


class TOD(H5IO):
    def __init__(self):
        pass
        # self.bs_altaz_az = None
        # self.bs_altaz_el = None
        # self.bs_altaz_dk = 0.

def beam_charc(fp, mmp, dks):
    """
    do the beam characterization process.
    For each dk angle and each detector, we will plot the beam pattern.
    :param fp: focal plane object
    :param mmp: mission model parameters object
    :param dks: deck angle values
    :return: save the beam pattern into different file
    """
    for dk in dks:
        print('dk={}'.format(dk))
        # load the tod data
        tod = TOD()
        filename = mmp.toddir+'FF_Calibration_dk{}_{}_faster.h5'.format(dk, mmp.srcopt)
        tod.loadh5(filename, opt=['da', 'db', 'bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk'])
        # bs2px_with_dk(tod, fp, pltopt=False)
        # calculate the trajectory of each detector pixel
        bs2px4dadb(tod, fp)
        for i in range(len(tod.da)):
            # for each detector, we will make a map
            # beam_mapmaking(tod, mmp, i, dk)
            beam_mapmaking_on_detector(tod, mmp, fp, i, dk)

def beam_parameter_fit(fp, mmp):
    for i in range(fp.lenfp):
        # load the beam map
        bm_data = BM()
        filename = mmp.bmdir+'Output/FinalBM_{}_{}Hz.h5'.format(i, mmp.srcopt)
        bm_data.loadh5(filename=filename, opt=['map', 'bm_x', 'bm_y'])
        # fit the parameters
        pv, pcov = bmparams_fit(bm_data)
        A = pv[0]
        sigma = pv[1]
        # print(pcov)
        print(sigma*180/np.pi*60)
        print(pv)
        bm_in_obj = beam_map_func(None, fp, mmp, i, sopt=True)
        bm_in = np.reshape(bm_in_obj.Ba, (mmp.ny, mmp.nx))
        pltopt = 1
        if pltopt:
            obsmap = np.reshape(bm_data.map, (mmp.ny, mmp.nx))
            fitmap = np.reshape(bm_func(bm_data, *pv), (mmp.ny, mmp.nx))
            # bmplot(mmp, obsmap, fitmap)
            bmplot(mmp, bm_in, fitmap)

def bmparams_fit(bm):
    # we have [A, sigma, xc, yc, p, c] parameters to fit
    A0 = 1
    sigma0 = 5. / 60 * np.pi / 180
    xc0 = 0
    yc0 = 0
    p0 = 0
    c0 = 0
    factor_tmp = np.pi / 180
    p0 = [A0, sigma0, xc0, yc0, p0, c0]
    lower = np.array([1e-5, 1. / 60 * factor_tmp, -1 * factor_tmp, -1 * factor_tmp, -1, -1])
    upper = np.array([1, 6. / 60 * factor_tmp, 1 * factor_tmp, 1 * factor_tmp, 3, 1])
    # p0 = [A0, sigma0]
    # lower = np.array([1e-2, 1./60*factor_tmp])
    # upper = np.array([10, 6./60*factor_tmp])
    pv, pcov = curve_fit(bm_func, bm, bm.map.flatten(), p0=p0, bounds=(lower, upper), check_finite=True)
    return pv, pcov

def beam_combination(fp, mmp, dks):
    for i in range(fp.lenfp):
        j = 0
        for dk in dks:
            bm_data = BM()
            filename = mmp.bmdir + 'Output/bm{}_dk{}_{}.h5'.format(i, dk, mmp.srcopt)
            bm_data.loadh5(filename=filename, opt=['map', 'bm_x', 'bm_y'])
            if j==0:
                bm = bm_data.map
            else:
                bm += bm_data.map
            j += 1
            # rotate the map
            # from scipy import ndimage
            # bmi = ndimage.rotate(bm_data.map, dk*180/np.pi, reshape=False)
        bm = bm/np.max(bm)
        bm_obj = BM()
        bm_obj.map = bm
        bm_obj.bm_x = bm_data.bm_x
        bm_obj.bm_y = bm_data.bm_y
        # save the beam pattern
        bm_obj.saveh5(filename=mmp.bmdir+'Output/FinalBM_{}_{}Hz.h5'.format(i, mmp.srcopt), opt=['map', 'bm_x', 'bm_y'])
        # cmap = plt.get_cmap('jet')
        # plt.imshow(bm, extent=[np.min(bm_data.bm_x)*180/np.pi, np.max(bm_data.bm_x)*180/np.pi, np.min(bm_data.bm_y)*180/np.pi, np.max(bm_data.bm_y)*180/np.pi, ], cmap=cmap)
        # plt.show()
        # compare with input beam
        bm_in_obj = beam_map_func(None, fp, mmp, i, sopt=True)
        bm_in = np.reshape(bm_in_obj.Ba, (mmp.ny, mmp.nx))
        bmplot(mmp, bm_in, bm)



if __name__=='__main__':
    # load the mission model parameters object
    mmp = Mission_Model_Parameters(src_center_opt=False)
    # load the focal plane object
    fp = FocalPlane()
    # set the dk angles, this angle should same to the tod simulation program
    dks = np.array([0, 45, 90, 135, 180])
    # dks = [0]
    # beam_charc(fp, mmp, dks)

    # beam_combination(fp, mmp, dks)
    beam_parameter_fit(fp, mmp)