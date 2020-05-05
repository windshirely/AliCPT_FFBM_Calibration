import numpy as np
from ctypes import *
import matplotlib.pyplot as plt
from .h5_util import H5IO
from .bs2px_func import src2bmcoord
from util.filter_func import filter_tod_v2, filter_tod

class BM(H5IO):
    def __init__(self):
        pass

def beam_mapmaking(tod, mmp, ch_idx, dk):
    """
    this function is to map making for each detector
    :return:
    """
    npix = mmp.ny * mmp.nx
    m = np.zeros(npix, dtype=np.float64)
    mw = np.zeros(npix, dtype=np.float64)
    mf = np.zeros(npix, dtype=np.float64)
    az = tod.px_altaz_az[ch_idx]
    el = tod.px_altaz_el[ch_idx]
    # plt.plot(az[0:2000]*180/np.pi)
    # plt.show()
    print(mmp.xmax*180/np.pi)
    ind = np.logical_and(az < mmp.xmax-mmp.xstep, np.logical_and(az > mmp.xmin+mmp.xstep, np.logical_and(el < mmp.ymax-mmp.ystep, el > mmp.ymin+mmp.ystep)))
    az1 = az[ind]
    el1 = el[ind]
    len_t = len(az1)
    print(len_t)
    w = np.ones(len_t, dtype=np.float64)
    lib = cdll.LoadLibrary('./clib/map_making.so')
    da = tod.da[ch_idx, ind]
    az1 = np.ascontiguousarray(az1, dtype=np.float64)
    el1 = np.ascontiguousarray(el1, dtype=np.float64)
    da = np.ascontiguousarray(da, dtype=np.float64)
    w = np.ascontiguousarray(w, dtype=np.float64)
    m = np.ascontiguousarray(m, dtype=np.float64)
    mw = np.ascontiguousarray(mw, dtype=np.float64)
    func = lib.rectang_mapmaking
    func(c_void_p(az1.ctypes.data), c_void_p(el1.ctypes.data), c_void_p(da.ctypes.data), c_void_p(w.ctypes.data), c_void_p(m.ctypes.data), c_void_p(mw.ctypes.data), \
            c_float(mmp.xmin), c_float(mmp.ymin), c_int(mmp.nx), c_int(mmp.ny), c_float(mmp.xstep),
            c_float(mmp.ystep), c_int(len_t)
    )
    idxt = mw!=0.
    mf[idxt] = m[idxt]/mw[idxt]
    mf = mf.reshape((mmp.ny, mmp.nx))
    # normalize the beam map using the max value
    mf = mf/np.max(mf)
    bm = BM()
    # bm.map = mf.flatten()
    bm.map = mf
    x, y = np.meshgrid(mmp.xrange, mmp.yrange)
    # bm.bm_x = x.flatten()-mmp.az_ref
    # bm.bm_y = y.flatten()-mmp.el_ref
    bm.bm_x = x-mmp.az_ref
    bm.bm_y = y-mmp.el_ref
    filename = mmp.bmdir+'Output/bm{}_dk{}_const.h5'.format(ch_idx, dk)
    bm.saveh5(filename=filename, opt=['map', 'bm_x', 'bm_y'])
    cmap = plt.get_cmap('jet')
    plt.imshow(mf, extent=[(mmp.xmin)*180/np.pi, (mmp.xmax)*180/np.pi, (mmp.ymin)*180/np.pi, (mmp.ymax)*180/np.pi, ], cmap=cmap)
    plt.show()
    # rotate the map
    # use the x_p ang y_p to rotate the map
    # from scipy import ndimage
    # mf_p = ndimage.rotate(mf, 45, reshape=False)
    # plt.imshow(mf_p, extent=[(mmp.xmin - mmp.az_ref) * 180 / np.pi, (mmp.xmax - mmp.az_ref) * 180 / np.pi,
    #                        (mmp.ymin - mmp.el_ref) * 180 / np.pi, (mmp.ymax - mmp.el_ref) * 180 / np.pi, ], cmap=cmap)
    # plt.show()
    return

def beam_mapmaking_on_detector(tod, mmp, fp, ch_idx, dk4filename):
    """
    make the beam map from the tod data into detector beam coordinate.
    :param tod: tod object
    :param mmp: mission model parameters object
    :param fp: focal plane object
    :param ch_idx: the detector channel index
    :param dk4filename: dk angle related to the filename
    :return: save the beam pattern into different h5 file
    """
    # numter of pixel in the beam pattern
    npix = mmp.ny * mmp.nx
    # initialize the beam sum map, map weight, and final map
    m = np.zeros(npix, dtype=np.float64)
    mw = np.zeros(npix, dtype=np.float64)
    mf = np.zeros(npix, dtype=np.float64)
    az = tod.px_altaz_az[ch_idx]
    el = tod.px_altaz_el[ch_idx]
    dk = tod.px_dk[ch_idx]
    da = np.copy(tod.da[ch_idx])
    # plot the tod data to check
    plt.plot(tod.da[ch_idx, int(4.5e6):int(5e6)])
    plt.show()
    # first remove the ground emission
    da = filter_tod(da, mmp)
    plt.plot(da[int(4.5e6):int(5e6)])
    plt.show()
    # we need demodulate the tod
    f = mmp.srcopt
    step = int(mmp.f_sample / f)
    idx0 = np.arange(len(da))
    idx1 = idx0 % step
    idx4demodul = idx1 < (int(step / 2))
    da = da[idx4demodul]
    az = az[idx4demodul]
    el = el[idx4demodul]
    dk = dk[idx4demodul]
    scale = 7
    ind = np.logical_and(az-mmp.az_ref < mmp.xmax-scale*mmp.xstep, np.logical_and(az-mmp.az_ref > mmp.xmin+scale*mmp.xstep, np.logical_and(el-mmp.el_ref < mmp.ymax-scale*mmp.ystep, el-mmp.el_ref > mmp.ymin+scale*mmp.ystep)))
    az1 = az[ind]
    el1 = el[ind]
    dk1 = dk[ind]
    # convert the src al el to beam coordinate
    theta = (fp.pxa_theta[ch_idx]+fp.pxb_theta[ch_idx])/2
    az_bmcoord, el_bmcoord = src2bmcoord(mmp, az1, el1, dk1, theta)
    len_t = len(az_bmcoord)
    w = np.ones(len_t)
    lib = cdll.LoadLibrary('./clib/map_making.so')
    # da0 = tod.da[ch_idx, ind]
    da = da[ind]
    az_bmcoord = np.ascontiguousarray(az_bmcoord, dtype=np.float64)
    el_bmcoord = np.ascontiguousarray(el_bmcoord, dtype=np.float64)
    da = np.ascontiguousarray(da, dtype=np.float64)
    # plt.scatter(el_bmcoord * 180 / np.pi, da)
    # plt.show()
    w = np.ascontiguousarray(w, dtype=np.float64)
    m = np.ascontiguousarray(m, dtype=np.float64)
    mw = np.ascontiguousarray(mw, dtype=np.float64)
    func = lib.rectang_mapmaking
    func(c_void_p(az_bmcoord.ctypes.data), c_void_p(el_bmcoord.ctypes.data), c_void_p(da.ctypes.data), c_void_p(w.ctypes.data), c_void_p(m.ctypes.data), c_void_p(mw.ctypes.data), \
            c_double(mmp.xmin), c_double(mmp.ymin), c_int(mmp.nx), c_int(mmp.ny), c_double(mmp.xstep),
            c_double(mmp.ystep), c_int(len_t)
    )
    idxt = mw>0.5
    mf[idxt] = m[idxt]/mw[idxt]
    mf = mf.reshape((mmp.ny, mmp.nx))
    # normalize the beam map using the max value
    mf = mf/np.max(mf)
    bm = BM()
    bm.map = mf
    x, y = np.meshgrid(mmp.xrange, mmp.yrange)
    bm.bm_x = x
    bm.bm_y = y
    filename = mmp.bmdir+'Output/bm{}_dk{}_{}.h5'.format(ch_idx, dk4filename, mmp.srcopt)
    bm.saveh5(filename=filename, opt=['map', 'bm_x', 'bm_y'])
    cmap = plt.get_cmap('jet')
    if 1:
        plt.imshow(mf, extent=[(mmp.xmin)*180/np.pi, (mmp.xmax)*180/np.pi, (mmp.ymin)*180/np.pi, (mmp.ymax)*180/np.pi, ], cmap=cmap)
        plt.show()
    return
