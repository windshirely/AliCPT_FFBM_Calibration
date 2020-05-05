import numpy as np
from ctypes import *
import matplotlib.pyplot as plt

# @time_consuming
def bs2px(tod, fp, pltopt=False):
    """
    transform the boresight trajectory to pixel trajectory with the focal plane structure.
    here we use the c language to make it faster.
    and in altaz we do not calculate the psi angle, because in altaz we only do the ground subtraction only depend on the azimuth.
    :param tod: contans the boresight trajectory
    :param fp: focal plane structure object
    :param pltopt: whether we want to show the result for each pixel
    :return: tod object contains pixel trajectory
    """
    # for altaz we will not calculate the psi angle
    bs_keys = ['bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk']
    px_keys = ['pxa_altaz_az', 'pxa_altaz_el', 'pxb_altaz_az', 'pxb_altaz_el']
    for k in bs_keys:
        setattr(tod, k, np.ascontiguousarray(getattr(tod, k), dtype=np.float64))
    lent = len(getattr(tod, bs_keys[0]))
    # convert the focal plane data into
    for k in fp.__dict__.keys():
        if k!='lenfp':
            setattr(fp, k, np.ascontiguousarray(getattr(fp, k), dtype=np.float64))
    # initialize the px pointing
    for k in px_keys:
        setattr(tod, k, np.ascontiguousarray(np.zeros((fp.lenfp, lent)), dtype=np.float64))
    # lib = cdll.LoadLibrary('{}/bs2px.so'.format(os.path.abspath(os.path.dirname(__file__))))
    lib = cdll.LoadLibrary('./clib/bs2px.so')
    func = lib.bs2px_altaz
    func(c_void_p(tod.bs_altaz_az.ctypes.data), c_void_p(tod.bs_altaz_el.ctypes.data), c_void_p(tod.bs_altaz_dk.ctypes.data), \
         c_void_p(fp.pxa_r.ctypes.data), c_void_p(fp.pxa_theta.ctypes.data), \
         c_void_p(fp.pxb_r.ctypes.data), c_void_p(fp.pxb_theta.ctypes.data), \
         c_void_p(tod.pxa_altaz_az.ctypes.data), c_void_p(tod.pxa_altaz_el.ctypes.data), \
         c_void_p(tod.pxb_altaz_az.ctypes.data), c_void_p(tod.pxb_altaz_el.ctypes.data), \
         c_int(lent), c_int(fp.lenfp))
    tod.px_altaz_az = (tod.pxa_altaz_az + tod.pxb_altaz_az) / 2
    tod.px_altaz_el = (tod.pxa_altaz_el + tod.pxb_altaz_el) / 2
    tod.px_altaz_dk = tod.bs_altaz_dk
    if pltopt:
        plt.plot(tod.px_altaz_az[0]*180/np.pi)
        # plt.plot(tod.px_altaz_el[0] * 180 / np.pi)
        # plt.plot(tod.px_altaz_el[1] * 180 / np.pi)
        # plt.plot(tod.px_altaz_el[582]*180/np.pi)
        plt.xlim([0,4000])
        plt.show()
        # print(np.shape(tod.pxa_altaz_az))
    return


def bs2px_with_dk(tod, fp, pltopt=False):
    """
    transform the boresight trajectory to pixel trajectory with the focal plane structure.
    here we use the c language to make it faster.
    and in altaz we do not calculate the psi angle, because in altaz we only do the ground subtraction only depend on the azimuth.
    :param tod: contans the boresight trajectory
    :param fp: focal plane structure object
    :param pltopt: whether we want to show the result for each pixel
    :return: tod object contains pixel trajectory
    """
    # for altaz we will not calculate the psi angle
    bs_keys = ['bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk']
    px_keys = ['pxa_altaz_az', 'pxa_altaz_el', 'pxb_altaz_az', 'pxb_altaz_el', 'px_dk']
    for k in bs_keys:
        setattr(tod, k, np.ascontiguousarray(getattr(tod, k), dtype=np.float64))
    lent = len(getattr(tod, bs_keys[0]))
    # convert the focal plane data into
    for k in fp.__dict__.keys():
        if k!='lenfp':
            setattr(fp, k, np.ascontiguousarray(getattr(fp, k), dtype=np.float64))
    # initialize the px pointing
    for k in px_keys:
        setattr(tod, k, np.ascontiguousarray(np.zeros((fp.lenfp, lent)), dtype=np.float64))
    # lib = cdll.LoadLibrary('{}/bs2px.so'.format(os.path.abspath(os.path.dirname(__file__))))
    lib = cdll.LoadLibrary('./clib/bs2px.so')
    func = lib.bs2px_with_dk
    func(c_void_p(tod.bs_altaz_az.ctypes.data), c_void_p(tod.bs_altaz_el.ctypes.data), c_void_p(tod.bs_altaz_dk.ctypes.data), \
         c_void_p(fp.pxa_r.ctypes.data), c_void_p(fp.pxa_theta.ctypes.data), \
         c_void_p(fp.pxb_r.ctypes.data), c_void_p(fp.pxb_theta.ctypes.data), \
         c_void_p(tod.pxa_altaz_az.ctypes.data), c_void_p(tod.pxa_altaz_el.ctypes.data), \
         c_void_p(tod.pxb_altaz_az.ctypes.data), c_void_p(tod.pxb_altaz_el.ctypes.data), c_void_p(tod.px_dk.ctypes.data),\
         c_int(lent), c_int(fp.lenfp))
    tod.px_altaz_az = (tod.pxa_altaz_az + tod.pxb_altaz_az) / 2
    tod.px_altaz_el = (tod.pxa_altaz_el + tod.pxb_altaz_el) / 2
    if pltopt:
        plt.plot(tod.px_altaz_az[0]*180/np.pi)
        # plt.plot(tod.px_altaz_el[0] * 180 / np.pi)
        # plt.plot(tod.px_altaz_el[1] * 180 / np.pi)
        # plt.plot(tod.px_altaz_el[582]*180/np.pi)
        plt.xlim([0,4000])
        plt.show()
        # print(np.shape(tod.pxa_altaz_az))
    return

def bs2px4dadb(tod, fp):
    """
    calculate the boresight to pixel in a python version, this version is memory consumming....
    :param tod: tod object
    :param fp: focal plane object
    :return: modified tod including the position for each detector pixel
    """
    lent = len(tod.bs_altaz_dk)
    # for detector a
    zbp_a = np.tile(tod.bs_altaz_dk, (fp.lenfp, 1)) + (np.tile(fp.pxa_theta, (lent, 1))).T - np.pi / 2
    bp_a = (np.tile(fp.pxa_r, (lent, 1))).T
    bz_a = np.pi/2-tod.bs_altaz_el
    # calculaion
    pz_a, pzb_a, bpz_a = bs2px_util(zbp_a, bp_a, bz_a)
    # set back to az, el and dk
    tod.pxa_altaz_el = np.pi / 2 - pz_a
    tod.pxa_altaz_az = tod.bs_altaz_az - pzb_a
    tod.pxa_dk = np.pi-bpz_a
    # ==================================================
    zbp_b = np.tile(tod.bs_altaz_dk, (fp.lenfp, 1)) + (np.tile(fp.pxb_theta, (lent, 1))).T - np.pi / 2
    bp_b = (np.tile(fp.pxb_r, (lent, 1))).T
    bz_b = np.pi / 2 - tod.bs_altaz_el
    pz_b, pzb_b, bpz_b = bs2px_util(zbp_b, bp_b, bz_b)
    tod.pxb_altaz_el = np.pi / 2 - pz_b
    tod.pxb_altaz_az = tod.bs_altaz_az - pzb_b
    tod.pxb_dk = np.pi - bpz_b
    tod.px_altaz_az = (tod.pxa_altaz_az+tod.pxb_altaz_az)/2.
    tod.px_altaz_el = (tod.pxa_altaz_el+tod.pxb_altaz_el)/2.
    tod.px_dk = (tod.pxa_dk+tod.pxb_dk)/2.
    return

def bs2px_util(zbp, bp, bz):
    """
    do the calculation in spherical triangle
    :param zbp: angle
    :param bp: spherical radia angle
    :param bz: spherical radia angle
    :return: three angles
    """
    cos_pz = np.cos(bp) * np.cos(bz) + np.sin(bp) * np.sin(bz) * np.cos(zbp)
    pz = np.arccos(cos_pz)
    sin_pzb = np.sin(zbp) / np.sin(pz) * np.sin(bp)
    cos_pzb = (np.cos(bp) - np.cos(bz) * cos_pz) / (np.sin(bz) * np.sin(pz))
    pzb = np.arctan2(sin_pzb, cos_pzb)
    sin_bpz = np.sin(zbp) / np.sin(pz) * np.sin(bz)
    cos_bpz = (np.cos(bz) - np.cos(bp) * np.cos(pz)) / (np.sin(bp) * np.sin(pz))
    bpz = np.arctan2(sin_bpz, cos_bpz)
    return pz, pzb, bpz

def bs2px_bm_v2(tod, bm):
    """
    We also can use the python version, but this is very memory consuming, and it is very slow.
    :param tod: tod object
    :param bm: beam object
    :return: modified tod object.
    """
    lent = len(tod.bs_altaz_dk)
    len_bm = np.size(bm.r)
    zbp = np.tile(tod.bs_altaz_dk, (len_bm, 1)) + (np.tile(bm.theta, (lent, 1))).T - np.pi / 2
    print(np.shape(zbp))
    input()
    bp = (np.tile(bm.r, (lent, 1))).T
    bz = np.pi / 2 - tod.bs_altaz_el
    pz, pzb, bpz = bs2px_util(zbp, bp, bz)
    tod.bm_az = tod.bs_altaz_az - pzb
    tod.bm_el = np.pi / 2 - pz
    return

def src2bmcoord(mmp, px_az, px_el, px_dk, theta):
    """
    this fucntion is to transform the source in the sky to the position is the beam coordinate.
    We use the spherical triangle formula.
    Positive angle direction is clockwise.
    :param mmp: mission model parmeters object
    :param px_az: the az position of detector in sky coordinate
    :param px_el: the el position of detector in sky coordinate
    :param px_dk: the dk position of detector in sky coordinate
    :param theta: the theta of the detector in focal plane
    :return: az and el in beam coordiante
    """
    # in szp spheriacal triangle
    szp = mmp.az_ref-px_az
    zp = np.pi/2-px_el
    sz = np.pi/2-mmp.el_ref
    cos_sp = np.cos(sz)*np.cos(zp)+np.sin(sz)*np.sin(zp)*np.cos(szp)
    sp = np.arccos(cos_sp)
    sin_spz = np.sin(szp)/np.sin(sp)*np.sin(sz)
    cos_spz = (np.cos(sz)-np.cos(zp)*np.cos(sp))/(np.sin(zp)*np.sin(sp))
    spz = np.arctan2(sin_spz, cos_spz)
    # in spf spheriacal triangle
    spf = theta-px_dk-spz
    sin_sf = np.sin(spf)*np.sin(sp)
    sf = np.arcsin(sin_sf)
    el_bmcoord = sf
    cos_fp = np.cos(sp)/np.cos(sf)
    sin_fp = (np.cos(sf)-np.cos(sp)*cos_fp)/(np.sin(sp)*np.cos(spf))
    az_bmcoord = np.arctan2(sin_fp, cos_fp)
    return az_bmcoord, el_bmcoord




def bs2px_bm(tod, bm, pltopt=False):
    """
    transform the boresight trajectory to pixel trajectory with the focal plane structure.
    here we use the c language to make it faster.
    and in altaz we do not calculate the psi angle, because in altaz we only do the ground subtraction only depend on the azimuth.
    :param tod: contans the boresight trajectory
    :param fp: focal plane structure object
    :param pltopt: whether we want to show the result for each pixel
    :return: tod object contains pixel trajectory
    """
    bm.length = len(bm.r)
    # for altaz we will not calculate the psi angle
    bs_keys = ['bs_altaz_az', 'bs_altaz_el', 'bs_altaz_dk']
    px_keys = ['bm_az', 'bm_el',]
    for k in bs_keys:
        setattr(tod, k, np.ascontiguousarray(getattr(tod, k), dtype=np.float64))
    lent = len(getattr(tod, bs_keys[0]))
    # convert the focal plane data into
    bm.r = np.ascontiguousarray(bm.r, dtype=np.float64)
    bm.theta = np.ascontiguousarray(bm.theta, dtype=np.float64)
    # initialize the px pointing
    for k in px_keys:
        setattr(tod, k, np.ascontiguousarray(np.zeros((bm.length, lent)), dtype=np.float64))
    # lib = cdll.LoadLibrary('{}/bs2px.so'.format(os.path.abspath(os.path.dirname(__file__))))
    lib = cdll.LoadLibrary('./clib/bs2px.so')
    func = lib.bs2px_bm
    func(c_void_p(tod.bs_altaz_az.ctypes.data), c_void_p(tod.bs_altaz_el.ctypes.data), c_void_p(tod.bs_altaz_dk.ctypes.data), \
         c_void_p(bm.r.ctypes.data), c_void_p(bm.theta.ctypes.data), \
         c_void_p(tod.bm_az.ctypes.data), c_void_p(tod.bm_el.ctypes.data), \
         c_int(lent), c_int(bm.length))
    if pltopt:
        plt.plot(tod.bs_altaz_el*180/np.pi)
        plt.plot(tod.bm_el[0]*180/np.pi)
        plt.plot(tod.bm_el[-1] * 180 / np.pi)
        # plt.xlim([0,4000])
        plt.show()
        # print(np.shape(tod.pxa_altaz_az))
    return

