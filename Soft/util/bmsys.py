import numpy as np
from .h5_util import H5IO
from .beam_characterization import BM
import matplotlib.pyplot as plt


class bmobj(H5IO):
    def __init__(self):
        pass


class BeamSysParams(object):
    def __init__(self, fp, mmp, sopt):
        lenfp = fp.lenfp
        self.dx = np.ones(lenfp)*mmp.sigma
        self.dy = np.zeros(lenfp)
        self.sigma = np.ones(lenfp) * mmp.sigma
        self.delta_g = np.ones(lenfp) * 0.
        self.delta_sigma = 0. * self.sigma
        self.delta_p = np.ones(lenfp) * 1.
        self.delta_c = np.ones(lenfp) * 0.
        self.nbin = mmp.nbin_bm
        return


def bm_load(filename, mmp):
    bm = BM()
    bm_map0 = np.loadtxt(filename, encoding='utf-16', skiprows=7).T
    # construct the x,y grid
    x0 = np.arange(-0.5, 0.5 + 0.01, 0.02) * np.pi / 180
    y0 = x0
    x0, y0 = np.meshgrid(x0, y0)
    # bm.bm_x, bm.bm_y = np.meshgrid(x, y)
    # do the interpolation
    x = mmp.xrange
    y = mmp.yrange
    # x, y = np.meshgrid(x, y)
    from scipy import interpolate
    f = interpolate.interp2d(x0, y0, bm_map0, kind='cubic')
    # bm_map = f(x.flatten(), y.flatten())
    bm_map = f(x, y)
    bm.Ba = bm_map/np.max(bm_map)
    bm.Bb = bm.Ba
    bm.bm_x = x
    bm.bm_y = y
    # plt.imshow(bm.Ba)
    # plt.show()
    return bm


def beam_map_func(bmfile, fp, mmp, chind, sopt=False):

    if bmfile is None:
        # if the beam map filename is None, we will genrate the beam pattern
        bmsysp = BeamSysParams(fp, mmp, sopt)
        dx = bmsysp.dx[chind]
        dy = bmsysp.dy[chind]
        dg = bmsysp.delta_g[chind]
        sigma = bmsysp.sigma[chind]
        delta_sigma = bmsysp.delta_sigma[chind]
        dp = bmsysp.delta_p[chind]
        dc = bmsysp.delta_c[chind]
        xc = 0
        yc = 0
        xca = xc-1./2*dx
        xcb = xc+1./2*dx
        yca = yc-1./2*dy
        ycb = yc+1./2*dy

        nbin = bmsysp.nbin
        ga = 1. / 2 + 1. / 2 * dg;
        gb = 1. / 2 - 1. / 2 * dg
        sigma_a = sigma + delta_sigma / 2.;
        sigma_b = sigma - delta_sigma / 2.
        pa = dp / 2.;
        pb = -dp / 2.
        ca = dc / 2.;
        cb = -dc / 2.
        x = mmp.xrange
        y = mmp.yrange
        x, y = np.meshgrid(x, y)

        Ba = 1. / (2 * np.pi * sigma_a ** 2) * (np.exp(-1. / (2 * sigma_a ** 2) * (
                    (x - (xca)) ** 2 / (1 + pa) + (y - (yca)) ** 2 / (1 - pa) - 2 * ca * (x - (xca)) * (y - (yca)))))
        Bb = 1. / (2 * np.pi * sigma_b ** 2) * (np.exp(-1. / (2 * sigma_b ** 2) * (
                (x - (xcb)) ** 2 / (1 + pb) + (y - (ycb)) ** 2 / (1 - pb) - 2 * cb * (x - (xcb)) * (y - (ycb)))))
        # Ba = np.exp(-1. / (2 * sigma_a ** 2) * ((x - (xca)) ** 2 + (y - (yca)) ** 2))
        # Bb = np.exp(-1. / (2 * sigma_b ** 2) * ((x - (xcb)) ** 2 + (y - (ycb)) ** 2))
        # Ba = Ba/np.sum(Ba)*ga; Bb = Bb/np.sum(Bb)*gb
        Bamax = Ba.max()
        Bbmax = Bb.max()
        Ba = Ba / Bamax
        Bb = Bb / Bbmax
    # Ba = Ba+np.random.normal(0, Bamax/50., np.shape(Ba))
    # Bb = Bb+np.random.normal(0, Bbmax/50., np.shape(Bb))
    else:
        # else, we will load the beam map
        bm = bm_load(bmfile, mmp)
        return bm
    pltopt = 0
    if pltopt:
        cmap = plt.get_cmap('jet')
        scale = 180 / np.pi
        plt.imshow(Ba, extent=[mmp.xmin*scale, mmp.xmax*scale, mmp.ymin*scale, mmp.ymax*scale], origin='lower', cmap=cmap)
        plt.show()
    # beamplot(x-xc, y-yc, Ba, Bb)
    beam_map = bmobj()
    beam_map.x = x
    beam_map.y = y
    # this r and theta are related to the boresight
    beam_map.Ba = Ba
    beam_map.Bb = Bb
    if sopt:
        filename = mmp.bmdir + 'Input/bm{}.h5'.format(chind)
        beam_map.saveh5(filename=filename, opt=['Ba'])
    return beam_map


def bm_func(bm_data, *args):
    bm_x = bm_data.bm_x.flatten()
    bm_y = bm_data.bm_y.flatten()
    A = args[0]
    # print(A)
    sigma = args[1]
    xc = args[2]
    yc = args[3]
    p = args[4]
    c = args[5]
    # xc = 0
    # yc = 0
    # p = 0
    # c = 0
    # B = A/(np.sqrt(2*np.pi)*sigma)*np.exp(-1. / (2 * sigma ** 2) * ((bm_x) ** 2 + (bm_y) ** 2))
    # B = A/(2*np.pi*sigma**2)*(np.exp(-1. / (2 * sigma ** 2) * ((bm_x - (xc)) ** 2 / (1 + p) + (bm_y - (yc)) ** 2 / (1 - p))) \
    # + np.exp(-1. / (2 * sigma ** 2) * ((bm_x - (xc)) ** 2 + (bm_y - (yc)) ** 2 - 2 * c * (bm_x - (xc)) * (bm_y - (yc)))))
    B = A / (2 * np.pi * sigma ** 2) * (np.exp(-1. / (2 * sigma ** 2) * (
                (bm_x - (xc)) ** 2 / (1 + p) + (bm_y - (yc)) ** 2 / (1 - p) - 2 * c * (bm_x - (xc)) * (bm_y - (yc)))))
    return B


def beamplot(x, y, Ba, Bb):
    import matplotlib.pyplot as plt
    vmax = max(Ba.max(), Bb.max())
    cmap = plt.get_cmap('jet')
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
    ax1.set_title('detector a')
    # ax1.imshow(Ba, extent=[x.min(), x.max(), y.min(), y.max()])
    # ax1.imshow(Ba, extent=[x.min()*180/np.pi, x.max()*180/np.pi, y.min()*180/np.pi, y.max()*180/np.pi], vmin=-vmax*0.8, vmax=vmax*0.8, cmap=cmap)
    ax1.imshow(Ba, extent=[x.min() * 180 / np.pi, x.max() * 180 / np.pi, y.min() * 180 / np.pi, y.max() * 180 / np.pi],
               cmap=cmap)
    ax2.set_title('detector b')
    # ax2.imshow(Bb, extent=[x.min(), x.max(), y.min(), y.max()])
    ax2.imshow(Bb, extent=[x.min() * 180 / np.pi, x.max() * 180 / np.pi, y.min() * 180 / np.pi, y.max() * 180 / np.pi],
               vmin=-vmax * 0.8, vmax=vmax * 0.8, cmap=cmap)
    ax3.set_title('diff')
    # ax3.imshow(Bb, extent=[x.min(), x.max(), y.min(), y.max()])
    im = ax3.imshow(Ba - Bb,
                    extent=[x.min() * 180 / np.pi, x.max() * 180 / np.pi, y.min() * 180 / np.pi, y.max() * 180 / np.pi],
                    vmin=-vmax * 0.8, vmax=vmax * 0.8, cmap=cmap)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.19, 0.015, 0.6])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()
    return


def bmplot(mmp, Ba, Bb):
    import matplotlib.pyplot as plt
    vmax = max(Ba.max(), Bb.max())
    scale = 0.1
    cmap = plt.get_cmap('jet')
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
    ax1.set_title('Input Beam')
    ax1.imshow(Ba, origin='lower',
               extent=[mmp.xmin * 180 / np.pi, mmp.xmax * 180 / np.pi, mmp.ymin * 180 / np.pi, mmp.ymax * 180 / np.pi],
               vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    # ax2.set_title('Fitted Beam')
    ax2.set_title('Fit Beam')
    ax2.imshow(Bb, origin='lower',
               extent=[mmp.xmin * 180 / np.pi, mmp.xmax * 180 / np.pi, mmp.ymin * 180 / np.pi, mmp.ymax * 180 / np.pi],
               vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    # ax2.imshow(Bb, extent=[x.min()*180/np.pi, x.max()*180/np.pi, y.min()*180/np.pi, y.max()*180/np.pi], vmin=-vmax*0.8, vmax=vmax*0.8, cmap=cmap)
    ax3.set_title('diff')
    # ax3.imshow(Bb, extent=[x.min(), x.max(), y.min(), y.max()])
    im = ax3.imshow(Ba - Bb, origin='lower',
                    extent=[mmp.xmin * 180 / np.pi, mmp.xmax * 180 / np.pi, mmp.ymin * 180 / np.pi,
                            mmp.ymax * 180 / np.pi],
                    vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.19, 0.015, 0.6])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()
    return
