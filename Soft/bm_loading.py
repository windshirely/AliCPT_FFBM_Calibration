import numpy as np
import matplotlib.pyplot as plt
from util.beam_characterization import BM
from ffbm_analysis import bmparams_fit
from util.bmsys import bm_func, bmplot

def bm_load(filename):
    bm = BM()
    bm.map = np.loadtxt(filename, encoding='utf-16', skiprows=7).T
    # print(np.shape(bm.map))
    # construct the x,y grid
    x = np.arange(-0.5, 0.5+0.01, 0.02)*np.pi/180
    y = x
    bm.bm_x, bm.bm_y = np.meshgrid(x, y)
    # plt.imshow(bm.map, extent=[-0.5, 0.5, -0.5, 0.5])
    # plt.show()
    return bm


if __name__=='__main__':
    # filename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_17.4355.txt'
    # filename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_8.9243.txt'
    # filename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_13.2542.txt'
    # filename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_0.0000.txt'
    filename = '../Data/Map/Input/BM_Map/POPBeammap_0.0000_4.4894.txt'
    bm = bm_load(filename)
    inmap = bm.map
    pv, pcov = bmparams_fit(bm)
    # print(pv)
    print('fwhm=', pv[1]*180/np.pi*60*(2*np.sqrt(2*np.log(2))))
    print('xc=', pv[2]*180/np.pi)
    print('yc=', pv[3] * 180 / np.pi)
    # print(pcov)
    fitmap = np.reshape(bm_func(bm, *pv), (51, 51))
    # print(fitmap)
    # input()
    vmax = max(inmap.max(), fitmap.max())
    scale = 0.1
    cmap = plt.get_cmap('jet')
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
    ax1.set_title('Input Beam')
    # ax1.imshow(Ba, extent=[x.min(), x.max(), y.min(), y.max()])
    # ax1.imshow(Ba, extent=[x.min()*180/np.pi, x.max()*180/np.pi, y.min()*180/np.pi, y.max()*180/np.pi], vmin=-vmax*0.8, vmax=vmax*0.8, cmap=cmap)
    ax1.imshow(inmap,
               extent=[-0.5, 0.5, -0.5, 0.5], origin='lower',
               vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    # ax2.set_title('Fitted Beam')
    ax2.set_title('Fit Beam')
    ax2.imshow(fitmap,
               extent=[-0.5, 0.5, -0.5, 0.5], origin='lower',
               vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    # ax2.imshow(Bb, extent=[x.min()*180/np.pi, x.max()*180/np.pi, y.min()*180/np.pi, y.max()*180/np.pi], vmin=-vmax*0.8, vmax=vmax*0.8, cmap=cmap)
    ax3.set_title('diff')
    # ax3.imshow(Bb, extent=[x.min(), x.max(), y.min(), y.max()])
    im = ax3.imshow(inmap-fitmap,
                    extent=[-0.5, 0.5, -0.5, 0.5], origin='lower',
                    vmin=-vmax * scale, vmax=vmax * scale, cmap=cmap)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.19, 0.015, 0.6])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()