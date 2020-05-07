import numpy as np

class Scan_Model(object):
    """
    Contain all the parameters for scan strategy.
    """
    def __init__(self):
        ###########################################################################
        # Parameters for the target parameters
        ###########################################################################
        # number of scanset per day
        self.n_scanset = 4
        # the reference point of az, el as well as dk
        # option 1, we fix the ref point in the encoder coord
        self.az_ref = 22*np.pi/180
        self.el_ref = 8*np.pi/180
        self.dk = np.array([0., 45.]) * np.pi / 180
        ###########################################################################
        # Paramters for scan file generation
        ###########################################################################
        # ==============================================================================
        self.el_step = 0.02*np.pi/180
        self.n_steps = 400
        # the parameters for ces
        # the ces will split into a start part, several full scan parts, and a end part.
        self.acc_ces_v = 4.
        self.v_ces = 2.
        self.tacc_ces = self.v_ces / self.acc_ces_v
        self.tc1_ces = 7. # const velocity for peak left same of elnod
        self.tc2_ces = self.tc1_ces + 1. / 2 * self.tacc_ces # const velocity for peak right same of elnod
        self.n_fullscans = 3 # number of full scans
        ####################################################################################
        # other parameters
        self.p_rms = 1./60*np.pi/180 # root mean square for each encoder axis about 0.2arcmin.

class PinkNoiseParams(object):
    """
    this is for pink noise plus white noise, when f<fknee belong to pink noise, f>knee belong to white noise.
    """
    def __init__(self):
        self.wnl = 300./np.sqrt(2) # white noise level NET
        self.fknee_s = 1. # knee frequency sum
        self.alpha_s = 1.5 # index for sum
        self.fknee_d = 0.5 # knee frequency diff
        self.alpha_d = 0.5 # index for diff
        return


class Mission_Model_Parameters(Scan_Model, PinkNoiseParams):
    """
    It constains all the parameter for mission model, this class will contain all the configure parameters except the fp structure.
    """
    def __init__(self, src_center_opt=True):
        Scan_Model.__init__(self)
        PinkNoiseParams.__init__(self)
        # super(Scan_Model, self).__init__(self)
        # print(self.tacc_elnod)
        # super(PinkNoiseParams, self).__init__()
        self.mapoutopt = 'Healpix'
        self.f_sample = 200
        self.toddir = '../Data/TOD/'
        self.bmdir = '../Data/Map/'
        # self.srcopt = 'const' # describe the emission property of the source, can be 'const' or some Hz
        self.srcopt = 20
        # for healpix map nside, if we just want to use camb to simute map as our input map, this parameter will be used
        self.nside_in = 1024
        # map input
        self.map_in = '../map_input/PlanckMap/HFI_SkyMap_353_2048_R2.02_full.fits'
        # for site location
        self.lon = 80.
        self.lat = 32.
        self.height = 5250
        self.D = 0.72 #meter
        # frequncy
        self.freq = 95
        self.fwhm = 1.22*(3e8/(self.freq*1e9))/(self.D)
        self.sigma = self.fwhm/(2*np.sqrt(2*np.log(2)))
        # print(self.sigma*180/np.pi*60)
        # input()
        # other parameters for simulation
        self.bmsys_opt = False
        self.ndays = 1
        self.t_atm_zenith = 10 # KCMB
        # chunck the tod into different file
        self.chunk_size = 6 * 60  # 6 minutes
        self.nbin_bm = 50
        self.nx = self.nbin_bm
        self.ny = self.nbin_bm
        # number of sigma that the beam pattern range
        n_sigma = 12
        self.n_sigma = n_sigma
        self.xstep = n_sigma*self.sigma/self.nx
        self.ystep = n_sigma*self.sigma/self.ny
        if src_center_opt:
            self.xrange = np.arange(-n_sigma/2*self.sigma, n_sigma/2*self.sigma, self.xstep)+self.az_ref
            self.yrange = np.arange(-n_sigma/2 * self.sigma, n_sigma/2 * self.sigma, self.ystep)+self.el_ref
        else:
            self.xrange = np.arange(-n_sigma/2*self.sigma, n_sigma/2*self.sigma, self.xstep)
            self.yrange = np.arange(-n_sigma/2 * self.sigma, n_sigma/2 * self.sigma, self.ystep)
        self.xmin = np.min(self.xrange)
        self.xmax = np.max(self.xrange)
        self.ymin = np.min(self.yrange)
        self.ymax = np.max(self.yrange)
        return

class FocalPlane(object):
    """
    this class will be used to construct (simulate) the focal place sturcture,
    """
    def __init__(self):
        self.__a = 0.52
        self.__b = 0.8
        self.__fp_structure()
        self.__generate_fp()
        self.index = np.arange(len(self.px_r))

        self.pxa_r = self.px_r
        self.pxa_theta = self.px_theta
        self.pxa_chi = self.px_chisky
        self.pxb_r = self.px_r
        self.pxb_theta = self.px_theta
        self.pxb_chi = self.px_chisky+np.pi/2
        self.lenfp = len(self.pxa_r)

    def __fp_structure(self):
        # edge means the distance between the edge detectors to the edge;
        # dist mean the distance between the near two detectors.
        # distance between two line parallel to the edge
        a = self.__a
        b = self.__b
        AB = np.sqrt(3) * a + np.sqrt(3) / 3. * a + 11 * b
        BE = AB * np.sin(np.pi / 3.)
        AF = 2 * a
        DG = 2 * np.sqrt(3) / 3. * a
        c = np.sqrt(3) / 2. * b

        pixx = np.array([])
        pixy = np.array([])
        # for tile1
        tile1x = np.array([])
        tile1y = np.array([])
        # for tile1_1
        tile1_1x = np.array([])
        tile1_1y = np.array([])
        for i in range(12):
            x = a + i * c
            dd = a / (np.sqrt(3) / 2.)
            y0 = dd - np.sqrt(3) / 3 * x
            y = y0 + np.arange(12) * b
            x = np.ones(12) * x
            tile1_1x = np.append(tile1_1x, x)
            tile1_1y = np.append(tile1_1y, y)
        tile1x = np.append(tile1x, tile1_1x)
        tile1y = np.append(tile1y, tile1_1y)
        # for tile1_2
        tile1_2x = -tile1_1x
        tile1_2y = tile1_1y
        tile1x = np.append(tile1x, tile1_2x)
        tile1y = np.append(tile1y, tile1_2y)
        # for tile1_3
        tile1_3x = np.array([])
        tile1_3y = np.array([])
        # for middle line
        tile1_3x0 = np.zeros(12)
        tile1_3y0 = -DG - np.arange(12) * b
        tile1_3x = np.append(tile1_3x, tile1_3x0)
        tile1_3y = np.append(tile1_3y, tile1_3y0)
        # for right
        tile1_3xr = np.array([])
        tile1_3yr = np.array([])
        i = 0
        for j in np.arange(12, 0, -1):
            x = np.arange(1, j) * c
            y = tile1_3y0[i] - np.sqrt(3) / 3. * x
            i = i + 1
            tile1_3xr = np.append(tile1_3xr, x)
            tile1_3yr = np.append(tile1_3yr, y)
        tile1_3x = np.append(tile1_3x, tile1_3xr)
        tile1_3y = np.append(tile1_3y, tile1_3yr)
        # for left
        tile1_3xl = -tile1_3xr
        tile1_3yl = tile1_3yr
        tile1_3x = np.append(tile1_3x, tile1_3xl)
        tile1_3y = np.append(tile1_3y, tile1_3yl)

        # combine to tile1
        tile1x = np.append(tile1x, tile1_3x)
        tile1y = np.append(tile1y, tile1_3y)

        # combine into pixx and pixy
        pixx = np.append(pixx, tile1x)
        pixy = np.append(pixy, tile1y)

        # tile2 form tile1
        tile2x = tile1x + 2 * BE
        tile2y = tile1y
        # combine into pixx and pixy
        pixx = np.append(pixx, tile2x)
        pixy = np.append(pixy, tile2y)

        # tile3 from tile1
        tile3x = tile1x + BE
        tile3y = tile1y + 3. / 2 * AB
        # combine into pixx and pixy
        pixx = np.append(pixx, tile3x)
        pixy = np.append(pixy, tile3y)

        # tile4 from tile1
        tile4x = tile1x - BE
        tile4y = tile1y + 3. / 2 * AB
        # combine into pixx and pixy
        self.__pixx = np.append(pixx, tile4x)
        self.__pixy = np.append(pixy, tile4y)

    def __generate_fp(self):
        rfp, thetafp = cart2pol(self.__pixx, self.__pixy)
        rfov = 10.4 * np.pi / 180
        ratio = rfov / max(rfp)
        rsky = ratio * rfp
        thetasky = thetafp
        chisky = np.pi / 2 - thetasky  # relative to the radia direction
        self.px_r = rsky[[0,100]]
        self.px_theta = thetasky[[0,100]]
        self.px_chisky = chisky[[0,100]]
        print(self.px_r*180/np.pi, self.px_theta*180/np.pi)

class Beam_Sys_Parameters(object):
    def __init__(self, fp):
        lenfp = fp.lenfp
        self.sigma = np.ones(lenfp) * 30. / 60 * np.pi / 180 / (2 * np.sqrt(2 * np.log(2)))
        self.delta_x = 0.2 * np.pi / 180 * np.ones(lenfp)
        self.delta_y = 0. * self.sigma
        self.delta_g = np.ones(lenfp) * 0.
        self.delta_sigma = 0. * self.sigma
        self.delta_p = np.ones(lenfp) * 0.
        self.delta_c = np.ones(lenfp) * 0.
        self.nbin = 100


def cart2pol(x, y):
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    # here the theta generated are ccw, we neeed transform it to cw which can be used as input of reckon
    theta = np.pi * 2 - theta
    return r, theta

def pol2cart(r, theta):
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x, y

if __name__=='__main__':
    fp = FocalPlane()
    print(fp.fp)