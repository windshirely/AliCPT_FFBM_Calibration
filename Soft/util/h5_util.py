import h5py

class H5IO(object):
    def __init__(self):
        pass

    def saveh5(self, filename='', opt=''):
        """
        save the data object into hdf5 file
        :param data_tuple:
        :param filename:
        :param datatype:
        :return:
        """
        f = h5py.File(filename, 'w')
        if opt=='':
            keys = self.__dict__.keys()
        elif isinstance(opt, list):
            keys = opt
        else:
            keys = [opt]
        # save all the attributes into file
        for k in keys:
            # data = getattr(self, k)
            # if k!='t_mjd':
            #     # this is only for tod data
            #     import numpy as np
            #     data = data.astype(np.float32)
            f.create_dataset(k, data=getattr(self, k))
        f.close()

    # @time_consuming
    def loadh5(self, filename='', opt=''):
        """
        load the hdf5 file into a data object.
        :param filename:
        :param self:
        :param opt:
        :return:
        """
        f = h5py.File(filename, 'r')
        if opt=='':
            keys = f.keys()
        elif isinstance(opt, list):
            keys = opt
        else:
            keys = [opt]
        for k in keys:
            setattr(self, k, f[k][()])
        return self