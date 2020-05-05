from util.bmsys import beam_map_func
import numpy as np
from util.calibration_model import Mission_Model_Parameters, FocalPlane
import matplotlib.pyplot as plt

mmp = Mission_Model_Parameters()
fp = FocalPlane()
bm = beam_map_func('', fp, mmp, 1)
lenx = int(np.sqrt(len(bm.Ba)))
Ba = np.reshape(bm.Ba, (int(np.sqrt(len(bm.Ba))), int(np.sqrt(len(bm.Ba)))))
cmap = plt.get_cmap('jet')
scale = 180/np.pi
# plt.imshow(Ba, cmap=cmap)
# plt.show()
x0 = np.reshape(bm.x, (lenx, lenx))
y0 = np.reshape(bm.y, (lenx, lenx))
fig, axs = plt.subplots(1, 2, )
axs[0].contourf(x0, y0, Ba)
axs[0].imshow(Ba)
axs[0].set_aspect('equal')
# plt.contourf(x0, y0, Ba)
# plt.show()
# rotate the beam
# az1 = mmp.az_ref+(az1-mmp.az_ref)*np.cos(dk1)+(el1-mmp.el_ref)*np.sin(dk1)
# el1 = mmp.el_ref+(az1-mmp.az_ref)*np.sin(dk1)+(el1-mmp.el_ref)*np.cos(dk1)
dk = np.pi/3
x1 = bm.xc+(bm.x-bm.xc)*np.cos(dk)+(bm.y-bm.yc)*np.sin(dk)
y1 = bm.yc-(bm.x-bm.xc)*np.sin(dk)+(bm.y-bm.yc)*np.cos(dk)
# use the x1 and y1 to remap the beam

# x1 = x0*np.cos(dk)+y0*np.sin(dk)
# y1 = -x0*np.sin(dk)+y0*np.cos(dk)
x1 = np.reshape(x1, (lenx, lenx))
y1 = np.reshape(y1, (lenx, lenx))
# axs[1].contourf(x1, y1, Ba)
axs[1].imshow(Ba)
axs[1].set_aspect('equal')
# plt.contourf(x1, y1, Ba)
plt.show()
