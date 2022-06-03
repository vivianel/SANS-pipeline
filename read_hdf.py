# -*- coding: utf-8 -*-
"""
Spyder Editor

conda activate spyder-env

This is a temporary script file.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

file_open = 'C:/Users/lutzbueno_v/Documents/Beamtimes/2021-II/0940_DeAngelis/Data/sans2021n080413.hdf'

f = h5py.File(file_open, 'r')
img = np.array(f['entry1/data1/counts'])
img1 = np.where(img==0, 1e-4, img)
imgplot = plt.imshow(np.log(img1), clim=(0, 8), origin='lower')
imgplot.set_cmap('jet')
plt.colorbar()